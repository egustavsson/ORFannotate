import os
import re
import shutil
import logging
import gffutils
import subprocess
import pandas as pd
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from rich.text import Text
from rich.progress import (
    TaskProgressColumn,
    TimeElapsedColumn,
)


def _keep_tx_and_exons(feat):
    """Filter GTF features to only keep transcript and exon rows."""
    return feat if feat.featuretype in {"transcript", "exon"} else False

def _load_transcript_sequences(transcript_fasta: Path):
    """Load transcript sequences into a dictionary keyed by transcript ID."""
    return SeqIO.to_dict(SeqIO.parse(str(transcript_fasta), "fasta"))

def _count_unique_transcripts(path):
    cmd = f"awk '$3 == \"transcript\" || $3 == \"mRNA\" || $3 == \"lnc_RNA\" || $3 == \"pseudogenic_transcript\" || $3 == \"miRNA\" || $3 == \"mRNA\" || $3 == \"snRNA\" || $3 == \"transcript\" || $3 == \"unconfirmed_transcript\" || $3 == \"snoRNA\" || $3 == \"scRNA\" || $3 == \"rRNA\" || $3 == \"processed_transcript\" || $3 == \"tRNA\" {{print}}' {path} | wc -l"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return int(result.stdout.strip())

def _validate_inputs(gtf_path: str, genome_fasta: str):
    """Check that required input files exist and are valid."""
    if not os.path.isfile(gtf_path):
        raise FileNotFoundError(f"GTF file not found or invalid: {gtf_path}")
    if not os.path.isfile(genome_fasta):
        raise FileNotFoundError(f"FASTA file not found or invalid: {genome_fasta}")

    has_transcripts = has_exons = False
    with open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if "\ttranscript\t" in line or "\tmRNA\t" in line:
                has_transcripts = True
            elif "\texon\t" in line:
                has_exons = True
            if has_transcripts and has_exons:
                break
    if not has_transcripts or not has_exons:
        raise ValueError("GTF file must contain both 'transcript' and 'exon' features.")

def _map_junctions_to_tx(junctions_genomic, exons):
    """Map genomic splice junction coordinates to transcript coordinates."""
    offsets = []
    cum_len = 0
    for exon in exons:
        offsets.append((exon.start, exon.end, cum_len))
        cum_len += exon.end - exon.start + 1

    mapped = []
    for donor, acceptor in junctions_genomic:
        donor_tx = acceptor_tx = None
        for start, end, offset in offsets:
            if start <= donor <= end:
                donor_tx = offset + (donor - start + 1)
            if start <= acceptor <= end:
                acceptor_tx = offset + (acceptor - start + 1)
        if donor_tx and acceptor_tx:
            mapped.append((donor_tx, acceptor_tx))
    return mapped

def _create_db(gtf_path, only_exons=True):
    """Create a DB using GFF utils."""
    annotated_db = gffutils.create_db(
            str(gtf_path) ,
            dbfn=":memory:",
            force=True,
            keep_order=True,
            disable_infer_transcripts=True,
            disable_infer_genes=True,
            sort_attribute_values=True,
            merge_strategy="create_unique",
            transform=_keep_tx_and_exons if only_exons else None,
            pragmas={
                "journal_mode": "OFF",
                "synchronous": "OFF",
                "temp_store": "MEMORY",
                "cache_size": -1000000
            },
            id_spec=None
        )
    return annotated_db

def _print_top_database_info(db):
    """
    Print the first 5 rows of the gffutils database for debugging.
    :param db: Database object created by gffutils
    """

    conn = db.conn  # connection gffutils uses internally
    
    try:
        rows = conn.execute("SELECT * FROM features LIMIT 10").fetchall()
    except Exception as e:
        print("Error querying database:", e)
        raise
    if not rows:
        raise Exception("Database is empty or has no entries in 'features' table.")
    else:
        df = pd.DataFrame([dict(row) for row in rows])
        
    print(df)

def _make_slim_gtf(in_gtf, only_exons=True, essential_keys=None):
    """
    Handles both GTF and GFF3 inputs.
    1. Converts GFF3 to GFF format, namely GFF3 key=value to GTF key "value".
    2. Keeps only essential attributes (gene_id, transcript_id by default, but you can pass more keys on the 'essential_kesy' param).
    3. Preserves all coordinates, strand, seqid, etc.
    4. Returns the path to the slim GTF (tmp/slim.gtf).
    
    :param in_gtf: GTF or GFF3 file
    :param essential_keys: currently only gene_id and transcript_id (and optionally gene_name, etc.)
    """   

    if essential_keys is None:
        essential_keys = {"gene_id", "transcript_id"}

    path = Path('./tmp')
    path.mkdir(parents=True, exist_ok=True)

    # Convert GFF3 -> GTF if needed
    in_path = Path(in_gtf)
    if in_path.suffix.lower() in {".gff3", ".gff"}:
        converted_gtf = path / "converted.gtf"

        subprocess.run(["gffread", str(in_path), "-T", "-o", str(converted_gtf)], check=True,stderr=subprocess.DEVNULL)
        converted_gtf_clean = path / "converted.cleaned.gtf"

        with open(converted_gtf, "rt", encoding="utf-8") as rin, open(converted_gtf_clean, "wt", encoding="utf-8") as rout:
            for line in rin:
                if line.startswith("#") or not line.strip():
                    rout.write(line)
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    rout.write(line)
                    continue
                attrs = cols[8]
                parsed = _parse_gtf_attributes(attrs)  # will already strip gene:/transcript:
                # Rebuild attributes: prefer gene_id then transcript_id then other keys
                order = []
                for k in ("gene_id", "transcript_id"):
                    if k in parsed and parsed[k]:
                        order.append(f'{k} "{parsed[k]}"')
                # append any others (optional)
                for k, v in parsed.items():
                    if k not in {"gene_id", "transcript_id"} and v:
                        order.append(f'{k} "{v}"')
                cols[8] = "; ".join(order) + ";" if order else "."
                rout.write("\t".join(cols) + "\n")
        in_gtf_to_use = converted_gtf_clean

    else:
        in_gtf_to_use = in_path

    slim_gtf_path = path / "slim.gtf"

    with open(in_gtf_to_use) as inp, open(slim_gtf_path, "w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            ftype = fields[2]
            if only_exons and ftype not in {"transcript", "exon"}:
                continue  # skip everything else

            attrs_dict = _parse_gtf_attributes(fields[8])

            new_attrs = []
            for key in essential_keys:
                if key in attrs_dict:
                    new_attrs.append(f'{key} "{attrs_dict[key]}"')
            fields[8] = "; ".join(new_attrs) + ";"

            out.write("\t".join(fields) + "\n")

    return slim_gtf_path

def _parse_gtf_attributes(attr_str):
    """Parse GTF attributes into a dict, stripping 'gene:' and 'transcript:' prefixes."""
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue

        if "=" in part:
            # handle GFF3 style
            key, val = part.split("=", 1)
            val = val.strip().strip('"')   # <-- FIX: remove quotes if present

        elif '"' in part:
            # handle GTF style
            key, val = part.split('"', 1)
            val = val.split('"')[0]

        else:
            continue

        key = key.strip()
        val = val.strip()

        # remove gene: and transcript: prefixes
        val = re.sub(r'^(?:gene:|transcript:)', '', val)

        attrs[key] = val

    return attrs
    
def _cleanup_tmp_folder():
    """
    Removes the local ./tmp folder with all its contents
    """
    tmp_dir = Path("./tmp")

    try:
        shutil.rmtree(tmp_dir)
    except FileNotFoundError:
        pass

# Logging filter: only keep selected messages on console
class _SelectiveConsoleFilter(logging.Filter):
    ALLOWED_TERMS = [
        "Starting ORFannotate",
        "Found",
        "Step 1:",
        "Step 2:",
        "Step 3:",
        "Step 4:",
        "Step 5:",
        "Step 6:",
        "Step 7:",
        "Processed",
        "ORFannotate completed"
    ]
    def filter(self, record):
        return any(term in record.getMessage() for term in self.ALLOWED_TERMS)
    
# Logging setup
def _setup_logging(output_dir):
    log_path = os.path.join(output_dir, "ORFannotate.log")

    file_formatter = logging.Formatter(
        fmt='[%(asctime)s] [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_formatter = logging.Formatter('[%(levelname)s] %(message)s')

    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(file_formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(console_formatter)
    stream_handler.addFilter(_SelectiveConsoleFilter())

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers = []
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger

## Define custom progress bar columns
class _ColoredTimeElapsedColumn(TimeElapsedColumn):
    """A custom TimeElapsedColumn with colored text - by default bright_yellow."""
    def __init__(self, style="bright_yellow"):
        super().__init__()
        self._style = style

    def render(self, task):
        base = super().render(task)
        return Text(str(base), style=self._style)
    
class _ColoredTaskProgressColumn(TaskProgressColumn):
    """A custom TaskProgressColumn with colored text - by default bright_magenta."""
    def __init__(self, style="bright_magenta"):
        super().__init__()
        self._style = style

    def render(self, task):
        base = super().render(task)
        return Text(str(base), style=self._style)