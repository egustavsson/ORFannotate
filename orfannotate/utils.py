import os
import re
import logging
import gffutils
from pathlib import Path
from turtle import pd
from Bio import SeqIO
from rich.text import Text
from rich.progress import (
    TaskProgressColumn,
    TimeElapsedColumn,
)

TX_RE = re.compile(r'transcript_id\s+"([^"]+)"')


def _keep_tx_and_exons(feat):
    """Filter GTF features to only keep transcript and exon rows."""
    return feat if feat.featuretype in {"transcript", "exon"} else False


def _load_transcript_sequences(transcript_fasta: Path):
    """Load transcript sequences into a dictionary keyed by transcript ID."""
    return SeqIO.to_dict(SeqIO.parse(str(transcript_fasta), "fasta"))


# Collect splice junctions into memory
def _collect_junctions(gtf_db):
    junctions = {}
    for tx in gtf_db.features_of_type("transcript"):
        exons = list(gtf_db.children(tx, featuretype="exon", order_by="start"))
        if len(exons) > 1:
            junctions[tx.id] = [(exons[i].end, exons[i+1].start) for i in range(len(exons)-1)]
    return junctions


def _count_unique_transcripts(gtf_path: str) -> int:
    """Count unique transcript IDs in a GTF file."""
    ids = set()
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#") or "\ttranscript\t" not in line:
                continue
            m = TX_RE.search(line)
            if m:
                ids.add(m.group(1))
    return len(ids)

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
            if "\ttranscript\t" in line:
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
            gtf_path,
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
                "chache_size": -1000000
            },
            id_spec=None
        )

    return annotated_db

def _print_top_database_info(db):
    
    print("Connecting....")
    conn = db.conn  # connection gffutils uses internally
    
    print("Connection stablished")
    
    rows=conn.execute("SELECT * FROM features LIMIT 5").fetchall()
    df = pd.DataFrame([dict(row) for row in rows])
    print(df.head(100))

   

def _make_slim_gtf(in_gtf, essential_keys=None):
    """
    Create a slim GTF containing only required GTF attributes:
    - gene_id
    - transcript_id
    - exon_number (optional)
    preserving all seqid, start, end, strand columns.
    """
    if essential_keys is None:
        essential_keys = {"gene_id", "transcript_id", "exon_number"}

    path = Path('./tmp')
    path.mkdir(parents=True, exist_ok=True)
    with open(in_gtf) as inp, open(Path(path / "slim.gtf"), "w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            attrs = fields[8]

            new_attrs = []
            for item in attrs.split(";"):
                item = item.strip()
                if not item:
                    continue
                key = item.split(" ")[0]
                if key in essential_keys:
                    new_attrs.append(item)

            fields[8] = "; ".join(new_attrs) + ";"
            out.write("\t".join(fields) + "\n")
    
def _cleanup_slim_gtf(slim_file):
    """Removes the slim GTF file and its directory if empty."""
    slim_file = Path(slim_file)
    folder = slim_file.parent

    # Delete file
    if slim_file.exists():
        slim_file.unlink()

    # Delete folder ONLY if empty
    try:
        folder.rmdir()
    except OSError:
        pass  # folder not empty → ignore

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