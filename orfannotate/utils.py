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
from rich.progress import TaskProgressColumn, TimeElapsedColumn

# Common GFF3/GTF transcript-like feature types.
# This set is used to normalise transcript-like features to "transcript" during slimming,
# and also to count transcript-like rows for progress reporting and input validation.
TRANSCRIPT_LIKE_FEATURES = {
    "transcript",
    "mRNA",
    "lnc_RNA",
    "pseudogenic_transcript",
    "miRNA",
    "snRNA",
    "unconfirmed_transcript",
    "snoRNA",
    "scRNA",
    "rRNA",
    "processed_transcript",
    "tRNA",
    "ncRNA",
}

# CPAT uppercases IDs  but some GTFs have lowercase transcript_ids, so we normalise to uppercase for lookups.
def _normalize_tid(tid):
    if tid is None:
        return None
    return str(tid).upper()


def _check_seqid_compatibility(db, genome, max_missing_examples: int = 20):
    """
    Check whether all seqids (contig/chromosome names) referenced in the annotation
    exist in the FASTA. If any are missing, raise an informative error listing
    examples of missing seqids.
    """
    logger = logging.getLogger(__name__)

    annotation_seqids = {row[0] for row in db.execute("SELECT DISTINCT seqid FROM features")}
    fasta_seqids = set(genome.keys())

    if not annotation_seqids:
        msg = "No seqids were found in the annotation database. Is the annotation empty or malformed?"
        logger.error(msg)
        raise ValueError(msg)

    missing = annotation_seqids - fasta_seqids
    if missing:
        example_fasta = sorted(list(fasta_seqids))[:10]

        msg = (
            "Sequence ID mismatch between annotation and FASTA.\n"
            f"Missing annotation seqids not found in FASTA: {missing}\n"
            f"FASTA seqids: {fasta_seqids}\n"
            f"Total missing seqids: {len(missing)} (out of {len(annotation_seqids)} in the GTF annotation).\n"
            "One or more contigs present in the annotation were not found in the FASTA. "
            "This usually means the files come from different genome builds or use different "
            "naming conventions (for example '3' vs 'chr3'). Please ensure both files originate "
            "from the same source and release."
        )

        logger.error(msg)
        raise ValueError(msg)


def _keep_tx_and_exons(feat):
    """Filter features to only keep transcript and exon rows (after normalisation)."""
    return feat if feat.featuretype in {"transcript", "exon"} else False


def _load_transcript_sequences(transcript_fasta: Path):
    """Load transcript sequences into a dictionary keyed by transcript ID."""
    return SeqIO.to_dict(SeqIO.parse(str(transcript_fasta), "fasta"))


def _count_unique_transcripts(path):
    """Count transcript-like feature rows (used for progress reporting)."""
    count = 0
    with open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 2 and fields[2] in TRANSCRIPT_LIKE_FEATURES:
                count += 1
    return count


def _validate_inputs(gtf_path: str, genome_fasta: str):
    """Check that required input files exist and contain transcript-like + exon features."""
    if not os.path.isfile(gtf_path):
        raise FileNotFoundError(f"GTF/GFF file not found or invalid: {gtf_path}")
    if not os.path.isfile(genome_fasta):
        raise FileNotFoundError(f"FASTA file not found or invalid: {genome_fasta}")

    has_transcripts = False
    has_exons = False

    with open(gtf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue

            ftype = fields[2]
            if ftype in TRANSCRIPT_LIKE_FEATURES:
                has_transcripts = True
            elif ftype == "exon":
                has_exons = True

            if has_transcripts and has_exons:
                break

    if not has_transcripts or not has_exons:
        raise ValueError(
            "Annotation must contain exon features and transcript-like features "
            "(for example 'transcript' or 'mRNA')."
        )


def _map_junctions_to_tx(junctions_genomic, exons):
    """Map genomic splice junction coordinates to transcript coordinates."""
    offsets = []
    cum_len = 0
    for exon in exons:
        offsets.append((exon.start, exon.end, cum_len))
        cum_len += exon.end - exon.start + 1

    mapped = []
    for donor, acceptor in junctions_genomic:
        donor_tx = None
        acceptor_tx = None
        for start, end, offset in offsets:
            if start <= donor <= end:
                donor_tx = offset + (donor - start + 1)
            if start <= acceptor <= end:
                acceptor_tx = offset + (acceptor - start + 1)
        if donor_tx is not None and acceptor_tx is not None:
            mapped.append((donor_tx, acceptor_tx))
    return mapped


def _create_db(gtf_path, only_exons=True):
    """Create a gffutils DB (in-memory)."""
    annotated_db = gffutils.create_db(
        str(gtf_path),
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
            "cache_size": -1000000,
        },
        id_spec=None,
    )
    return annotated_db


def _print_top_database_info(db):
    """
    Print the first 10 rows of the gffutils database for debugging.
    :param db: Database object created by gffutils
    """
    conn = db.conn  # connection gffutils uses internally

    try:
        rows = conn.execute("SELECT * FROM features LIMIT 10").fetchall()
    except Exception as e:
        print("Error querying database:", e)
        raise

    if not rows:
        raise RuntimeError("Database is empty or has no entries in 'features' table.")

    df = pd.DataFrame([dict(row) for row in rows])
    print(df)


def _make_slim_gtf(in_gtf, only_exons=True, essential_keys=None):
    """
    Handles both GTF and GFF3 inputs.
    1. Converts GFF3/GFF to GTF (via gffread -T).
    2. Rewrites attributes into a minimal GTF-like key "value"; format.
    3. Normalises transcript-like feature types (e.g. mRNA/ncRNA) to "transcript".
    4. Keeps only transcript/exon rows (if only_exons=True).
    5. Returns the path to the slim GTF (./tmp/slim.gtf).

    :param in_gtf: GTF or GFF3 file
    :param essential_keys: keys to keep (default: gene_id, transcript_id)
    """
    if essential_keys is None:
        essential_keys = {"gene_id", "transcript_id"}

    tmp_dir = Path("./tmp")
    tmp_dir.mkdir(parents=True, exist_ok=True)

    in_path = Path(in_gtf)

    # Convert GFF3/GFF -> GTF if needed
    if in_path.suffix.lower() in {".gff3", ".gff"}:
        converted_gtf = tmp_dir / "converted.gtf"

        res = subprocess.run(
            ["gffread", str(in_path), "-T", "-o", str(converted_gtf)],
            capture_output=True,
            text=True,
        )
        if res.returncode != 0:
            raise subprocess.CalledProcessError(
                res.returncode, res.args, output=res.stdout, stderr=res.stderr
            )
        if res.stderr:
            logging.getLogger(__name__).warning("gffread warnings:\n%s", res.stderr.strip())

        converted_gtf_clean = tmp_dir / "converted.cleaned.gtf"

        with open(converted_gtf, "rt", encoding="utf-8") as rin, open(
            converted_gtf_clean, "wt", encoding="utf-8"
        ) as rout:
            for line in rin:
                if line.startswith("#") or not line.strip():
                    rout.write(line)
                    continue

                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    rout.write(line)
                    continue

                parsed = _parse_gtf_attributes(cols[8])

                # Rebuild attributes: prefer gene_id then transcript_id then other keys
                order = []
                for k in ("gene_id", "transcript_id"):
                    if k in parsed and parsed[k]:
                        order.append(f'{k} "{parsed[k]}"')

                for k, v in parsed.items():
                    if k not in {"gene_id", "transcript_id"} and v:
                        order.append(f'{k} "{v}"')

                cols[8] = ("; ".join(order) + ";") if order else "."
                rout.write("\t".join(cols) + "\n")

        in_gtf_to_use = converted_gtf_clean
    else:
        in_gtf_to_use = in_path

    slim_gtf_path = tmp_dir / "slim.gtf"

    with open(in_gtf_to_use, "rt", encoding="utf-8") as inp, open(
        slim_gtf_path, "wt", encoding="utf-8"
    ) as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            ftype = fields[2]

            # Normalise transcript-like feature types to "transcript" for downstream compatibility
            if ftype in TRANSCRIPT_LIKE_FEATURES:
                fields[2] = "transcript"
                ftype = "transcript"

            if only_exons and ftype not in {"transcript", "exon"}:
                continue  # skip everything else

            attrs_dict = _parse_gtf_attributes(fields[8])

            new_attrs = []
            for key in essential_keys:
                if key in attrs_dict and attrs_dict[key]:
                    new_attrs.append(f'{key} "{attrs_dict[key]}"')

            fields[8] = ("; ".join(new_attrs) + ";") if new_attrs else "."
            out.write("\t".join(fields) + "\n")

    return slim_gtf_path

def _parse_gtf_attributes(attr_str):
    """Parse GTF/GFF3 attributes into a dict, stripping 'gene:' and 'transcript:' prefixes."""
    attrs = {}

    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue

        if "=" in part:
            # GFF3 style: key=value
            key, val = part.split("=", 1)
            val = val.strip().strip('"')
        elif '"' in part:
            # GTF style: key "value"
            key, rest = part.split('"', 1)
            val = rest.split('"')[0]
        else:
            continue

        key = key.strip()
        val = val.strip()

        # remove gene: and transcript: prefixes
        val = re.sub(r"^(?:gene:|transcript:)", "", val)

        attrs[key] = val

    return attrs


def _cleanup_tmp_folder():
    """Removes the local ./tmp folder with all its contents."""
    tmp_dir = Path("./tmp")
    try:
        shutil.rmtree(tmp_dir)
    except FileNotFoundError:
        pass


class _SelectiveConsoleFilter(logging.Filter):
    # Logging filter: only keep selected messages on console
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
        "ORFannotate completed",
    ]

    def filter(self, record):
        return any(term in record.getMessage() for term in self.ALLOWED_TERMS)


def _setup_logging(output_dir):
    log_path = os.path.join(output_dir, "ORFannotate.log")

    file_formatter = logging.Formatter(
        fmt="[%(asctime)s] [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    console_formatter = logging.Formatter("[%(levelname)s] %(message)s")

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


# Define custom progress bar columns
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


def _build_tx_lookup(db):
    lut = {}
    for ftype in TRANSCRIPT_LIKE_FEATURES:
        for tx in db.features_of_type(ftype):
            tid_full = tx.attributes["transcript_id"][0]
            tid_base = tid_full.split(".")[0]

            lut[tid_full] = tx
            lut[tid_base] = tx
            lut[_normalize_tid(tid_full)] = tx
            lut[_normalize_tid(tid_base)] = tx
    return lut


def _extract_tid_from_attrs(attr_field: str):
    
    # GTF style
    m = re.search(r'\btranscript_id\s+"([^"]+)"', attr_field)
    if m:
        return m.group(1)

    # # Very simple GFF3-style handling
    # m = re.search(r'(?:^|;)\s*Parent=([^;]+)', attr_field)
    # if m:
    #     print(m)
    #     parent = m.group(1)
    #     # strip optional transcript: prefix
    #     return re.sub(r'^(?:transcript:)?', '', parent)

    return None