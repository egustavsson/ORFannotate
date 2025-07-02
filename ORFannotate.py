import os
import re
import io
import logging
import gffutils
import argparse

# Internal modules
from orfannotate.orf_prediction import run_cpat
from orfannotate.transcript_extraction import extract_transcripts_from_gtf
from orfannotate.orf_filter import get_best_orfs_by_cpat, build_cds_features
from orfannotate.gtf_annotation import annotate_gtf_with_cds
from orfannotate.summarise import generate_summary

# logging info. Only high-level written to terminal
class SelectiveConsoleFilter(logging.Filter):
    ALLOWED_TERMS = [
        "Starting ORFannotate",
        "Found",
        "Step 1:",
        "Step 2:",
        "Step 3:",
        "Step 4:",
        "Step 5:",
        "Step 6:",
        "ORFannotate completed"
    ]

    def filter(self, record):
        return any(term in record.getMessage() for term in self.ALLOWED_TERMS)

# Help functions

def _setup_logging(output_dir):
    log_path = os.path.join(output_dir, "ORFannotate.log")

    # Detailed file format
    file_formatter = logging.Formatter(
        fmt='[%(asctime)s] [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Minimal terminal format
    console_formatter = logging.Formatter('[%(levelname)s] %(message)s')

    # File handler – everything goes here
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(file_formatter)

    # Terminal handler – filtered output
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(console_formatter)
    stream_handler.addFilter(SelectiveConsoleFilter())

    # Set up the logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers = []
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger

# count number of transcripts in GTF
TX_RE = re.compile(r'transcript_id\s+"([^"]+)"')

def _count_unique_transcripts(gtf_path):
    ids = set()
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#") or "\ttranscript\t" not in line:
                continue
            m = TX_RE.search(line)
            if m:
                ids.add(m.group(1))
    return len(ids)

# remove CDS features if present in input GTF
def _keep_tx_and_exons(feat):
    return feat if feat.featuretype in {"transcript", "exon"} else False

# Main function

def main():
    parser = argparse.ArgumentParser(
        description="ORFannotate – predict coding ORFs, annotate GTF, and generate summaries."
    )
    parser.add_argument("--gtf", required=True, help="Input GTF or GFF file with transcript and exon features")
    parser.add_argument("--fa", required=True, help="Reference genome in FASTA format")
    parser.add_argument("--outdir", required=True, help="Directory to write all outputs")
    parser.add_argument("--coding-cutoff", type=float, default=0.364,
                        help="CPAT cutoff for classifying coding vs noncoding transcripts (default: 0.364 for human)")
    args = parser.parse_args()

    gtf_path = args.gtf
    genome_fasta = args.fa
    output_dir = args.outdir
    coding_cutoff = args.coding_cutoff

    os.makedirs(output_dir, exist_ok=True)
    logger = _setup_logging(output_dir)

    logger.info("Starting ORFannotate")
    num_transcripts = _count_unique_transcripts(gtf_path)
    logger.info(f"Found {num_transcripts:,} unique transcripts in the GTF")

    # Building the GTF database once that will be used throughout
    logger.info("Step 1: Building GTF database in memory...")
    gtf_db = gffutils.create_db(
        gtf_path,
        dbfn=":memory:",
        force=True,
        keep_order=True,
        disable_infer_transcripts=True,
        disable_infer_genes=True,
        merge_strategy="create_unique",
        sort_attribute_values=True,
        transform=_keep_tx_and_exons,
        pragmas={"journal_mode": "OFF",
                 "synchronous": "OFF",
                 "temp_store": "MEMORY"},
        id_spec={"transcript": "transcript_id",
                 "exon": "transcript_id",
                 "CDS": "transcript_id"}
    )

    # extract transcript FASTA sequences using the db
    logger.info("Step 2: Extracting transcript sequences...")
    transcript_fasta = os.path.join(output_dir, "transcripts.fa")
    extract_transcripts_from_gtf(gtf_db, genome_fasta, transcript_fasta)

    # ORF prediction using CPAT
    logger.info("Step 3: Predicting and scoring ORFs...")
    hexamer_path = os.path.abspath("data/Human_Hexamer.tsv")
    logit_model_path = os.path.abspath("data/Human_logitModel.RData")
    run_cpat(transcript_fasta, output_dir, hexamer_path, logit_model_path)

    # parse ORF result to select the best per transcript
    logger.info("Step 4: Parsing ORF results..." )
    cpat_results = os.path.join(output_dir, "cpat.ORF_prob.best.tsv")
    debug_output = os.path.join(output_dir, "cpat_debug.tsv")
    best_orfs = get_best_orfs_by_cpat(cpat_results, debug_output_path=debug_output)

    coding_orfs = {
        tid: info for tid, info in best_orfs.items() if info["coding_prob"] >= coding_cutoff
    }
    logger.info(
        f"Selected {len(coding_orfs):,} transcripts classified as coding (cutoff = {coding_cutoff})",
    )

    # add CDS features to GTF
    logger.info("Step 5: Annotating GTF with CDS features...")

    logger.info(f"Adding CDS to {len(coding_orfs):,} coding transcripts")

    cds_features = build_cds_features(gtf_db, coding_orfs)
    annotated_gtf = os.path.join(output_dir, "ORFannotate_annotated.gtf")
    annotate_gtf_with_cds(gtf_path, cds_features, annotated_gtf)
    logger.info(f"Annotated GTF written to {annotated_gtf}")

    # build a lightweight DB from the annotated GTF so that generate_summary() sees the newly added CDS.
    logger.info("Building database with CDS features for summary generation")
    annotated_db = gffutils.create_db(
        annotated_gtf,
        dbfn=":memory:",
        force=True,
        keep_order=True,
        disable_infer_transcripts=True,
        disable_infer_genes=True,
        merge_strategy="create_unique",
        sort_attribute_values=True,
        pragmas={"journal_mode": "OFF",
                 "synchronous": "OFF",
                 "temp_store": "MEMORY"},
        id_spec={"transcript": "transcript_id",
                 "exon": "transcript_id",
                 "CDS": "transcript_id"},
    )

    # generate summary and write table
    logger.info("Step 6: Annotating and generating final summary TSV...")
    summary_tsv = os.path.join(output_dir, "ORFannotate_summary.tsv")
    generate_summary(
        best_orfs,
        transcript_fasta,
        annotated_db,
        summary_tsv,
        coding_cutoff=coding_cutoff,
    )

    # Close databases
    try:
        gtf_db.conn.close()
        annotated_db.conn.close()
    except Exception:
        pass

    logger.info(f"Summary written to {summary_tsv}")
    logger.info("ORFannotate completed successfully.")

if __name__ == "__main__":
    main()
