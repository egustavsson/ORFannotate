import os
import re
import logging
import gffutils
import argparse
import json
from datetime import datetime
from pathlib import Path

# Internal modules
from orfannotate.orf_prediction import run_cpat
from orfannotate.transcript_extraction import extract_transcripts_from_gtf
from orfannotate.orf_filter import get_best_orfs_by_cpat
from orfannotate.gtf_annotation import build_cds_features, annotate_gtf_with_cds
from orfannotate.summarise import generate_summary
from orfannotate.utils import count_unique_transcripts, validate_inputs, keep_tx_and_exons

# version
from orfannotate import __version__

# Path relative to ORFannotate for Snakemake, Nextflow workflows, etc.
MODULE_DIR = Path(__file__).resolve().parent
DATA_DIR = MODULE_DIR / "data"

# Logging filter: only keep selected messages on console
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
    stream_handler.addFilter(SelectiveConsoleFilter())

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers = []
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger

# Collect splice junctions into memory
def _collect_junctions(gtf_db):
    junctions = {}
    for tx in gtf_db.features_of_type("transcript"):
        exons = list(gtf_db.children(tx, featuretype="exon", order_by="start"))
        if len(exons) > 1:
            junctions[tx.id] = [(exons[i].end, exons[i+1].start) for i in range(len(exons)-1)]
    return junctions

# Main function
def main():
    parser = argparse.ArgumentParser(
        description="ORFannotate - predict coding ORFs, annotate GTF, and generate summaries."
    )
    parser.add_argument("--gtf", required=True, help="Input GTF or GFF file with transcript and exon features")
    parser.add_argument("--fa", required=True, help="Reference genome in FASTA format")
    parser.add_argument("--outdir", required=True, help="Directory to write all outputs")
    parser.add_argument("--species", choices=["human", "mouse", "fly", "zebrafish"], default="human",
                        help="Species to use for CPAT model (default: human)")
    parser.add_argument("--coding-cutoff", type=float, default=None,
                        help="CPAT cutoff for classifying coding vs noncoding transcripts")
    parser.add_argument(
    "--version",
    action="version",
    version=f"ORFannotate {__version__}")
    args = parser.parse_args()

    gtf_path = args.gtf
    genome_fasta = args.fa
    output_dir = args.outdir
    species = args.species.lower()

    validate_inputs(gtf_path, genome_fasta)
    os.makedirs(output_dir, exist_ok=True)
    logger = _setup_logging(output_dir)

    logger.info("Starting ORFannotate")
    num_transcripts = count_unique_transcripts(gtf_path)
    logger.info(f"Found {num_transcripts:,} unique transcripts in the GTF")

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
        transform=keep_tx_and_exons,
        pragmas={
            "journal_mode": "OFF",
            "synchronous": "OFF",
            "temp_store": "MEMORY"
        },
        id_spec={
            "transcript": "transcript_id",
            "exon": "transcript_id",
            "CDS": "transcript_id"
        }
    )

    logger.info("Step 2: Extracting transcript sequences...")
    transcript_fasta = os.path.join(output_dir, "transcripts.fa")
    extract_transcripts_from_gtf(gtf_db, genome_fasta, transcript_fasta)

    # Validate transcript ID consistency
    from Bio import SeqIO
    fasta_tx_ids = {rec.id for rec in SeqIO.parse(transcript_fasta, "fasta")}
    gtf_tx_ids = {t.id for t in gtf_db.features_of_type("transcript")}
    missing = gtf_tx_ids - fasta_tx_ids
    if missing:
        logger.warning(f"{len(missing)} transcripts from GTF missing in extracted FASTA: e.g. {sorted(list(missing))[:5]}")

    logger.info("Step 3: Predicting and scoring ORFs...")
    
    # Load species config from JSON
    config_path = DATA_DIR / "config.json"
    with open(config_path) as f:
        config = json.load(f)

    if species not in config["species_models"]:
        logger.error(f"Species '{species}' not found in config.json")
        raise ValueError(f"Species '{species}' not supported. Check config.json.")

    species_data = config["species_models"][species]
    hexamer_path = str(DATA_DIR / species_data["hexamer"])
    logit_model_path = str(DATA_DIR / species_data["model"])
    coding_cutoff = args.coding_cutoff if args.coding_cutoff is not None else species_data["cutoff"]

    logger.info(f"Using CPAT model directory: {DATA_DIR}")
    logger.info(f"Selected species: {species} | Cutoff: {coding_cutoff}")
    logger.info(f"Using coding cutoff: {coding_cutoff}")

    cpat_dir = os.path.join(output_dir, "CPAT")
    os.makedirs(cpat_dir, exist_ok=True)
    run_cpat(transcript_fasta, cpat_dir, hexamer_path, logit_model_path)

    logger.info("Step 4: Parsing ORF results...")
    cpat_results = os.path.join(cpat_dir, "cpat.ORF_prob.best.tsv")
    debug_output = os.path.join(cpat_dir, "cpat_debug.tsv")
    best_orfs = get_best_orfs_by_cpat(cpat_results, debug_output_path=debug_output)

    coding_orfs = {tid: info for tid, info in best_orfs.items() if info["coding_prob"] >= coding_cutoff}
    logger.info(f"Selected {len(coding_orfs):,} transcripts classified as coding (cutoff = {coding_cutoff})")

    if not coding_orfs:
        logger.warning("No coding transcripts predicted above the cutoff. Downstream outputs may be empty or incomplete.")

    logger.info("Step 5: Annotating GTF with CDS features...")
    logger.info(f"Adding CDS to {len(coding_orfs):,} coding transcripts")

    cds_features = build_cds_features(gtf_db, coding_orfs)
    annotated_gtf = os.path.join(output_dir, "ORFannotate_annotated.gtf")
    annotate_gtf_with_cds(gtf_path, cds_features, annotated_gtf)
    logger.info(f"Annotated GTF written to {annotated_gtf}")

    # Cleanup first db
    try:
        gtf_db.conn.close()
        del gtf_db
        logger.info("Closed and removed first GTF database to free memory")
    except Exception:
        logger.warning("Could not close first GTF database early")

    logger.info("Step 6: Annotating transcripts and generating final summary TSV...")
    annotated_db = gffutils.create_db(
        annotated_gtf,
        dbfn=":memory:",
        force=True,
        keep_order=True,
        disable_infer_transcripts=True,
        disable_infer_genes=True,
        merge_strategy="create_unique",
        sort_attribute_values=True,
        pragmas={
            "journal_mode": "OFF",
            "synchronous": "OFF",
            "temp_store": "MEMORY"
        },
        id_spec={
            "transcript": "transcript_id",
            "exon": "transcript_id",
            "CDS": "transcript_id"
        }
    )

    junctions_mem = _collect_junctions(annotated_db)
    summary_tsv = os.path.join(output_dir, "ORFannotate_summary.tsv")
    generate_summary(best_orfs, transcript_fasta, annotated_db, summary_tsv, coding_cutoff=coding_cutoff, junctions_by_tx=junctions_mem)

    try:
        annotated_db.conn.close()
    except Exception:
        pass

    logger.info(f"Summary written to {summary_tsv}")

    params = {
        "timestamp": datetime.now().isoformat(),
        "gtf": os.path.abspath(gtf_path),
        "fasta": os.path.abspath(genome_fasta),
        "outdir": os.path.abspath(output_dir),
        "species": species,
        "coding_cutoff": coding_cutoff
    }
    with open(os.path.join(output_dir, "parameters.json"), "w") as f:
        json.dump(params, f, indent=2)

    logger.info("Parameters saved to parameters.json")
    logger.info("ORFannotate completed successfully.")

if __name__ == "__main__":
    main()
