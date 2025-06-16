import os
import subprocess
import logging
import gffutils
import argparse

# Suppress warnings unless explicitly needed
logging.basicConfig(level=logging.ERROR)

# Internal modules
from orfannotate.transcript_extraction import extract_transcripts_from_gtf
from orfannotate.orf_filter import get_best_orfs_by_cpat, build_cds_features
from orfannotate.gtf_annotation import annotate_gtf_with_cds
from orfannotate.summarise import generate_summary


def run_cpat(transcript_fasta, output_dir):
    """
    Run CPAT to predict and score ORFs on transcript sequences.
    """
    output_prefix = os.path.join(output_dir, "cpat")
    hexamer_path = os.path.abspath("data/Human_Hexamer.tsv")
    logit_model_path = os.path.abspath("data/Human_logitModel.RData")

    cpat_cmd = [
        "cpat.py",
        "-g", transcript_fasta,
        "-x", hexamer_path,
        "-d", logit_model_path,
        "--top-orf=10",
        "--min-orf=75",
        "-o", output_prefix
    ]

    cpat_log_path = os.path.join(output_dir, "CPAT_run_info.log")
    with open(cpat_log_path, "w") as log_file:
        subprocess.run(
            cpat_cmd,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=True,
            cwd=output_dir
        )


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

    # Unpack arguments for clarity
    gtf_path = args.gtf
    genome_fasta = args.fa
    output_dir = args.outdir
    coding_cutoff = args.coding_cutoff

    os.makedirs(output_dir, exist_ok=True)

    print("[Step 1] Extracting transcript sequences...")
    transcript_fasta = os.path.join(output_dir, "transcripts.fa")
    extract_transcripts_from_gtf(gtf_path, genome_fasta, transcript_fasta)

    print("[Step 2] Predicting and scoring ORFs...")
    run_cpat(transcript_fasta, output_dir)

    print("[Step 3] Annotating GTF with CDS features...")
    db_path = os.path.join(output_dir, "gtf.db")
    gtf_db = gffutils.create_db(
        gtf_path,
        dbfn=db_path,
        force=True,
        keep_order=True,
        disable_infer_transcripts=True,
        disable_infer_genes=True,
        merge_strategy="merge",
        sort_attribute_values=True
    )

    print("[Step 4] Parsing ORF results...")
    cpat_results = os.path.join(output_dir, "cpat.ORF_prob.best.tsv")
    debug_output = os.path.join(output_dir, "cpat_debug.tsv")
    best_orfs = get_best_orfs_by_cpat(cpat_results, debug_output_path=debug_output)

    coding_orfs = {
        tid: info for tid, info in best_orfs.items()
        if info["coding_prob"] >= coding_cutoff
    }

    cds_features = list(build_cds_features(gtf_db, coding_orfs))
    annotated_gtf = os.path.join(output_dir, "ORFannotate_annotated.gtf")
    annotate_gtf_with_cds(gtf_path, cds_features, annotated_gtf)

    print("[Step 5] Generating final summary TSV...")
    summary_tsv = os.path.join(output_dir, "ORFannotate_summary.tsv")
    generate_summary(best_orfs, transcript_fasta, annotated_gtf, summary_tsv, coding_cutoff=coding_cutoff)

    print("ORFannotate completed successfully.")


if __name__ == "__main__":
    main()
