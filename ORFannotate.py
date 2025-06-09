import os
import sys
import subprocess
import shutil
import gffutils
import logging

# Configure logging level to suppress warnings
logging.basicConfig(level=logging.ERROR)

from orfannotate.orf_filter import get_best_orfs_by_cpat, build_cds_features
from orfannotate.gtf_annotation import extract_transcripts_from_gtf, annotate_gtf_with_cds
from orfannotate.summarise import generate_summary


def run_cpat(transcript_fasta, output_dir):
    output_prefix = os.path.join(output_dir, "cpat")

    # Use absolute paths for model files
    hexamer_path = os.path.abspath("data/Human_Hexamer.tsv")
    logit_model_path = os.path.abspath("data/Human_logitModel.RData")

    cpat_cmd = [
        "cpat.py",
        "-g", transcript_fasta,
        "-x", hexamer_path,
        "-d", logit_model_path,
        "--top-orf=5",
        "--min-orf=75",
        "-o", output_prefix
    ]

    cpat_log = os.path.join(output_dir, "CPAT_run_info.log")
    with open(cpat_log, "w") as log_fh:
        subprocess.run(
            cpat_cmd,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            check=True,
            cwd=output_dir  # ensures CPAT_run_info.log is written here
        )

def main():
    gtf_path = sys.argv[1]
    genome_fa = sys.argv[2]
    output_dir = sys.argv[3]

    os.makedirs(output_dir, exist_ok=True)

    print("[Step 1] Extracting transcript sequences...")
    transcript_fa = os.path.join(output_dir, "transcripts.fa")
    extract_transcripts_from_gtf(gtf_path, genome_fa, transcript_fa)

    print("[Step 2] Running CPAT to predict and score ORFs...")
    run_cpat(transcript_fa, output_dir)

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

    print("[Step 4] Parsing CPAT results...")
    cpat_scores_path = os.path.join(output_dir, "cpat.ORF_prob.best.tsv")
    cpat_debug = os.path.join(output_dir, "cpat_debug.tsv")
    best_orfs = get_best_orfs_by_cpat(cpat_scores_path, debug_output_path=cpat_debug)

    cds_features = list(build_cds_features(gtf_db, best_orfs))


    annotated_gtf = os.path.join(output_dir, "annotated.gtf")
    annotate_gtf_with_cds(gtf_path, cds_features, annotated_gtf)

    print("[Step 5] Generating final summary TSV...")
    summary_path = os.path.join(output_dir, "orf_summary.tsv")

    generate_summary(
    best_orfs,
    transcript_fa,
    gtf_path,
    summary_path
)

    print("ORFannotate completed successfully.")

if __name__ == "__main__":
    main()
