import os
import argparse
import subprocess
from orfannotate.io_utils import extract_transcripts
from orfannotate.orf_parser import parse_orfs
from orfannotate.nmd import detect_nmd
from orfannotate.gtf_annotation import annotate_gtf
from orfannotate.summarise import generate_summary

def main():
    parser = argparse.ArgumentParser(description="Run ORFannotate pipeline.")
    parser.add_argument("gtf", help="Input GTF or GFF file")
    parser.add_argument("genome_fasta", help="Reference genome FASTA")
    parser.add_argument("out_dir", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    transcript_fasta = os.path.join(args.out_dir, "transcripts.fa")
    annotated_gtf = os.path.join(args.out_dir, "annotated.gtf")
    output_tsv = os.path.join(args.out_dir, "transcript_summary.tsv")
    orf_fasta = os.path.join(args.out_dir, "predicted_orfs.fasta")

    print("[Step 1] Extracting transcript sequences...")
    extract_transcripts(args.gtf, args.genome_fasta, transcript_fasta)

    print("[Step 2] Running ORFipy...")
    subprocess.run([
    "orfipy", transcript_fasta,
    "--outfmt", "fasta",
    "--min", "30",
    "--strand", "both",
    "--out", orf_fasta
], check=True)

    print("[Step 3] Annotating GTF with CDS features...")
    annotate_gtf(args.gtf, orf_fasta, annotated_gtf)

    print("[Step 4] Parsing ORF results...")
    orf_info = parse_orfs(orf_fasta)

    print("[Step 5] Generating final summary TSV...")
    generate_summary(transcript_fasta, orf_info, args.gtf, output_tsv)

    print("Pipeline completed successfully.")

if __name__ == "__main__":
    main()