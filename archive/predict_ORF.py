# This script predicts open reading frames and detects nonsense-mediated decay
# We use orfipy ORF predictions - https://github.com/urmi-21/orfipy

import argparse
import subprocess
import os

def run_orfipy(fasta_file, output_dir, min_length=30, strand='both'):
    """
    Runs ORF prediction using orfipy for each entry in a FASTA file.
    
    :param fasta_file: Path to the input FASTA file.
    :param output_dir: Directory where the output files will be saved.
    :param min_length: Minimum length of ORFs to be predicted.
    :param strand: Strand selection ('plus', 'minus', or 'both').
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_file = os.path.join(output_dir, "predicted_orfs.fasta")
    
    cmd = [
        "orfipy", fasta_file,
        "--outfmt", "fasta",
        "--min", str(min_length),
        "--strand", strand,
        "--out", output_file
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"ORF prediction completed. Results saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running orfipy: {e}")

def detect_nmd(isoform_hit, rec):
    """
    NMD detection function.
    
    :param isoform_hit: ORF hit containing CDS genomic coordinates.
    :param rec: Record containing junction information.
    """
    if len(rec.junctions) > 0:
        if rec.strand == '+':
            dist_to_last_junc = isoform_hit.CDS_genomic_end - rec.junctions[-1][0]
        else: # - strand
            dist_to_last_junc = rec.junctions[0][1] - isoform_hit.CDS_genomic_end
        isoform_hit.is_NMD = "TRUE" if dist_to_last_junc < -50 else "FALSE"

def process_orfs(fasta_file, output_dir):
    """
    Process predicted ORFs and run NMD detection.
    
    :param fasta_file: Path to the input FASTA file.
    :param output_dir: Directory where the output files will be saved.
    """
    predicted_orfs = os.path.join(output_dir, "predicted_orfs.fasta")
    
    # Placeholder: Read ORFs and their records
    orfs = []  # Replace with actual ORF parsing logic
    records = {}  # Replace with actual record fetching logic
    
    for orf in orfs:
        record = records.get(orf.id)
        if record:
            detect_nmd(orf, record)
            print(f"ORF {orf.id}: NMD status = {orf.is_NMD}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ORF prediction using orfipy on a FASTA file.")
    parser.add_argument("fasta_file", type=str, help="Path to the input FASTA file")
    parser.add_argument("output_dir", type=str, help="Path to the output directory")
    parser.add_argument("--min_length", type=int, default=30, help="Minimum ORF length (default: 30)")
    parser.add_argument("--strand", type=str, choices=['plus', 'minus', 'both'], default='both', help="Strand to analyze (default: both)")
    
    args = parser.parse_args()
    
    run_orfipy(args.fasta_file, args.output_dir, args.min_length, args.strand)
    process_orfs(args.fasta_file, args.output_dir)
