import subprocess

def extract_transcripts(gtf, genome_fasta, output_fasta):
    cmd = ["gffread", gtf, "-g", genome_fasta, "-w", output_fasta]
    subprocess.run(cmd, check=True)