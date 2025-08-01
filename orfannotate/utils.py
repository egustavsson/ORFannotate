import os
import re
from pathlib import Path
from Bio import SeqIO

TX_RE = re.compile(r'transcript_id\s+"([^"]+)"')

def count_unique_transcripts(gtf_path: str) -> int:
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

def validate_inputs(gtf_path: str, genome_fasta: str):
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

def load_transcript_sequences(transcript_fasta: Path):
    """Load transcript sequences into a dictionary keyed by transcript ID."""
    return SeqIO.to_dict(SeqIO.parse(str(transcript_fasta), "fasta"))

def write_junctions(gtf_db, output_path):
    """Write junctions (donor, acceptor positions) per transcript."""
    with open(output_path, "w") as out:
        out.write("transcript_id\tdonor\tacceptor\n")
        for tx in gtf_db.features_of_type("transcript"):
            exons = list(gtf_db.children(tx, featuretype="exon", order_by="start"))
            if len(exons) > 1:
                for i in range(len(exons) - 1):
                    donor = exons[i].end
                    acceptor = exons[i + 1].start
                    out.write(f"{tx.id}\t{donor}\t{acceptor}\n")

def keep_tx_and_exons(feat):
    """Filter GTF features to only keep transcript and exon rows."""
    return feat if feat.featuretype in {"transcript", "exon"} else False

def map_junctions_to_tx(junctions_genomic, exons):
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
