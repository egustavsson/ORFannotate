from __future__ import annotations

import csv
import logging
import os
from pathlib import Path
from typing import Dict, Any, List, Tuple

import pandas as pd
from orfannotate.nmd import predict_nmd

from Bio import SeqIO
from Bio.Seq import Seq
import gffutils
from gffutils.exceptions import FeatureNotFoundError



LOGGER = logging.getLogger(__name__)
logging.basicConfig(format="%(levelname)s:%(name)s:%(message)s", level=logging.INFO)

# Helper functions

def _load_transcript_sequences(transcript_fasta: Path):
    """Return a dict keyed by transcript ID → SeqRecord."""
    return SeqIO.to_dict(SeqIO.parse(str(transcript_fasta), "fasta"))


def _get_gtf_db(gtf_path: Path):
    """Create (if necessary) and open a **gffutils** DB from the supplied GTF."""
    db_path = gtf_path.with_suffix(".db")
    if not db_path.exists():
        LOGGER.info("Building gffutils cache … this is a one‑time cost ∼")
        gffutils.create_db(
            str(gtf_path),
            dbfn=str(db_path),
            force=True,
            keep_order=True,
            disable_infer_transcripts=True,
            disable_infer_genes=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )
    return gffutils.FeatureDB(str(db_path))


def _extract_orf_seqs(tx_record, start: int, end: int) -> Tuple[str, str]:
    """Return (nt_seq, aa_seq) for the 1‑based inclusive ORF coords."""
    nt = tx_record.seq[start - 1 : end]
    return str(nt), str(nt.translate(to_stop=True))

# Function to count junctions for each ORF
def count_exon_junctions(db, transcript, orf_start, orf_end):
    """Count exon–exon junctions spanned by the ORF."""
    exons = sorted(db.children(transcript, featuretype="exon"), key=lambda x: x.start)
    junctions = []
    for i in range(len(exons) - 1):
        donor = exons[i].end
        acceptor = exons[i+1].start
        if donor < orf_end and acceptor > orf_start:
            junctions.append((donor, acceptor))
    return len(junctions)

# Main function to generate the summary TSV
def generate_summary(best_orfs, transcript_fa, gtf_path, output_path, coding_cutoff=0.364):
   
    db = gffutils.create_db(
        gtf_path,
        dbfn=":memory:",
        force=True,
        keep_order=True,
        disable_infer_transcripts=True,
        disable_infer_genes=True,
        merge_strategy="merge",
        sort_attribute_values=True,
    )

    transcript_seqs = SeqIO.to_dict(SeqIO.parse(transcript_fa, "fasta"))
    summary = []

    for tid, orf_data in best_orfs.items():
        coding_prob = orf_data["coding_prob"]
        coding_class = "coding" if coding_prob >= coding_cutoff else "noncoding"

        if tid not in transcript_seqs:
            continue

        try:
            tx = db[tid]
        except Exception:
            try:
                tx = db[tid.split(".")[0]]
            except Exception:
                continue

        chrom = tx.chrom
        strand = tx.strand
        gene_id = tx.attributes.get("gene_id", ["NA"])[0]

        orf_start = orf_data["start"]
        orf_end = orf_data["end"]
        frame = orf_data["frame"]
        orf_nt_len = orf_end - orf_start + 1
        orf_seq = str(transcript_seqs[tid].seq[orf_start - 1 : orf_end])
        orf_aa_len = len(orf_seq) // 3
        protein_seq = str(transcript_seqs[tid].seq[orf_start - 1 : orf_end].translate(to_stop=True))

        cds_seq = orf_seq
        
        exons = list(db.children(tx, featuretype="exon", order_by="start"))
        if strand == '-':
            exons = exons[::-1]
        
        junctions = [
            (exons[i].end, exons[i + 1].start)
            for i in range(len(exons) - 1)
        ]
        
        junction_count = len(junctions)
        
        if coding_class == "noncoding":
            nmd_flag = "FALSE"
        else:
            nmd_flag = predict_nmd(orf_end, junctions, strand)

        summary.append({
            "transcript_id": tid,
            "gene_id": gene_id,
            "chrom": chrom,
            "strand": strand,
            "orf_start": orf_start,
            "orf_end": orf_end,
            "frame": frame,
            "orf_nt_len": orf_nt_len,
            "orf_aa_len": orf_aa_len,
            "coding_prob": coding_prob,
            "coding_class": coding_class,
            "junction_count": junction_count,
            "NMD_sensitive": nmd_flag,
            "cds_sequence": cds_seq,
            "protein_sequence": protein_seq
        })

    df = pd.DataFrame(summary)
    df.to_csv(output_path, sep="\t", index=False)


# CLI helper
if __name__ == "__main__":
    import json
    import argparse

    p = argparse.ArgumentParser(description="Write ORF summary TSV")
    p.add_argument("--best-orfs-json", required=True, help="Path to best_orfs JSON")
    p.add_argument("--transcripts-fa", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", required=True)
    p.add_argument(
        "--coding-cutoff",
        type=float,
        default=0.364,
        help="CPAT cutoff for classifying coding vs noncoding (default: 0.364)"
    )
    args = p.parse_args()

    best_orfs = json.load(open(args.best_orfs_json))
    generate_summary(
        best_orfs,
        args.transcripts_fa,
        args.gtf,
        args.out,
        coding_cutoff=args.coding_cutoff
    )

