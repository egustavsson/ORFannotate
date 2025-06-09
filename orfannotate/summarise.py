from __future__ import annotations

import csv
import logging
import os
from pathlib import Path
from typing import Dict, Any, List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
import gffutils
from gffutils.exceptions import FeatureNotFoundError

from orfannotate.nmd import predict_nmd


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

# Main function

def generate_summary(
    best_orfs: Dict[str, Dict[str, Any]],
    transcript_fasta: str | os.PathLike,
    gtf_path: str | os.PathLike,
    output_tsv: str | os.PathLike,
):
    """Write a richly annotated TSV to *output_tsv*.

    Parameters
    ----------
    best_orfs
        Mapping returned by :pyfunc:`orf_filter.get_best_orfs_by_cpat`.
    transcript_fasta
        FASTA with one record per transcript (output of
        :pyfunc:`gtf_annotation.extract_transcripts_from_gtf`).
    gtf_path
        The original gene annotation GTF.
    output_tsv
        Destination path for the summary.
    """

    tx_seqs = _load_transcript_sequences(Path(transcript_fasta))
    db = _get_gtf_db(Path(gtf_path))

    header: List[str] = [
        "transcript_id",
        "gene_id",
        "gene_name",
        "chrom",
        "strand",
        "orf_start",
        "orf_end",
        "frame",
        "orf_len_nt",
        "orf_len_aa",
        "coding_prob",
        "junction_count",
        "predicted_nmd",
        "orf_nt_seq",
        "orf_aa_seq",
    ]

    missing_in_gtf = 0
    missing_in_fa = 0

    with open(output_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)

        for tid, info in best_orfs.items():
            
            # Resolve transcript feature in GTF DB (robust to version suffix)
            
            try:
                tx = db[tid]
            except FeatureNotFoundError:
                try:
                    tx = db[tid.split(".")[0]]
                except FeatureNotFoundError:
                    LOGGER.warning("Transcript %s not found in GTF – writing NA" % tid)
                    missing_in_gtf += 1
                    tx = None

            if tx is not None:
                gene_id = tx.attributes.get("gene_id", ["NA"])[0]
                gene_name = tx.attributes.get("gene_name", ["NA"])[0]
                chrom = tx.chrom
                strand = tx.strand
            else:
                gene_id = gene_name = chrom = strand = "NA"

            # Sequences & NMD flag
            
            orf_len_nt = int(info["length"])
            orf_len_aa = orf_len_nt // 3
            junctions = info.get("junctions", [])
            nmd_flag = predict_nmd(info["end"], junctions, strand)

            # Fetch nt + aa sequences; handle absent FASTA entries gracefully.
            try:
                nt_seq, aa_seq = _extract_orf_seqs(
                    tx_seqs[tid], int(info["start"]), int(info["end"])
                )
            except KeyError:
                LOGGER.warning("Transcript %s not found in FASTA – empty sequence" % tid)
                missing_in_fa += 1
                nt_seq = aa_seq = "NA"

            # Write TSV line
       
            writer.writerow(
                [
                    tid,
                    gene_id,
                    gene_name,
                    chrom,
                    strand,
                    info["start"],
                    info["end"],
                    info["frame"],
                    orf_len_nt,
                    orf_len_aa,
                    info["coding_prob"],
                    len(junctions),
                    nmd_flag,
                    nt_seq,
                    aa_seq,
                ]
            )

    if missing_in_gtf or missing_in_fa:
        LOGGER.info(
            "Summary complete – %d transcripts absent from GTF, %d absent from FASTA",
            missing_in_gtf,
            missing_in_fa,
        )

    return output_tsv


# CLI helper

if __name__ == "__main__":
    import json
    import argparse

    p = argparse.ArgumentParser(description="Write ORF summary TSV")
    p.add_argument("--best-orfs-json", required=True, help="Path to best_orfs JSON")
    p.add_argument("--transcripts-fa", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    best_orfs = json.load(open(args.best_orfs_json))
    generate_summary(best_orfs, args.transcripts_fa, args.gtf, args.out)
