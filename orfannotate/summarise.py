from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Tuple
import tempfile
import json
import argparse

import pandas as pd
from Bio import SeqIO
import gffutils

from orfannotate.nmd import predict_nmd

LOGGER = logging.getLogger(__name__)
logging.basicConfig(format="%(levelname)s:%(name)s:%(message)s", level=logging.INFO)

# ---------------------------- helpers ---------------------------- #

def _load_transcripts(fa: Path):
    """Return transcript‑id → SeqRecord dict."""
    return SeqIO.to_dict(SeqIO.parse(str(fa), "fasta"))

# ------------------------- main routine -------------------------- #

def generate_summary(best_orfs: dict, transcript_fa: str | Path, gtf_path: str | Path,
                     output_path: str | Path, coding_cutoff: float = 0.364) -> None:
    """Write per‑transcript summary TSV including NMD prediction."""

    # Build temporary gffutils DB that keeps CDS rows
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmp:
        gffutils.create_db(
            str(gtf_path),
            dbfn=tmp.name,
            force=True,
            keep_order=True,
            disable_infer_transcripts=True,
            disable_infer_genes=True,
            merge_strategy="merge",
            sort_attribute_values=True,
            id_spec={
                "transcript": "transcript_id",
                "exon":        "transcript_id",
                "CDS":         "transcript_id",
            },
        )
        db = gffutils.FeatureDB(tmp.name)

    tx_seqs = _load_transcripts(Path(transcript_fa))
    rows: list[dict] = []

    for tid, orf in best_orfs.items():
        if tid not in tx_seqs:
            continue  # no sequence

        # fetch transcript feature (try version‑stripped id as fallback)
        try:
            tx = db[tid]
        except KeyError:
            try:
                tx = db[tid.split(".")[0]]
            except KeyError:
                continue

        strand = tx.strand

        # CDS stop in genomic coords
        cds_feats = list(db.children(tx, featuretype="CDS", order_by="start"))
        orf_end_gen: int | None = None
        if cds_feats:
            orf_end_gen = max(f.end for f in cds_feats) if strand == "+" else min(f.start for f in cds_feats)

        coding_prob = orf["coding_prob"]
        coding_class = "coding" if (coding_prob >= coding_cutoff and orf_end_gen is not None) else "noncoding"

        # build junction list in transcript order
        exons = list(db.children(tx, featuretype="exon", order_by="start"))
        if strand == "-":
            exons.reverse()
        junctions: List[Tuple[int, int]] = [(exons[i].end, exons[i + 1].start) for i in range(len(exons) - 1)]

        # distance stop → last EJ
        if orf_end_gen is not None and junctions:
            last_j = junctions[-1][0] if strand == "+" else junctions[0][1]
            dist_stop_last = last_j - orf_end_gen if strand == "+" else orf_end_gen - last_j
        else:
            dist_stop_last = "NA"

        # NMD call
        nmd_flag = (
            predict_nmd(orf_end_gen, junctions, strand)
            if coding_class == "coding" else "FALSE"
        )

        # sequences (transcript‑relative slice)
        nt_seq = str(tx_seqs[tid].seq[orf["start"] - 1 : orf["end"]])
        aa_seq = str(tx_seqs[tid].seq[orf["start"] - 1 : orf["end"]].translate(to_stop=True))

        rows.append({
            "transcript_id": tid,
            "gene_id": tx.attributes.get("gene_id", ["NA"])[0],
            "chrom": tx.chrom,
            "strand": strand,
            "orf_start": orf["start"],
            "orf_end": orf_end_gen if orf_end_gen is not None else orf["end"],
            "orf_nt_len": orf["end"] - orf["start"] + 1,
            "orf_aa_len": (orf["end"] - orf["start"] + 1) // 3,
            "coding_prob": coding_prob,
            "coding_class": coding_class,
            "junction_count": len(junctions),
            "stop_to_last_EJ": dist_stop_last,
            "NMD_sensitive": nmd_flag,
            "cds_sequence": nt_seq if coding_class == "coding" else "NA",
            "protein_sequence": aa_seq if coding_class == "coding" else "NA",
        })

    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":

    cli = argparse.ArgumentParser(description="Write ORF summary TSV")
    cli.add_argument("--best-orfs-json", required=True)
    cli.add_argument("--transcripts-fa", required=True)
    cli.add_argument("--gtf", required=True)
    cli.add_argument("--out", required=True)
    cli.add_argument("--coding-cutoff", type=float, default=0.364)
    args = cli.parse_args()

    best_orfs = json.load(open(args.best_orfs_json))
    generate_summary(best_orfs, args.transcripts_fa, args.gtf, args.out, args.coding_cutoff)
