from __future__ import annotations

import logging
import json
import tempfile
from pathlib import Path
from typing import List, Tuple
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from Bio import SeqIO
import gffutils

from orfannotate.nmd import predict_nmd

LOGGER = logging.getLogger(__name__)
logging.basicConfig(format="%(levelname)s:%(name)s:%(message)s", level=logging.INFO)


# Help functions
def _load_transcript_sequences(transcript_fasta: Path):
    
    return SeqIO.to_dict(SeqIO.parse(str(transcript_fasta), "fasta"))


# Main routine
def generate_summary(best_orfs, transcript_fa, gtf_path, output_path, coding_cutoff=0.364):

    # build a temp gffutils DB that links CDS to transcript_id
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmpdb:
        gffutils.create_db(
            str(gtf_path),
            dbfn=tmpdb.name,
            force=True,
            keep_order=True,
            disable_infer_transcripts=True,
            disable_infer_genes=True,
            merge_strategy="create_unique",
            sort_attribute_values=True,
            id_spec={
                "transcript": "transcript_id",
                "exon": "transcript_id",
                "CDS": "transcript_id",
            },
        )
        db = gffutils.FeatureDB(tmpdb.name)

    transcript_seqs = _load_transcript_sequences(Path(transcript_fa))
    summary = []

    # iterate over EVERY transcript that has a sequence
    all_tx_ids = {t.id for t in db.features_of_type("transcript")}
    all_tx_ids &= set(transcript_seqs.keys())

    for tid in sorted(all_tx_ids):
        has_orf = tid in best_orfs
        orf_data = best_orfs.get(tid, None)

        try:
            tx = db[tid]
        except KeyError:
            tx = db.get(tid.split(".")[0], None)
            if tx is None:
                continue

        chrom = tx.chrom
        strand = tx.strand
        gene_id = tx.attributes.get("gene_id", ["NA"])[0]

        # ORF fields
        if has_orf:
            orf_start = orf_data["start"]
            orf_end_tx = orf_data["end"]
            frame = orf_data["frame"]
            coding_prob = orf_data["coding_prob"]
            nt_seq = str(transcript_seqs[tid].seq[orf_start - 1 : orf_end_tx])
            aa_seq = str(transcript_seqs[tid].seq[orf_start - 1 : orf_end_tx].translate(to_stop=True))
            orf_nt_len = orf_end_tx - orf_start + 1
            orf_aa_len = len(aa_seq)
        else:
            orf_start = orf_end_tx = frame = "NA"
            coding_prob = "NA"
            nt_seq = aa_seq = "NA"
            orf_nt_len = orf_aa_len = "NA"

        # CDS stop in genomic coords
        cds_feats = list(db.children(tx, featuretype="CDS", order_by="start"))
        orf_end_gen = None
        if cds_feats:
            orf_end_gen = max(f.end for f in cds_feats) if strand == "+" else min(f.start for f in cds_feats)

        # junction list
        exons = list(db.children(tx, featuretype="exon", order_by="start"))
        if strand == '-':
            exons = exons[::-1]

        junctions: List[Tuple[int, int]] = [
            (exons[i].end, exons[i + 1].start) for i in range(len(exons) - 1)
        ]
        junction_count = len(junctions)

        # distance stop → last junction
        if has_orf and orf_end_gen is not None and junctions:
            last_j = junctions[-1][0] if strand == '+' else junctions[0][1]
            stop_to_last_ej = last_j - orf_end_gen if strand == '+' else orf_end_gen - last_j
        else:
            stop_to_last_ej = "NA"

        # coding / NMD flags
        coding_class = (
            "coding"
            if (has_orf and isinstance(coding_prob, (float, int)) and coding_prob >= coding_cutoff and orf_end_gen)
            else "noncoding"
        )
        nmd_flag = predict_nmd(orf_end_gen, junctions, strand) if coding_class == "coding" else "FALSE"

        summary.append({
            "transcript_id": tid,
            "gene_id": gene_id,
            "chrom": chrom,
            "strand": strand,
            "has_orf": "TRUE" if has_orf else "FALSE",
            "orf_start": orf_start,
            "orf_end": orf_end_gen if orf_end_gen is not None else orf_end_tx,
            "orf_nt_len": orf_nt_len,
            "orf_aa_len": orf_aa_len,
            "coding_prob": coding_prob,
            "coding_class": coding_class,
            "junction_count": junction_count,
            "stop_to_last_EJ": stop_to_last_ej,
            "NMD_sensitive": nmd_flag,
            "cds_sequence": nt_seq if coding_class == "coding" else "NA",
            "protein_sequence": aa_seq if coding_class == "coding" else "NA",
        })

    pd.DataFrame(summary).to_csv(output_path, sep="\t", index=False)


# CLI
if __name__ == "__main__":
    import argparse

    cli = argparse.ArgumentParser(description="Write ORF summary TSV")
    cli.add_argument("--best-orfs-json", required=True)
    cli.add_argument("--transcripts-fa", required=True)
    cli.add_argument("--gtf", required=True)
    cli.add_argument("--out", required=True)
    cli.add_argument("--coding-cutoff", type=float, default=0.364)
    args = cli.parse_args()

    best_orfs_dict = json.load(open(args.best_orfs_json))
    generate_summary(
        best_orfs_dict,
        args.transcripts_fa,
        args.gtf,
        args.out,
        args.coding_cutoff
    )
