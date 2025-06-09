import pandas as pd
import logging
import re
from collections import defaultdict

logger = logging.getLogger(__name__)


def _build_tx_lookup(db):
    """
    Create a dictionary that lets us fetch a transcript feature in
    O(1) time by either its full ID (including version) or the ID
    without the version suffix.
    """
    lut = {}
    for tx in db.features_of_type("transcript"):
        tid_full = tx.attributes["transcript_id"][0]
        lut[tid_full] = tx            # ENST00000619216.1
        lut[tid_full.split(".")[0]] = tx  # ENST00000619216
    return lut

def get_best_orfs_by_cpat(cpat_best_path, all_orfs_df=None, debug_output_path=None):
    df = pd.read_csv(cpat_best_path, sep='\t')

    # Strip ORF suffix to group by transcript
    df["base_id"] = df["ID"].str.replace(r'_ORF_\d+$', '', regex=True)

    # Keep only highest-scoring ORF per base transcript ID
    best_per_transcript = df.sort_values('Coding_prob', ascending=False).groupby('base_id').first()

    best_orfs = {}
    for base_id, row in best_per_transcript.iterrows():
        best_orfs[base_id] = {
            "start": row["ORF_start"],
            "end": row["ORF_end"],
            "strand": row["ORF_strand"],
            "frame": row["ORF_frame"],
            "length": row["ORF"],
            "coding_prob": row["Coding_prob"]
        }

    # Optional: write debug output
    if isinstance(all_orfs_df, pd.DataFrame) and debug_output_path:
        all_orfs_df.to_csv(debug_output_path, sep='\t', index=False)
    else:
        logging.warning("Skipping debug output: all_orfs_df is not a DataFrame")

    return best_orfs


def build_cds_features(gtf_db, best_orfs):
    """
    Convert the best CPAT ORF on each transcript to one CDS record per overlapping exon.
    """
    tx_lookup = _build_tx_lookup(gtf_db)
    cds_features = []

    for tid, orf in best_orfs.items():
        tx = tx_lookup.get(tid) or tx_lookup.get(tid.split('.')[0])
        if tx is None:
            logger.warning(f"[CDS] No transcript found for {tid}")
            continue

        exons = list(gtf_db.children(tx, featuretype="exon", order_by="start"))
        if tx.strand == "-":
            exons.reverse()

        cds_start_tr, cds_end_tr = sorted((int(orf["start"]), int(orf["end"])))
        transcript_pos = 1

        for exon in exons:
            exon_len = exon.end - exon.start + 1
            exon_tr_start = transcript_pos
            exon_tr_end   = transcript_pos + exon_len - 1

            if exon_tr_end < cds_start_tr or exon_tr_start > cds_end_tr:
                transcript_pos += exon_len
                continue

            cds_exon_start_tr = max(exon_tr_start, cds_start_tr)
            cds_exon_end_tr   = min(exon_tr_end,   cds_end_tr)

            if tx.strand == "+":
                cds_start_gen = exon.start + (cds_exon_start_tr - exon_tr_start)
                cds_end_gen   = exon.start + (cds_exon_end_tr   - exon_tr_start)
            else:
                cds_end_gen   = exon.end   - (cds_exon_start_tr - exon_tr_start)
                cds_start_gen = exon.end   - (cds_exon_end_tr   - exon_tr_start)

            frame = (cds_exon_start_tr - cds_start_tr) % 3

            cds_features.append({
                "seqid":  exon.seqid,
                "source": "ORFannotate",
                "feature": "CDS",
                "start":  min(cds_start_gen, cds_end_gen),
                "end":    max(cds_start_gen, cds_end_gen),
                "score":  ".",
                "strand": tx.strand,
                "frame":  str(frame),
                "attributes": {
                    "gene_id":       tx.attributes.get("gene_id", [""])[0],
                    "transcript_id": tid,
                    "gene_name":     tx.attributes.get("gene_name", [""])[0],
                    "ref_gene_id":   tx.attributes.get("ref_gene_id", [""])[0],
                },
            })
            transcript_pos += exon_len

    logger.debug("Built %d CDS records", len(cds_features))
    return cds_features
