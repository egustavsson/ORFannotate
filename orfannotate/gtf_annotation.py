import gffutils
import logging
from collections import defaultdict

logger = logging.getLogger(__name__)

# Create a dictionary that lets us fetch a transcript feature in
# O(1) time by either its full ID (including version) or the ID
# without the version suffix.
def _build_tx_lookup(db):
    
    lut = {}
    for tx in db.features_of_type("transcript"):
        tid_full = tx.attributes["transcript_id"][0]
        lut[tid_full] = tx
        lut[tid_full.split(".")[0]] = tx
    return lut

# Convert the best CPAT ORF on each transcript to one CDS record per overlapping exon.
def build_cds_features(gtf_db, best_orfs):
    """
    Converts best CPAT ORFs to CDS features aligned to exon structures in a GTF database.

    Parameters:
    - gtf_db: gffutils FeatureDB created from input GTF
    - best_orfs: dict mapping transcript_id to ORF info (start, end, strand, etc.)

    Returns:
    - List of CDS feature dicts, suitable for annotation in GTF.
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
            exon_tr_end = transcript_pos + exon_len - 1

            if exon_tr_end < cds_start_tr or exon_tr_start > cds_end_tr:
                transcript_pos += exon_len
                continue

            cds_exon_start_tr = max(exon_tr_start, cds_start_tr)
            cds_exon_end_tr = min(exon_tr_end, cds_end_tr)

            if tx.strand == "+":
                cds_start_gen = exon.start + (cds_exon_start_tr - exon_tr_start)
                cds_end_gen = exon.start + (cds_exon_end_tr - exon_tr_start)
            else:
                cds_end_gen = exon.end - (cds_exon_start_tr - exon_tr_start)
                cds_start_gen = exon.end - (cds_exon_end_tr - exon_tr_start)

            frame = (cds_exon_start_tr - cds_start_tr) % 3

            cds_features.append({
                "seqid": exon.seqid,
                "source": "ORFannotate",
                "feature": "CDS",
                "start": min(cds_start_gen, cds_end_gen),
                "end": max(cds_start_gen, cds_end_gen),
                "score": ".",
                "strand": tx.strand,
                "frame": str(frame),
                "attributes": {
                    "gene_id": tx.attributes.get("gene_id", [""])[0],
                    "transcript_id": tid,
                    "gene_name": tx.attributes.get("gene_name", [""])[0],
                    "ref_gene_id": tx.attributes.get("ref_gene_id", [""])[0],
                },
            })
            transcript_pos += exon_len

    logger.debug("Built %d CDS records", len(cds_features))
    return cds_features



def annotate_gtf_with_cds(gtf_path, cds_features, output_path):
    cds_by_transcript = defaultdict(list)
    for feat in cds_features:
        cds_by_transcript[feat['attributes']['transcript_id']].append(feat)

    written_cds_for = set()

    with open(gtf_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            stripped = line.strip()
            if stripped == "" or stripped.startswith("#"):
                fout.write(line)
                continue

            parts = stripped.split("\t")
            if len(parts) != 9:
                fout.write(line)
                continue

            feature_type = parts[2]
            attributes = parts[8]
            tid = None
            for attr in attributes.split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    tid = attr.split(" ")[1].replace('"', '')
                    break

            fout.write(line)

            if feature_type == "transcript" and tid in cds_by_transcript and tid not in written_cds_for:
                for cds in sorted(cds_by_transcript[tid], key=lambda x: x['start']):
                    cds_attrs = [
                        f'gene_id "{cds["attributes"]["gene_id"]}"',
                        f'transcript_id "{cds["attributes"]["transcript_id"]}"'
                    ]
                    if cds["attributes"].get("gene_name"):
                        cds_attrs.append(f'gene_name "{cds["attributes"]["gene_name"]}"')
                    if cds["attributes"].get("ref_gene_id"):
                        cds_attrs.append(f'ref_gene_id "{cds["attributes"]["ref_gene_id"]}"')

                    cds_line = [
                        cds['seqid'], cds['source'], cds['feature'],
                        str(cds['start']), str(cds['end']), cds['score'],
                        cds['strand'], cds['frame'],
                        "; ".join(cds_attrs) + ";"
                    ]
                    fout.write("\t".join(cds_line) + "\n")

                written_cds_for.add(tid)
