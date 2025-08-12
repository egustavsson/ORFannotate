import re
import gffutils
import logging
from collections import defaultdict, OrderedDict

logger = logging.getLogger(__name__)

# Help functions

def _build_tx_lookup(db):
    
    lut = {}
    for tx in db.features_of_type("transcript"):
        tid_full = tx.attributes["transcript_id"][0]
        lut[tid_full] = tx
        lut[tid_full.split(".")[0]] = tx
    return lut


def _extract_tid_from_attrs(attr_field: str):
    
    # GTF style
    m = re.search(r'\btranscript_id\s+"([^"]+)"', attr_field)
    if m:
        return m.group(1)

    # Very simple GFF3-style handling
    m = re.search(r'(?:^|;)\s*Parent=([^;]+)', attr_field)
    if m:
        parent = m.group(1)
        # strip optional transcript: prefix
        return re.sub(r'^(?:transcript:)?', '', parent)

    return None

# Main functions

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

            # No overlap with CDS in transcript coordinates
            if exon_tr_end < cds_start_tr or exon_tr_start > cds_end_tr:
                transcript_pos += exon_len
                continue

            # Overlap in transcript coords
            cds_exon_start_tr = max(exon_tr_start, cds_start_tr)
            cds_exon_end_tr = min(exon_tr_end, cds_end_tr)

            if tx.strand == "+":
                cds_start_gen = exon.start + (cds_exon_start_tr - exon_tr_start)
                cds_end_gen = exon.start + (cds_exon_end_tr - exon_tr_start)
            else:
                cds_end_gen = exon.end - (cds_exon_start_tr - exon_tr_start)
                cds_start_gen = exon.end - (cds_exon_end_tr - exon_tr_start)

            cds_features.append({
                "seqid": exon.seqid,
                "source": "ORFannotate",
                "feature": "CDS",
                "start": min(cds_start_gen, cds_end_gen),
                "end": max(cds_start_gen, cds_end_gen),
                "score": ".",
                "strand": tx.strand,
                "frame": ".",
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
    """
    Write a cleaned GTF that contains, per transcript:
      - transcript line (preserving original line)
      - all exon lines (sorted by genomic start)
      - all ORFannotate-predicted CDS lines (sorted by genomic start)

    Everything else from the original (start_codon, stop_codon, UTR, CDS from
    other sources, etc.) is removed. Transcripts with no predicted CDS are
    kept as transcript + exon only.
    """
    # Group predicted CDS by transcript_id
    cds_by_tid = defaultdict(list)
    for feat in cds_features:
        cds_by_tid[feat["attributes"]["transcript_id"]].append(feat)

    # Collect transcript + exon lines from original GTF
    transcript_lines = OrderedDict()
    exons_by_tid = defaultdict(list)

    with open(gtf_path, "rt", encoding="utf-8") as fin:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            ftype = parts[2]
            attrs = parts[8]
            tid = _extract_tid_from_attrs(attrs)
            if not tid:
                continue

            if ftype == "transcript":
                # keep first occurrence to preserve input order
                if tid not in transcript_lines:
                    transcript_lines[tid] = line
            elif ftype == "exon":
                exons_by_tid[tid].append(parts)

    # Write cleaned, sorted output
    with open(output_path, "wt", encoding="utf-8") as fout:
        for tid, t_line in transcript_lines.items():
            # 1) transcript
            fout.write(t_line)

            # 2) exons sorted by genomic start (col 4)
            exons = exons_by_tid.get(tid, [])
            for exon in sorted(exons, key=lambda p: int(p[3])):
                fout.write("\t".join(exon) + "\n")

            # 3) predicted CDS (if any), sorted by start
            if tid in cds_by_tid:
                for cds in sorted(cds_by_tid[tid], key=lambda x: int(x["start"])):
                    a = cds["attributes"]
                    attr_items = [
                        f'gene_id "{a.get("gene_id","")}"',
                        f'transcript_id "{a.get("transcript_id","")}"'
                    ]
                    if a.get("gene_name"):
                        attr_items.append(f'gene_name "{a["gene_name"]}"')
                    if a.get("ref_gene_id"):
                        attr_items.append(f'ref_gene_id "{a["ref_gene_id"]}"')

                    cds_line = [
                        cds["seqid"],
                        cds.get("source", "ORFannotate"),
                        "CDS",
                        str(cds["start"]),
                        str(cds["end"]),
                        cds.get("score", "."),
                        cds["strand"],
                        cds.get("frame", "."),
                        "; ".join(attr_items) + ";"
                    ]
                    fout.write("\t".join(cds_line) + "\n")
