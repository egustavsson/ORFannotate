import re
import gffutils
import logging
from pathlib import Path
from collections import defaultdict, OrderedDict

from orfannotate.utils import _build_tx_lookup, _extract_tid_from_attrs, TRANSCRIPT_LIKE_FEATURES, _normalize_tid

logger = logging.getLogger(__name__)


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
        tid_base = tid.split(".")[0]
        tx = (
            tx_lookup.get(tid)
            or tx_lookup.get(tid_base)
            or tx_lookup.get(_normalize_tid(tid))
            or tx_lookup.get(_normalize_tid(tid_base))
        )
        
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
            
            # Overlap in transcript coords
            cds_exon_start_tr = max(exon_tr_start, cds_start_tr)
            cds_exon_end_tr = min(exon_tr_end, cds_end_tr)

            if tx.strand == "+":
                cds_start_gen = exon.start + (cds_exon_start_tr - exon_tr_start)
                cds_end_gen = exon.start + (cds_exon_end_tr - exon_tr_start)
            else:
                cds_start_gen = exon.end - (cds_exon_end_tr - exon_tr_start)
                cds_end_gen = exon.end - (cds_exon_start_tr - exon_tr_start)
                

            common_attrs = {
                "gene_id":      tx.attributes.get("gene_id", [""])[0],
                "transcript_id": tx.attributes.get("transcript_id", [tid])[0],
                "gene_name":    tx.attributes.get("gene_name", [""])[0],
                "ref_gene_id":  tx.attributes.get("ref_gene_id", [""])[0],
            }
            
           
            # Full UTR exon & CDS definition
            if exon_tr_end < cds_start_tr:
                feature = "five_prime_utr"
                start_coord = min(exon.start, exon.end)
                end_coord = max(exon.start, exon.end)
            elif exon_tr_start > cds_end_tr:
                feature = "three_prime_utr"
                start_coord = min(exon.start, exon.end)
                end_coord = max(exon.start, exon.end)
            else:
                feature = "CDS"
                start_coord = min(cds_start_gen, cds_end_gen)
                end_coord   = max(cds_start_gen, cds_end_gen)

                # 5' partial UTR (upstream of CDS start within this exon)
                if tx.strand == "+" and exon.start < start_coord:
                    cds_features.append({
                        "seqid": exon.seqid, "source": "ORFannotate",
                        "feature": "five_prime_utr",
                        "start": exon.start, "end": start_coord - 1,
                        "score": ".", "strand": tx.strand, "frame": ".",
                        "attributes": common_attrs,
                    })
                elif tx.strand == "-" and exon.end > end_coord:
                    cds_features.append({
                        "seqid": exon.seqid, "source": "ORFannotate",
                        "feature": "five_prime_utr",
                        "start": end_coord + 1, "end": exon.end,
                        "score": ".", "strand": tx.strand, "frame": ".",
                        "attributes": common_attrs,
                    })

                # 3' partial UTR (downstream of CDS end within this exon)
                if tx.strand == "+" and exon.end > end_coord:
                    cds_features.append({
                        "seqid": exon.seqid, "source": "ORFannotate",
                        "feature": "three_prime_utr",
                        "start": end_coord + 1, "end": exon.end,
                        "score": ".", "strand": tx.strand, "frame": ".",
                        "attributes": common_attrs,
                    })
                elif tx.strand == "-" and exon.start < start_coord:
                    cds_features.append({
                        "seqid": exon.seqid, "source": "ORFannotate",
                        "feature": "three_prime_utr",
                        "start": exon.start, "end": start_coord - 1,
                        "score": ".", "strand": tx.strand, "frame": ".",
                        "attributes": common_attrs,
                    })

            
            cds_features.append({
                    "seqid": exon.seqid,
                    "source": "ORFannotate",
                    "feature": feature,
                    "start": start_coord,
                    "end": end_coord,
                    "score": ".",
                    "strand": tx.strand,
                    "frame": ".",
                    "attributes": common_attrs,
                })
            transcript_pos += exon_len

    logger.debug("Built %d CDS records", len(cds_features))
    cds_features.sort(key=lambda f: (f["start"], f["end"]))

    return cds_features



def annotate_gtf_with_cds(gtf_path, cds_features, output_path, restrict_to_tids=None):
    """
    Write a cleaned GTF that contains, per transcript:
      - transcript line (from the original file)
      - exons + ORFannotate-predicted CDS interleaved by genomic start
    """
    # Group predicted CDS by transcript_id
    cds_by_tid = defaultdict(list)
    for feat in cds_features:
        cds_by_tid[feat["attributes"]["transcript_id"]].append(feat)
    
    # Collect transcript + exon lines from original GTF (preserve transcript order)
    transcript_lines = OrderedDict()
    exons_by_tid = defaultdict(list)
    
    if Path(gtf_path).suffix.lower() in {".gff3", ".gff"}:
        gtf_path = "tmp/converted.cleaned.gtf"

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

            if ftype in TRANSCRIPT_LIKE_FEATURES:
                
                if tid not in transcript_lines:
                    transcript_lines[tid] = line
            elif ftype == "exon":
                exons_by_tid[tid].append(parts)

    # Write cleaned, interleaved output
    with open(output_path, "wt", encoding="utf-8") as fout:
        for tid, t_line in transcript_lines.items():

            # If a filter set is provided, skip transcripts not in it
            if restrict_to_tids is not None and tid not in restrict_to_tids:
                continue
            
            # 1) transcript header
            fout.write(t_line)

            # 2) build a combined list of exon + predicted CDS rows, then sort by start
            combined_rows = []

            # add exons from original (tag so we can break ties consistently)
            for exon in exons_by_tid.get(tid, []):
                combined_rows.append((
                    int(exon[3]),
                    0,
                    "\t".join(exon) + "\n"
                ))
            
            # add predicted CDS (converted to GTF row strings)
            for cds in cds_by_tid.get(tid, []):
                
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
                    str(cds["feature"]),
                    str(cds["start"]),
                    str(cds["end"]),
                    cds.get("score", "."),
                    cds["strand"],
                    cds.get("frame", "."),
                    "; ".join(attr_items) + ";"
                ]
                combined_rows.append((
                    int(cds["start"]),
                    1,
                    "\t".join(cds_line) + "\n"
                ))
          
            # sort and write
            for _, _, row in sorted(combined_rows, key=lambda x: (x[0], x[1])):
                fout.write(row)
