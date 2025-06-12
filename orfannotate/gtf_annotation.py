import gffutils
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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
