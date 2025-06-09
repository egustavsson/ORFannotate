import gffutils
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def extract_transcripts_from_gtf(gtf_path, genome_fa, out_fasta):
    db_path = gtf_path + ".db"
    gffutils.create_db(
        gtf_path,
        dbfn=db_path,
        force=True,
        keep_order=True,
        disable_infer_transcripts=True,
        disable_infer_genes=True,
        merge_strategy="merge",
        sort_attribute_values=True
    )
    db = gffutils.FeatureDB(db_path)
    genome = SeqIO.to_dict(SeqIO.parse(genome_fa, "fasta"))

    transcript_seqs = []
    for transcript in db.features_of_type("transcript"):
        exons = list(db.children(transcript, featuretype="exon", order_by="start"))
        seq = ""
        for exon in exons:
            chrom = exon.chrom
            start = exon.start - 1
            end = exon.end
            exon_seq = genome[chrom].seq[start:end]
            seq += exon_seq
        if transcript.strand == "-":
            seq = seq.reverse_complement()
        transcript_seqs.append(SeqRecord(Seq(seq), id=transcript.id, description=""))

    SeqIO.write(transcript_seqs, out_fasta, "fasta")

def annotate_gtf_with_cds(gtf_path, cds_features, output_path):
    from collections import defaultdict

    # bucket CDS by transcript
    cds_by_transcript = defaultdict(list)
    for feat in cds_features:
        cds_by_transcript[feat['attributes']['transcript_id']].append(feat)

    written_cds_for = set()                      # ← new

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
            attributes  = parts[8]
            tid = None
            for attr in attributes.split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    tid = attr.split(" ")[1].replace('"', '')
                    break

            fout.write(line) 

            if (
                feature_type == "transcript"      # pick your trigger
                and tid in cds_by_transcript
                and tid not in written_cds_for
            ):
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

if __name__ == "__main__":
    annotate_gtf_with_cds(
        gtf_path="annotated.gtf",
        cds_features=[],
        output_path="annotated_with_cds.gtf"
    )
