import gffutils
from Bio import SeqIO

def parse_orfipy_output(orf_file):
    orf_data = {}
    for record in SeqIO.parse(orf_file, "fasta"):
        fields = record.description.split()
        transcript_id = fields[0]
        start, end = int(fields[1]), int(fields[2])

        if transcript_id not in orf_data:
            orf_data[transcript_id] = []
        orf_data[transcript_id].append((start, end))
    return orf_data

def annotate_gtf(gtf_file, orf_file, output_gtf):
    orf_data = parse_orfipy_output(orf_file)
    db = gffutils.create_db(gtf_file, dbfn="transcripts.db", force=True, keep_order=True)

    with open(output_gtf, "w") as out_gtf:
        for feature in db.all_features():
            transcript_id = feature.attributes.get("transcript_id", [""])[0]
            out_gtf.write(str(feature) + "\n")

            if feature.featuretype == "exon" and transcript_id in orf_data:
                for start, end in orf_data[transcript_id]:
                    if start >= feature.start and end <= feature.end:
                        cds_feature = feature.copy()
                        cds_feature.featuretype = "CDS"
                        cds_feature.start = start
                        cds_feature.end = end
                        cds_feature.attributes["ID"] = [f"CDS_{transcript_id}"]
                        out_gtf.write(str(cds_feature) + "\n")