import argparse
import gffutils
from Bio import SeqIO

def parse_orfipy_output(orf_file):
    """
    Parses ORFipy output to extract ORF coordinates.
    :param orf_file: Path to the predicted ORFs FASTA file.
    :return: Dictionary with transcript IDs as keys and ORF coordinates as values.
    """
    orf_data = {}
    for record in SeqIO.parse(orf_file, "fasta"):
        fields = record.description.split()
        transcript_id = fields[0]  # ORF identifier (linked to a transcript)
        start, end = int(fields[1]), int(fields[2])  # ORF coordinates

        if transcript_id not in orf_data:
            orf_data[transcript_id] = []
        orf_data[transcript_id].append((start, end))
    
    return orf_data

def annotate_gtf_with_cds(gtf_file, orf_file, output_gtf):
    """
    Annotates a GTF file with predicted CDS regions from ORFipy.
    :param gtf_file: Path to the input GTF file.
    :param orf_file: Path to the predicted ORFs FASTA file.
    :param output_gtf: Path to the output GTF file with CDS annotations.
    """
    # Load ORF predictions
    orf_data = parse_orfipy_output(orf_file)
    
    # Create GFF database
    db = gffutils.create_db(gtf_file, dbfn="transcripts.db", force=True, keep_order=True)

    with open(output_gtf, "w") as out_gtf:
        for feature in db.all_features():
            transcript_id = feature.attributes.get("transcript_id", [""])[0]
            
            # Keep exons unchanged
            out_gtf.write(str(feature) + "\n")

            # Add CDS features if found
            if feature.featuretype == "exon" and transcript_id in orf_data:
                for start, end in orf_data[transcript_id]:
                    if start >= feature.start and end <= feature.end:
                        # Create a new CDS feature
                        cds_feature = feature.copy()
                        cds_feature.featuretype = "CDS"
                        cds_feature.start = start
                        cds_feature.end = end
                        cds_feature.attributes["ID"] = [f"CDS_{transcript_id}"]
                        out_gtf.write(str(cds_feature) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate a GTF file with CDS features based on ORF predictions.")
    parser.add_argument("gtf_file", type=str, help="Path to the input GTF file")
    parser.add_argument("orf_file", type=str, help="Path to the predicted ORFs FASTA file")
    parser.add_argument("output_gtf", type=str, help="Path to the output annotated GTF file")

    args = parser.parse_args()
    annotate_gtf_with_cds(args.gtf_file, args.orf_file, args.output_gtf)
