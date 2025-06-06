from Bio import SeqIO
from Bio.Seq import Seq

def parse_orfs(orf_fasta):
    orf_data = {}
    for record in SeqIO.parse(orf_fasta, "fasta"):
        desc = record.description.split()
        transcript_id = desc[0]
        start, end = int(desc[1]), int(desc[2])
        strand = desc[3] if len(desc) > 3 else '+'
        seq = str(record.seq)
        aa_seq = str(Seq(seq).translate(to_stop=True))

        if transcript_id not in orf_data or (end - start) > orf_data[transcript_id]["length"]:
            orf_data[transcript_id] = {
                "start": start,
                "end": end,
                "length": end - start,
                "strand": strand,
                "nt_seq": seq,
                "aa_seq": aa_seq,
                "predicted_NMD": "NA"
            }
    return orf_data