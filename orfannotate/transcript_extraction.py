from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffutils

def extract_transcripts_from_gtf(gtf_path, genome_fa, out_fasta):
    db_path = gtf_path + ".db"
    gffutils.create_db(
        gtf_path,
        dbfn=db_path,
        force=True,
        keep_order=True,
        disable_infer_transcripts=True,
        disable_infer_genes=True,
        merge_strategy="create_unique",
        sort_attribute_values=True
    )
    db = gffutils.FeatureDB(db_path)
    genome = SeqIO.to_dict(SeqIO.parse(genome_fa, "fasta"))

    transcript_seqs = []
    for transcript in db.features_of_type("transcript"):
        exons = list(db.children(transcript, featuretype="exon", order_by="start"))
        seq = Seq("".join(str(genome[exon.chrom].seq[exon.start - 1:exon.end]) for exon in exons))
        if transcript.strand == "-":
            seq = seq.reverse_complement()
        transcript_seqs.append(SeqRecord(seq, id=transcript.id, description=""))

    SeqIO.write(transcript_seqs, out_fasta, "fasta")
