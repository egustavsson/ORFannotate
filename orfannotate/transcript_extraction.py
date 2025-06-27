from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pyfaidx import Fasta
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
        sort_attribute_values=True,
        pragmas={"journal_mode": "OFF", "synchronous": "OFF"}
    )

    db = gffutils.FeatureDB(db_path)
    genome = Fasta(genome_fa, as_raw=True)

    transcript_seqs = []
    for transcript in db.features_of_type("transcript"):
        exons = list(db.children(transcript, featuretype="exon", order_by="start"))

        # Extract exon sequences using pyfaidx
        seq = "".join(genome[exon.chrom][exon.start - 1:exon.end] for exon in exons)
        if transcript.strand == "-":
            seq = str(Seq(seq).reverse_complement())

        record = SeqRecord(Seq(seq), id=transcript.id, description="")
        transcript_seqs.append(record)

    SeqIO.write(transcript_seqs, out_fasta, "fasta")
