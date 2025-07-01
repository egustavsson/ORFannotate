from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pyfaidx import Fasta
import gffutils
import os


def extract_transcripts_from_gtf(gtf_db_or_path, genome_fa, out_fasta):
 
    if isinstance(gtf_db_or_path, gffutils.FeatureDB):
        db = gtf_db_or_path

    else:
        
        if os.path.splitext(gtf_db_or_path)[1] == ".db":
            db_path = gtf_db_or_path
        else:
            db_path = gtf_db_or_path + ".db"

            if not os.path.exists(db_path):
                
                gffutils.create_db(
                    gtf_db_or_path,
                    dbfn=db_path,
                    force=True,
                    keep_order=True,
                    disable_infer_transcripts=True,
                    disable_infer_genes=True,
                    merge_strategy="create_unique",
                    sort_attribute_values=True,
                    pragmas={
                        "journal_mode": "OFF",
                        "synchronous": "OFF",
                        "temp_store": "MEMORY",
                    },
                    id_spec={
                        "transcript": "transcript_id",
                        "exon": "transcript_id",
                        "CDS": "transcript_id",
                    },
                )

        db = gffutils.FeatureDB(db_path)

    genome = Fasta(genome_fa, as_raw=True)
    records = []

    for tx in db.features_of_type("transcript"):
        exons = list(db.children(tx, featuretype="exon", order_by="start"))
        seq = "".join(genome[exon.chrom][exon.start - 1 : exon.end] for exon in exons)

        if tx.strand == "-":
            seq = str(Seq(seq).reverse_complement())

        records.append(SeqRecord(Seq(seq), id=tx.id, description=""))

    SeqIO.write(records, out_fasta, "fasta")
