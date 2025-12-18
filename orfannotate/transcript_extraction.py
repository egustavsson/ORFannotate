from turtle import pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pyfaidx import Fasta
import gffutils 
import os

from orfannotate.utils import _create_db


def extract_transcripts_from_gtf(gtf_db_or_path, genome_fa, out_fasta):
 
    if isinstance(gtf_db_or_path, gffutils.FeatureDB):
        db = gtf_db_or_path
    else:
        if os.path.splitext(gtf_db_or_path)[1] == ".db":
            db_path = gtf_db_or_path
        else:
            db_path = gtf_db_or_path + ".db"
            if not os.path.exists(db_path):                
                # Create DB from GTF using GFFUtils
                _create_db(gtf_db_or_path)
        db = gffutils.FeatureDB(db_path)

    genome = Fasta(genome_fa, as_raw=True)

    tx_ids = [tx.id for tx in db.features_of_type("transcript")]

    exon_map = _fetch_exons_bulk(db, tx_ids)
    
    # ---- Stream FASTA writing ----
    with open(out_fasta, "w") as out:
        for tx_id, exons in exon_map.items():
            
            # Build concatenated exon sequence
            records = []
            for start, end, strand, seqid in exons:
                part = genome[seqid][start - 1:end]
                records.append(part)

            seq = "".join(records)
           
            # Reverse complement 
            if exons[0][2] == "-":
                seq = str(Seq(seq).reverse_complement())

            # Write one record
            rec = SeqRecord(Seq(seq), id=tx_id, description="")
            SeqIO.write(rec, out, "fasta")


def _fetch_exons_bulk(db, tx_ids):
    tx_str = ",".join(f'"{t}"' for t in tx_ids)

    query = f"""
        SELECT
            COALESCE(
                json_extract(attributes, '$.transcript_id[0]'),
                json_extract(attributes, '$.Parent[0]')
            ) AS transcript_id,
            seqid, start, end, strand
        FROM features
        WHERE featuretype = 'exon'
        AND COALESCE(
                json_extract(attributes, '$.transcript_id[0]'),
                json_extract(attributes, '$.Parent[0]')
            ) IN ({tx_str})
        ORDER BY transcript_id, start;
    """
    
    exon_map = {}
    for row in db.execute(query):
        tx, seqid, start, end, strand = row
        exon_map.setdefault(tx, []).append((start, end, strand, seqid))
    
    return exon_map


