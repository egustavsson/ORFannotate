import logging
import json
import tempfile
import argparse
import os
import datetime
from pathlib import Path
from typing import List, Tuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffutils
from rich.progress import Progress

# Internal modules
from orfannotate.nmd import predict_nmd
from orfannotate.kozak import score_kozak
from orfannotate.utils import _load_transcript_sequences, _map_junctions_to_tx, _create_db, _normalize_tid

LOGGER = logging.getLogger(__name__)
logging.basicConfig(format="%(levelname)s:%(name)s:%(message)s", level=logging.INFO)

# Main function

def generate_summary(best_orfs, transcript_fa, gtf_db_or_path, output_path, output_dir, progress, task, STEP_UNITS, coding_cutoff=0.364):
    output_path = Path(output_path)
    output_dir = output_path.parent
   
    # if junctions_all is None:
    #     raise ValueError("junctions_by_tx dictionary must be provided for summarisation")

    # Load GTF database
    if isinstance(gtf_db_or_path, gffutils.FeatureDB):
        db = gtf_db_or_path
    else:
        db =_create_db(str(Path(gtf_db_or_path)), only_exons=False)
    
    # Load transcript sequences from FASTA
    transcript_seqs = _load_transcript_sequences(Path(transcript_fa))
    
    summary = []
    cds_records = []
    protein_records = []
    utr5_records = []
    utr3_records = []

    fasta_ids_norm = {_normalize_tid(tid) for tid in transcript_seqs.keys()}
    
    all_tx_ids = {
        tid for tid in (t.id for t in db.features_of_type("transcript"))
        if _normalize_tid(tid) in fasta_ids_norm
    }
    
    children = db.children
    
    inc=(STEP_UNITS*3)/len(all_tx_ids)

    for tid in sorted(all_tx_ids): 
        
        tid_base = tid.split(".")[0]
        tid_norm = _normalize_tid(tid)
        tid_base_norm = _normalize_tid(tid_base)
        
        orf_data = (
            best_orfs.get(tid)
            or best_orfs.get(tid_base)
            or best_orfs.get(tid_norm)
            or best_orfs.get(tid_base_norm)
        )
        
        has_orf = orf_data is not None
        
        try:
            tx = db[tid]
        except KeyError:
            tx = db.get(tid.split(".")[0], None)
            if tx is None:
                continue
        
        chrom = tx.chrom
        strand = tx.strand
        tx_start = tx.start
        tx_end = tx.end
        gene_id = tx.attributes.get("gene_id", ["NA"])[0]

        if has_orf:
            orf_start = orf_data["start"]
            orf_end_tx = orf_data["end"]
            coding_prob = orf_data["coding_prob"]
        else:
            orf_start = orf_end_tx = "NA"
            coding_prob = "NA"
    
        # Get all CDS exons
        cds_feats = list(children(tx, featuretype="CDS", order_by="start"))
        
        orf_end_gen = None
        if cds_feats:
            if strand == "+":
                orf_end_gen = max(f.end for f in cds_feats)
            else:
                orf_end_gen = min(f.start for f in cds_feats)
       
        coding_class = (
            "coding"
            if (has_orf and isinstance(coding_prob, (float, int)) and coding_prob >= coding_cutoff and orf_end_gen)
            else "noncoding"
        )
        
        seq_record = (
            transcript_seqs.get(tid)
            or transcript_seqs.get(tid_base)
            or transcript_seqs.get(tid_norm)
            or transcript_seqs.get(tid_base_norm)
        )
        
        if seq_record is None:
            continue
        
        full_seq = str(seq_record.seq)
        
        seq_len = len(full_seq)

        if coding_class == "coding":
            nt_seq = full_seq[int(orf_start) - 1:orf_end_tx]
            aa_seq = str(Seq(nt_seq).translate(to_stop=True))
            orf_nt_len = orf_end_tx - orf_start + 1
            orf_aa_len = len(aa_seq)
            kozak_strength, kozak_seq = score_kozak(full_seq, orf_start)
        else:
            nt_seq = aa_seq = "NA"
            orf_nt_len = orf_aa_len = "NA"
            kozak_strength = kozak_seq = "NA"

        
        if has_orf and orf_end_gen is not None and cds_feats:
            last_j = cds_feats[-1].end if strand == '+' else cds_feats[0].start
            stop_to_last_ej = last_j - orf_end_gen if strand == '+' else orf_end_gen - last_j
        else:
            stop_to_last_ej = "NA"

        nmd_flag = predict_nmd(orf_end_gen, cds_feats, strand) if coding_class == "coding" else "FALSE"
        

        # UTR junction count 
        if coding_class == "coding":
            utr5_exons = list(children(tx, featuretype="five_prime_utr"))
            utr3_exons = list(children(tx, featuretype="three_prime_utr"))
            utr5_junctions = len(utr5_exons)-1 if len(utr5_exons) > 1 else 0
            utr3_junctions = len(utr3_exons)-1 if len(utr3_exons) > 1 else 0
            cds_junctions = len(cds_feats)-1 if len(cds_feats) > 1 else 0
            total_junctions = (utr5_junctions + cds_junctions + utr3_junctions)
        else:
            utr5_junctions = cds_junctions = utr3_junctions = "NA"
            total_junctions = len(list(children(tx, featuretype="exon")))
        
        
        # UTR length calculation
        utr5_seq = utr3_seq = "NA"
        utr5_len = utr3_len = "NA"
        if coding_class == "coding":
            if orf_start > 1:
                utr5_seq = full_seq[:orf_start - 1]
                utr5_len = len(utr5_seq)
            if orf_end_tx < seq_len:
                utr3_seq = full_seq[orf_end_tx:]
                utr3_len = len(utr3_seq)

            desc = f"gene_id={gene_id};coding_prob={coding_prob}"
            cds_records.append(SeqRecord(Seq(nt_seq), id=tid, description=desc))
            protein_records.append(SeqRecord(Seq(aa_seq), id=tid, description=desc))

            if utr5_len != "NA":
                utr5_records.append(SeqRecord(Seq(utr5_seq), id=tid, description=desc))
            if utr3_len != "NA":
                utr3_records.append(SeqRecord(Seq(utr3_seq), id=tid, description=desc))
        
        
        
        summary.append({
            "transcript_id": tid,
            "gene_id": gene_id,
            "chrom": chrom,
            "strand": strand,
            "transcript_start": tx_start,
            "transcript_end": tx_end,
            "has_orf": "TRUE" if has_orf else "FALSE",
            "orf_nt_len": orf_nt_len,
            "orf_aa_len": orf_aa_len,
            "utr5_nt_len": utr5_len,
            "utr3_nt_len": utr3_len,
            "coding_prob": coding_prob,
            "coding_class": coding_class,
            "total_junctions": total_junctions,
            "utr5_junctions": utr5_junctions,
            "cds_junctions": cds_junctions,
            "utr3_junctions": utr3_junctions,
            "stop_to_last_EJ": stop_to_last_ej,
            "NMD_sensitive": nmd_flag,
            "kozak_strength": kozak_strength,
            "kozak_sequence": kozak_seq
        })
        progress.update(task, advance=inc)

    pd.DataFrame(summary).to_csv(output_path, sep="\t", index=False)

    SeqIO.write(cds_records, output_dir / "cds.fa", "fasta")
    SeqIO.write(protein_records, output_dir / "protein.fa", "fasta")
    SeqIO.write(utr5_records, output_dir / "utr5.fa", "fasta")
    SeqIO.write(utr3_records, output_dir / "utr3.fa", "fasta")

    return utr5_records

if __name__ == "__main__":

    cli = argparse.ArgumentParser(description="Write ORF summary TSV")
    cli.add_argument("--best-orfs-json", required=True)
    cli.add_argument("--transcripts-fa", required=True)
    cli.add_argument("--gtf", required=True)
    cli.add_argument("--out", required=True)
    cli.add_argument("--coding-cutoff", type=float, default=0.364)
    args = cli.parse_args()

    best_orfs_dict = json.load(open(args.best_orfs_json))
    raise ValueError("This script expects in-memory junctions when run directly. Use ORFannotate.py instead.")
