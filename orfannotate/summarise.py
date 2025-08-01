import logging
import json
import tempfile
import argparse
from pathlib import Path
from typing import List, Tuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import gffutils

# Internal modules
from orfannotate.nmd import predict_nmd
from orfannotate.kozak import score_kozak

LOGGER = logging.getLogger(__name__)
logging.basicConfig(format="%(levelname)s:%(name)s:%(message)s", level=logging.INFO)

# Help functions

def _load_transcript_sequences(transcript_fasta: Path):
    return SeqIO.to_dict(SeqIO.parse(str(transcript_fasta), "fasta"))

def _map_junctions_to_tx(junctions_genomic, exons):
    offsets = []
    cum_len = 0
    for exon in exons:
        offsets.append((exon.start, exon.end, cum_len))
        cum_len += exon.end - exon.start + 1

    mapped = []
    for donor, acceptor in junctions_genomic:
        donor_tx = acceptor_tx = None
        for start, end, offset in offsets:
            if start <= donor <= end:
                donor_tx = offset + (donor - start + 1)
            if start <= acceptor <= end:
                acceptor_tx = offset + (acceptor - start + 1)
        # Handle boundary cases
        if donor_tx is None:
            for start, end, offset in offsets:
                if donor == end + 1:
                    donor_tx = offset + (end - start + 1)
        if acceptor_tx is None:
            for start, end, offset in offsets:
                if acceptor == start - 1:
                    acceptor_tx = offset + 1
        if donor_tx and acceptor_tx:
            mapped.append((donor_tx, acceptor_tx))
    return mapped

# Main function

def generate_summary(best_orfs, transcript_fa, gtf_db_or_path, output_path, coding_cutoff=0.364, junctions_by_tx=None):
    output_path = Path(output_path)
    output_dir = output_path.parent

    if junctions_by_tx is None:
        raise ValueError("junctions_by_tx dictionary must be provided for summarisation")

    # Load GTF database
    if isinstance(gtf_db_or_path, gffutils.FeatureDB):
        db = gtf_db_or_path
    else:
        gtf_db_or_path = Path(gtf_db_or_path)
        if gtf_db_or_path.suffix == ".db" and gtf_db_or_path.exists():
            db = gffutils.FeatureDB(str(gtf_db_or_path))
        else:
            with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmpdb:
                gffutils.create_db(
                    str(gtf_db_or_path),
                    dbfn=tmpdb.name,
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
                db = gffutils.FeatureDB(tmpdb.name)

    # Load transcript sequences
    transcript_seqs = _load_transcript_sequences(Path(transcript_fa))

    summary = []
    cds_records = []
    protein_records = []
    utr5_records = []
    utr3_records = []

    all_tx_ids = {t.id for t in db.features_of_type("transcript")}
    all_tx_ids &= set(transcript_seqs.keys())

    children = db.children

    for tid in tqdm(sorted(all_tx_ids), desc="Processing transcripts"):

        has_orf = tid in best_orfs
        orf_data = best_orfs.get(tid, None)

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

        full_seq = str(transcript_seqs[tid].seq)
        seq_len = len(full_seq)

        if coding_class == "coding":
            nt_seq = full_seq[orf_start - 1:orf_end_tx]
            aa_seq = str(Seq(nt_seq).translate(to_stop=True))
            orf_nt_len = orf_end_tx - orf_start + 1
            orf_aa_len = len(aa_seq)
            kozak_strength, kozak_seq = score_kozak(full_seq, orf_start)
        else:
            nt_seq = aa_seq = "NA"
            orf_nt_len = orf_aa_len = "NA"
            kozak_strength = kozak_seq = "NA"

        exons = list(children(tx, featuretype="exon", order_by="start"))
        if strand == '-':
            exons = exons[::-1]

        # Get junctions
        junctions_genomic = junctions_by_tx.get(tid, [])
        total_junctions = len(junctions_genomic)

        # Map to transcript coordinates
        junctions_tx = _map_junctions_to_tx(junctions_genomic, exons)

        # UTR classification
        if coding_class == "coding" and isinstance(orf_start, int) and isinstance(orf_end_tx, int):
            utr5_junctions = cds_junctions = utr3_junctions = 0
            for donor_tx, acceptor_tx in junctions_tx:
                if strand == '+':
                    if donor_tx <= orf_start - 1 and acceptor_tx <= orf_start - 1:
                        utr5_junctions += 1
                    elif donor_tx >= orf_end_tx + 1 and acceptor_tx >= orf_end_tx + 1:
                        utr3_junctions += 1
                    else:
                        cds_junctions += 1
                else:  # negative strand: 5' UTR is after ORF end
                    if donor_tx >= orf_end_tx + 1 and acceptor_tx >= orf_end_tx + 1:
                        utr5_junctions += 1
                    elif donor_tx <= orf_start - 1 and acceptor_tx <= orf_start - 1:
                        utr3_junctions += 1
                    else:
                        cds_junctions += 1
        else:
            utr5_junctions = cds_junctions = utr3_junctions = "NA"

        if has_orf and orf_end_gen is not None and junctions_genomic:
            last_j = junctions_genomic[-1][0] if strand == '+' else junctions_genomic[0][1]
            stop_to_last_ej = last_j - orf_end_gen if strand == '+' else orf_end_gen - last_j
        else:
            stop_to_last_ej = "NA"

        nmd_flag = predict_nmd(orf_end_gen, junctions_genomic, strand) if coding_class == "coding" else "FALSE"

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
            "kozak_sequence": kozak_seq,
        })

    pd.DataFrame(summary).to_csv(output_path, sep="\t", index=False)

    SeqIO.write(cds_records, output_dir / "cds.fa", "fasta")
    SeqIO.write(protein_records, output_dir / "protein.fa", "fasta")
    SeqIO.write(utr5_records, output_dir / "utr5.fa", "fasta")
    SeqIO.write(utr3_records, output_dir / "utr3.fa", "fasta")


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
