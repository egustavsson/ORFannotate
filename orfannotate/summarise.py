from __future__ import annotations
import logging
import json
import tempfile
from pathlib import Path
from typing import List, Tuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffutils

# Internal modules
from orfannotate.nmd import predict_nmd
from orfannotate.kozak import score_kozak

LOGGER = logging.getLogger(__name__)
logging.basicConfig(format="%(levelname)s:%(name)s:%(message)s", level=logging.INFO)

# Help functions
def _load_transcript_sequences(transcript_fasta: Path):
    return SeqIO.to_dict(SeqIO.parse(str(transcript_fasta), "fasta"))

# Main function
def generate_summary(best_orfs, transcript_fa, gtf_db_or_path, output_path, coding_cutoff=0.364):
    output_path = Path(output_path)
    output_dir = output_path.parent

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

    transcript_seqs = _load_transcript_sequences(Path(transcript_fa))
    summary = []
    cds_records = []
    protein_records = []
    utr5_records = []
    utr3_records = []

    all_tx_ids = {t.id for t in db.features_of_type("transcript")}
    all_tx_ids &= set(transcript_seqs.keys())

    for tid in sorted(all_tx_ids):
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
            frame = orf_data["frame"]
            coding_prob = orf_data["coding_prob"]
        else:
            orf_start = orf_end_tx = frame = "NA"
            coding_prob = "NA"

        cds_feats = list(db.children(tx, featuretype="CDS", order_by="start"))
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

        if coding_class == "coding":
            full_seq = str(transcript_seqs[tid].seq)
            nt_seq = full_seq[orf_start - 1:orf_end_tx]
            aa_seq = str(Seq(nt_seq).translate(to_stop=True))
            orf_nt_len = orf_end_tx - orf_start + 1
            orf_aa_len = len(aa_seq)
            kozak_strength, kozak_seq = score_kozak(full_seq, orf_start)
        else:
            nt_seq = aa_seq = "NA"
            orf_nt_len = orf_aa_len = "NA"
            kozak_strength = kozak_seq = "NA"
            kozak_strength = kozak_seq = "NA"

        cds_feats = list(db.children(tx, featuretype="CDS", order_by="start"))
        orf_end_gen = None
        if cds_feats:
            if strand == "+":
                orf_end_gen = max(f.end for f in cds_feats)
            else:
                orf_end_gen = min(f.start for f in cds_feats)

        exons = list(db.children(tx, featuretype="exon", order_by="start"))
        if strand == '-':
            exons = exons[::-1]

        junctions: List[Tuple[int, int]] = [
            (exons[i].end, exons[i + 1].start) for i in range(len(exons) - 1)
        ]
        junction_count = len(junctions)

        if has_orf and orf_end_gen is not None and junctions:
            last_j = junctions[-1][0] if strand == '+' else junctions[0][1]
            stop_to_last_ej = last_j - orf_end_gen if strand == '+' else orf_end_gen - last_j
        else:
            stop_to_last_ej = "NA"

        coding_class = (
            "coding"
            if (has_orf and isinstance(coding_prob, (float, int)) and coding_prob >= coding_cutoff and orf_end_gen)
            else "noncoding"
        )
        nmd_flag = predict_nmd(orf_end_gen, junctions, strand) if coding_class == "coding" else "FALSE"

        utr5_seq = utr3_seq = "NA"
        utr5_len = utr3_len = "NA"
        if coding_class == "coding":
            if orf_start > 1:
                utr5_seq = full_seq[:orf_start - 1]
                utr5_len = len(utr5_seq)
            if orf_end_tx < len(full_seq):
                utr3_seq = full_seq[orf_end_tx:]
                utr3_len = len(utr3_seq)

            desc = f"gene_id={gene_id};coding_prob={coding_prob}"
            cds_records.append(SeqRecord(Seq(nt_seq), id=tid, description=desc))
            protein_records.append(SeqRecord(Seq(aa_seq), id=tid, description=desc))

            if utr5_seq != "NA":
                utr5_records.append(SeqRecord(Seq(utr5_seq), id=tid, description=desc))
            if utr3_seq != "NA":
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
            "coding_prob": coding_prob,
            "coding_class": coding_class,
            "junction_count": junction_count,
            "stop_to_last_EJ": stop_to_last_ej,
            "NMD_sensitive": nmd_flag,
            "kozak_strength": kozak_strength,
            "kozak_sequence": kozak_seq,
            "utr5_nt_len": utr5_len,
            "utr3_nt_len": utr3_len,
        })

    pd.DataFrame(summary).to_csv(output_path, sep="\t", index=False)

    SeqIO.write(cds_records, output_dir / "cds.fa", "fasta")
    SeqIO.write(protein_records, output_dir / "protein.fa", "fasta")
    SeqIO.write(utr5_records, output_dir / "utr5.fa", "fasta")
    SeqIO.write(utr3_records, output_dir / "utr3.fa", "fasta")


if __name__ == "__main__":
    import argparse

    cli = argparse.ArgumentParser(description="Write ORF summary TSV")
    cli.add_argument("--best-orfs-json", required=True)
    cli.add_argument("--transcripts-fa", required=True)
    cli.add_argument("--gtf", required=True)
    cli.add_argument("--out", required=True)
    cli.add_argument("--coding-cutoff", type=float, default=0.364)
    args = cli.parse_args()

    best_orfs_dict = json.load(open(args.best_orfs_json))
    generate_summary(
        best_orfs_dict,
        args.transcripts_fa,
        args.gtf,
        args.out,
        args.coding_cutoff
    )
