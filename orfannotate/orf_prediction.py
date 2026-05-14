import os
import subprocess
import logging
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from gffutils.exceptions import FeatureNotFoundError



from orfannotate.orf_filter import get_best_orfs_by_cpat
from orfannotate.utils import _normalize_tid
from orfannotate.nmd import predict_nmd
from orfannotate.kozak import score_kozak
from orfannotate.gtf_annotation import build_cds_features, annotate_gtf_with_cds

logger = logging.getLogger(__name__)

def run_cpat(transcript_fasta, output_dir, hexamer_path, logit_model_path, top_orf, min_length_orf, force=False):
    output_prefix = os.path.join(output_dir, "cpat")
    output_file = output_prefix + ".ORF_prob.best.tsv"

    cpat_cmd = [
        "cpat.py",
        "-g", os.path.abspath(transcript_fasta),
        "-x", os.path.abspath(hexamer_path),
        "-d", os.path.abspath(logit_model_path),
        "--top-orf", str(top_orf),
        "--min-orf", str(min_length_orf),
        "--best-orf=p",
        "--log-file", os.path.join(output_dir, "CPAT.log"),
        "-o", output_prefix
    ]

    logger.info("Running CPAT...")
    with open(os.path.join(output_dir, "CPAT_run_info.log"), "w") as log_file:
        subprocess.run(
            cpat_cmd,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=True
        )

    # Clean up unneeded cpat.r
    cpat_r = os.path.join(output_dir, "cpat.r")
    if os.path.exists(cpat_r):
        os.remove(cpat_r)

    return output_file

def predict_uorf(output_dir, hexamer_path, logit_model_path, coding_cutoff, top_orf, min_length_uorf, gtf_db, gtf_path):
    """
    Run CPAT per-transcript only in 5'UTR sequences:
      - Ensure only non-overlapping uORFs above a certain protein cutoff are selected
      - Incorporate selected uORFs to the 'summary.tsv' file
      - Generate GTF file with selected uORFs
    """

    # Canonical ORF 
    canonical_orf_cpat_dir = os.path.join(output_dir, "CPAT")
    canonical_orf_cpat_df = pd.read_csv(os.path.join(canonical_orf_cpat_dir, "cpat.ORF_prob.best.tsv"), sep="\t")

    # Create directory to store CPAT uORF results
    uorf_cpat_dir = os.path.join(output_dir, "CPAT_uORF")
    os.makedirs(uorf_cpat_dir, exist_ok=True)

    # Run CPAT in current 5'UTR sequences alone
    utr5_fasta_path = os.path.join(output_dir, "utr5.fa")
    run_cpat(utr5_fasta_path, uorf_cpat_dir, hexamer_path, logit_model_path, top_orf, min_length_orf = min_length_uorf)
 
    # Get best uORF 
    uorf_cpat_results = os.path.join(uorf_cpat_dir, "cpat.ORF_prob.best.tsv")
    uorf_debug_output = os.path.join(uorf_cpat_dir, "cpat_debug.tsv")
    
    # Extract best uORF values per transcript
    best_uorf = get_best_orfs_by_cpat(uorf_cpat_results, debug_output_path=uorf_debug_output)

    # Lookup dict for quick canonical tx_id and ORF_start access
    canonical_orf_cpat_df["norm_seq_ID"] = canonical_orf_cpat_df["seq_ID"].map(_normalize_tid)
    
    canon_orf_start = (
        canonical_orf_cpat_df
        .set_index("norm_seq_ID")["ORF_start"]
        .to_dict()
    )
    
    # Filter uORFs so that only TRUE uORFs remain (end < canonical ORF start)
    # Initialize selected_uorfs for ALL transcripts
    selected_uorfs = {}

    all_transcripts = set(
        canonical_orf_cpat_df["norm_seq_ID"]
    ).union(set(best_uorf.keys()))

    # Load once before the loop
    transcript_seqs = {
        _normalize_tid(rec.id): rec
        for rec in SeqIO.parse(utr5_fasta_path, "fasta")
    }

    for tx in all_transcripts:

        selected_uorfs[tx] = {
            "coding_prob": None,
            "nmd_flag": "NA",
            "kozak_strength": "NA",
            "kozak_seq": "NA"
        }

        uorf = best_uorf.get(tx)
        orf_start = canon_orf_start.get(tx)

        # Skip if no canonical ORF
        if orf_start is None or uorf is None:
            continue

        # Apply filtering
        if uorf["end"] < orf_start and uorf["coding_prob"] > coding_cutoff:

            selected_uorfs[tx]["coding_prob"] = uorf["coding_prob"]
            

            # ---- KOZAK prediction for uORF ----
            try:
                seq_record = transcript_seqs.get(_normalize_tid(tx))
                
                if seq_record is not None:
                    # Get the uORF lengt
                    uorf_start = int(uorf["start"])  
                    selected_uorfs[tx]["uorf_length"] = abs(int(uorf["start"]) - int(uorf["end"]) )

                    # Get the uORF sequence
                    full_seq = str(seq_record.seq)
                    kozak_strength, kozak_seq = score_kozak(full_seq, uorf_start)
                    
                    selected_uorfs[tx]["kozak_strength"] = kozak_strength
                    selected_uorfs[tx]["kozak_seq"] = kozak_seq

                else:
                    selected_uorfs[tx]["kozak_strength"] = None
                    selected_uorfs[tx]["kozak_seq"] = None

            except Exception as e:
                selected_uorfs[tx]["kozak_strength"] = None
                selected_uorfs[tx]["kozak_seq"] = None
          
        
                
    # ---- GTF -------------
    filtered_best_uorf = {
        tx: uorf for tx, uorf in best_uorf.items()
        if selected_uorfs.get(tx, {}).get("coding_prob") is not None
    }
    cds_features = build_cds_features(gtf_db, filtered_best_uorf)

    # Group CDS features by transcript for NMD prediction
    cds_by_tx = defaultdict(list)
    for f in cds_features:
        if f["feature"] == "CDS":
            cds_by_tx[f["attributes"]["transcript_id"]].append(f)

    # ---- NMD analysis ----
    for tx, uorf_cds_only in cds_by_tx.items():
        tx_norm = _normalize_tid(tx)  # ensure consistent key
        try:
            db_tx = gtf_db[tx]
            strand = db_tx.strand

            exon_features = [
                f for f in gtf_db.children(db_tx, featuretype="exon", order_by="start")
                if f.attributes.get("transcript_id", [None])[0] == tx
            ]

            if strand == "+":
                first_cds_start = uorf_cds_only[0]["start"]
            else:
                first_cds_start = uorf_cds_only[-1]["end"]

            selected_uorfs[tx_norm]["nmd_flag"] = predict_nmd(
                first_cds_start, exon_features, strand
            )

        except (KeyError, FeatureNotFoundError):
            print(f"Transcript {tx} not found in GTF database for NMD prediction. Setting NMD flag to 'NA'.")
            pass

    annotated_gtf = os.path.join(output_dir, "uORFannotate_annotated.gtf")
    annotate_gtf_with_cds(gtf_path, cds_features, annotated_gtf)

    
    

    # Load final summary.tsv file
    summary_tsv_path = os.path.join(output_dir, "ORFannotate_summary.tsv")
    # Set column types for the final tsv summary
    orf_summary = pd.read_csv(
        summary_tsv_path,
        sep='\t',
        dtype={
            "chrom": str,
            "gene_id": str,
            "transcript_id": str,
            "kozak_strength": str,
            "kozak_sequence": str,
        }
    )

    # Update the final 'summary.tsv': add the new two uORF columns
    # Add has_uORF flag (based on coding_prob)
    orf_summary["norm_tid"] = orf_summary["transcript_id"].map(_normalize_tid)

    orf_summary["has_uORF"] = orf_summary["norm_tid"].map(
        lambda x: selected_uorfs.get(x, {}).get("coding_prob") is not None
    )

    orf_summary["coding_prob_best_uORF"] = orf_summary["norm_tid"].map(
        lambda x: selected_uorfs.get(x, {}).get("coding_prob")
    )

    orf_summary["length_best_uORF"] = orf_summary["norm_tid"].map(
        lambda x: selected_uorfs.get(x, {}).get("uorf_length")
    )

    orf_summary["nmd_sensitivity_best_uORF"] = orf_summary["norm_tid"].map(
        lambda x: selected_uorfs.get(x, {}).get("nmd_flag", "NA")
    )
    orf_summary["kozak_strength_best_uORF"] = orf_summary["norm_tid"].map(
        lambda x: selected_uorfs.get(x, {}).get("kozak_strength")
    )

    orf_summary["kozak_seq_best_uORF"] = orf_summary["norm_tid"].map(
        lambda x: selected_uorfs.get(x, {}).get("kozak_seq")
    )

    orf_summary.drop(columns=["norm_tid"], inplace=True)
   
    # Save results
    orf_summary.to_csv(summary_tsv_path,
                       sep="\t",
                       index=False,
                       na_rep="NA"
                       )
    