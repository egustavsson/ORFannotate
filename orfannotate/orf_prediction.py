import os
import subprocess
import logging
import pandas as pd

from orfannotate.orf_filter import get_best_orfs_by_cpat
from orfannotate.utils import _normalize_tid

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

def predict_uorf(output_dir, hexamer_path, logit_model_path, coding_cutoff, top_orf, min_length_uorf):
    """
    Run CPAT per-transcript only in 5'UTR sequences:
      - Ensure only non-overlapping uORFs above a certain protein cutoff are selected
      - Incorporate selected uORFs to the 'summary.tsv' file
    """
    # Canonical ORF 
    canonical_orf_cpat_dir = os.path.join(output_dir, "CPAT")
    canonical_orf_cpat_df = pd.read_csv(os.path.join(canonical_orf_cpat_dir, "cpat.ORF_prob.best.tsv"), sep="\t")

    # Create directory to store CPAT uORF results
    uorf_cpat_dir = os.path.join(output_dir, "CPAT_uORF")
    os.makedirs(uorf_cpat_dir, exist_ok=True)

    # Run CPAT in current 5'UTR sequences alone
    run_cpat(os.path.join(output_dir, "utr5.fa"), uorf_cpat_dir, hexamer_path, logit_model_path, top_orf, min_length_orf = min_length_uorf)
 
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
    selected_uorfs = {}
    for tx, uorf in best_uorf.items():

        orf_start = canon_orf_start.get(_normalize_tid(tx))

        if orf_start is None:
            # no canonical ORF found > skip
            selected_uorfs[tx] = None
            continue
   
        # Only keep uORFs ending before the canonical ORF start and with a coding probability above the cutoff
        if uorf["end"] < orf_start and uorf["coding_prob"] > coding_cutoff:
            selected_uorfs[tx] = uorf["coding_prob"]
        else:
            selected_uorfs[tx] = None

    # Load final summary.tsv file
    summary_tsv_path = os.path.join(output_dir, "ORFannotate_summary.tsv")
    orf_summary = pd.read_csv(summary_tsv_path, sep='\t')    

    # Update the final 'summary.tsv': add the new two uORF columns
    orf_summary["has_uORF"] = orf_summary["transcript_id"].map(
        lambda x: selected_uorfs.get(_normalize_tid(x)) is not None
    )
    orf_summary["coding_prob_best_uORF"] = orf_summary["transcript_id"].map(
        lambda x: selected_uorfs.get(_normalize_tid(x))
    )
   
    # Save results
    orf_summary.to_csv(summary_tsv_path,
                       sep="\t",
                       index=False,
                       na_rep="NA"
                       )