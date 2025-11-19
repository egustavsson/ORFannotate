import pandas as pd
import logging
import re
from collections import defaultdict

logger = logging.getLogger(__name__)

def get_best_orfs_by_cpat(cpat_best_path, coding_cutoff, all_orfs_df=None, debug_output_path=None):
    df = pd.read_csv(cpat_best_path, sep='\t')
    
    # Strip ORF suffix to group by transcript
    df["base_id"] = df["ID"].str.replace(r'_ORF_\d+$', '', regex=True)
    
    # Keep only highest-scoring ORF per base transcript ID
    best_per_transcript = df.sort_values('Coding_prob', ascending=False).groupby('base_id').first()

    best_orfs = {}
    for base_id, row in best_per_transcript.iterrows():  

        best_orfs[base_id] = {
            "start": row["ORF_start"],
            "end": row["ORF_end"],
            "strand": row["ORF_strand"],
            "frame": row["ORF_frame"],
            "length": row["ORF"],
            "coding_prob": row["Coding_prob"]
        }
        ## SG ---------------------------------
        # Filter the df by the current 'base_id' (this will get all ORFs including the best one)
        group_df = df[df["base_id"] == base_id]
        # uORF defined by:
        #  - Coding_prob >= threshold
        #  - ORF_start < best ORF starting position
        #  - ORF_end < best ORF starting position
        uORF_entry = group_df[
            (group_df["Coding_prob"] >= coding_cutoff)
            & (group_df["ORF_start"] < row["ORF_start"])
            & (group_df["ORF_end"] < row["ORF_start"])
        ]
        if not uORF_entry.empty:
            # Mark that an upstream uORF exists
            best_orfs[base_id]["has_uORF"] = True
            # Get coding_prob of the uORF with the highest Coding_prob
            best_orfs[base_id]["coding_prob_uORF"] = uORF_entry["Coding_prob"].max()
        else:
            # No uORFs for this transcript
            best_orfs[base_id]["has_uORF"] = False
            best_orfs[base_id]["coding_prob_uORF"] = None    # or np.nan
        ## --------------------------------------
   

    # Optional: write debug output
    if isinstance(all_orfs_df, pd.DataFrame) and debug_output_path:
        all_orfs_df.to_csv(debug_output_path, sep='\t', index=False)
    else:
        logging.warning("Skipping debug output: all_orfs_df is not a DataFrame")

    return best_orfs
