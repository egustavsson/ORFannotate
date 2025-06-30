import os
import subprocess
import logging

logger = logging.getLogger(__name__)

def run_cpat(transcript_fasta, output_dir, hexamer_path, logit_model_path, force=False):
    output_prefix = os.path.join(output_dir, "cpat")
    output_file = output_prefix + ".ORF_prob.best.tsv"

    if os.path.exists(output_file) and not force:
        logger.info("CPAT output exists, skipping CPAT run")
        return output_file

    cpat_cmd = [
        "cpat.py",
        "-g", os.path.abspath(transcript_fasta),
        "-x", os.path.abspath(hexamer_path),
        "-d", os.path.abspath(logit_model_path),
        "--top-orf=10",
        "--min-orf=75",
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
