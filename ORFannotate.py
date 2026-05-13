import os
import re
import logging
import gffutils
import argparse
import json
from datetime import datetime
from pathlib import Path
from Bio import SeqIO

# Progress bar
from rich.text import Text
from rich.progress import (
    BarColumn,
    Progress
)


# Internal modules
from orfannotate.orf_prediction import run_cpat, predict_uorf
from orfannotate.transcript_extraction import extract_transcripts_from_gtf
from orfannotate.orf_filter import get_best_orfs_by_cpat
from orfannotate.gtf_annotation import build_cds_features, annotate_gtf_with_cds
from orfannotate.summarise import generate_summary
from orfannotate.utils import _cleanup_tmp_folder, _count_unique_transcripts, _make_slim_gtf, _print_top_database_info, _validate_inputs, _create_db, _setup_logging,_ColoredTaskProgressColumn,_ColoredTimeElapsedColumn #_collect_junctions,_collect_junctions_batch, 

# version
from orfannotate import __version__


# Path relative to ORFannotate for Snakemake, Nextflow workflows, etc.
MODULE_DIR = Path(__file__).resolve().parent
DATA_DIR = MODULE_DIR / "data"


    
## Define global progress bar parameters
TOTAL_UNITS = 1000          # Global bar total
STEP_UNITS = TOTAL_UNITS//12 # = 200 units per step
progress_columns = (
    "[progress.description]{task.description}",
    BarColumn(complete_style="white", finished_style="white", pulse_style="white"),
    _ColoredTaskProgressColumn("white"),   # recolored %
    "Elapsed:",
    _ColoredTimeElapsedColumn("white"),     # recolored elapsed
)


    
# Main function
def main():
    with Progress(*progress_columns) as progress:
        
        parser = argparse.ArgumentParser(
            description="ORFannotate - predict coding ORFs, annotate GTF, and generate summaries."
        )
        parser.add_argument("--gtf", required=True, help="Input GTF or GFF file with transcript and exon features")
        parser.add_argument("--fa", required=True, help="Reference genome in FASTA format")
        parser.add_argument("--outdir", required=True, help="Directory to write all outputs")
        parser.add_argument("--species", choices=["human", "mouse", "fly", "zebrafish"], default="human",
                            help="Species to use for CPAT model (default: human)")
        parser.add_argument("--coding-cutoff", type=float, default=None,
                            help="CPAT cutoff for classifying coding vs noncoding transcripts")
        parser.add_argument("--version", action="version", version=f"ORFannotate {__version__}")
        args = parser.parse_args()

        gtf_path = args.gtf
        genome_fasta = args.fa
        output_dir = args.outdir
        species = args.species.lower()

        _validate_inputs(gtf_path, genome_fasta)
        os.makedirs(output_dir, exist_ok=True)

        logger = _setup_logging(output_dir)
        logger.info("Starting ORFannotate")
        

        num_transcripts = _count_unique_transcripts(gtf_path)
        logger.info(f"Found {num_transcripts:,} unique transcripts in the GTF")


        task = progress.add_task("[white]ORFannotate progress...", total=TOTAL_UNITS)
        

        progress.update(task, advance=STEP_UNITS)
        logger.info("Step 1: Building GTF database in memory...")

        # Create slim GTF with only transcript and exon features
        slim_gtf_path = _make_slim_gtf(gtf_path)

        # Create database from GTF ussing GFFUtils
        gtf_db = _create_db(slim_gtf_path, only_exons=True)

        progress.update(task, advance=STEP_UNITS)


        logger.info("Step 2: Extracting transcript sequences...")
        transcript_fasta = os.path.join(output_dir, "transcripts.fa")
        extract_transcripts_from_gtf(gtf_db, genome_fasta, transcript_fasta)
        progress.update(task, advance=STEP_UNITS)

        # Validate transcript ID consistency
        fasta_tx_ids = {rec.id for rec in SeqIO.parse(transcript_fasta, "fasta")}
        progress.update(task, advance=STEP_UNITS)
        
        # Collect all transcript IDs
        gtf_tx_ids = {
            row[0]
            for row in gtf_db.execute("SELECT id FROM features WHERE featuretype='transcript'")
        }
        missing = gtf_tx_ids - fasta_tx_ids
        if missing:
            logger.warning(f"{len(missing)} transcripts from GTF missing in extracted FASTA: e.g. {sorted(list(missing))[:5]}")

        progress.update(task, advance=STEP_UNITS)
        logger.info("Step 3: Predicting and scoring ORFs...")

        # Load species config from JSON
        config_path = DATA_DIR / "config.json"
        with open(config_path) as f:
            config = json.load(f)

        if species not in config["species_models"]:
            logger.error(f"Species '{species}' not found in config.json")
            raise ValueError(f"Species '{species}' not supported. Check config.json.")

        # Extract species-specific config parameters
        species_data = config["species_models"][species]
        top_orf = species_data["top_orf"]
        min_length_orf = species_data["min_length_orf"]
        hexamer_path = str(DATA_DIR / species_data["hexamer"])
        logit_model_path = str(DATA_DIR / species_data["model"])
        coding_cutoff = args.coding_cutoff if args.coding_cutoff is not None else species_data["cutoff"]
        
        logger.info(f"Using CPAT model directory: {DATA_DIR}")
        logger.info(f"Selected species: {species} | Cutoff: {coding_cutoff}")
        logger.info(f"Using coding cutoff: {coding_cutoff}")

        cpat_dir = os.path.join(output_dir, "CPAT")
        os.makedirs(cpat_dir, exist_ok=True)

        # Run CPAT for canonical ORF
        run_cpat(transcript_fasta, cpat_dir, hexamer_path, logit_model_path, top_orf, min_length_orf)

        progress.update(task, advance=STEP_UNITS)
        logger.info("Step 4: Parsing ORF results...")
        
        cpat_results = os.path.join(cpat_dir, "cpat.ORF_prob.tsv")
        debug_output = os.path.join(cpat_dir, "cpat_debug.tsv")
        best_orfs = get_best_orfs_by_cpat(cpat_results, debug_output_path=debug_output)
        
        coding_orfs = {tid: info for tid, info in best_orfs.items() if info["coding_prob"] >= coding_cutoff}
        logger.info(f"Selected {len(coding_orfs):,} transcripts classified as coding (cutoff = {coding_cutoff})")
  
        if not coding_orfs:
            logger.warning("No coding transcripts predicted above the cutoff. Downstream outputs may be empty or incomplete.")

        progress.update(task, advance=STEP_UNITS)
        logger.info("Step 5: Annotating GTF with CDS features...")


        logger.info(f"Adding CDS to {len(coding_orfs):,} coding transcripts")
        cds_features = build_cds_features(gtf_db, coding_orfs)
        
    
        # Cleanup first db
        try:
            gtf_db.conn.close()
            del gtf_db
            logger.info("Closed and removed first GTF database to free memory")
        except Exception:
            logger.warning("Could not close first GTF database early")
            
    
        annotated_gtf = os.path.join(output_dir, "ORFannotate_annotated.gtf")
        annotate_gtf_with_cds(gtf_path, cds_features, annotated_gtf)
        logger.info(f"Annotated GTF written to {annotated_gtf}")


        progress.update(task, advance=STEP_UNITS)
        logger.info("Step 6: Annotating transcripts and generating final summary TSV...")


        # Slim the GTF for summary annotation
        slim_annotated_path = _make_slim_gtf(annotated_gtf, only_exons=False)
        slim_annotated_db = _create_db(slim_annotated_path, only_exons=False)


        summary_tsv = os.path.join(output_dir, "ORFannotate_summary.tsv")
        utr5_records = generate_summary(best_orfs = best_orfs, transcript_fa = transcript_fasta, gtf_db_or_path = slim_annotated_db, output_path = summary_tsv, progress=progress, task=task, STEP_UNITS=STEP_UNITS, coding_cutoff=coding_cutoff, output_dir=output_dir)
        
        _cleanup_tmp_folder()

        logger.info("Step 7: Predicting and scoring uORFs...")
        min_length_uorf = species_data["min_length_uorf"]

        if not utr5_records:
            logger.warning("No 5'UTR sequences found. Skipping uORF prediction step.")
        elif len(utr5_records) < 2:
            logger.warning(f"Only {len(utr5_records)} 5'UTR sequences found. Too few for reliable uORF prediction. Skipping.")
        else:
            # Check that at least one 5'UTR sequence is above the minimum length threshold
            all_short = all(len(rec.seq) < ((min_length_uorf*3)+10) for rec in utr5_records)
            if all_short:
                logger.warning("All 5'UTR sequences are shorter than the minimum '{min_length_uorf}' amino acids uORF length configured. Eg, 'min_length_uorf=100', requires at least 300 nt plus start and stop codon constraints. Skipping uORF prediction step.")
            else:
                predict_uorf(output_dir, hexamer_path, logit_model_path, coding_cutoff, top_orf, min_length_uorf, slim_annotated_db, gtf_path)
        
        progress.update(task, advance=STEP_UNITS)

        try:
            gtf_db.conn.close()
        except Exception:
            pass

        logger.info(f"Summary written to {summary_tsv}")
        
        params = {
            "timestamp": datetime.now().isoformat(),
            "gtf": os.path.abspath(gtf_path),
            "fasta": os.path.abspath(genome_fasta),
            "outdir": os.path.abspath(output_dir),
            "species": species,
            "coding_cutoff": coding_cutoff
        }
        with open(os.path.join(output_dir, "parameters.json"), "w") as f:
            json.dump(params, f, indent=2)
        
        logger.info("Parameters saved to parameters.json")
        logger.info("ORFannotate completed successfully.")

        progress.stop()

if __name__ == "__main__":
    main()
