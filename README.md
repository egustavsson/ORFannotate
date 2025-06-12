# ORFannotate

`ORFannotate` is a modular pipeline for predicting open reading frames (ORFs), annotating coding sequences (CDS) in GTF/GFF files, and summarizing transcript coding potential, including a simple heuristic for nonsense-mediated decay (NMD).For a comprehensive isoform-level quality assessment, consider using `SQANTI3` in parallel.

---

## Features
- Extract transcript sequences from a GTF/GFF and genome FASTA
- Predict and score ORFs using [`CPAT`](https://cpat.readthedocs.io/en/latest/#introduction)
- Annotate GTF/GFF files with CDS features (only for coding transcripts)
- Identify likely NMD targets (simple heuristic)
- Generate a rich tab-separated summary of ORF and coding properties, including:
  - Transcript and gene ID, strand, chromosome
  - ORF start/end, frame, coding probability
  - ORF/CDS length (nt/aa), junction count
  - Predicted NMD flag
  - Nucleotide and protein sequence
  - `coding_class` classification (coding / noncoding)

---

## Installation

1. Clone the repository:
```
git clone https://github.com/egustavsson/ORFannotate.git
cd ORFannotate
```

2. Create the environment:
```
conda env create -f ORFannotate.conda_env.yml
conda activate ORFannotate
```

> Note: All required tools are installed automatically when creating the environment.

---
## Input

The pipeline requires:

| File           | Format       | Notes                                                                 |
|----------------|--------------|-----------------------------------------------------------------------|
| GTF/GFF file   | `.gtf` or `.gff` | Must contain `transcript` and `exon` features. Transcript IDs should be unique and consistent. |
| Genome FASTA   | `.fa`, `.fasta` | Reference genome matching the GTF. Must have transcript chromosome names matching those in the annotation. |

> Both `.gtf` and `.gff` formats are supported as long as feature labels are compatible (`transcript`, `exon`, etc.).

Ensure that your GTF includes transcript-level features (not just exons), or the pipeline may skip some records or produce incomplete outputs.

---

## Usage

To see available arguments:
```bash
python ORFannotate.py --help
```
This will print:

```
usage: ORFannotate.py [-h] --gtf GTF --fa FA --outdir OUTDIR [--coding-cutoff CODING_CUTOFF]

ORFannotate – predict coding ORFs, annotate GTF, and generate summaries.

required arguments:
  --gtf             Input GTF or GFF file with transcript and exon features
  --fa              Reference genome in FASTA format
  --outdir          Directory to write all outputs

optional arguments:
  --coding-cutoff   CPAT probability threshold to classify coding vs noncoding (default: 0.364)
```
To run the full pipeline:
```bash
python ORFannotate.py --gtf annotations.gtf --fa genome.fa --outdir output/
```

**Example:**
```bash
python ORFannotate.py transcripts.gtf genome.fa output/
```

After a successful run, the following files will be saved in <output_dir>:

| **File**                      | **Description**                                 |
| ----------------------------- | ----------------------------------------------- |
| `transcripts.fa`              | FASTA file of transcript sequences              |
| `cpat.ORF_prob.best.tsv`      | CPAT output for the best ORF per transcript     |
| `cpat_debug.tsv`              | (Optional) full CPAT-scored ORFs                |
| `ORFannotate_annotated.gtf`   | GTF with CDS features added to transcripts      |
| `ORFannotate_summary.tsv`     | Final summary table with ORF/NMD annotations    |
| `CPAT_run_info.log`           | CPAT runtime log (captured in output directory) |


> Only transcripts classified as coding (by CPAT probability ≥ 0.364 by default) are annotated with CDS in the output GTF. You can override the cutoff using the --coding-cutoff argument.

## Directory Structure
```
ORFannotate/
├── ORFannotate.py                # Main script (entry point)
├── ORFannotate.conda_env.yml     # Conda environment file
├── orfannotate/                  # Modular Python package
│   ├── __init__.py
|   ├── transcript_extraction.py  # Handles extraction of transcript sequences
│   ├── gtf_annotation.py         # GTF handling and CDS annotation
│   ├── nmd.py                    # NMD prediction logic
│   ├── orf_filter.py             # CPAT result parsing and filtering
│   ├── summarise.py              # Final summary writer
├── LICENSE
└── README.md

```

## NMD Prediction
Nonsense-mediated decay (NMD) is predicted using a simple rule:
If the stop codon lies >50 nt upstream of the final exon–exon junction, the transcript is flagged as a likely NMD target.

This conservative approach is fast and works well for general transcriptome-level analyses, but may not capture all context-dependent cases.

## License
This project is licensed under the GPLv3 License. See `LICENSE` for details.

## Acknowledgements
- [CPAT](https://github.com/urmi-21/orfipy)
- [gffutils](https://github.com/daler/gffutils)
- [Biopython](https://biopython.org/)
