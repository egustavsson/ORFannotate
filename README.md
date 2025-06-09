# ORFannotate

`ORFannotate` is a modular pipeline for predicting open reading frames (ORFs), annotating coding sequences (CDS) in GTF files, and summarizing transcript coding status including potential nonsense-mediated decay (NMD).

---

## Features
- Extract transcript sequences from a GTF and genome FASTA
- Predict ORFs using [`orfipy`](https://github.com/urmi-21/orfipy)
- Annotate GTF files with CDS features
- Identify likely NMD targets (simple heuristic)
- Output detailed TSV with:
  - transcript ID, gene, chrom, strand
  - coding status
  - ORF length, CDS length
  - predicted NMD
  - transcript, CDS, and amino acid sequences

---

## Installation

1. Clone the repository:
```bash
git clone https://github.com/egustavsson/ORFannotate.git
cd ORFannotate
```

2. Create the environment:
```bash
conda env create -f ORFannotate.conda_env.yml
conda activate ORFannotate
```

> Note: All required tools including `orfipy` and `gffread` are installed automatically from the `bioconda` channel when creating the environment.

---

## Usage

```bash
python ORFannotate.py <input.gtf> <genome.fasta> <output_dir>
```

**Example:**
```bash
python ORFannotate.py annotations.gtf genome.fa output/
```

This will produce:
- `output/transcripts.fa`: transcript sequences
- `output/predicted_orfs.fasta`: ORF predictions from orfipy
- `output/annotated.gtf`: GTF annotated with CDS features
- `output/transcript_summary.tsv`: full summary table

---

## 📁 Directory Structure
```bash
ORFannotate/
├── ORFannotate.py          # Main script
├── environment.yml         # Conda environment
├── orfannotate/            # Modular Python package
│   ├── __init__.py
│   ├── nmd.py
│   ├── orf_parser.py
│   ├── gtf_annotation.py
│   ├── summarise.py
│   └── io_utils.py
├── README.md
├── LICENSE
└── tests/                  # Optional test cases
```

---

## NMD Prediction
NMD is predicted based on the distance between the ORF stop codon and the last exon junction. Transcripts with a stop codon >50 nt upstream of the final junction are flagged as likely NMD targets.

---

## License
This project is licensed under the MIT License. See `LICENSE` for details.

---

## Acknowledgements
- [orfipy](https://github.com/urmi-21/orfipy)
- [gffutils](https://github.com/daler/gffutils)
- [Biopython](https://biopython.org/)
