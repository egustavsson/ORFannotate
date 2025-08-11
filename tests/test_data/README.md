# Minimal test dataset

This minimal dataset verifies that ORFannotate is installed and functioning correctly.
It runs on a small, realistic subset of the genome (mitochondrial contig only) and completes in just a few seconds using minimal memory.

## Included files 

- `toy.gtf` – a minimal GTF file with a few transcripts and exons.
- `toy.fa` – a matching reference FASTA containing only the sequences required by `toy.gtf`.

Running this test takes only a few seconds and completes the full ORFannotate workflow, including sequence reconstruction, ORF prediction, CDS annotation, and summary table generation.

## Test command

```
conda activate ORFannotate

python ORFannotate.py \
  --gtf tests/test_data/toy.gtf \
  --fa tests/test_data/toy.fa \
  --species human \
  --outdir output/
```

This will create:

- `output/transcripts.fa` – reconstructed transcript sequences.
- `output/CPAT/` – CPAT results.
- `output/ORFannotate_annotated.gtf` – GTF with CDS features inserted.
- `output/ORFannotate_summary.tsv` – transcript-level summary table.
- `output/ORFannotate.log` – run log.