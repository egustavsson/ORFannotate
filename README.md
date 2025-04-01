# ORFannotate

<!-- badges: start -->
![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

The `ORFannotate` tool is designed to predict open reading frames (ORFs) and assess nonsense-mediated decay (NMD) from a given GTF/GFF file containing transcript annotations. It processes transcript structures to identify coding sequences (CDS), determine potential translation start and stop sites, and infer whether transcripts are likely to be subject to NMD.

The tool outputs ORF and NMD predictions and updates the input GTF/GFF file with annotated coding exons (CDS) for downstream analyses. ORFannotate provides an efficient way to refine transcript annotations with coding information, enhancing downstream functional analyses of RNA sequencing data.