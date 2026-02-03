# raw-data

This directory contains the full preprocessing workflows used to generate
the example datasets distributed with the OmniCorr package.

These materials are provided to ensure transparency and reproducibility,
as requested by reviewers, and are not executed during package installation.

## Contents

### transcriptomics/
- Raw count data (`Transcriptomics_count.txt.gz`)
- Sample metadata
- Gene annotations
- `Transcriptomics_preprocessing.Rmd`: normalization, filtering, and preparation
  of transcriptomics data prior to WGCNA and OmniCorr analysis.

### metagenomics/
- Raw sequencing and feature tables (FASTA, BIOM, taxonomy)
- `Metageomics_preprocessing.Rmd`: preprocessing and feature summarization
  used to generate metagenomics inputs for OmniCorr.

### metatranscriptomics/
- Intermediate R objects (`.rds`)
- `Metatranscriptomics_preprocessing.Rmd`: preprocessing and formatting steps
  used prior to OmniCorr analysis.

The final processed objects derived from these workflows are saved as `.rda`
files in the package `data/` directory and are used in the README and manuscript
examples.
