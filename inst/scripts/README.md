# Reproducibility of OmniCorr Example Datasets

This directory contains R scripts used to generate and validate the example datasets distributed with the OmniCorr package.

## Source data

The source datasets are maintained separately from the R package to avoid including large intermediate files in the package distribution.

A source-data ZIP archive has been prepared for deposition on Zenodo. The archive is not yet publicly available.

Until a permanent Zenodo URL is available, the archive can be supplied locally using the `OMNICORR_SOURCE_ARCHIVE` environment variable.

Run all scripts from the root directory of the OmniCorr repository.

## 1. Prepare the source data

Set the environment variable to the location of the source ZIP archive:

```bash
export OMNICORR_SOURCE_ARCHIVE="/path/to/OmniCorr_source_data.zip"
Rscript inst/scripts/download_source_data.R
```

The script extracts the archive into `raw-data/zenodo_source/` and verifies the presence of 17 expected source files.

## 2. Generate the example datasets

The following scripts reproduce the example datasets or their corresponding processed representations.

| Dataset | Script | Input | Reproducibility |
|---|---|---|---|
| Cattle | `generate_cattle_data.R` | Cattle source `.RData` files | Reproduced from supplied source data |
| Salmon metagenomics | `generate_salmon_metagenomics.R` | BIOM table, taxonomy, metadata, sequences and phylogenetic tree | Reproduced from supplied source data |
| Salmon metatranscriptomics | `generate_salmon_metatranscriptomics.R` | Saved DESeq2 object | Reproduced from supplied intermediate |
| Salmon transcriptomics | `generate_salmon_transcriptomics.R` | Saved liver WGCNA results | Reproduced from supplied intermediate |
| Salmon metadata | `generate_salmon_metadata.R` | `sample_info.csv` | Reproduced from supplied source metadata |

Execute each script from the repository root:

```bash
Rscript inst/scripts/generate_cattle_data.R
Rscript inst/scripts/generate_salmon_metagenomics.R
Rscript inst/scripts/generate_salmon_metatranscriptomics.R
Rscript inst/scripts/generate_salmon_transcriptomics.R
Rscript inst/scripts/generate_salmon_metadata.R
```

The scripts write reconstructed outputs to locations under `raw-data/`, leaving the packaged datasets unchanged.

## 3. R package dependencies

The workflows use the following R packages:

- `WGCNA`
- `DESeq2`
- `phyloseq`
- `Biostrings`
- `stringr`
- `decontam`
- `vegan`
- `dplyr`
- `tibble`

These dependencies are required for the corresponding generation scripts, not necessarily for every OmniCorr workflow.

## 4. Validation

The generation scripts have been checked against the corresponding datasets distributed with OmniCorr.

The validated comparisons include:

- Cattle module eigengenes: numerical agreement within tolerance `1e-6`.
- Salmon metagenomics: numerical agreement within tolerance `1e-10`.
- Salmon metatranscriptomics: numerical agreement with the packaged dataset.
- Salmon transcriptomics: numerical agreement within tolerance `1e-10`.
- Salmon metadata: numerical agreement within tolerance `1e-10`.

The comparisons were performed using the reconstructed outputs and the existing packaged datasets.

## 5. Salmon transcriptomics reproducibility limitation

The packaged salmon transcriptomics dataset contains 36 tanks and 76 module eigengenes.

`generate_salmon_transcriptomics.R` reproduces this dataset by averaging module eigengenes from the supplied `stage2results_Transcriptomics_Liver.rds` intermediate, which contains 144 samples and 77 eigengenes, including the grey module (`ME0`).

The resulting 36 × 76 dataset matches the packaged `Transcriptomics.rda` numerically.

A separate exploratory script, `reconstruct_salmon_transcriptomics_wgcna.R`, reconstructs a 142-sample transcriptomics analysis from count data using DESeq2 and WGCNA.

Although the reconstructed gene-to-module assignments match the corresponding original 142-sample analysis, the eigengene values differ. This exploratory reconstruction does not reproduce the packaged 36 × 76 transcriptomics dataset.

Therefore, the upstream generation of the liver WGCNA intermediate has not yet been independently reproduced from the original count data.

## 6. Data provenance and availability

The source archive contains both source datasets and previously generated analysis intermediates. It should not be interpreted as containing raw sequencing data for all omics layers.

The Zenodo DOI and permanent download URL will be added after publication.

Redistribution permissions and licensing must be confirmed with the relevant data owners before the archive is made public.
