# =============================================================================
# Reproduce the salmon metatranscriptomics dataset used in OmniCorr
# =============================================================================
#
# Input:
#   raw-data/zenodo_source/salmon/metatranscriptomics/metatranscriptomics.rds
#
# Output:
#   Metatranscriptomics: 36 samples x 11 module eigengenes
#
# Run this script from the OmniCorr package root directory.
#
# =============================================================================
# Software environment used for successful reproduction:
# R version: 4.6.2 (2026-10-29)
# Platform: aarch64-apple-darwin23 (macOS)
#
# Main package versions:
# DESeq2: 1.52.0
# WGCNA: 1.74
# dplyr: 1.2.1
# tibble: 3.3.1
#
# The reproduced dataset was compared against
# data/Metatranscriptomics.rda using all.equal().
# Result: TRUE (numerical values match within tolerance).
#
# Note: WGCNA results may depend on software versions
# and random initialization.

library(DESeq2)
library(WGCNA)
library(dplyr)
library(tibble)

# 1. Load raw data ------------------------------------------------------------

dds <- readRDS(
  "raw-data/zenodo_source/salmon/metatranscriptomics/metatranscriptomics.rds"
)

# 2. Preprocess the data ------------------------------------------------------

# Remove T0 samples
dds <- dds[, dds$TimePoint != "T0"]

# Retain genes with at least 5 counts in at least 5 samples
keep <- rowSums(counts(dds) >= 5) >= 5
dds <- dds[keep, ]

# Variance-stabilizing transformation
metatranscriptomics <- t(assay(vst(dds)))

stopifnot(
  nrow(metatranscriptomics) == 139,
  ncol(metatranscriptomics) == 73128
)

# 3. Check data quality -------------------------------------------------------

gsg <- goodSamplesGenes(
  metatranscriptomics,
  verbose = 3
)

stopifnot(gsg$allOK)

# 4. Select WGCNA soft-thresholding power ------------------------------------

powers <- c(1:10, seq(12, 20, by = 2))

sft <- pickSoftThreshold(
  metatranscriptomics,
  powerVector = powers,
  verbose = 1,
  networkType = "signed",
  corFn = "bicor"
)

fit <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]

idx <- which(fit > 0.90)

if (length(idx) > 0) {
  st <- sft$fitIndices[min(idx), 1]
} else {
  idx <- which(fit > 0.80)

  if (length(idx) > 0) {
    st <- sft$fitIndices[min(idx), 1]
  } else {
    st <- sft$fitIndices[which.max(fit), 1]
  }
}

message("Selected WGCNA power: ", st)

stopifnot(st == 9)

# 5. Construct WGCNA modules -------------------------------------------------

modules.omics.MetaTrans <- blockwiseModules(
  metatranscriptomics,
  power = st,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  minModuleSize = 21,
  nThreads = 10,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  saveTOMFileBase = "HOST_TOM",
  verbose = 3,
  maxBlockSize = 8000
)

# 6. Extract module eigengenes ------------------------------------------------

MEs <- as.data.frame(modules.omics.MetaTrans$MEs)

# 7. Remove unassigned module and average replicates -------------------------

Metatranscriptomics <- MEs %>%
  rownames_to_column("sampleName") %>%
  select(-ME0) %>%
  mutate(
    sampleName = sub("_HGm[0-9]{1,2}$", "", sampleName),
    sampleName = paste0(
      sub("_.*$", "", sampleName),
      "_",
      toupper(sub("^.*?_([^_]+)_.*$", "\\1", sampleName)),
      "_",
      sub("^.*_[^_]+_", "", sampleName)
    )
  ) %>%
  group_by(sampleName) %>%
  summarise(
    across(everything(), mean),
    .groups = "drop"
  ) %>%
  column_to_rownames("sampleName")

# 8. Verify the final dataset -------------------------------------------------

stopifnot(
  nrow(Metatranscriptomics) == 36,
  ncol(Metatranscriptomics) == 11
)

message("Final dimensions: ",
        paste(dim(Metatranscriptomics), collapse = " x "))

# 9. Compare with the packaged OmniCorr dataset -----------------------------

# Preserve the newly generated dataset
Metatranscriptomics_reproduced <- Metatranscriptomics

# Load the packaged dataset into a separate environment
reference_env <- new.env(parent = baseenv())

load(
  "data/Metatranscriptomics.rda",
  envir = reference_env
)

Metatranscriptomics_reference <- reference_env$Metatranscriptomics

# Verify that sample and module names match
stopifnot(
  setequal(
    rownames(Metatranscriptomics_reproduced),
    rownames(Metatranscriptomics_reference)
  ),
  setequal(
    colnames(Metatranscriptomics_reproduced),
    colnames(Metatranscriptomics_reference)
  )
)

# Align rows and columns
Metatranscriptomics_reproduced <-
  Metatranscriptomics_reproduced[
    rownames(Metatranscriptomics_reference),
    colnames(Metatranscriptomics_reference),
    drop = FALSE
  ]

# Verify numerical reproducibility
stopifnot(
  isTRUE(all.equal(
    as.matrix(Metatranscriptomics_reproduced),
    as.matrix(Metatranscriptomics_reference)
  ))
)

message("Verification successful: reproduced data match packaged data.")

# 10. Save the reproduced dataset ---------------------------------------------

dir.create(
  "raw-data/metatranscriptomics",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save a separate copy without overwriting the packaged dataset.
saveRDS(
  Metatranscriptomics,
  file = "raw-data/metatranscriptomics/Metatranscriptomics_reproduced.rds"
)

message("Salmon metatranscriptomics dataset successfully generated.")

# 11. Record software environment --------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = "raw-data/metatranscriptomics/sessionInfo.txt"
)
