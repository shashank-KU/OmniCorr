# =============================================================================
# Reconstruct salmon transcriptomics WGCNA modules from count data
# =============================================================================
#
# Input:
#   raw-data/zenodo_source/salmon/transcriptomics/
#     Transcriptomics_count.txt.gz
#     sample_info.csv
#
# Purpose:
#   Reconstruct the original 142-sample transcriptomics WGCNA analysis.
#
# IMPORTANT:
#   The reconstructed gene-to-module assignments match the original
#   142-sample analysis, but the module eigengene values differ.
#
#   This script does NOT reproduce the packaged OmniCorr Transcriptomics
#   dataset (36 tanks x 76 modules).
#
#   The packaged dataset is reproduced separately by:
#     inst/scripts/generate_salmon_transcriptomics.R
#
# Run from the OmniCorr package root directory.
# =============================================================================

library(DESeq2)
library(WGCNA)
library(dplyr)
library(tibble)

set.seed(13118)

# 1. Load source data ---------------------------------------------------------

source_dir <- "raw-data/zenodo_source/salmon/transcriptomics"

counts_raw <- read.table(
  file.path(source_dir, "Transcriptomics_count.txt.gz"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

sample_info <- read.csv(
  file.path(source_dir, "sample_info.csv"),
  row.names = 1
)

# 2. Reconstruct and match sample identifiers --------------------------------

sample_info$ID_New <- paste0(
  sample_info$Time, "_",
  sample_info$New_Diet, "_",
  sample_info$Tank_number, "_0",
  sample_info$Sample_Number
)

rownames(sample_info) <- sample_info$ID_New

stopifnot(
  !anyDuplicated(sample_info$ID_New),
  all(sample_info$ID_New %in% colnames(counts_raw))
)

counts_mat <- as.matrix(
  counts_raw[, sample_info$ID_New, drop = FALSE]
)

stopifnot(
  identical(colnames(counts_mat), rownames(sample_info)),
  ncol(counts_mat) == 142
)

message("Matched samples: ", ncol(counts_mat))

# 3. Apply original gene filtering -------------------------------------------

keep <- rowSums(counts_mat) >= 50

counts_filtered <- round(
  counts_mat[keep, , drop = FALSE]
)

stopifnot(
  nrow(counts_filtered) == 36878
)

message("Genes retained: ", nrow(counts_filtered))

# 4. DESeq2 preprocessing -----------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_info,
  design = ~ 1
)

dds$group <- factor(
  paste0(dds$Time, dds$New_Diet)
)

design(dds) <- ~ group - 1

dds <- DESeq(
  dds,
  parallel = TRUE
)

omics_data_host <- t(assay(vst(dds)))

stopifnot(
  nrow(omics_data_host) == 142,
  ncol(omics_data_host) == 36878
)

# 5. Check data quality -------------------------------------------------------

gsg <- goodSamplesGenes(
  omics_data_host,
  verbose = 3
)

stopifnot(gsg$allOK)

# 6. Select WGCNA soft-thresholding power ------------------------------------

powers <- c(1:10, seq(12, 20, by = 2))

sft <- pickSoftThreshold(
  omics_data_host,
  powerVector = powers,
  verbose = 1,
  networkType = "signed",
  corFn = "bicor"
)

fit <- -sign(sft$fitIndices[, 3]) *
  sft$fitIndices[, 2]

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

# 7. Construct WGCNA modules -------------------------------------------------

set.seed(13118)

modules <- blockwiseModules(
  omics_data_host,
  power = st,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  minModuleSize = 20,
  nThreads = 10,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  saveTOMFileBase = "HOST_TOM",
  verbose = 3,
  maxBlockSize = 8000
)

# Save reconstructed WGCNA object for reproducibility checks
dir.create(
  "raw-data/transcriptomics",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  modules,
  "raw-data/transcriptomics/WGCNA_modules_reconstructed.rds"
)

# 8. Extract and aggregate module eigengenes ---------------------------------

MEs <- as.data.frame(modules$MEs)

MEs <- MEs[
  , colnames(MEs) != "ME0",
  drop = FALSE
]

tank_ids <- sub(
  "_[0-9]{1,2}$",
  "",
  rownames(MEs)
)

aggregated <- MEs %>%
  mutate(common = tank_ids) %>%
  group_by(common) %>%
  summarise(
    across(everything(), mean),
    .groups = "drop"
  )

Transcriptomics_reproduced <- as.matrix(
  aggregated[, -1, drop = FALSE]
)

rownames(Transcriptomics_reproduced) <- aggregated$common

message(
  "Final dimensions: ",
  paste(dim(Transcriptomics_reproduced), collapse = " x ")
)


# 9. Compare with original 142-sample transcriptomics analysis ----------------

original <- readRDS(
  "raw-data/transcriptomics/stage2results_Transcriptomics.rds"
)

reference_MEs <- as.data.frame(original$modules$MEs)

reference_MEs <- reference_MEs[
  , colnames(reference_MEs) != "ME0",
  drop = FALSE
]

reference_tanks <- sub(
  "_[0-9]{1,2}$",
  "",
  rownames(reference_MEs)
)

reference_aggregated <- reference_MEs %>%
  mutate(common = reference_tanks) %>%
  group_by(common) %>%
  summarise(
    across(everything(), mean),
    .groups = "drop"
  )

reference <- as.matrix(
  reference_aggregated[, -1, drop = FALSE]
)

rownames(reference) <- reference_aggregated$common

cat("Reconstructed dimensions:\n")
print(dim(Transcriptomics_reproduced))

cat("\nOriginal transcriptomics dimensions:\n")
print(dim(reference))

cat("\nMatching tank identifiers:\n")
print(setequal(
  rownames(Transcriptomics_reproduced),
  rownames(reference)
))

cat("\nMatching module names:\n")
print(setequal(
  colnames(Transcriptomics_reproduced),
  colnames(reference)
))

if (
  setequal(
    rownames(Transcriptomics_reproduced),
    rownames(reference)
  ) &&
  setequal(
    colnames(Transcriptomics_reproduced),
    colnames(reference)
  )
) {
  aligned <- Transcriptomics_reproduced[
    rownames(reference),
    colnames(reference),
    drop = FALSE
  ]
  
  comparison <- all.equal(
    unname(aligned),
    unname(reference),
    tolerance = 1e-6
  )
  
  cat("\nNumerical comparison:\n")
  print(comparison)
  
  if (!isTRUE(comparison)) {
    warning(
      "Reconstructed eigengenes differ from the original ",
      "142-sample transcriptomics analysis. ",
      "Exact numerical reproducibility has not been established."
    )
  } else {
    message("Original 142-sample eigengenes reproduced numerically.")
  }
} else {
  stop("Reconstructed module or tank identifiers do not match.")
}

# 10. Save reproduced dataset -------------------------------------------------

output_dir <- "raw-data/transcriptomics"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  Transcriptomics_reproduced,
  file.path(output_dir, "Transcriptomics_WGCNA_reconstructed.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

message(
  "Exploratory salmon transcriptomics WGCNA reconstruction completed. ",
  "See the numerical comparison and any warnings above."
)
