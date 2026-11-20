# =============================================================================
# Reproduce the salmon transcriptomics dataset used in OmniCorr
# =============================================================================
#
# Input:
#   raw-data/zenodo_source/salmon/transcriptomics/
#     stage2results_Transcriptomics_Liver.rds
#
# Output:
#   raw-data/transcriptomics/Transcriptomics_reproduced.rds
#
# The input contains WGCNA module eigengenes from 144 salmon samples.
# This script removes the grey module and averages eigengenes by tank.
#
# IMPORTANT:
# This reproduces the packaged dataset from an existing WGCNA intermediate.
# It does not reproduce the upstream liver WGCNA analysis from raw counts.
#
# Run from the OmniCorr package root directory.
# =============================================================================

library(dplyr)

# 1. Load the original liver WGCNA intermediate ------------------------------

input_file <- paste0(
  "raw-data/zenodo_source/salmon/transcriptomics/",
  "stage2results_Transcriptomics_Liver.rds"
)

stopifnot(file.exists(input_file))

original <- readRDS(input_file)

MEs <- as.data.frame(original$modules$MEs)

stopifnot(
  nrow(MEs) == 144,
  ncol(MEs) == 77,
  "ME0" %in% colnames(MEs)
)

# 2. Remove the grey module --------------------------------------------------

MEs <- MEs[
  , colnames(MEs) != "ME0",
  drop = FALSE
]

# 3. Extract tank identifiers ------------------------------------------------

tank_ids <- sub(
  "_[0-9]{1,2}$",
  "",
  rownames(MEs)
)

stopifnot(
  !anyNA(tank_ids),
  all(tank_ids != rownames(MEs))
)

# 4. Average module eigengenes by tank ---------------------------------------

aggregated <- MEs %>%
  mutate(Tank = tank_ids) %>%
  group_by(Tank) %>%
  summarise(
    across(everything(), mean),
    .groups = "drop"
  )

Transcriptomics_reproduced <- as.matrix(
  aggregated[, -1, drop = FALSE]
)

rownames(Transcriptomics_reproduced) <- aggregated$Tank

stopifnot(
  identical(dim(Transcriptomics_reproduced), c(36L, 76L))
)

message(
  "Reproduced dimensions: ",
  paste(dim(Transcriptomics_reproduced), collapse = " x ")
)

# 5. Validate against the packaged dataset ----------------------------------

reference_env <- new.env()

load(
  "data/Transcriptomics.rda",
  envir = reference_env
)

reference <- as.matrix(reference_env$Transcriptomics)

stopifnot(
  setequal(
    rownames(Transcriptomics_reproduced),
    rownames(reference)
  ),
  setequal(
    colnames(Transcriptomics_reproduced),
    colnames(reference)
  )
)

aligned <- Transcriptomics_reproduced[
  rownames(reference),
  colnames(reference),
  drop = FALSE
]

comparison <- all.equal(
  unname(aligned),
  unname(reference),
  tolerance = 1e-10
)

print(comparison)

if (!isTRUE(comparison)) {
  stop("Reproduced dataset differs from packaged Transcriptomics.")
}

message("SUCCESS: Packaged Transcriptomics dataset reproduced.")

# 6. Save reproduced dataset -------------------------------------------------

output_dir <- "raw-data/transcriptomics"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  Transcriptomics_reproduced,
  file.path(output_dir, "Transcriptomics_reproduced.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

message("Salmon transcriptomics dataset successfully generated.")
