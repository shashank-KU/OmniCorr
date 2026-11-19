# ============================================================
# OmniCorr - Reproducible cattle example data generation
# ============================================================
#
# Run from the root of the OmniCorr repository.
#
# Source data are obtained using:
# inst/scripts/download_source_data.R
#
# ============================================================

library(WGCNA)

# 1. Source datasets

source_dir <- Sys.getenv(
  "OMNICORR_SOURCE_DIR",
  unset = "raw-data/zenodo_source/cattle"
)

if (!dir.exists(source_dir)) {
  stop(
    "Cattle source directory not found: ", source_dir,
    "\nRun inst/scripts/download_source_data.R first."
  )
}

# 2. WGCNA parameters for each layer

layer_settings <- list(

  Archbac = list(
    file = "archbac.RData",
    object = "archbac",
    power = 9,
    minModuleSize = 20
  ),

  Ciliates = list(
    file = "ciliates.RData",
    object = "ciliates",
    power = 9,
    minModuleSize = 20
  ),

  Fungi = list(
    file = "fungi.RData",
    object = "fungi",
    power = 14,
    minModuleSize = 5
  ),

  Metazoa = list(
    file = "metazoa.RData",
    object = "metazoa",
    power = 10,
    minModuleSize = 5
  ),

  Liver = list(
    file = "metL.RData",
    object = "metL",
    power = 14,
    minModuleSize = 5
  ),

  Wall = list(
    file = "metW.RData",
    object = "metW",
    power = 28,
    minModuleSize = 5
  )
)

# 3. Check source files

source_files <- vapply(
  layer_settings,
  function(x) file.path(source_dir, x$file),
  character(1)
)

missing_files <- source_files[!file.exists(source_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing cattle source files:\n",
    paste(missing_files, collapse = "\n")
  )
}

# 4. Reproduce WGCNA for each omics layer

modules <- lapply(names(layer_settings), function(layer) {

  settings <- layer_settings[[layer]]

  input_env <- new.env()

  load(
    file.path(source_dir, settings$file),
    envir = input_env
  )

  if (!exists(settings$object, envir = input_env)) {
    stop(
      "Object '", settings$object,
      "' not found in ", settings$file
    )
  }

  datExpr <- t(input_env[[settings$object]])

  set.seed(13118)

  result <- WGCNA::blockwiseModules(

    datExpr = datExpr,

    power = settings$power,

    maxBlockSize = 8000,

    networkType = "signed",

    corType = "bicor",

    maxPOutliers = 0.05,

    replaceMissingAdjacencies = TRUE,

    pamStage = FALSE,

    deepSplit = 4,

    minModuleSize = settings$minModuleSize,

    minKMEtoStay = 0.5,

    minCoreKME = 0.5,

    minCoreKMESize = 2,

    reassignThreshold = 0,

    mergeCutHeight = 0.2,

    verbose = 0

  )

  message(
    layer, ": ",
    ncol(result$MEs),
    " module eigengenes generated"
  )

  result$MEs

})

names(modules) <- names(layer_settings)

# 5. Identify shared samples

common_samples <- Reduce(
  intersect,
  lapply(modules, rownames)
)

stopifnot(length(common_samples) == 20)

# 6. Align samples

cattle_data <- lapply(
  modules,
  function(x) {
    as.data.frame(
      x[common_samples, , drop = FALSE]
    )
  }
)

# 7. Compare with original packaged dataset

original_file <- system.file(
  "extdata",
  "cattle_module_eigengenes.rds",
  package = "OmniCorr"
)

if (original_file == "") {
  original_file <- "inst/extdata/cattle_module_eigengenes.rds"
}

if (!file.exists(original_file)) {
  stop("Original cattle dataset not found.")
}

original_data <- readRDS(original_file)

comparison <- all.equal(
  original_data,
  cattle_data,
  tolerance = 1e-6
)

cat("\nComparison with original cattle dataset:\n")

print(comparison)

stopifnot(isTRUE(comparison))

message(
  "\nSUCCESS: Cattle dataset reproduced exactly."
)
