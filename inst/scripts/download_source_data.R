# ============================================================
# OmniCorr - Source data download and extraction
# ============================================================
#
# Temporary local archive version.
# Replace the local archive path with the permanent Zenodo
# download URL after publication.
#
# Run from the OmniCorr repository root.
#
# ============================================================

# 1. Locate the source archive

archive_file <- Sys.getenv(
  "OMNICORR_SOURCE_ARCHIVE",
  unset = ""
)

if (!nzchar(archive_file)) {
  stop(
    "Please define the OMNICORR_SOURCE_ARCHIVE ",
    "environment variable with the path to the source ZIP archive."
  )
}

archive_file <- path.expand(archive_file)

if (!file.exists(archive_file)) {
  stop(
    "Source data ZIP archive not found: ",
    archive_file
  )
}

# 2. Define extraction directory

source_dir <- "raw-data/zenodo_source"

dir.create(
  source_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. Define expected source files

expected_files <- c(
  "cattle/archbac.RData",
  "cattle/ciliates.RData",
  "cattle/fungi.RData",
  "cattle/metazoa.RData",
  "cattle/metL.RData",
  "cattle/metW.RData",
  "cattle/metadata.RData",
  "salmon/metagenomics/feature-table_taxonomy.biom",
  "salmon/metagenomics/taxonomy.tsv",
  "salmon/metagenomics/metadata.txt",
  "salmon/metagenomics/tree.nwk",
  "salmon/metagenomics/dna-sequences.fasta",
  "salmon/metatranscriptomics/metatranscriptomics.rds",
  "salmon/transcriptomics/Transcriptomics_count.txt.gz",
  "salmon/transcriptomics/sample_info.csv",
  "salmon/transcriptomics/Salmo_salar-GCA_905237065.2_gene_annotations.tsv",
  "salmon/transcriptomics/stage2results_Transcriptomics_Liver.rds"
)

# 4. Check whether extraction is necessary

missing_files <- expected_files[
  !file.exists(file.path(source_dir, expected_files))
]

if (length(missing_files) > 0) {

  # Extract archive
  utils::unzip(
    zipfile = archive_file,
    exdir = source_dir
  )

}

# 5. Verify extracted files

missing_files <- expected_files[
  !file.exists(file.path(source_dir, expected_files))
]

if (length(missing_files) > 0) {

  stop(
    "The following source files are missing:\n",
    paste(missing_files, collapse = "\n")
  )

}

message(
  "All ", length(expected_files),
  " source files are available."
)

message(
  "Source data directory: ",
  normalizePath(source_dir)
)

# 6. Display extracted files

print(
  list.files(
    source_dir,
    recursive = TRUE
  )
)

# ============================================================
# Source data preparation completed
# ============================================================
