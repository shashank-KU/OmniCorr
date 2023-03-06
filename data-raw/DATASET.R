## code to prepare `DATASET` dataset goes here
Metagenomics <- readRDS("data-raw/wgcna_metagenomics_eigengenes.rds")
Transcriptomics <- readRDS("data-raw/wgcna_transcriptomics_Liver_eigengenes.rds")
usethis::use_data(Transcriptomics, Metagenomics, overwrite = TRUE)
usethis::use_data(Metagenomics, overwrite = TRUE)
