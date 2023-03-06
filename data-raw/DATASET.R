## code to prepare `DATASET` dataset goes here
Metagenomics <- readRDS("data-raw/wgcna_metagenomics_eigengenes.rds")
Transcriptomics <- readRDS("data-raw/wgcna_transcriptomics_Liver_eigengenes.rds")
Metatranscriptomics <- readRDS("data-raw/Metatranscriptomics.rds")
usethis::use_data(Metagenomics,Transcriptomics, Metatranscriptomics, overwrite = TRUE)
