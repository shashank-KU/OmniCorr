# ============================================================
# OmniCorr - Reproducible salmon metagenomics data generation
# ============================================================
#
# Reconstruct the packaged Metagenomics dataset directly
# from the original Zenodo source data.
#
# Run from the root of the OmniCorr repository.
#
# Source data:
# raw-data/zenodo_source/salmon/metagenomics/
# raw-data/zenodo_source/salmon/transcriptomics/
#
# Original workflow:
# Metageomics_preprocessing.Rmd
#
# ============================================================


# 1. Load required packages
# ============================================================

library(phyloseq)
library(Biostrings)
library(stringr)
library(decontam)
library(vegan)
library(WGCNA)
library(dplyr)


# 2. Define source directories
# ============================================================

source_dir <- file.path(
  "raw-data",
  "zenodo_source",
  "salmon",
  "metagenomics"
)

sample_info_file <- file.path(
  "raw-data",
  "zenodo_source",
  "salmon",
  "transcriptomics",
  "sample_info.csv"
)

required_files <- c(
  "feature-table_taxonomy.biom",
  "tree.nwk",
  "dna-sequences.fasta",
  "metadata.txt"
)

missing_files <- required_files[
  !file.exists(file.path(source_dir, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "Missing metagenomics source files:\n",
    paste(missing_files, collapse = "\n")
  )
}

if (!file.exists(sample_info_file)) {
  stop("Salmon sample metadata file not found: ", sample_info_file)
}


# 3. Import original datasets
# ============================================================

raw <- import_biom(
  file.path(source_dir, "feature-table_taxonomy.biom")
)

tree <- read_tree(
  file.path(source_dir, "tree.nwk")
)

refseq <- Biostrings::readDNAStringSet(
  file.path(source_dir, "dna-sequences.fasta")
)

dat <- read.table(
  file.path(source_dir, "metadata.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

all <- merge_phyloseq(
  raw,
  sample_data(dat),
  tree,
  refseq
)

cat("Imported samples:", nsamples(all), "\n")
cat("Imported taxa:", ntaxa(all), "\n")


# 4. Clean taxonomy
# ============================================================

tax <- as.data.frame(
  tax_table(all),
  stringsAsFactors = FALSE
)

tax <- tax[, 1:7]

colnames(tax) <- c(
  "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species"
)

tax.clean <- data.frame(
  row.names = rownames(tax),
  Kingdom = str_replace(tax[, 1], "d__", ""),
  Phylum  = str_replace(tax[, 2], "p__", ""),
  Class   = str_replace(tax[, 3], "c__", ""),
  Order   = str_replace(tax[, 4], "o__", ""),
  Family  = str_replace(tax[, 5], "f__", ""),
  Genus   = str_replace(tax[, 6], "g__", ""),
  Species = str_replace(tax[, 7], "s__", ""),
  stringsAsFactors = FALSE
)

tax.clean[is.na(tax.clean)] <- ""


# Class

bad <- c(
  "Ambiguous_taxa", "uncultured", "Subgroup_21"
)

tax.clean[tax.clean$Class %in% bad, 3:7] <- ""


# Order

bad <- c(
  "0319-6G20", "1-20", "11-24", "ADurb.Bin180",
  "D8A-2", "Group_1.1c", "JGI_0000069-P22",
  "Marine_Group_II", "Pla3_lineage", "Run-SP154",
  "Ambiguous_taxa", "uncultured", "UBA10353_marine_group",
  "Subgroup_17", "SAR86_clade", "SAR11_clade",
  "SAR202_clade"
)

tax.clean[tax.clean$Order %in% bad, 4:7] <- ""


# Family

bad <- c(
  "Ambiguous_taxa", "11-24", "67-14", "uncultured",
  "SAR116_clade", "Run-SP154", "Marine_Group_II",
  "env.OPS_17", "S085", "S-70",
  "NS9_marine_group", "Mitochondria"
)

tax.clean[tax.clean$Family %in% bad, 5:7] <- ""


# Genus

bad <- c(
  "Ambiguous_taxa", "Unknown_Family", "uncultured",
  "Subgroup_10", "1174-901-12", "67-14"
)

tax.clean[tax.clean$Genus %in% bad, 6:7] <- ""


# Species

bad <- c(
  "Ambiguous_taxa", "marine_metagenome", "low_GC",
  "wastewater_metagenome", "unidentified",
  "uncultured_synthetic", "uncultured_organism"
)

tax.clean[tax.clean$Species %in% bad, 6:7] <- ""


# Standardize taxonomy names

for (i in seq_len(ncol(tax.clean))) {

  tax.clean[, i] <- str_replace_all(
    tax.clean[, i], "[.]", ""
  )

  tax.clean[, i] <- str_replace_all(
    tax.clean[, i], "[(]", ""
  )

  tax.clean[, i] <- str_replace_all(
    tax.clean[, i], "[)]", ""
  )

  tax.clean[, i] <- str_replace_all(
    tax.clean[, i], "-", "_"
  )

  tax.clean[, i] <- str_replace_all(
    tax.clean[, i], " ", "_"
  )
}


# Fill missing taxonomic ranks

for (i in seq_len(nrow(tax.clean))) {

  if (tax.clean[i, 2] == "") {

    kingdom <- paste("Kingdom_", tax.clean[i, 1], sep = "")
    tax.clean[i, 2:7] <- kingdom

  } else if (tax.clean[i, 3] == "") {

    phylum <- paste("Phylum_", tax.clean[i, 2], sep = "")
    tax.clean[i, 3:7] <- phylum

  } else if (tax.clean[i, 4] == "") {

    class <- paste("Class_", tax.clean[i, 3], sep = "")
    tax.clean[i, 4:7] <- class

  } else if (tax.clean[i, 5] == "") {

    order <- paste("Order_", tax.clean[i, 4], sep = "")
    tax.clean[i, 5:7] <- order

  } else if (tax.clean[i, 6] == "") {

    family <- paste("Family_", tax.clean[i, 5], sep = "")
    tax.clean[i, 6:7] <- family

  } else if (tax.clean[i, 7] == "") {

    tax.clean$Species[i] <- paste(
      "Genus_", tax.clean$Genus[i], sep = "_"
    )
  }
}

stopifnot(!any(tax.clean == ""))

tax_table(all) <- as.matrix(tax.clean)


# Remove unassigned taxa and chloroplasts

all.clean <- subset_taxa(
  all,
  Kingdom != "Unassigned"
)

all.clean <- subset_taxa(
  all.clean,
  Order != "Chloroplast"
)

all.clean <- prune_taxa(
  taxa_sums(all.clean) > 0,
  all.clean
)

cat("Taxa after taxonomy cleaning:", ntaxa(all.clean), "\n")


# 5. Identify and remove contaminants
# ============================================================

df <- as.data.frame(sample_data(all.clean))

df$Sample_or_Control <- ifelse(
  df$New_Diet %in% c("ext-ctrl"),
  "Control_Sample",
  "True_Sample"
)

sample_data(all.clean) <- df

sample_data(all.clean)$is.neg <-
  sample_data(all.clean)$Sample_or_Control == "Control_Sample"

contamdf.prev05 <- isContaminant(
  all.clean,
  method = "prevalence",
  neg = "is.neg",
  threshold = 0.5
)

cat(
  "Identified contaminants:",
  sum(contamdf.prev05$contaminant, na.rm = TRUE),
  "\n"
)

all.noncontam <- prune_taxa(
  !contamdf.prev05$contaminant,
  all.clean
)


# Retain true samples only

psdata <- subset_samples(
  all.noncontam,
  Sample_or_Control == "True_Sample"
)

psdata <- prune_taxa(
  taxa_sums(psdata) > 0,
  psdata
)

cat("Samples after control removal:", nsamples(psdata), "\n")


# 6. Exclude T0 samples
# ============================================================

psdata.p <- subset_samples(
  psdata,
  samplingTime != "T0"
)

psdata.p <- prune_taxa(
  taxa_sums(psdata.p) > 0,
  psdata.p
)

cat("Samples after T0 removal:", nsamples(psdata.p), "\n")


# 7. Feature filtering and relative abundance
# ============================================================

psdata.r1 <- filter_taxa(
  psdata.p,
  function(x) sum(x > 4) > (0.01 * length(x)),
  TRUE
)

psdata.r1 <- transform_sample_counts(
  psdata.r1,
  function(x) x / sum(x)
)

cat("Samples after feature filtering:", nsamples(psdata.r1), "\n")
cat("Features after filtering:", ntaxa(psdata.r1), "\n")

stopifnot(
  nsamples(psdata.r1) == 144,
  ntaxa(psdata.r1) == 1041
)


# 8. Match microbiome samples with sample metadata
# ============================================================

sample_info <- read.csv(
  sample_info_file,
  row.names = 1
)

gut_16SrRNA_mapping_file <- data.frame(
  sample_data(psdata.r1)
)

gut_16SrRNA_mapping_file$ID_New <- paste(
  gut_16SrRNA_mapping_file$samplingTime, "_",
  gut_16SrRNA_mapping_file$New_Diet, "_",
  gut_16SrRNA_mapping_file$samplingTank, "_",
  "0",
  gut_16SrRNA_mapping_file$sampleNumber,
  sep = ""
)


# Check sample identifier uniqueness

stopifnot(
  !anyDuplicated(gut_16SrRNA_mapping_file$ID_New),
  !anyDuplicated(sample_info$ID_New)
)


# Match metadata to microbiome samples

host_gut_mapping_file <- sample_info

host_gut_mapping_file <- host_gut_mapping_file[
  host_gut_mapping_file$ID_New %in%
    gut_16SrRNA_mapping_file$ID_New,
]

gut_16SrRNA_mapping_file <- gut_16SrRNA_mapping_file[
  gut_16SrRNA_mapping_file$ID_New %in%
    host_gut_mapping_file$ID_New,
]


# Preserve original microbiome sample names

gut_16SrRNA_mapping_file$rownames <-
  rownames(gut_16SrRNA_mapping_file)

rownames(gut_16SrRNA_mapping_file) <-
  gut_16SrRNA_mapping_file$ID_New


# Extract abundance matrix

omics_data <- otu_table(psdata.r1)

ID_New_name_list <- gut_16SrRNA_mapping_file$rownames

gut_16SrRNA_OTU_table <- omics_data[
  ,
  colnames(omics_data) %in% ID_New_name_list
]


# Check sample alignment

stopifnot(
  all(
    colnames(gut_16SrRNA_OTU_table) ==
      gut_16SrRNA_mapping_file$rownames
  )
)


# Rename columns using matched identifiers

colnames(gut_16SrRNA_OTU_table) <-
  gut_16SrRNA_mapping_file$ID_New

cat("Matched samples:", ncol(gut_16SrRNA_OTU_table), "\n")
cat("Matched features:", nrow(gut_16SrRNA_OTU_table), "\n")

stopifnot(
  ncol(gut_16SrRNA_OTU_table) == 142,
  nrow(gut_16SrRNA_OTU_table) == 1041
)


# 9. Hellinger transformation
# ============================================================

omics_data <- t(gut_16SrRNA_OTU_table)

HellingerData <- vegan::decostand(
  omics_data,
  method = "hellinger"
)

omics_data <- HellingerData

stopifnot(
  nrow(omics_data) == 142,
  ncol(omics_data) == 1041,
  all(is.finite(omics_data))
)


# 10. Soft-threshold selection
# ============================================================

set.seed(13118)

powers <- c(1:10, seq(12, 20, 2))

sft <- suppressWarnings(
  WGCNA::pickSoftThreshold(
    omics_data,
    powerVector = powers,
    verbose = 1,
    networkType = "signed",
    corFnc = "bicor"
  )
)

signed_R2 <- -sign(sft$fitIndices[, 3]) *
  sft$fitIndices[, 2]

idx <- which(signed_R2 > 0.90)

if (length(idx) > 0) {

  idx <- idx[1]

} else {

  idx <- which(signed_R2 > 0.80)

  if (length(idx) > 0) {

    idx <- idx[1]

  } else {

    idx <- which.max(signed_R2)

  }
}

st <- sft$fitIndices[idx, 1]

cat("Selected soft-thresholding power:", st, "\n")

stopifnot(st == 9)


# 11. WGCNA module detection
# ============================================================

allowWGCNAThreads()


# Store TOM files outside the repository

tom_dir <- tempfile("OmniCorr_salmon_TOM_")

dir.create(tom_dir)

modules.omics.Y <- WGCNA::blockwiseModules(
  omics_data,
  power = st,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  minModuleSize = 12,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(tom_dir, "TOM"),
  verbose = 3,
  maxBlockSize = 8000,
  nThreads = 10
)

moduleLabels <- modules.omics.Y$colors

moduleColors <- labels2colors(
  modules.omics.Y$colors
)


# Document sparse-data bicor behavior:
# Features with zero MAD trigger an automatic Pearson
# correlation fallback in bicor. This behavior is retained
# to reproduce the original analysis.


cat("\nReproduced module distribution:\n")
print(table(moduleLabels))

stopifnot(
  length(moduleLabels) == 1041,
  identical(
    as.integer(table(factor(moduleLabels, levels = 0:7))),
    c(945L, 18L, 15L, 13L, 13L, 13L, 12L, 12L)
  )
)


# 12. Extract module eigengenes
# ============================================================

MEs <- modules.omics.Y$MEs


# Exclude the grey module

MEs <- MEs[
  ,
  colnames(MEs) != "ME0",
  drop = FALSE
]

cat("\nModule eigengene dimensions:\n")
print(dim(MEs))

stopifnot(ncol(MEs) == 7)


# 13. Aggregate module eigengenes by tank
# ============================================================

tank_ids <- sub(
  "_[0-9]{1,2}$",
  "",
  rownames(MEs)
)

aggregated <- as.data.frame(MEs) %>%
  mutate(common = tank_ids) %>%
  group_by(common) %>%
  summarise(
    across(everything(), mean),
    .groups = "drop"
  )

Metagenomics_reproduced <- as.matrix(
  aggregated[, -1, drop = FALSE]
)

rownames(Metagenomics_reproduced) <- aggregated$common

cat("\nFinal Metagenomics dimensions:\n")
print(dim(Metagenomics_reproduced))


# 14. Compare with the packaged dataset
# ============================================================

data(
  "Metagenomics",
  package = "OmniCorr"
)

reference <- as.matrix(
  Metagenomics[
    rownames(Metagenomics_reproduced),
    colnames(Metagenomics_reproduced),
    drop = FALSE
  ]
)

comparison <- all.equal(
  reference,
  Metagenomics_reproduced,
  tolerance = 1e-10,
  check.attributes = FALSE
)

cat("\nComparison with packaged Metagenomics dataset:\n")
print(comparison)

stopifnot(
  identical(
    dim(reference),
    dim(Metagenomics_reproduced)
  ),
  isTRUE(comparison)
)

message(
  "\nSUCCESS: Metagenomics dataset reproduced exactly."
)
