# Reproduce the salmon metadata dataset.
# Run from the root of the OmniCorr repository.

library(dplyr)

sample_info <- read.csv(
  "raw-data/zenodo_source/salmon/transcriptomics/sample_info.csv",
  row.names = 1
)

# Phenotypic variables retained in the packaged dataset
variables <- c(
  "BW", "Length", "Sex", "Gut_struc", "Gut_light",
  "Anus_bleed", "Gut_bleed", "Fat_Score", "Heart_Sco",
  "softliver", "patchy_liver", "Liver_Sco", "HSI", "CSI"
)

# Calculate group means, preserving missing values
metadata <- sample_info %>%
  mutate(Group = paste(Time, New_Diet, Tank_number, sep = "_")) %>%
  group_by(Group) %>%
  summarise(
    across(all_of(variables), ~ mean(.x)),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(metadata) <- metadata$Group
metadata$Group <- NULL

# Save the reproduced metadata separately
output_dir <- "raw-data/metadata"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  metadata,
  file.path(output_dir, "metadata_reproduced.rds")
)
