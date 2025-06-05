# Load libraries
library(here)
library(tidyverse)

# Load raw species abundance data
species_raw <- read.table(here("data/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt"),
                          header = TRUE, sep = "\t", row.names = NULL)

# Load final merged dataset
data <- readRDS(here("scripts/GITHUB/output/df_NG_ID_continuous_species_filtered.rds"))

# Rename first column to NG_ID
colnames(species_raw)[1] <- "NG_ID"

# Filter to matching NG_IDs
species_filtered <- species_raw %>%
  filter(NG_ID %in% data$NG_ID)

# Identify species-level columns (ending in .s__)
species_cols <- grep("\\.s__.*$", colnames(species_filtered), value = TRUE)

# Calculate richness (number of species with abundance > 0)
species_filtered <- species_filtered %>%
  mutate(richness = rowSums(across(all_of(species_cols), ~ . > 0)))

# Merge richness into main dataset
data <- data %>%
  left_join(species_filtered %>% dplyr::select(NG_ID, richness), by = "NG_ID")

# Reorder columns to place richness after shannon
data <- data %>%
  relocate(richness, .after = shannon)

# Save final dataset
saveRDS(data, file = "output/NG_ID_continuous_species_richness.rds")
