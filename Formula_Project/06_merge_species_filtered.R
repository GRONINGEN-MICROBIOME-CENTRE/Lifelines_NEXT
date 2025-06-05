# Load libraries
library(here)
library(tidyverse)
library(readr)

# Load cleaned nutrient + metadata dataset from Script 5
data2 <- readRDS(here("scripts/GITHUB/output/data2.rds"))

# Load filtered CLR-transformed species abundance data
species_filtered <- read_tsv(here("data/clr_trans_with_prev_and_detection_filters_10april2025.tsv"))

# Optional: Check number of species columns
cat("Number of species columns in filtered data:", ncol(species_filtered), "\n")

# Merge datasets by NG_ID
merged_data <- data2 %>%
  inner_join(species_filtered, by = "NG_ID")

# Save the final merged dataset
saveRDS(merged_data, file = "output/df_NG_ID_continuous_species_filtered.rds")
