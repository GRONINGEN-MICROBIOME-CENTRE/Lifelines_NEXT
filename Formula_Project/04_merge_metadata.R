# Load libraries
library(dplyr)
library(tidyverse)
library(here)

# Load nutrient data
NG_ID_continuous <- readRDS("output/df_NG_ID_continuous_median_long.rds")

# Load metadata
metadata <- read_delim(here("data/LLNEXT_metadata_15_04_2024.txt"))

# Merge datasets by NG_ID
merged_nutrient_microbiome <- NG_ID_continuous %>%
  left_join(metadata, by = "NG_ID")

# Save merged dataset
saveRDS(merged_nutrient_microbiome, file = "output/df_NG_ID_continuous_metadata.rds")
