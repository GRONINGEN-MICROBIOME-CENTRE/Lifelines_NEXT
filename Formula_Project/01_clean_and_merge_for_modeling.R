# 01_clean_and_merge_for_modeling.R

# Load libraries
library(here)
library(tidyverse)
library(dplyr)

# Load final dataset
data <- readRDS(here("scripts/GITHUB/output/NG_ID_continuous_species_richness.rds"))

# Function to clean and rename columns
rename_columns <- function(col_names) {
  sapply(col_names, function(name) {
    prefix <- if (grepl("^Continuous_", name)) "continuous_" else if (grepl("^median group_", name)) "median_group_" else return(name)
    component <- sub("^.*?_", "", name)
    component <- sub(" \\(.*", "", component)
    unit <- if (grepl("mg/100mL|mg/100ml|µg/100mL|µg/100ml", name)) "_mg"
    else if (grepl("gr/100mL|gr/100ml", name)) "_g" else ""
    paste0(prefix, gsub(" ", "_", component), unit)
  })
}

# Apply renaming
names(data) <- rename_columns(names(data))

# Load phenotype data
data_pheno <- read_delim(here("data/masterfile_cross_sectional_2023_11_15.txt"))

# Select and rename relevant column
data_pheno2 <- data_pheno %>%
  dplyr::select(next_id_infant, birth_deliverybirthcard_mode_binary) %>%
  rename(NEXT_ID = next_id_infant)

# Merge phenotype into main data
data2 <- data %>%
  as_tibble() %>%
  left_join(data_pheno2, by = "NEXT_ID")

# Save cleaned dataset
saveRDS(data2, file = here("scripts/GITHUB/output/NG_ID_shannon_richness_continuous_median_species_birth_mode.rds"))
