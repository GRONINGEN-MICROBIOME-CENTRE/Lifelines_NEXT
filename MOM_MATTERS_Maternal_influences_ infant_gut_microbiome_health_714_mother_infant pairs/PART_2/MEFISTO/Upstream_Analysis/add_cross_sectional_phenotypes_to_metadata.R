# Purpose: This script creates the metadata that is used in downstream MEFISTO analysis.
# N.B. some babies may have duplicate mother IDs because twins.

# opening line ----
print("loading and merging meta data and cross sectional phenotypes...")

#load meta data ----
meta_data_path <- here("data_meta/LLNEXT_metadata_family_structure_28_08_2023.txt")

meta_data <- read.delim(meta_data_path, header = TRUE)

#load and tidy cross-sectional phenotypes ----
data_cross <- read.delim2(file = here("data_phenotypes/masterfile_cross_sectional_2023_09_29.txt"), stringsAsFactors = TRUE)

data_cross <-
data_cross %>%
  rename(NEXT_ID = next_id_infant)


# combine cross-sectional phenotypes with meta_data ----
meta_data <-  
  meta_data %>%
  filter(Type == "infant") %>%
  left_join(., data_cross, by = c("NEXT_ID", "FAMILY")) 

#remove unneeded objects and/or functions from global environment
rm(data_cross)

# closing line ----
print("loading and merging meta data and cross sectional phenotypes COMPLETED.")

