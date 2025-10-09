# Purpose: This script loads the sequencing (i.e., METAPHLAN data) and also creates a feature matrix of that data, which can be used in downstream analysis.

# opening line ----
print("loading bacterial relative abundance data...")

# load abundance data ---
data_path <- here("data_sequencing/LLNEXT_metaphlan_4_CLR_transformed_fil_all_levels_infants_20_07_2023.txt")

#rows are samples
#columns are species found
data <- read.delim(data_path, header = TRUE)

#samples IDs are not yet assigned to tibble
#add them
data$NG_ID <- rownames(data)

#make data into a tibble
data <- as_tibble(data) 

# remove columns with zero variance----

#N.B. uncomment if using untransformed data

# calc column variances
#col_variances <- data %>% summarize(across(-last_col(), var, na.rm=TRUE))

# get columns with zero variance
#zero_var_cols <- colnames(col_variances)[col_variances == 0]

# create new data tibble with zero variances features that are filtered
#data <- data %>% select(-one_of(zero_var_cols))

# create feature matrix from data
feature_matrix <-
  as_tibble(colnames(data)) %>% 
  rename(taxon = value) %>%
  filter(taxon != "NG_ID") %>% 
  mutate(taxonID = NA) %>%
  mutate(taxonID = purrr::map_chr(taxonID, ~ stri_rand_strings(1, length = 7, pattern = "[A-Za-z0-9]"))) %>% #assign unique ID to each taxon 
  relocate(taxonID)

# check unique codes to make sure there are no duplicates
potential_duplicates <- feature_matrix %>%
  filter(duplicated(taxonID) | duplicated(taxonID, fromLast = TRUE))

if (any(duplicated(potential_duplicates$taxonID))) {
  stop("Duplicated taxonID values found. Script stopped.")
}

# save feature matrix
feature_matrix_file_name <- paste0(here("feature_matrix_"), current_date, ".rda")
save(feature_matrix, file = feature_matrix_file_name)

# closing line ----
print("loading bacterial relative abundance data COMPLETED.")
