# Purpose: Create the tibble that will be input to MEFISTO function, which will trian the model.

# opening line ----
print("creating tibble for MOFA object...")

# create bacterial data ----
data_bacteria <- 
  data %>%
  pivot_longer(col = -c(NG_ID)) %>%
  rename(taxon = name) %>% #nrow = # nrow = 30195680
  full_join(., feature_matrix, by = "taxon") %>% # nrow = 30195680
  full_join(.,meta_data, by = "NG_ID") %>% # nrow = 30195680
  rename(group = NEXT_ID,
         sample = NG_ID,
         feature = taxonID,
         value = value) %>% 
  mutate(time = if_else(Timepoint_categorical == "W2", 0.5, as.numeric(str_extract(Timepoint_categorical, "\\d+")))) %>%
  mutate(view = "microbiome",
         feature = as.character(feature)) %>%
  #select(sample, group, FAMILY, exact_age_months_at_collection, shannon, feature, view, value) %>%
  filter(!is.na(group)) %>%
  select(!group) %>%
  relocate(time, feature, value, view) %>%
  mutate(view = as.factor(view)) %>%
  mutate(feature = as.factor(as.character(feature))) %>%
  mutate(sample = as.factor(sample))


# sanity check

# how many samples are there?
# how many individuals are th
#n_samples <-
#data_bacteria %>%
#  reframe(n_samples = unique(sample)) %>%
#  nrow()

#n_individuals <-
#  data_bacteria %>%
#  reframe(n_individuals = unique(group)) %>%
#  nrow()

#print(paste("This run is analyzing", as.character(n_samples), "different samples.",
#            "Amongst these samples are", as.character(n_individuals), "individuals."))
  
# do all samples contain the same number of features?
tib <-
data_bacteria %>%
  group_by(sample) %>%
  reframe(n_features_per_sample = n()) 

are_all_values_the_same <- all(length(unique(tib$n_features_per_sample)) == 1)
  
if (are_all_values_the_same == FALSE) {
  stop("ERROR: The number of features is not equal across all samples.")
}

# put all data views together ----
data_for_MOFA <- data_bacteria

# save final tibble ----

# make file name
file_name <- paste0(here("data_for_MOFA_"), current_date, ".rda")
#file_name <- paste0(here("data_for_MOFA_"), current_date, ".rda")

# save tibble for later use, if needed
save(data_for_MOFA, file = file_name)

# remove all unneeded object to free memory

all_objects <- ls()

object_to_keep <- c("data_for_MOFA", "max_age", "feature_matrix","current_date")

objects_to_remove <- setdiff(all_objects, object_to_keep)

rm(list = objects_to_remove)

# closing line ----
print("creating tibble for MOFA object COMPLETED.")
