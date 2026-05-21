############################################################
# Prep data and prep splits/folds for mefisto stability analysis
# Factors: 1–15, Folds: 1–5
############################################################

# load libraries

library(here)
library(tidyverse)
library(MOFA2)

# load meta and helper data----

sample_IDs <- readRDS(here("sample_info/NG_IDs_infants_with_cazyme_data.rds"))

meta_data_path <- here("data_meta/LLNEXT_metadata_family_structure_28_08_2023.txt")
data_meta <- read.delim(meta_data_path, header = TRUE)

data_meta2 <- data_meta %>% 
  filter(Type == "infant") %>% 
  select(NG_ID, NEXT_ID)

data_meta3 <- as_tibble(data_meta) %>% 
  filter(Type == "infant") %>% 
  select(NG_ID, Timepoint_categorical)

# check
all(sample_IDs %in% data_meta2$NG_ID)

# load microbiome data----

data_microbiome <- read.delim(
 file = here("data_sequencing/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
)


data_microbiome <- data_microbiome[rownames(data_microbiome) %in% sample_IDs,]

data_microbiome2 <- data_microbiome %>%
  rownames_to_column("NG_ID") %>%
  as_tibble() %>%
  pivot_longer(-NG_ID, names_to = "feature", values_to = "value") %>%
  left_join(data_meta2, by = "NG_ID") %>% 
  relocate(NEXT_ID) %>%
  mutate(view = "microbiome") %>%
  left_join(data_meta3, by = "NG_ID")


# add timepoint and fill missing combinations----

data_microbiome3 <- data_microbiome2 %>%
  mutate(time = Timepoint_categorical) %>%
  mutate(
    time_numeric = case_when(
      time == "W2"  ~ 0.5,
      time == "M1"  ~ 1,
      time == "M2"  ~ 2,
      time == "M3"  ~ 3,
      time == "M6"  ~ 6,
      time == "M9"  ~ 9,
      time == "M12" ~ 12,
      TRUE          ~ NA_real_
    )
  )

all_groups   <- unique(data_microbiome3$NEXT_ID)
all_times    <- unique(data_microbiome3$time)
all_features <- unique(data_microbiome3$feature)

full_grid <- expand_grid(
  NEXT_ID = all_groups,
  time = all_times,
  feature = all_features
)

data_filled_mic <- full_grid %>%
  left_join(data_microbiome3, by = c("NEXT_ID", "time", "feature")) %>%
  mutate(
    sample = paste(NEXT_ID,time,sep="_"),
    view = "microbiome"
  ) %>%
  mutate(group = NEXT_ID) %>%
  select(NG_ID, view, sample, time, feature, value, NEXT_ID)


# shape data for MEFISTO----
data_for_MEFISTO <-
  bind_rows(data_filled_mic) %>% 
  #print(n = 300)
  relocate(view) %>%
  #select(-NG_ID) %>%
  mutate(
    time = case_when(
      time == "W2"  ~ 0.5,
      time == "M1"  ~ 1,
      time == "M2"  ~ 2,
      time == "M3"  ~ 3,
      time == "M6"  ~ 6,
      time == "M9"  ~ 9,
      time == "M12" ~ 12,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(is.na(NG_ID) | NG_ID != "NG_IDX") %>% # duplicates
  select(-NG_ID)


# compute splits
set.seed(123)

individuals <- unique(data_for_MEFISTO$NEXT_ID)
n_ind <- length(individuals)
n_boot <- 15         
prop <- 0.80

bootstrap_sets <- vector("list", n_boot)

for (b in 1:n_boot) {
  boot_ids <- sample(individuals,
                     size = round(prop * n_ind),
                     replace = FALSE)   # no replacement
  bootstrap_sets[[b]] <- boot_ids
}

saveRDS(bootstrap_sets, "bootstrap_sets_80pct_n20.rds")

str(bootstrap_sets)


