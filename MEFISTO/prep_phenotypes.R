source(here("run_10FAC_BAC_NOGRP_14:59:36_26_01_2024_copy/Lifelines_NEXT_gut_microbiome_taxa_cross_sectional_infants.R"))

head(phenos)
names(phenos)

pheno_parity <- names(phenos)[19] #parity

#examine phenos if needed
as_tibble(phenos) %>%
  select(pheno_parity) %>% print(n = nrow(.))
  ggplot() +
  geom_histogram(aes_string(x= pheno_parity))

phenos <-
as_tibble(phenos) %>% 
mutate(mother_birthcard_parity = case_when( # categorize parity
    mother_birthcard_parity < 0 ~ "First",
    mother_birthcard_parity >= 0 & mother_birthcard_parity < 1 ~ "Second",
    mother_birthcard_parity >= 1 ~ "Third+",
    TRUE ~ NA
  )) %>%
  mutate(mother_birthcard_parity = factor(mother_birthcard_parity, levels = c("First", "Second", "Third+"))) %>%
  mutate(across(where(is.factor), as.character)) %>% 
  mutate(across(where(is.character), ~na_if(.x, "<NA>"))) %>% #replace <NA> values with NA
  select(-c(birth_delivery_mode_complex, birth_delivery_mode_simple, infant_ffq_ever_never_breastfed)) %>% # remove overlapping phenotypes
  mutate(NG_ID_short = substr(NG_ID, 1, nchar(NG_ID) - 4))
# names(phenos)

# inspect if needed
# test[1:5,1:7]
# test[1:5,8:14]
# test[1:5,15:17]
