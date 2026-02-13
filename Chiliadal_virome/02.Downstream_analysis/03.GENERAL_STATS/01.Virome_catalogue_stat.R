setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

#############################################################
# Descriptive statistics on NEXT virome catalogue
# Staton on novels
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
get_iqr <- function(vec) {
  
  SMR <- summary(vec)
  
  print(paste0(round(SMR[3], 2), ' (', round(SMR[2], 2), ' - ', round(SMR[5], 2), ')')) 
  
}
#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(purrr)
library(tidyverse)
#############################################################
# 2. Load Input Data
#############################################################
ETOF <- read.table('06.CLEAN_DATA/VLP_MGS_ETOF_full_rep.txt', sep='\t', header = T)

ETOF_only <- ETOF %>%
  filter(grepl('NEXT_', ETOF$New_CID)) %>%
  mutate(sample = gsub('_.*', '', Original_CID)) %>%
  mutate( method = ifelse(grepl('NEXT_V', New_CID), 'VLP', 'MGS') )

smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)
smeta$Timepoint_new <- factor(smeta$Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)
#############################################################
# 3.1 Analysis: redundant & hq genomes per method
#############################################################
genomic_stat <- ETOF_only %>%
  group_by(method) %>%
  summarise(N_total = n()) %>%
  left_join(ETOF_only %>%
              filter(miuvig_quality == "High-quality") %>%
              group_by(method) %>%
              summarise(N_HQ = n()) ) %>%
  left_join(ETOF_only %>%
              filter(miuvig_quality == "High-quality" & POST_CHV_length >= 10000) %>%
              group_by(method) %>%
              summarise(N_HQ10 = n())) %>%
  pivot_longer(cols = !method) %>%
  pivot_wider(names_from = method,
              values_from = value) %>%
  rename(feature = name) %>%
  mutate(across(where(is.numeric), as.character ))

sample_stat <- ETOF_only %>%
  group_by(sample) %>%
  summarise(N_discovered = n()) %>%
  left_join(ETOF_only %>%
              filter(miuvig_quality == "High-quality") %>%
              group_by(sample) %>%
              summarise(N_discovered_HQ = n())) %>%
  left_join(ETOF_only %>%
              filter(miuvig_quality == "High-quality" & POST_CHV_length >= 10000) %>%
              group_by(sample) %>%
              summarise(N_discovered_HQ10 = n()))

smeta <- smeta %>%
  left_join(sample_stat, by = c("Sequencing_ID" = "sample"))

disc_staton <- map_dfr(c('MGS', 'VLP'), function(metavirome){
  
  lengths <- ETOF_only %>%
    filter(miuvig_quality == "High-quality" & method == metavirome) %>%
    pull(POST_CHV_length) %>%
    get_iqr()
  
  lengths_df <- data.frame(median = lengths,
                           method = metavirome,
                           feature = "Genome_length")
  
  metrics_df <- map_dfr(c("N_discovered", "N_discovered_HQ", "N_discovered_HQ10"), function(metric) {
    
    med_and_iqr <- smeta %>%
    filter(seq_type == metavirome) %>%
      pull(metric) %>%
      get_iqr()
    
    data.frame(median = med_and_iqr,
               method = metavirome,
               feature = metric)
    
  })
  
  rbind(lengths_df, metrics_df)
})

disc_staton <- genomic_stat %>%
  bind_rows(disc_staton %>%
              pivot_wider(names_from = method, 
                          values_from = median) %>%
              mutate(feature = if_else(feature != "Genome_length", paste0(feature, "_per_sample"), feature)))

writexl::write_xlsx(disc_staton, '07.RESULTS/Virus_discovery_summary_stat.xlsx')
# tidying up
rm(list=c("ETOF_only", "genomic_stat", "sample_stat"))


