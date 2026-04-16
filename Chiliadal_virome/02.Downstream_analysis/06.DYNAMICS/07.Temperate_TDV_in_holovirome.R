

setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# the plan is to move all PPV holovirome shit here
# please make a MVP, not whole-blown analysis, I am begging you
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################

#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(MetBrewer)
#############################################################
# 2. Load Input Data
#############################################################
# metadatas:
smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)

smeta <- smeta %>%
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) %>% # treating 2nd pregnancy as a separate family:
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)) %>%
  group_by(NEXT_ID, seq_type) %>%
  mutate(n_tp = n()) %>%
  ungroup() %>%
  filter(n_tp >= 3)

# abundance tables etc:
# VLP:
VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)
VLP <- VLP[,colnames(VLP) %in% smeta$Sequencing_ID]
VLP <- VLP[rowSums(VLP) > 0,]

# HLV:
holovirome <- read.table('06.CLEAN_DATA/Intermediate/Holovirome_RPKM_1110samples_120997vOTUs.txt', sep='\t', header=T)
holovirome <- holovirome[,colnames(holovirome) %in% smeta$Universal_ID]
holovirome <- holovirome[rownames(holovirome) %in% rownames(VLP),]
holovirome <- holovirome[rowSums(holovirome) > 0,]

# virus metadata
ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)
#############################################################
# 3.0 Analysis: VLP TDV in HLV
#############################################################
# first, infant-wise vOTU detection in VLP:

VLP_detection <- VLP %>%
  rownames_to_column("vOTU") %>%
  pivot_longer(-vOTU, names_to = "Sequencing_ID", values_to = "RPKM") %>%
  filter(RPKM != 0) %>% # essential, otherwise some merging will run out of memory
  left_join(smeta %>% select(Sequencing_ID, NEXT_ID, Universal_ID, Timepoint_new), by = "Sequencing_ID") %>%
  mutate(VLP_detected = T) %>%
  group_by(NEXT_ID, vOTU) %>%
  summarise(n_tp_VLP = n_distinct(Timepoint_new), .groups = "drop") %>%
  mutate(is_persistent_VLP = ifelse(n_tp_VLP >= 3, T, F),
         VLP_ever_detected = T) %>%
  filter(vOTU %in% ETOF_vOTUr$New_CID[ETOF_vOTUr$lifestyle == "Temperate"]) # keeping only temperates

# holovirome long:
HLV_long <- holovirome %>%
  rownames_to_column("vOTU") %>%
  pivot_longer(-vOTU, names_to = "Universal_ID", values_to = "RPKM") %>%
  filter(RPKM != 0) %>% # essential, otherwise some merging will run out of memory
  left_join(smeta %>% select(Sequencing_ID, NEXT_ID, Universal_ID, Timepoint_new), by = "Universal_ID")

# persistent temperate regardles of metavirome:
HLV_persistent_lookup <- HLV_long %>%
  group_by(NEXT_ID, vOTU) %>%
  summarise(n_tps_HLV = n_distinct(Timepoint_new), .groups = "drop") %>%
  filter(n_tps_HLV >= 3) %>%
  mutate(is_persistent_HLV = TRUE) %>%
  filter(vOTU %in% ETOF_vOTUr$New_CID[ETOF_vOTUr$lifestyle == "Temperate"]) %>%
  left_join(VLP_detection) %>%
  mutate(n_tp_VLP = replace_na(n_tp_VLP, 0),
         VLP_ever_detected = replace_na(VLP_ever_detected, F),
         is_persistent_VLP = replace_na(is_persistent_VLP, F)) %>% 
  filter(VLP_ever_detected == T & is_persistent_VLP == F) # necessary and sufficient conditions for TDV in VLP!

PPV_HLV <- HLV_persistent_lookup %>%
  filter(is_persistent_VLP != T) %>%
  group_by(NEXT_ID) %>%
  summarise(PPV_size_HLV = n_distinct(vOTU), .groups = "drop") %>%
  right_join(smeta %>% group_by(NEXT_ID) %>% summarise(n_tp = n_distinct(Timepoint_new), .groups = "drop")) %>%
  mutate(PPV_size_HLV = replace_na(PPV_size_HLV, 0))

testik <- map_dfr(c("PPV_size_HLV"), function(variable) {
  
  F1 <- as.formula(paste0(variable, " ~ n_tp"))
  
  saver <- t.test(F1, data = PPV_HLV)
  
  data.frame(tested = variable,
             p_value = saver$p.value)
  
} ) 

# no need to rarefy, because it does not differ between 3 and 4 timepoints (I guess because we take a subset of a subset)

summary(PPV_HLV$PPV_size_HLV)



