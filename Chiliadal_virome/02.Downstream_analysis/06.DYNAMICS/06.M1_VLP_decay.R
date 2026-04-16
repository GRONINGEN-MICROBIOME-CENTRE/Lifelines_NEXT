setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# Here we work specifically on M1-colonizers decay, 
# because I did not want to mess 05.PPV script
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
# longitudinal vOTU sharing
bagira <- function(df, diag_incl, meta, merging_ID){
  
  binary_matrix <- as.matrix(df > 0) + 0 
  
  shared_species_matrix <- t(binary_matrix) %*% binary_matrix
  
  upper_tri_idx <- which(upper.tri(shared_species_matrix, diag = diag_incl), arr.ind = T)
  
  sharing <- data.frame(
    sample_1 = rownames(shared_species_matrix)[upper_tri_idx[, 1]],
    sample_2 = colnames(shared_species_matrix)[upper_tri_idx[, 2]],
    N_shared = shared_species_matrix[upper_tri_idx]) %>%
    left_join(meta %>% select(Sequencing_ID,
                              Universal_ID,
                              Modified_NEXT_ID_without_preg_number,
                              FAMILY, 
                              FAMILYupd,
                              Type, 
                              seq_type,
                              Timepoint_new), by = c("sample_1" = merging_ID )) %>%
    left_join(meta %>% select(Sequencing_ID,
                              Universal_ID,
                              Modified_NEXT_ID_without_preg_number,
                              FAMILY,
                              FAMILYupd,
                              Type, seq_type,
                              Timepoint_new), by = c("sample_2" = merging_ID), suffix = c("_1", "_2")) %>%
    select(sample_1, sample_2, N_shared, sort(colnames(.))) %>%
    mutate(longitudinal = ifelse(FAMILYupd_1 == FAMILYupd_2 &
                                   Modified_NEXT_ID_without_preg_number_1 == Modified_NEXT_ID_without_preg_number_2, T, F))
  
  return(sharing)
}
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
  filter(Type == "K")

# abundance tables etc:
VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)
VLP <- VLP[,colnames(VLP) %in% smeta$Sequencing_ID]
VLP <- VLP[rowSums(VLP) > 0,]

# virus metadata
ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)
#############################################################
# 3.0 Analysis: M1 colonizers dynamics
#############################################################
sharing <- bagira(VLP, T, smeta, "Sequencing_ID")

baring <- sharing %>%
  filter(longitudinal == T) %>%
  mutate(
    temp_s1 = sample_1,
    sample_1 = ifelse(Timepoint_new_1 != "M1", sample_2, sample_1),
    sample_2 = ifelse(Timepoint_new_1 != "M1", temp_s1, sample_2)
  ) %>%
  select(-temp_s1) %>%
  select(sample_1, sample_2, N_shared, longitudinal) %>%
  left_join(smeta %>% select(Universal_ID, Sequencing_ID, vir_richness_cf, Timepoint_new) %>% rename(M1_richness = vir_richness_cf), by = c("sample_1" = "Sequencing_ID")) %>%
  filter(grepl('M1$', Universal_ID)) %>%
  select(-Universal_ID) %>%
  left_join(smeta, by = c("sample_2" = "Sequencing_ID"), suffix = c("", "_2")) %>%
  mutate(prop_pres = N_shared / M1_richness,
         prop_rich_by_M1 = N_shared / vir_richness_cf) %>%
  relocate(Timepoint_new_2)

M1_in_VLP <- baring %>%
  mutate(Timepoint_numeric = as.integer(as.numeric(Timepoint_new_2))) %>%
  ggplot(aes(x = Timepoint_new_2, y = prop_pres*100)) +
  geom_line(aes(group = NEXT_ID), alpha = 0.7, linewidth = 0.2, color = "#FF9797") +
  geom_smooth(aes(Timepoint_numeric, prop_pres*100), method = "loess",  span = 1, linewidth = 1, se = TRUE, alpha = 0.3, fill = "#360982", color = "#360982") +
  ggnewscale::new_scale_color() +
  stat_summary(aes(x = Timepoint_new_2, 
                   group = Timepoint_new_2), color = "#B72A67",
               fun.data = mean_cl_boot, 
               geom = "errorbar", 
               width = 0.25,
               size = 0.5) +
  stat_summary(aes(x = Timepoint_new_2,
                   group = Timepoint_new_2), color = "#B72A67",
               fun = mean, 
               geom = "point", 
               size = 1) +
  theme_minimal() + 
  labs(x = "Timepoint", y = "% of preserved M1 vOTUs") +
  ylim(0, 100) + 
  theme(axis.title = element_text(size = 9), 
        axis.text = element_text(size = 8))

summary(baring$prop_pres[baring$Timepoint_new_2 == "M3"])
sd(baring$prop_pres[baring$Timepoint_new_2 == "M3"])

ggsave("05.PLOTS/06.DYNAMICS/M1_decay_proportion.pdf",
       M1_in_VLP, "pdf", width = 6, height = 7, units = "cm", dpi = 300)


M1_in_VLP_abs <- baring %>%
  mutate(Timepoint_numeric = as.integer(as.numeric(Timepoint_new_2))) %>%
  ggplot(aes(x = Timepoint_new_2, y = N_shared)) +
  geom_line(aes(group = NEXT_ID), alpha = 0.7, linewidth = 0.2, color = "#FF9797") +
  geom_smooth(aes(Timepoint_numeric, N_shared), method = "loess",  span = 1, linewidth = 1, se = TRUE, alpha = 0.3, fill = "#360982", color = "#360982") +
  ggnewscale::new_scale_color() +
  stat_summary(aes(x = Timepoint_new_2,
                   group = Timepoint_new_2), color = "#B72A67",
               fun.data = mean_cl_boot,
               geom = "errorbar",
               width = 0.25,
               size = 0.5) +
  stat_summary(aes(x = Timepoint_new_2,
                   group = Timepoint_new_2), color = "#B72A67",
               fun = mean,
               geom = "point",
               size = 0.5) +
  theme_minimal() +
  labs(x = "Timepoint", y = "N M1-detected vOTUs") +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))

ggsave("05.PLOTS/06.DYNAMICS/M1_decay_N_vOTUs.pdf",
       M1_in_VLP_abs, "pdf", width = 6, height = 7, units = "cm", dpi = 300)

### decays with time? (decided to go with absolute numbers for this analysis; otherwise there is a p-value inflation for proportions)
summary(lmer(N_shared ~ Timepoint_new_2 + (1|NEXT_ID), data = baring, REML = F))$coefficients

