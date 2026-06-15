setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# Here we explore which phenotypes shape the virome
#############################################################

#############################################################
# 0. Used files source
#############################################################
phenos_select <- read.table("06.CLEAN_DATA/Intermediate/Phenotype_selection.txt", sep='\t', header = T)
#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(tidyverse)
library(MetBrewer)
library(vegan)
library(skimr)
library(lme4)
library(lmerTest)
#############################################################
# 2. Load Input Data
#############################################################
# metadatas:
smeta_w_phenos <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_phenotypes.txt', sep='\t', header=T)

smeta_w_phenos <- smeta_w_phenos %>%
  filter(Type == "K") %>%
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) %>% # treating 2nd pregnancy as a separate family:
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T))

metaphlan <- read.table("06.CLEAN_DATA/02.FINAL/MGS_Chiliadal_metaphlan_full_taxonomy_ver_01_07102025.txt", sep='\t', header=T)

# Longitudinal phenotypes:
L_phenos <- read.table('06.CLEAN_DATA/Phenotypes/masterfile_longitudinal_2023_09_29.txt', sep='\t', header=T)
L_phenos <- L_phenos[L_phenos$SAMPLE_ID %in% smeta_w_phenos$Universal_ID,]

# Cross-sectional phenotypes:
C_phenos <- read.table('06.CLEAN_DATA/Phenotypes/masterfile_cross_sectional_2023_11_15.txt', sep = '\t', header=T)
C_phenos <- C_phenos[C_phenos$FAMILY %in% smeta_w_phenos$FAMILY,]

long_smeta <- read.delim("06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt", sep='\t', header=T)
long_smeta <- long_smeta %>%
  filter(seq_type == "MGS" & Type == "K") %>%
  left_join(C_phenos %>% 
              select(all_of(phenos_select %>%
                              filter(type == "cross-sectional") %>%
                              pull(phenotype)
                            ), next_id_infant), 
            by = c("NEXT_ID" = "next_id_infant")) %>%
  left_join(L_phenos %>% 
              select(all_of(phenos_select %>%
                              filter(type == "longitudinal") %>%
                              pull(phenotype)), next_id_infant, timepoint), 
            by = c("NEXT_ID" = "next_id_infant", "Timepoint_new" = "timepoint"),
            suffix = c("_cross", "_long")) %>%
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12"), ordered = T))


metaphlan_genus <- metaphlan[grep('g__', row.names(metaphlan)),] 
metaphlan_genus <- metaphlan_genus[-grep('s__', row.names(metaphlan_genus)),] 

metaphlan_genus <- metaphlan_genus[rowSums(metaphlan_genus) > 0,]

#############################################################
# 3.1 Analysis: 
#############################################################
beta <- as.data.frame(t(metaphlan_genus))

gamma <- beta %>%
  select(grep('Bacteroides', colnames(beta))) %>%
  rename("Bacteroides" = 1) %>%
  rownames_to_column(var = "Sequencing_ID") %>%
  left_join(long_smeta)

Bacteroides_abundance <- gamma %>%
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)) %>%
  ggplot(aes(Timepoint_new, Bacteroides)) +
  ggrastr::rasterise(geom_jitter(color = "#b12625", position = position_jitterdodge(), size=0.1, alpha=1), dpi = 300)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, fill = "#b12625") +
  theme_minimal() +
  theme(strip.background = element_rect(NA),
        legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm")) + 
  labs(y = "Bacteroides abundance", x="Timepoint")

ggsave("05.PLOTS/06.DYNAMICS/Bacteroides_genus_abundance.pdf",
       Bacteroides_abundance, "pdf", width=6, height=7, units="cm", dpi = 300)

summary(lmer(log(Bacteroides + 0.0003) ~ Timepoint_new + birth_deliverybirthcard_mode_binary + (1|NEXT_ID), data = gamma))

length(unique(long_smeta$Sequencing_ID[!is.na(gamma$birth_deliverybirthcard_mode_binary)]))


Bacteroides_DM <- gamma %>%
  right_join(smeta_w_phenos %>% select(Universal_ID, birth_deliverybirthcard_mode_binary)) %>%
  filter(!is.na(birth_deliverybirthcard_mode_binary)) %>%
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12"), ordered = T)) %>%
  ggplot(aes(Timepoint_new, Bacteroides, fill = birth_deliverybirthcard_mode_binary, color = birth_deliverybirthcard_mode_binary)) +
  ggrastr::rasterise(geom_jitter(position = position_jitterdodge(), size=0.1, alpha=1), dpi = 300)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color = "black") +
  theme_minimal() +
  theme(strip.background = element_rect(NA),
        legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm")) + 
  scale_fill_manual(values = met.brewer("Morgenstern")[c(2,5)], labels = c("C-section", "Vaginal delivery")) +
  scale_color_manual(values = met.brewer("Morgenstern")[c(2,5)], labels = c("C-section", "Vaginal delivery")) +
  labs(y = "Bacteroides abundance", x="Timepoint", fill = "Delivery mode", color = "Delivery mode")

ggsave("05.PLOTS/07.Health_outcomes/Bacteroides_genus_abundance_DM.pdf",
       Bacteroides_DM, "pdf", width=7, height=7, units="cm", dpi = 300)


bacshannon_DM <- gamma %>%
  right_join(smeta_w_phenos %>% select(Universal_ID, infant_ffq_feeding_mode_simple )) %>%
  filter(!is.na(infant_ffq_feeding_mode_simple) & Timepoint_new != "M12") %>%
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12"), ordered = T)) %>%
  ggplot(aes(Timepoint_new, bacShannon, fill = infant_ffq_feeding_mode_simple, color = infant_ffq_feeding_mode_simple)) +
  ggrastr::rasterise(geom_jitter(position = position_jitterdodge(), size=0.1, alpha=1), dpi = 300)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color = "black") +
  theme_minimal() +
  theme(strip.background = element_rect(NA),
        legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm")) + 
  scale_fill_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  scale_color_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  labs(y = "Bacterial alpha diversity", x="Timepoint", fill = "Feeding mode", color = "Feeding mode")

ggsave("05.PLOTS/07.Health_outcomes/BacShannon_FM.pdf",
       bacshannon_DM, "pdf", width=6, height=7, units="cm", dpi = 300)

summary(lmer(bacShannon ~ Timepoint_new + infant_ffq_feeding_mode_simple + (1|NEXT_ID), data = long_smeta))

length(unique(long_smeta$Sequencing_ID[!is.na(long_smeta$infant_ffq_feeding_mode_simple)]))

