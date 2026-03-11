setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# MGS recapitulation of VLP composition
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
# virome_recapitulation(MGS, VLP, smeta) quantifies vOTU 
# overlap between MGS and VLP metaviromes from the same
# fecal sample (denoted by Universal_ID)
#
# INPUTS:
# - MGS & VLP: RPKM tables (vOTUs = rows, samples = columns).
# - smeta: metadata containing Sequencing_ID, Universal_ID, 
# and seq_type
#
# !! PAY ATTENTION !! :
# 1. FIXED NAMING: this function is very rigid for variable names
# 2. OUTPUT STANDARDIZATION: Regardless of input order, 
# 'sample_1' perceived as MGS ID and 'sample_2' - VLP ID.
virome_recapitulation <- function(MGS, VLP, smeta){
  
  DF_common <- merge(MGS, VLP, by = "row.names", all = T)
  row.names(DF_common) <- DF_common$Row.names
  DF_common$Row.names <- NULL
  DF_common[is.na(DF_common)] <- 0
  
  binary_matrix <- as.matrix(DF_common > 0) + 0
  
  shared_species_matrix <-  t(binary_matrix) %*% binary_matrix
  
  upper_tri_idx <- which(upper.tri(shared_species_matrix, diag = FALSE), arr.ind = TRUE)
  
  smeta$new_richness <- colSums(DF_common >0)[match(smeta$Sequencing_ID, colnames(DF_common))] # richness for MGS and VLP samples when respective counter-method is depleted
  
  smeta$new_temperate_richness <- colSums(DF_common[row.names(DF_common) %in% temperate,] >0)[match(smeta$Sequencing_ID, colnames(DF_common))]
  smeta$new_virulent_richness <- colSums(DF_common[row.names(DF_common) %in% virulent,] >0)[match(smeta$Sequencing_ID, colnames(DF_common))]
  smeta$new_ssDNA_richness <- colSums(DF_common[row.names(DF_common) %in% ssDNA,] >0)[match(smeta$Sequencing_ID, colnames(DF_common))]
  smeta$new_RNA_richness <- colSums(DF_common[row.names(DF_common) %in% RNA,] >0)[match(smeta$Sequencing_ID, colnames(DF_common))]
  smeta$new_dsDNA_richness <- colSums(DF_common[row.names(DF_common) %in% dsDNA,] >0)[match(smeta$Sequencing_ID, colnames(DF_common))]
  
  # actual calculation
  sharing <- data.frame(
    sample_1 = rownames(shared_species_matrix)[upper_tri_idx[, 1]],
    sample_2 = colnames(shared_species_matrix)[upper_tri_idx[, 2]],
    N_shared = shared_species_matrix[upper_tri_idx]) %>%
    left_join(smeta %>% select(Sequencing_ID,
                               Universal_ID,
                               seq_type), by = c("sample_1" = "Sequencing_ID" )) %>%
    left_join(smeta %>% select(Sequencing_ID,
                               Universal_ID,
                               seq_type), by = c("sample_2" = "Sequencing_ID"), suffix = c("_1", "_2")) %>%
    mutate(same_feces = ifelse(Universal_ID_1 == Universal_ID_2, T, F)) %>%
    filter(same_feces == T) %>%
    mutate(
      temp_s1 = sample_1,
      sample_1 = ifelse(seq_type_1 != "MGS", sample_2, sample_1),
      sample_2 = ifelse(seq_type_1 != "MGS", temp_s1, sample_2)
    ) %>%
    select(-temp_s1) %>%
    select(sample_1, sample_2, N_shared, same_feces) %>%
    left_join(smeta %>% select(Sequencing_ID, new_richness, 
                               new_temperate_richness, 
                               new_virulent_richness, 
                               new_ssDNA_richness,
                               new_RNA_richness, 
                               new_dsDNA_richness,
                               Timepoint_new, Type) %>% rename_with(~paste0(., "_VLP"), -Sequencing_ID), by = c("sample_2" = "Sequencing_ID")) %>%
    mutate(perc_VLP_in_MGS = N_shared / new_richness_VLP * 100) %>%
    left_join(smeta %>% select(Sequencing_ID, new_richness) %>% rename_with(~paste0(., "_MGS"), -Sequencing_ID), by = c("sample_1" = "Sequencing_ID")) %>%
    mutate(perc_MGS_in_VLP = N_shared / new_richness_MGS * 100)
  
  return(sharing)
  
}
#############################################################
# 1. Loading libraries
#############################################################

library(dplyr)
library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
#############################################################
# 2. Load Input Data
#############################################################

esmeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)
esmeta$Timepoint_new <- factor(esmeta$Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)

esmeta <- esmeta %>%
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) # making the 2nd pregnancy as separate fams

ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

temperate <- ETOF_vOTUr$New_CID[ETOF_vOTUr$lifestyle == "Temperate"]
virulent <- ETOF_vOTUr$New_CID[ETOF_vOTUr$lifestyle == "Virulent"]
ssDNA <- ETOF_vOTUr$New_CID[ETOF_vOTUr$genome == "ssDNA"]
RNA <- ETOF_vOTUr$New_CID[ETOF_vOTUr$genome == "RNA"]
dsDNA <- ETOF_vOTUr$New_CID[ETOF_vOTUr$genome == "dsDNA"]

# bacteriome:
metaphlan <- read.table('06.CLEAN_DATA/02.FINAL/MGS_Chiliadal_metaphlan_full_taxonomy_ver_01_07102025.txt', sep='\t', header=T)
metaphlan <- metaphlan[grep('t__SGB', row.names(metaphlan)),] # SGB-level
metaphlan <- metaphlan[rowSums(metaphlan) > 0,]

# MGS virome:
MGS <- read.table('06.CLEAN_DATA/02.FINAL/MGS_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)

MGS_all <- MGS[row.names(MGS) %in% ETOF_vOTUr$New_CID[ETOF_vOTUr$vOTU_cluster_type != "NEXT_VLP"],] # removing all vOTUs that would not be there without VLP

MGS_all_RA <- as.data.frame(t(MGS_all)/colSums(MGS_all)*100) # relative abundance

# VLP virome:
VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)

VLP_all <- VLP[row.names(VLP) %in% ETOF_vOTUr$New_CID[ETOF_vOTUr$vOTU_cluster_type != "NEXT_MGS"],] # removing all vOTUs that would not be there without MGS

VLP_all_RA <- as.data.frame(t(VLP_all)/colSums(VLP_all)*100)

#############################################################
# 3.1 Analysis: VLP virome recapitulation
#############################################################
sharing_all <- virome_recapitulation(MGS_all, VLP_all, esmeta)

# check the conclusion of the whole section

sharing_all_UPD <- sharing_all %>%
  mutate(new_perc_tmeperate = new_temperate_richness_VLP / new_richness_VLP*100,
         new_perc_lytic = new_virulent_richness_VLP / new_richness_VLP*100,
         lytic_index = new_virulent_richness_VLP / (new_temperate_richness_VLP +1),
         ssDNAtodsDNA = new_ssDNA_richness_VLP / new_dsDNA_richness_VLP,
         RNAtodsDNA = new_RNA_richness_VLP / new_dsDNA_richness_VLP,
         lytic_index_scaled = scale(lytic_index),
         ssDNAtodsDNA_scaled = scale(ssDNAtodsDNA),
         RNAtodsDNA_scaled = scale(RNAtodsDNA)) %>%
  left_join(esmeta %>% select(Sequencing_ID, NEXT_ID), by = c("sample_2" = "Sequencing_ID"))


stat_summary <-  map_dfr(c("perc_VLP_in_MGS", "perc_MGS_in_VLP"), function(perc){
    
   sharing_all_UPD %>%
      summarise(
        fraction = perc,
        mean = mean(!!sym(perc), na.rm = TRUE),
        sd   = sd(!!sym(perc), na.rm = TRUE),
        median = median(!!sym(perc), na.rm = TRUE)
      )
    
  })

stat_mm <- map_dfr(c("lytic_index_scaled", 
                     "ssDNAtodsDNA_scaled",
                     "RNAtodsDNA_scaled"), function(lfs_metric){
  
  formula <- as.formula(paste0("perc_VLP_in_MGS ~ ", lfs_metric, "+ Timepoint_new_VLP + (1|NEXT_ID)"))
  
  model <- lmer(formula,  REML = FALSE, data = sharing_all_UPD)
  
    summary(model)$coefficients %>%
    as.data.frame() %>%
      rownames_to_column(var = "rowname") %>%
      filter(rowname == lfs_metric) 
})

stat_mm$FDR <- p.adjust(stat_mm$`Pr(>|t|)`, "BH")

# distribution of VLP recapitulation and correlation to index
xAxisBoxPlot_v <- ggplot(sharing_all_UPD, aes(x="", y = perc_VLP_in_MGS)) +
  geom_boxplot(fill = "#132440", alpha = 0.6) + 
  theme_void() +
  coord_flip()

p_lfs_index <- ggplot(sharing_all_UPD, aes(x = perc_VLP_in_MGS, y = lytic_index)) +
  ggrastr::rasterise(geom_point(alpha = 0.8, size = 1, color = "#132440"), dpi = 300)+
  geom_smooth(method = "lm", se = TRUE, size = 1.2, color = "#BF092F", fill = "#BF092F") +
  theme_minimal() +
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=9)) +
  labs(x = "% VLP-detected vOTUs recapitulated in MGS", 
       y = "Virulent-to-temperate phage ratio") + 
  annotate(geom = "text", x = 70, y = 2.5, label = "beta = -6.6\np-value = 3.9e-39", size = 2.5)

p_lfs_index_w <- xAxisBoxPlot_v / p_lfs_index + 
  plot_layout(heights = c(1,8))

ggsave('05.PLOTS/05.VLP_MGS/VLP_in_MGS_detection_vs_lytic_index.pdf',
       p_lfs_index_w,  "pdf", width=7.5, height=8.5 , units="cm", dpi = 300)


efsize <- ggplot(stat_mm, aes(x = reorder(rowname, Estimate), y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick") +
  geom_pointrange(aes(ymin = Estimate - 1.96*`Std. Error`, ymax = Estimate + 1.96*`Std. Error`), size = 0.5, color = "#132440") +
  scale_x_discrete(labels = c("Virulent-to-temperate vOTU ratio", "ssDNA-to-dsDNA vOTU ratio", "RNA-to-dsDNA vOTU ratio")) +
  coord_flip() +
  labs(x = "VLP virome feature", y = "Standardized effect size") +
  theme_minimal() +
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=9))

ggsave('05.PLOTS/05.VLP_MGS/VLP_in_MGS_detection_drivers.pdf',
       efsize,  "pdf", width=10, height=3.5, units="cm", dpi = 300)

#############################################################
# 3.2 Analysis: virome independence from bacteriome
#############################################################
# mantel test
# ordering:

colnames(metaphlan) <- esmeta$Universal_ID[match(colnames(metaphlan), esmeta$Sequencing_ID)]
metaphlan <- t(metaphlan)

row.names(VLP_all_RA) <- esmeta$Universal_ID[match(row.names(VLP_all_RA), esmeta$Sequencing_ID)]

row.names(MGS_all_RA) <- esmeta$Universal_ID[match(row.names(MGS_all_RA), esmeta$Sequencing_ID)]

VLP_all_RA <- VLP_all_RA[row.names(metaphlan),]
MGS_all_RA <- MGS_all_RA[row.names(metaphlan),]

dist_bac <- vegdist(metaphlan, method="bray")
dist_vlp <- vegdist(VLP_all_RA, method="bray")
dist_mgs <- vegdist(MGS_all_RA, method="bray", na.rm = T)

mantel_vlp <- mantel(dist_bac, dist_vlp, method="spearman", permutations=999, na.rm = T)
mantel_mgs <- mantel(dist_bac, dist_mgs, method="spearman", permutations=999, na.rm = T)



#### visualization:
library(ggplot2)

df_summary <- tibble(
  Bact = dist_bac[upper.tri(dist_bac)],
  MGS  = dist_mgs[upper.tri(dist_mgs)],
  VLP = dist_vlp[upper.tri(dist_vlp)]
)

# Calculate SD instead of SE for a more visible spread
plot_data <- df_summary %>%
  mutate(Bin = cut(Bact, breaks = seq(0, 1, by = 0.02), include.lowest = TRUE)) %>%
  pivot_longer(cols = c(MGS, VLP), names_to = "Fraction", values_to = "Viral_Dist") %>%
  group_by(Bin, Fraction) %>%
  summarise(
    Mean_Viral = mean(Viral_Dist, na.rm = TRUE),
    SD_Viral   = sd(Viral_Dist, na.rm = TRUE), # Standard Deviation shows the "cloud"
    .groups = 'drop'
  ) %>%
  filter(!is.na(Bin)) %>%
  mutate(Bact_Midpoint = as.numeric(sub(".*,([0-9.]*)\\]", "\\1", Bin)) - 0.01)

# Plot with SD ribbon
plotik <- ggplot(plot_data, aes(x = Bact_Midpoint, y = Mean_Viral, color = Fraction, fill = Fraction)) +
  # Use SD for the ribbon to see the actual range of the data
  geom_ribbon(aes(ymin = Mean_Viral - SD_Viral, ymax = Mean_Viral + SD_Viral), 
              alpha = 0.15, color = NA) + 
  geom_line(size = 1.2) +
  scale_color_manual(values = c("MGS" = "#E64B35FF", "VLP" = "#4DBBD5FF")) +
  scale_fill_manual(values = c("MGS" = "#E64B35FF", "VLP" = "#4DBBD5FF")) +
  #coord_cartesian(xlim = c(0.8, 1.0), ylim = c(0.8, 1.01)) + # Expanded Y-axis to see the ribbon
  labs(x = "Bacteriome dissimilarity", y = "Virome dissimilarity", color = "Metavirome type", fill = "Metavirome type") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text = element_text(size=7),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=7),
        legend.key.size=unit(0.7, "line"),
        legend.key.spacing.y = unit(1, 'pt')) +
  annotate("rect", xmin = 0.75, xmax = 0.95, ymin = 0.45, ymax = 0.6, 
           alpha = 0.2, fill = "#4DBBD5FF") +
  annotate(geom="text", x = 0.85, y=0.53, label = "r = 0.24\np-value < 0.001", size = 2.5) +
  annotate("rect", xmin = 0.75, xmax = 0.95, ymin = 0.3, ymax = 0.45, 
           alpha = 0.2, fill = "#E64B35FF") +
  annotate(geom="text", x = 0.85, y=0.38, label = "r = 0.88\np-value < 0.001", size = 2.5)
  
  

ggsave('05.PLOTS/05.VLP_MGS/mantel_results.pdf',
       plotik,  "pdf", width=13, height=10, units="cm", dpi = 300)
