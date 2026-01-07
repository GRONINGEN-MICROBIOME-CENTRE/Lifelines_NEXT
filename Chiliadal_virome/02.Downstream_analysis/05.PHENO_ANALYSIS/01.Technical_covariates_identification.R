setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

#############################################################
# Identification of intercorrelation between virome metrics
# Identification of technical covariates for the downstream
# analyses
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
library(purrr)
library(tidyverse)
library(broom)
library(skimr)
library(ggplot2)
#############################################################
# 2. Load Input Data
#############################################################
metadata <- read.delim('06.CLEAN_DATA/Intermediate/esmeta_12Dec25_enriched_try.txt', sep='\t', header=T, check.names = F)

metadata <- metadata %>%
  filter(Type %in% c('M', 'K')) %>%
  select(-temperate_RAb_ch) %>%
  rename("inv_ver" = "invertebrates, vertebrates")

# creating new Timepoint factor:
metadata$Timepoint_new <- factor(metadata$Timepoint_new, levels = c('M1', 'M3', 'M6', 'M12', 'Mother'), ordered = T)

tmp <- metadata %>%
  select(Universal_ID, bacterial_markers_alignment_rate, seq_type) %>%
  pivot_wider(id_cols = Universal_ID,
              names_from = seq_type,
              values_from = bacterial_markers_alignment_rate) %>%
  mutate(sc_enrichment_vlp = pmin(MGS / VLP, 100)) %>%
  mutate(sc_enrichment_mgs = pmin(MGS / median(MGS), 100))

metadata$sc_enrichment <- NA
metadata$sc_enrichment[metadata$seq_type == "VLP"] <- tmp$sc_enrichment_vlp[match(metadata$Universal_ID[metadata$seq_type == "VLP"], tmp$Universal_ID)]
metadata$sc_enrichment[metadata$seq_type == "MGS"] <- tmp$sc_enrichment_vlp[match(metadata$Universal_ID[metadata$seq_type == "MGS"], tmp$Universal_ID)]

# virome technical features:
tech_features <- c("isolation_method", "isolation_batch", "sequencing_batch", 
                   "dna_conc", "raw_reads", "human_reads",
                   "clean_reads", "reads_lost_QC", "bacterial_markers_alignment_rate",
                   "contigs_0_bp", "contigs_1000_bp", 
                   #"perc_aligned_d",
                   "perc_aligned_df", "perc_aligned_c", "perc_aligned_cf",
                   "sc_enrichment")

RPKM <- read.table("06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ab3kbp_1110_samples.txt", sep='\t', header=T)

#############################################################
# 3.1 Analysis: identify inter-correlation between virome metrices:
#############################################################

# VLP
virmetrics_vlp <- metadata %>%
  filter(seq_type == "VLP") %>%
  select(vir_richness_d:vertebrates)

metrics_vlp_rcorr <- Hmisc::rcorr(as.matrix(virmetrics_vlp))

png("05.PLOTS/07.Health_outcomes/VLP_virmetrics_intercorrelation.png", width = 27, height = 27, units = "cm", res = 300)
corrplot::corrplot(corr = metrics_vlp_rcorr[[1]],
                   p.mat = metrics_vlp_rcorr[[3]],
                   sig.level = 0.05,
                   insig='blank',
                   order = 'hclust')
dev.off()

# MGS
virmetrics_mgs <- metadata %>%
  filter(seq_type == "MGS") %>%
  select(vir_richness_d:vertebrates)

metrics_mgs_rcorr <- Hmisc::rcorr(as.matrix(virmetrics_mgs))

png("05.PLOTS/07.Health_outcomes/MGS_virmetrics_intercorrelation.png", width = 27, height = 27, units = "cm", res = 300)
corrplot::corrplot(metrics_mgs_cor)
dev.off()

#############################################################
# 3.2 Analysis: intercorrelation of technical features
#############################################################

numeric_tf <- metadata %>%
  filter(seq_type == "VLP") %>%
  select(all_of(tech_features)) %>%
  select(where(~ is.numeric(.x) && any(!is.na(.x))))

features_vlp_rcorr <- Hmisc::rcorr(as.matrix(numeric_tf))

png("05.PLOTS/07.Health_outcomes/VLP_tech_features_intercorrelation.png", width = 15, height = 15, units = "cm", res = 300)
corrplot::corrplot(corr = features_vlp_rcorr[[1]],
                   p.mat = features_vlp_rcorr[[3]],
                   sig.level = 0.05,
                   insig='blank',
                   order = 'hclust')
dev.off()

#############################################################
# 3.3 Analysis: explore how technical characteristics explain the variation of virome metrics:
#############################################################

virmetrics <- colnames(virmetrics_vlp)

results <- map_dfr(virmetrics, function(virvar) {
  
  map_dfr(tech_features, function(tech_feat) {
    
    formula <- as.formula(paste(virvar, "~ exact_age_days_at_collection + perc_aligned_d +", tech_feat))
    model <- lm(formula, data = metadata[metadata$seq_type == "VLP",])
    
    anova_res <- car::Anova(model, type = "II") %>%
      tidy() %>%
      filter(term == tech_feat)
    
    eff_res <- effectsize::eta_squared(model, partial = F) %>%
      as.data.frame() %>%
      filter(Parameter == tech_feat)
    
    data.frame(
      virvar_VLP = virvar,
      tech_feature = tech_feat,
      p_value = anova_res$p.value,
      eta_sq = eff_res$Eta2
    )
    
  })
})

res_plot <- results %>%
  mutate(
    neglog10_p = -log10(p_value + min(results[results$p_value!=0,]$p_value))
  ) %>%
  mutate(neglog10_p = ifelse(p_value >= 0.05, NA, neglog10_p))
  

png("05.PLOTS/07.Health_outcomes/VLP_virmetrics_tech_features_timepoint_aligned_reads.png", width = 17, height = 25, units = "cm", res = 300)
ggplot(res_plot, aes(x = tech_feature, y = virvar_VLP)) +
  geom_tile(color = "grey", aes(fill = eta_sq)) +
  geom_point(aes(size = neglog10_p), shape = 21, colour = "black", alpha = 1, fill = NA) +
  scale_fill_gradient2(name = "Eta sq",
                       low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  scale_size_continuous(name = "-log10(p)", range = c(0, 6)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
dev.off() 

# list to correct for besides age and percentage of aligned reads:
pupupu <- results[results$p_value < 0.05 & results$eta_sq > 0.1,]
pupupu <- pupupu[pupupu$virvar_VLP %in% c("vir_richness_cf", "vir_diversity", "temperate_richness"),]

#############################################################
# 3.4 Analysis: explore how technical characteristics explain the variation of virome composition:
#############################################################
dist <- vegan::vegdist(t(RPKM), method = "bray")
dm <- as.matrix(dist)


ad2_results <- map_dfr(tech_features, function(tech_feat) {
    
    formula <- as.formula(paste("dist ~ Timepoint_new +", tech_feat))
    model <- vegan::adonis2(formula, 
                            data=metadata[metadata$seq_type == "VLP",], 
                            permutations = 999)
    
    anova_res <- model %>%
      tidy() %>%
      filter(! term %in% c("Residual", "Total"))
    
    data.frame(
      tech_feature = tech_feat,
      p_value = anova_res$p.value,
      R2 = anova_res$R2
    )
    
  })

ad2_results <- ad2_results[order(ad2_results$R2),]

png("05.PLOTS/07.Health_outcomes/VLP_adonis2_tech_features_after_time.png", width = 13, height = 10, units = "cm", res = 300)
ggplot(ad2_results, aes(x = tech_feature, y = R2)) +
  geom_bar(stat = "identity", fill="#D25353") +
  scale_x_discrete(limits=ad2_results$tech_feature) + 
  coord_flip() +
  theme_bw()
dev.off()

#############################################################
# 4. OUTPUT
#############################################################
write.table(dm, "06.CLEAN_DATA/02.FINAL/BC_dist_VLP_only_RPKM_table_VLP_MGS_dec99ANI_ab3kbp_1110_samples.txt", sep='\t', quote=F)

