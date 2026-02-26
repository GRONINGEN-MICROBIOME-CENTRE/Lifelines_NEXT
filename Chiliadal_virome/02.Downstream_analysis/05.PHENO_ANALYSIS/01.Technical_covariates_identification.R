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
library(vegan)
#############################################################
# 2. Load Input Data
#############################################################
metadata <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T, check.names = F)

metadata <- metadata %>%
  rename("inv_ver" = "invertebrates, vertebrates") %>%
  mutate(Timepoint_new = factor(Timepoint_new, levels = c('M1', 'M3', 'M6', 'M12', 'Mother'), ordered = T))

# virome technical features:
tech_features <- c("isolation_method", "isolation_batch", "sequencing_batch", 
                   "dna_conc", "raw_reads", "human_reads",
                   "clean_reads", "reads_lost_QC", "bacterial_markers_alignment_rate",
                   "contigs_0_bp", "contigs_1000_bp", 
                   "perc_aligned_d", "perc_aligned_cf",
                   "N_shared_cs_d", "N_shared_cs_cf",
                   "sc_enrichment")

VLP <- read.table("06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt", sep='\t', header=T)

#############################################################
# 3.1 Analysis: identify inter-correlation between virome metrics:
#############################################################

# VLP
virmetrics_vlp <- metadata %>%
  filter(seq_type == "VLP" & Timepoint_new == "Mother") %>%
  select(vir_richness_d:ssDNA) %>%
  select(-c(sc_enrichment, genome_unclassified, N_shared_cs_d, N_shared_cs_cf))

metrics_vlp_rcorr <- Hmisc::rcorr(as.matrix(virmetrics_vlp))

png("05.PLOTS/07.Health_outcomes/Clean_VLP_virmetrics_intercorrelation.png", width = 27, height = 27, units = "cm", res = 300)
corrplot::corrplot(corr = metrics_vlp_rcorr[[1]],
                   p.mat = metrics_vlp_rcorr[[3]],
                   sig.level = 0.05,
                   insig='blank')
dev.off()

# MGS
virmetrics_mgs <- metadata %>%
  filter(seq_type == "MGS" & Type == "K") %>%
  filter(vir_richness_cf !=0) %>%
  select(vir_richness_d:ssDNA) %>%
  select(-all_of(c("PPV_fraction", "PPV_abundance", "sc_enrichment", "genome_unclassified")))

metrics_mgs_rcorr <- Hmisc::rcorr(as.matrix(virmetrics_mgs))


png("05.PLOTS/07.Health_outcomes/Clean_MGS_virmetrics_intercorrelation.png", width = 27, height = 27, units = "cm", res = 300)
corrplot::corrplot(corr = metrics_mgs_rcorr[[1]],
                   p.mat = metrics_mgs_rcorr[[3]],
                   sig.level = 0.05,
                   insig='blank')
dev.off()

#############################################################
# 3.2 Analysis: intercorrelation of technical features
#############################################################

numeric_tf <- metadata %>%
  filter(seq_type == "VLP") %>%
  select(all_of(tech_features)) %>%
  select(where(~ is.numeric(.x) && any(!is.na(.x))))

features_vlp_rcorr <- Hmisc::rcorr(as.matrix(numeric_tf))

png("05.PLOTS/07.Health_outcomes/Clean_VLP_tech_features_intercorrelation.png", width = 15, height = 15, units = "cm", res = 300)
corrplot::corrplot(corr = features_vlp_rcorr[[1]],
                   p.mat = features_vlp_rcorr[[3]],
                   sig.level = 0.05,
                   insig='blank')
dev.off()

#############################################################
# 3.3 Analysis: explore how technical covariates for virome metrics:
#############################################################
virmetrics <- colnames(virmetrics_vlp)

tech_cov_vlp <- map_dfr(unique(metadata$Timepoint_new), function(timepoint){
  
  map_dfr(tech_features, function(technica) {
    
    map_dfr(virmetrics, function(metrica){
      
      print(paste(timepoint, technica, metrica))
      
      formula <- as.formula(paste(metrica, "~", technica))
      
      dat <- metadata %>%
        filter(Timepoint_new == timepoint & seq_type == "VLP" & !is.na(!!sym(technica)) & !is.na(!!sym(metrica)))
      
      n_contrast <- dat %>%
        pull(technica) %>%
        unique() %>%
        length()
      
      if (n_contrast < 2) return(NULL)
      
      model <- lm(
        formula,
        data = metadata[metadata$Timepoint_new == timepoint,]
      )
      mod_anova <- anova(model)
      
      data.frame(at_tp = timepoint,
                 metrics_tested = metrica, 
                 technica_tested = technica,
                 RSq_adj = summary(model)$adj.r.squared,
                 p_value = mod_anova[["Pr(>F)"]][1])
      
    })
  })
  
})

tech_cov_vlp <- tech_cov_vlp %>%
  mutate(p_adj = p.adjust(p_value, "BH")) %>%
  arrange(metrics_tested, technica_tested, at_tp)

unique(tech_cov_vlp$metrics_tested)

saver_tech_cov_vlp <- tech_cov_vlp
writexl::write_xlsx(saver_tech_cov_vlp, "07.RESULTS/Timepoint-wise_VLP_lm_covariates_vs_metrics.xlsx")

summary_tech_cov_vlp <- tech_cov_vlp %>%
  mutate(p_adj = p.adjust(p_value, "BH")) %>%
  filter(p_adj < 0.05 & RSq_adj >= 0.1 & at_tp != "Mother") %>%
  group_by(metrics_tested, technica_tested) %>%
  summarize( times_significant = n(),
             avg_rsq = mean(RSq_adj),
             timepoints_active = paste(unique(at_tp), collapse = ", "),
             .groups = 'drop') %>%
  filter(times_significant >= 3)

tech_cov_vlp_plot <- tech_cov_vlp %>%
  filter(!technica_tested %in% c('dna_conc', 'human_reads','isolation_batch','
                                 isolation_method','raw_reads', 'N_shared_cs_d', 
                                 'contigs_0_bp', 'perc_aligned_d',
                                 'bacterial_markers_alignment_rate',
                                 'isolation_method')) %>% # based on summary_tech_cov_vlp results, prev visualization and\r colinearity
  filter(!metrics_tested %in% c('vir_richness_d', 'PPV_fraction', 'RNA')) %>%
  filter(p_adj < 0.05 & RSq_adj >= 0.1 & at_tp != "Mother") %>%
  mutate(
    neglog10_p = -log10(p_adj + min(tech_cov_vlp[tech_cov_vlp$p_adj!=0,]$p_adj))
  ) %>%
  mutate(technica_tested = fct_recode(technica_tested,
                                      "N contigs > 1kbp" = "contigs_1000_bp",
                                      "N shared to CTRLs" = "N_shared_cs_cf",
                                      "Reads aligned, %"  = "perc_aligned_cf",
                                      "Reads lost at QC, %"    = "reads_lost_QC",
                                      "Enrichment efficacy"          = "sc_enrichment",
                                      "Sequencing batch"          = "sequencing_batch")) %>%
  mutate(metrics_tested = fct_recode(metrics_tested,
                                      "Virus richness" = "vir_richness_cf",
                                      "Virus diversity" = "vir_diversity",
                                      "Temperate richness" = "temperate_richness",
                                      "Temperate abundance" = "temperate_RAb",
                                      "Temperate fraction" = "temperate_perc",
                                      "Temperate:lytic ratio" = "temp_to_lytic_ratio",
                                      "ssDNA abundance" = "ssDNA",
                                      "Lytic richness" = "lytic_richness",
                                      "Lytic abundance" = "lytic_RAb",
                                      "Lytic fraction" = "lytic_perc",
                                      "dsDNA abundance" = "dsDNA")) %>%
ggplot(aes(x = at_tp, y = metrics_tested)) +
  geom_tile(aes(fill = RSq_adj), color = "white") +
  facet_wrap(~technica_tested, ncol = 10) +
  geom_point(aes(size = neglog10_p), shape = 21, colour = "black", alpha = 1, fill = NA) +
  scale_fill_gradient2(name = "Eta sq",
                       low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  scale_size_continuous(name = "-log10(p)", range = c(0, 6)) + 
  theme_bw() +
  theme(strip.background = element_rect(NA),
        legend.position = "bottom",
        axis.text = element_text(size=7),
        strip.text = element_text(size = 8),
        axis.title = element_text(size=8)) + 
  labs(y = "Virome metric", x="Infant timepoint")

ggsave("05.PLOTS/07.Health_outcomes/Timepoint-wise_tech_virmetrics_correlation.pdf", tech_cov_vlp_plot, 
       "pdf", width=21, height=12, units="cm", dpi = 300)

#############################################################
# 3.4 Analysis: explore how technical characteristics explain the variation of virome composition:
#############################################################
infant_VLP <- VLP[,colnames(VLP) %in% VLP_meta$Sequencing_ID[VLP_meta$Type == "K"]]
infant_VLP <- infant_VLP[rowSums(infant_VLP)>0,]

VLP_rCLR <- as.data.frame( decostand(t(infant_VLP), method = "rclr"))
VLP_ait_dist <- vegdist(VLP_rCLR, method = "euclidean", parallel = 6)

dm <- as.data.frame(as.matrix(VLP_ait_dist))

VLP_meta <- metadata[metadata$seq_type == "VLP" & metadata$Type == "K",]
VLP_meta <- VLP_meta[match(colnames(dm), VLP_meta$Sequencing_ID),]

ad2_results <- map_dfr(tech_features, function(tech_feat) {
    
    formula <- as.formula(paste("VLP_ait_dist ~ Timepoint_new +", tech_feat))
    model <- vegan::adonis2(formula, 
                            data=VLP_meta,
                            strata = VLP_meta$NEXT_ID,
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

ad_inf <- ggplot(ad2_results, aes(x = tech_feature, y = R2)) +
  geom_bar(stat = "identity", fill="#D25353") +
  scale_x_discrete(limits=ad2_results$tech_feature) + 
  coord_flip() +
  theme_bw() +
  theme(axis.text = element_text(size=7),
        strip.text = element_text(size = 8),
        axis.title = element_text(size=8)) + 
  labs(y = "Technical feature", x="Adonis R-sqaured")

ggsave("05.PLOTS/07.Health_outcomes/VLP_adonis2_tech_features_after_time.pdf", ad_inf, 
       "pdf", width=10, height=7, units="cm", dpi = 300)

blob <- vegan::betadisper(VLP_ait_dist, group = VLP_meta$isolation_batch)

anova(blob)
summary(blob$group.distances)
#############################################################
# 4. OUTPUT
#############################################################

