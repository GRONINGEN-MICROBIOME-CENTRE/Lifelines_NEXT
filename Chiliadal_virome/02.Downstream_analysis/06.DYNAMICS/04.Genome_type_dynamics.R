setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# dynamics by genome
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
# original source: 01.Virome_metrics_dynamics.R
# primarily to calculate the difference between infant timepoints and maternal samples
mixed_model_tukey <- function(df, formula, variable){
  
  model <- lmer(as.formula(formula), data = df)
  
  formula2 <- paste0("~ ", variable)
  
  emm <- emmeans(model, as.formula(formula2))
  
  as.data.frame(contrast(emm, method = "pairwise", adjust = "tukey")) %>%
    mutate(tested_var = gsub(" ~.*", "", formula))
  
}
#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MetBrewer)
library(emmeans)
library(patchwork)
#############################################################
# 2. Load Input Data
#############################################################
# metadatas:
smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)

smeta <- smeta %>%
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) %>% # treating 2nd pregnancy as a separate family:
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T))

# abundance tables etc:
VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)

ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header = T)
ETOF_vOTUr <- ETOF_vOTUr %>%
  mutate(Host = ifelse(Host == "invertebrates, vertebrates", "invertebrates_vertebrates", Host))

#############################################################
# 3.0 Analysis: richness by genome type
#############################################################
richness_by_genome <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, genome)) %>%
  group_by(genome) %>%
  summarise(across(where(is.numeric), ~sum(. > 0, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(!genome) %>%
  pivot_wider(names_from = genome, values_from = value) %>%
  left_join(smeta %>% select(Sequencing_ID, Timepoint_new, NEXT_ID, vir_richness_cf, vir_diversity, bacShannon), by = c("name" = "Sequencing_ID")) %>%
  mutate(RNAtodsDNA = log((RNA + 1)/(dsDNA + 1)),
         ssDNAtodsDNA = log((ssDNA + 1)/(dsDNA + 1)))

genome_type_richness <- rbind(mixed_model_tukey(richness_by_genome, 
                                                "RNAtodsDNA ~ Timepoint_new + (1|NEXT_ID)",
                                                "Timepoint_new"),
                              mixed_model_tukey(richness_by_genome, 
                                                "ssDNAtodsDNA ~ Timepoint_new + (1|NEXT_ID)",
                                                "Timepoint_new"))

coda_genome_richness <- map_dfr(c('ssDNAtodsDNA', 'RNAtodsDNA'), function(prop_type) {
    
    formula <- as.formula(paste0(prop_type, " ~ ", " Timepoint_new + (1|NEXT_ID)"))
    
    model <- lmer(
      formula,
      REML = FALSE,
      data = richness_by_genome[richness_by_genome$Timepoint_new != "Mother",]
    )
    
    summary(model)$coefficients %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      filter(!rowname %in% c("(Intercept)", "Timepoint_new.Q", "Timepoint_new.C")) %>%
      mutate(genome = prop_type)
  })

coda_genome_richness$p_adjust <- p.adjust(coda_genome_richness$`Pr(>|t|)`, "BH")

rich_plot_data <- richness_by_genome %>%
  select(name, ssDNA, RNA, dsDNA, Timepoint_new, vir_richness_cf) %>%
  pivot_longer(cols = c('ssDNA', 'RNA', 'dsDNA'), names_to = "genome") %>%
  mutate(prop_rich = value/vir_richness_cf)

mom_means <- rich_plot_data %>%
  filter(Timepoint_new == "Mother") %>%
  group_by(genome) %>%
  summarise(mean_value = mean(prop_rich), .groups = "drop") %>%
  mutate(Type = "Mother")

inf_means <- rich_plot_data %>%
  filter(Timepoint_new != "Mother") %>%
  group_by(genome, Timepoint_new) %>%
  summarise(mean_value = mean(prop_rich), .groups = "drop")


rich_plot <- rich_plot_data %>%
  mutate(Type = "Infant") %>%
  filter(Timepoint_new != "Mother" & genome != "Unclassified") %>%
  ggplot(aes(Timepoint_new, prop_rich, group = genome, color = genome, fill = genome, linetype = Type)) +
  geom_hline(data = mom_means, aes(yintercept = mean_value, color = genome, linetype = Type), lwd = 0.5) +
  stat_summary(fun.data = "mean_se", geom = "ribbon", alpha = 0.2, color = NA) +
  stat_summary(fun = mean, geom = "line", lwd = 0.5) +
  #stat_summary(fun.data = "median_hilow", geom = "ribbon", alpha = 0.2, color = NA) +
  #stat_summary(fun = median, geom = "line", lwd = 1) +
  labs(x = "Timepoint", y = "Richness proportion", color = "Genome type", fill = "Genome type") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size=8)) +
  scale_color_manual(values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                                "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                                "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2])) +
  scale_fill_manual(values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                               "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                               "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2]))

#############################################################
# 3.0 Analysis: abundance by genome type
#############################################################

abundance_by_genome <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, genome)) %>%
  group_by(genome) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x, na.rm = TRUE))) %>%
  pivot_longer(!genome) %>%
  pivot_wider(names_from = genome, values_from = value) %>%
  left_join(smeta %>% select(Sequencing_ID, Timepoint_new, NEXT_ID, vir_richness_cf, vir_diversity, bacShannon), by = c("name" = "Sequencing_ID")) %>%
  mutate(RNAtodsDNA = log((RNA + 2.56e-06)/(dsDNA + 2.56e-06)), # adding pseudocount, I calculated it separately
         ssDNAtodsDNA = log((ssDNA + 2.56e-06)/(dsDNA + 2.56e-06)))


genome_type_abundance <- rbind(mixed_model_tukey(abundance_by_genome, 
                                                "RNAtodsDNA ~ Timepoint_new + (1|NEXT_ID)",
                                                "Timepoint_new"),
                              mixed_model_tukey(abundance_by_genome, 
                                                "ssDNAtodsDNA ~ Timepoint_new + (1|NEXT_ID)",
                                                "Timepoint_new"))

coda_genome_abundance <- map_dfr(c('ssDNAtodsDNA', 'RNAtodsDNA'), function(prop_type) {
  
  formula <- as.formula(paste0(prop_type, " ~ ", " Timepoint_new + (1|NEXT_ID)"))
  
  model <- lmer(
    formula,
    REML = FALSE,
    data = abundance_by_genome[abundance_by_genome$Timepoint_new != "Mother",]
  )
  
  summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(!rowname %in% c("(Intercept)", "Timepoint_new.Q", "Timepoint_new.C")) %>%
    mutate(genome = prop_type)
})

coda_genome_abundance$p_adjust <- p.adjust(coda_genome_abundance$`Pr(>|t|)`, "BH")

abu_plot_data <- abundance_by_genome %>%
  select(name, ssDNA, RNA, dsDNA, Timepoint_new, vir_richness_cf) %>%
  pivot_longer(cols = c('ssDNA', 'RNA', 'dsDNA'), names_to = "genome")

mom_means_abu <- abu_plot_data %>%
  filter(Timepoint_new == "Mother") %>%
  group_by(genome) %>%
  summarise(mean_value = mean(value), .groups = "drop") %>%
  mutate(Type = "Mother")

abu_plot <- abu_plot_data %>%
  mutate(Type = "Infant") %>%
  filter(Timepoint_new != "Mother" & genome != "Unclassified") %>%
  ggplot(aes(Timepoint_new, value, group = genome, color = genome, fill = genome, linetype = Type)) +
  geom_hline(data = mom_means_abu, aes(yintercept = mean_value, color = genome, linetype = Type), lwd = 0.5) +
  stat_summary(fun.data = "mean_se", geom = "ribbon", alpha = 0.2, color = NA) +
  stat_summary(fun = mean, geom = "line", lwd = 0.5) +
  #stat_summary(fun.data = "median_hilow", geom = "ribbon", alpha = 0.2, color = NA) +
  #stat_summary(fun = median, geom = "line", lwd = 1) +
  labs(x = "Timepoint", y = "Relative abundance", color = "Genome type", fill = "Genome type") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size=8)) +
  scale_color_manual(values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                                "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                                "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2])) +
  scale_fill_manual(values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                               "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                               "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2]))

# saving output:

## plot
better_together <- rich_plot + abu_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("05.PLOTS/06.DYNAMICS/VLP_richness_abundance_by_genome_lines.pdf",
       better_together, "pdf", width = 8, height = 9, units = "cm")

## dynamics in infants:
coda <- coda_genome_richness %>%
  mutate(Info_type = "Richness proportion") %>%
  bind_rows(coda_genome_abundance %>%
              mutate(Info_type = "Abundance"))

writexl::write_xlsx(coda, '07.RESULTS/Compare_genome_type_richness_abundance.xlsx')

## timepoint-wise comparison to mothers:
genome_type_tukeys <- genome_type_richness %>%
  mutate(Info_type = "Richness proportion") %>%
  bind_rows(genome_type_abundance %>%
              mutate(Info_type = "Abundance")) %>%
  #group_by(contrast) %>%
  mutate(p_adjusted = p.adjust(p.value, "BH")) %>%
  ungroup()

writexl::write_xlsx(genome_type_tukeys, '07.RESULTS/Compare_genome_type_richness_abundance_to_Mother.xlsx')

## numbers:

# contribution to maternal virome:
richness_by_genome %>%
  mutate(non_dsDNA = (RNA + ssDNA)/vir_richness_cf) %>%
  filter(Timepoint_new == "Mother") %>%
  pull(non_dsDNA) %>%
  sd()

# abu:
abundance_by_genome %>%
  mutate(non_dsDNA = (RNA)) %>%
  filter(Timepoint_new == "Mother") %>%
  pull(non_dsDNA) %>%
  mean()

# abu:
abundance_by_genome %>%
  mutate(non_dsDNA = (ssDNA)) %>%
  filter(Timepoint_new == "M12") %>%
  pull(non_dsDNA) %>%
  sd()
