setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# holovirome dynamics
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
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
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(MetBrewer)
library(emmeans)
library(patchwork)
#############################################################
# 2. Load Input Data
#############################################################
# abundance tables:
holovirome <- read.table('06.CLEAN_DATA/Intermediate/Holovirome_RPKM_1110samples_120997vOTUs.txt', sep='\t', header=T)

VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)
# 
metaphlan <- read.table('06.CLEAN_DATA/02.FINAL/MGS_Chiliadal_metaphlan_full_taxonomy_ver_01_07102025.txt', sep='\t', header=T)
metaphlan <- metaphlan[grepl('t__SGB', row.names(metaphlan)),]

# virus metadata
ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)
temperate_list <- ETOF_vOTUr$New_CID[ETOF_vOTUr$lifestyle == "Temperate"]
lytic_list <- ETOF_vOTUr$New_CID[ETOF_vOTUr$lifestyle == "Virulent"]

# sample metadata:
smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)

smeta <- smeta %>%
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) %>% # treating 2nd pregnancy as a separate family:
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T))

# holovirome metrics:
smeta$vir_richness_holo <- colSums(holovirome)[match(smeta$Universal_ID, colnames(holovirome))]

smeta$temperate_richness_holo <- colSums(holovirome[row.names(holovirome) %in% temperate_list,])[match(smeta$Universal_ID, colnames(holovirome))]
# 
# # richness of integrated only:
smeta$prophage_richness <- smeta$temperate_richness_holo - smeta$temperate_richness # not so confident in fairness of this metrics, will give it a thought later
# 
smeta$delta <- smeta$vir_richness_holo - smeta$vir_richness_cf
# 
# long smeta (including MGS samples as well)
long_smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)
long_smeta$Timepoint_new <- factor(long_smeta$Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)

long_smeta$bacRichness <- colSums(metaphlan > 0)[match(long_smeta$Sequencing_ID, colnames(metaphlan))]
# 
smeta <- smeta %>%
  left_join(long_smeta %>% filter(seq_type == "MGS") %>% select(Universal_ID, bacRichness))
# 
#############################################################
# 3.1 Analysis: holo vs VLP diversity/richness
#############################################################

# maybe holovirome grows faster / more substantially than VLP because of the temperate phage/prophage contribution
# maybe association w bacterial diversity/richness
# so more bacteria come but their phages do not really get out there or what?

holo_vs_vlp <- smeta %>%
  select(Timepoint_new, vir_richness_cf, vir_richness_holo) %>%
  pivot_longer(cols = !Timepoint_new) %>%
ggplot(aes(Timepoint_new, value, fill = name)) +
  geom_jitter(aes(color = name), position = position_jitterdodge(), size=0.1, alpha=1)+
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  theme_bw() +
  labs(x = "Timepoint", y = "N viruses", fill = "Metavirome\nfraction", color = "Metavirome\nfraction") +
  scale_fill_manual(labels = c("Active", "Holovirome"), values = c(met.brewer("Kandinsky")[2], "#547792")) +
  scale_color_manual(labels = c("Active", "Holovirome"), values = c(met.brewer("Kandinsky")[2], "#547792"))

ggsave('05.PLOTS/06.DYNAMICS/Richness_over_time_holo_VLP.pdf',
       holo_vs_vlp,  "pdf", width=10, height=8, units="cm", dpi = 300)

# COMPARING INFANT AND MATERNAL RICHNESS:
inf_mom_richness <- map_dfr(c("vir_richness_cf", "vir_richness_holo"), function(richness) {
  
  formula <- paste0(richness, "~ Timepoint_new + (1 | FAMILY/NEXT_ID)")
  model <- lmer(formula, 
                data = smeta)
  emm <- emmeans(model, ~ Timepoint_new)
  
  as.data.frame(contrast(emm, method = "pairwise", adjust = "tukey")) %>%
    mutate(richness = richness)
  
})

writexl::write_xlsx(inf_mom_richness, '07.RESULTS/Compare_infant_richness_to_Mother.xlsx')

# TEMPORAL RICHNESS DYMAMICS
rich_growth <- map_dfr(c("vir_richness_cf", "vir_richness_holo"), function(richness) {
  
  formula <- as.formula(paste0(richness, "~ exact_age_months_at_collection + (1|NEXT_ID)"))
  
  model <- lmer(
    formula,
    REML = FALSE,
    data = smeta[smeta$Type == "K",]
  )
  
  model_summary <- summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(rowname != "(Intercept)") %>%
    mutate(richness_type = richness)
})

rich_growth$p_adj <- p.adjust(rich_growth$`Pr(>|t|)`, "BH")

writexl::write_xlsx(rich_growth, '07.RESULTS/Test_infant_richness_over_time.xlsx')

# DOES HOLOVIROME GROW FASTER? (yes)
growth_diff_om <- smeta %>%
  filter(Type == "K") %>%
  select(Timepoint_new, exact_age_months_at_collection, NEXT_ID, vir_richness_cf, vir_richness_holo, bacShannon, bacRichness, prophage_richness) %>%
  pivot_longer(cols = !c(Timepoint_new, exact_age_months_at_collection, NEXT_ID, bacShannon, bacRichness, prophage_richness))

model <- lmer(value ~ exact_age_months_at_collection*name + (1|NEXT_ID), data = growth_diff_om)
growth_diff_om_res <- as.data.frame(summary(model)$coefficients) %>%
  rownames_to_column(var = "rowname")
writexl::write_xlsx(growth_diff_om_res, '07.RESULTS/Compare_holo_vlp_growth_infant_richness_over_time.xlsx')

# DOES HOLOVIROME GROW FASTER BECAUSE OF DIFF IN PROPHAGES?
model <- lmer(prophage_richness ~ exact_age_months_at_collection + (1|NEXT_ID), data = smeta[smeta$Type == "K",])
summary(model) # prophage_richness grows with time

# a separate temporary smeta for residue calculation (excluding infants w NA in exact age, otherwise problems downstream)
inf_smeta_temp <- smeta %>%
  filter(Type == "K" & !is.na(exact_age_months_at_collection))

# bacRichness is too much connceted to time:
cor(inf_smeta_temp$bacRichness,inf_smeta_temp$exact_age_months_at_collection) # 0.63
bac_resid <- residuals(lmer(bacRichness ~ exact_age_months_at_collection + (1 | NEXT_ID), data = inf_smeta_temp))

delta_growth_driver <- summary(lmer(delta ~ exact_age_months_at_collection + prophage_richness + bac_resid + (1 | NEXT_ID), data = inf_smeta_temp))$coefficients

delta_growth_driver <- as.data.frame(delta_growth_driver) %>%
  rownames_to_column(var = "rowname")

writexl::write_xlsx(delta_growth_driver, '07.RESULTS/Drivers_holo_vlp_growth_infant_richness_over_time.xlsx')

# housekeeping
rm(inf_mom_richness, rich_growth, growth_diff_om, model, growth_diff_om_res, delta_growth_driver, inf_smeta_temp)

# supporting graphs:

# increase in the number of integrated prophages (temperate_richness_holo - temperate_richness)
proph_growt <- smeta %>%
  ggplot(aes(Timepoint_new, prophage_richness)) +
  ggrastr::rasterise(geom_jitter(position = position_jitterdodge(), size=0.1, alpha=1, color = "firebrick"), dpi = 300) +
  geom_boxplot(width = 0.3, alpha = 0.3, fill = "firebrick", outlier.shape = NA) + 
  theme_minimal() +
  labs(x = "Timepoint", y = "N integrated prophages") +
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=8))

ggsave("05.PLOTS/06.DYNAMICS/Prophages_growth.pdf", proph_growt, "pdf", width = 10, height=10, dpi = 300, units = "cm")

# increase in the number of bacterial SGBs (from metaphlan)
bac_growt <- smeta %>%
  ggplot(aes(Timepoint_new, bacRichness)) +
  ggrastr::rasterise(geom_jitter(position = position_jitterdodge(), size=0.1, alpha=1, color = "firebrick"), dpi = 300) +
  geom_boxplot(width = 0.3, alpha = 0.3, fill = "firebrick", outlier.shape = NA) + 
  theme_minimal() +
  labs(x = "Timepoint", y = "Bacterial richness") +
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=8))

ggsave("05.PLOTS/06.DYNAMICS/BacRichness_growth.pdf", bac_growt, "pdf", width = 10, height=10, dpi = 300, units = "cm")

#############################################################
# 3.2 Analysis: eukaryotic vs prokaryotic over time
#############################################################

VLP_bin <- as.data.frame(VLP > 0) + 0

# simple host:
obos <- VLP_bin %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, Host_simple)) %>%
  group_by(Host_simple) %>%
  select(-New_CID) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  column_to_rownames(var = "Host_simple")

obos <- as.data.frame(t(obos)) %>%
  select(-Unknown)

smeta <- smeta %>%
  left_join(obos %>% rownames_to_column(var = "Sequencing_ID"), 
            by = "Sequencing_ID",suffix = c("", "_richness"))

# detailed host:
obos <- VLP_bin %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, Host)) %>%
  group_by(Host) %>%
  select(-New_CID) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  column_to_rownames(var = "Host")

obos <- as.data.frame(t(obos)) %>%
  select(-Unknown, -fungi, -`invertebrates, vertebrates`)

smeta <- smeta %>%
  left_join(obos %>% rownames_to_column(var = "Sequencing_ID"), 
            by = "Sequencing_ID",suffix = c("", "_richness"))

# percentage of eukaryotic viruses in the infant gut:
smeta <- smeta %>%
  mutate(euk_perc = Eukaryote_richness/(vir_richness_cf + 1)*100)

euk_perc_stat <- smeta %>%
  filter(Type %in% c("K", "M")) %>%
  group_by(Type) %>%
  summarise(
    mean = mean(euk_perc, na.rm = TRUE),
    median = median(euk_perc, na.rm = TRUE),
    sd = sd(euk_perc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.frame()


smeta <- smeta %>%
  mutate(euk_to_prok = log((Eukaryote_richness +1)/Prokaryote_richness))


euk_to_prok_richness_type <- mixed_model_tukey(smeta, 
                           "euk_to_prok ~ Timepoint_new + (1|NEXT_ID)",
                           "Timepoint_new")

writexl::write_xlsx(euk_to_prok_richness_type, '07.RESULTS/Compare_euk_to_prok_richness_ratio_between_mother_and_infant.xlsx')


# percentage of eukaryotic viruses in the infant gut:

diff_richness <- c("Eukaryote_richness", "Prokaryote_richness", "archaea_richness",
                   "bacteria_richness", "invertebrates_richness", "plants_richness",
                   "protists_richness", "vertebrates_richness")

diff_host_richness_dynamics <- map_dfr(diff_richness, function(richness){
            
            formula <- paste(richness, " ~ exact_age_months_at_collection + (1|NEXT_ID)")
            
            model <- lmer(formula, data = smeta %>% filter(Type == "K"))
            
            summary(model)$coefficients %>%
              as.data.frame() %>%
              rownames_to_column(var = "rowname") %>%
              filter(rowname != "(Intercept)") %>%
              mutate(richness = richness)
              
            
          })

diff_host_richness_dynamics <- diff_host_richness_dynamics %>%
  mutate(p_adj = p.adjust(`Pr(>|t|)`, "BH"))

writexl::write_xlsx(diff_host_richness_dynamics, '07.RESULTS/Infant_richness_by_host_over_time.xlsx')

# different richness over time (w Mother, boxplots):

euks <- smeta %>%
  select(all_of(diff_richness), Timepoint_new, Type, Sequencing_ID) %>%
  select(-Eukaryote_richness, -Prokaryote_richness) %>%
  pivot_longer(any_of(diff_richness)) %>%
  mutate(Host_type = ifelse(name %in% c("bacteria_richness", "archaea_richness"), "Prokaryote", "Eukaryote")) %>%
  mutate(name = gsub("_richness", "", name)) %>%
  filter(Host_type == "Eukaryote") %>%
  ggplot(aes(Timepoint_new, value)) +
  ggrastr::rasterise(geom_jitter(color = "firebrick", position = position_jitterdodge(), size=0.05, alpha=0.4), dpi = 300) +
  geom_boxplot(width = 0.5, alpha = 0.4, fill = "firebrick", outlier.shape = NA) +
  facet_wrap(~ name,  ncol = 1, nrow = 4, scale = "free_y", strip.position = 'right') +
  labs(x = "Timepoint", y = "Virus richness", title = "Eukaryotic hosts") +
  theme_bw() +
  theme(strip.background = element_rect(NA),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, hjust = 0.5))

# different richness over time (w Mother, boxplots):
proks <- smeta %>%
  select(all_of(diff_richness), Timepoint_new, Type, Sequencing_ID) %>%
  select(-Eukaryote_richness, -Prokaryote_richness) %>%
  pivot_longer(any_of(diff_richness)) %>%
  mutate(Host_type = ifelse(name %in% c("bacteria_richness", "archaea_richness"), "Prokaryote", "Eukaryote")) %>%
  mutate(name = gsub("_richness", "", name)) %>%
  filter(Host_type == "Prokaryote") %>%
  ggplot(aes(Timepoint_new, value)) +
  ggrastr::rasterise(geom_jitter(color = "firebrick", position = position_jitterdodge(), size=0.05, alpha=0.4), dpi = 300) +
  geom_boxplot(width = 0.5, alpha = 0.4, fill = "firebrick", outlier.shape = NA) +
  facet_wrap(~ name,  ncol = 1, nrow = 2, scale = "free_y", strip.position = 'right') +
  labs(x = "Timepoint", y = "Virus richness", title = "Prokaryotic hosts") +
  theme_bw() +
  theme(strip.background = element_rect(NA),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, hjust = 0.5))

prok_euks <- euks + proks

ggsave("05.PLOTS/06.DYNAMICS/Prok_euk_richness_growth.pdf",
       prok_euks, "pdf", width = 13, height=11, dpi = 300, units = "cm")

# different richness over time (only infants, lines):
prok_euks_lines <- smeta %>%
  filter(Type == "K") %>%
  select(all_of(diff_richness), Timepoint_new, Type, Sequencing_ID) %>%
  select(-Eukaryote_richness, -Prokaryote_richness) %>%
  pivot_longer(any_of(diff_richness)) %>%
  mutate(Host_type = ifelse(name %in% c("bacteria_richness", "archaea_richness"), "Prokaryote", "Eukaryote")) %>%
  mutate(name = gsub("_richness", "", name)) %>%
  ggplot(aes(Timepoint_new, log10(value + 1), group = name, color = name, fill = name)) +
  #stat_summary(fun.data = "mean_se", geom = "ribbon", alpha = 0.2, color = NA) +
  #stat_summary(fun = mean, geom = "line", lwd = 1) +
  stat_summary(fun.data = "median_hilow", geom = "ribbon", alpha = 0.2, color = NA) +
  stat_summary(fun = median, geom = "line", lwd = 1) +
  facet_wrap(~Host_type, nrow = 1, scales = "free") +
  labs(x = "Timepoint", y = "Virus richness, log10", color = "Virus host", fill = "Virus host") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(NA),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, hjust = 0.5)) +
  scale_color_manual(values = met.brewer("Austria")[c(1:5, 7)]) +
  scale_fill_manual(values = met.brewer("Austria")[c(1:5, 7)])

ggsave("05.PLOTS/06.DYNAMICS/Prok_euk_richness_growth_lines.pdf",
       prok_euks_lines, "pdf", width = 11, height=11, dpi = 300, units = "cm")
#############################################################
# 3.2 Analysis: phage by lifestyle in VLP
#############################################################
# does richness of phages of different lifestyle changes over time?
lifestyle_richness <- smeta %>%
  select(Timepoint_new, temperate_richness, lytic_richness) %>%
  rename(Temperate = temperate_richness, Virulent = lytic_richness) %>%
  pivot_longer(!Timepoint_new) %>%
  ggplot(aes(Timepoint_new, value, fill = name)) +
  ggrastr::rasterise(geom_jitter(aes(color = name), position = position_jitterdodge(), size=0.05, alpha=0.4), dpi = 300) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, linewidth = 0.2) +
  scale_color_manual(values = c("#DA1212", "#11468F")) + 
  scale_fill_manual(values = c("#DA1212", "#11468F")) +
  theme_minimal() + 
  labs(x = "Timepoint", y = "Phage richness in VLP", fill = "Lifestyle", color = "Lifestyle") +
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=8),
        legend.position = "bottom")

ggsave("05.PLOTS/06.DYNAMICS/Lifestyle_richness_by_time.pdf", lifestyle_richness, "pdf", width = 8, height=10, dpi = 300, units = "cm")

# does abundance of phages of different lifestyle changes over time?
lifestyle_abundnce <- smeta %>%
  filter(Type == "K") %>%
  select(Timepoint_new, temperate_RAb, lytic_RAb) %>%
  rename(Temperate = temperate_RAb, Virulent = lytic_RAb) %>%
  pivot_longer(!Timepoint_new) %>%
  ggplot(aes(x = Timepoint_new, y = value, color = name, group = name, fill = name)) +
  stat_summary(fun.data = "mean_se", geom = "ribbon", alpha = 0.2, color = NA) +
  stat_summary(fun = mean, geom = "line", lwd = 1) +
  scale_color_manual(values = c("#DA1212", "#11468F")) + 
  scale_fill_manual(values = c("#DA1212", "#11468F")) +
  theme_minimal() + 
  labs(x = "Timepoint", y = "Relative abundance of phages in VLP", fill = "Lifestyle", color = "Lifestyle") +
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=8),
        legend.position = "bottom")

ggsave("05.PLOTS/06.DYNAMICS/Lifestyle_abundnce_by_time.pdf", lifestyle_abundnce, "pdf", width = 10, height=10, dpi = 300, units = "cm")


# supp (boxplot): # does abundance of phages of different lifestyle changes over time?
lifestyle_abundnce_box <- smeta %>%
  select(Timepoint_new, temperate_RAb, lytic_RAb) %>%
  rename(Temperate = temperate_RAb, Virulent = lytic_RAb) %>%
  pivot_longer(!Timepoint_new) %>%
  ggplot(aes(x = Timepoint_new, y = value, fill = name)) +
  ggrastr::rasterise(geom_jitter(aes(color = name), position = position_jitterdodge(), size=0.05, alpha=0.4), dpi = 300) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  scale_color_manual(values = c("#DA1212", "#11468F")) + 
  scale_fill_manual(values = c("#DA1212", "#11468F")) +
  theme_minimal() + 
  labs(x = "Timepoint", y = "Relative abundance of phages in VLP", fill = "Lifestyle", color = "Lifestyle") +
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=8),
        legend.position = "right")

ggsave("05.PLOTS/06.DYNAMICS/Boxplot_Lifestyle_abundnce_by_time.pdf", lifestyle_abundnce_box, "pdf", width = 10, height=10, dpi = 300, units = "cm")

lifestyle_dynamics <- map_dfr(c("temperate_richness", "lytic_richness", "temperate_RAb", "lytic_RAb"), function(metric){
  
  formula <- paste(metric, "~ exact_age_months_at_collection + (1|NEXT_ID)")
  
  model <- lmer(formula, 
                data = smeta[smeta$Type == "K",])
  
  as.data.frame(summary(model)$coefficients) %>%
    rownames_to_column(var = "rowname") %>%
    filter(rowname != "(Intercept)") %>%
    mutate(metric = metric)
  
})

lifestyle_dynamics <- lifestyle_dynamics %>%
  mutate(p_adj = p.adjust(`Pr(>|t|)`, "BH"))

writexl::write_xlsx(lifestyle_dynamics, '07.RESULTS/Phage_lifestyle_dynamics.xlsx')

# does the infant gut resembles maternal gut in lifestyle?

smeta <- smeta %>%
  mutate(temp_to_lytic_rich_log = log((lytic_richness +1)/(temperate_richness + 1)),
         temp_to_lytic_ab_log = log((lytic_RAb+1)/(temperate_RAb + 1)))

temp_to_lytic_rich_log <- mixed_model_tukey(smeta, 
                                            "temp_to_lytic_rich_log ~ Timepoint_new + (1|NEXT_ID)",
                                            "Timepoint_new") 

temp_to_lytic_ab_log <- mixed_model_tukey(smeta, 
                                            "temp_to_lytic_ab_log ~ Timepoint_new + (1|NEXT_ID)",
                                            "Timepoint_new") 

