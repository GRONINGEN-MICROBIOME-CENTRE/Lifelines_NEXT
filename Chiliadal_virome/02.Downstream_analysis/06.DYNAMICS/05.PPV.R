setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# Here we work only with kids that have at least 3
# timepoints, i.e., 154 kids (536 samples)
# I am honestly running out of names for variables and plots
# One day, I'll change all the repeatedly used object names
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
# no comments
create_summary <- function(data) {
  data %>%
    summarise(across(c("ivr", "PPV_size", "PPV_perc"),
                     list(mean = mean,
                          median = median,
                          sd = sd,
                          q1 = ~quantile(., 0.25),
                          q3 = ~quantile(., 0.75)),
                     .names = "{.col}_{.fn}")) %>%
    pivot_longer(cols = everything(),
                 names_to = c("metric", "statistic"),
                 names_pattern = "(.+)_(mean|median|sd|q1|q3)$",
                 values_to = "value") %>%
    pivot_wider(names_from = statistic,
                values_from = value)
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
  group_by(NEXT_ID, seq_type) %>%
  mutate(n_tp = n()) %>%
  ungroup() %>%
  filter(n_tp >= 3)

# abundance tables etc:
VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)
VLP <- VLP[,colnames(VLP) %in% smeta$Sequencing_ID]
VLP <- VLP[rowSums(VLP) > 0,]

# virus metadata
ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

# shared to mom:
shared_to_mom <- read.table("06.CLEAN_DATA/Intermediate/Mom-shared_vOTUs_table.txt", sep='\t', header=T)
#############################################################
# 3.0 Analysis: PPV in VLP
#############################################################
VLP_long <- VLP %>%
  rownames_to_column("vOTU") %>%
  pivot_longer(-vOTU, names_to = "Sequencing_ID", values_to = "RPKM") %>%
  filter(RPKM != 0) %>% # essential, otherwise some merging will run out of memory
  left_join(smeta, by = "Sequencing_ID")

# individual virome repertoires
ivr <- VLP_long %>%
  group_by(NEXT_ID) %>%
  summarise(n_tps = n_distinct(Timepoint_new), 
            ivr = n_distinct(vOTU),
            .groups = "drop")

# identify persistent vOTUs
persistent_lookup <- VLP_long %>%
  group_by(NEXT_ID, vOTU) %>%
  summarise(n_tps = n_distinct(Timepoint_new), .groups = "drop") %>%
  relocate(n_tps) %>%
  filter(n_tps >= 3) %>%
  mutate(is_persistent = TRUE) 
  

# PPV size per infant
PPV_VLP <- persistent_lookup %>%
  group_by(NEXT_ID) %>%
  summarise(PPV_size = n_distinct(vOTU), .groups = "drop") %>%
  right_join(smeta %>% group_by(NEXT_ID) %>% summarise(n_tp = n_distinct(Timepoint_new), .groups = "drop")) %>%
  mutate(PPV_size = replace_na(PPV_size, 0)) %>%
  left_join(ivr %>% select(-n_tps)) %>% # to avoid confusion
  mutate(PPV_perc = PPV_size / ivr * 100)

testik <- map_dfr(c("ivr", "PPV_size", "PPV_perc"), function(variable) {
  
  F1 <- as.formula(paste0(variable, " ~ n_tp"))
  
  saver <- t.test(F1, data = PPV_VLP)
  
  data.frame(tested = variable,
             p_value = saver$p.value)
  
} ) # N of available timepoints significantly increases IVR (expected) and PPV (also kind of expected)

# summary:
babah <- PPV_VLP %>%
  create_summary() %>%
  mutate(n_infants = nrow(PPV_VLP)) %>%
  bind_rows(PPV_VLP %>%
              filter(n_tp == 3) %>%
              create_summary() %>%
              mutate(metric = paste0(metric, '_', 3),
                     n_infants = nrow(PPV_VLP %>%
                                        filter(n_tp == 3)) ) ) %>%
  bind_rows(PPV_VLP %>%
              filter(n_tp == 4) %>%
              create_summary() %>%
              mutate(metric = paste0(metric, '_', 4),
                     n_infants = nrow(PPV_VLP %>%
                                        filter(n_tp == 4))))
rm(ivr)
#############################################################
# 3.1 Analysis: PPV in VLP, rarefaction
#############################################################
# lucky kids with 4 timepoints
lucky <- smeta %>%
  group_by(NEXT_ID) %>%
  summarise(n_tp = n(), .groups = "drop") %>%
  filter(n_tp > 3) %>%
  pull(NEXT_ID)

combos <- combn(unique(as.character(smeta$Timepoint_new)), m = 3, simplify = FALSE)

ppv_rarefied <- map_dfr(combos, function(combo) {
  
  VLP_long %>%
    filter(NEXT_ID %in% lucky) %>%
    filter(Timepoint_new %in% combo) %>%
    group_by(NEXT_ID, vOTU) %>%
    summarise(n_tps = n_distinct(Timepoint_new), .groups = "drop") %>%
    filter(n_tps >= 3) %>%
    mutate(is_persistent = TRUE) %>%
    mutate(combos = paste0(combo, collapse = ",")) %>%
    group_by(NEXT_ID, combos) %>%
    summarise(PPV_size = n_distinct(vOTU), .groups = "drop") %>%
    right_join(smeta %>%
                 filter(NEXT_ID %in% lucky) %>%
                 group_by(NEXT_ID) %>%
                 summarise(n_tp = n())) %>%
    mutate(PPV_size = replace_na(PPV_size, 0),
           combos = replace_na(combos, paste0(combo, collapse = ",")))
    
  
})

ivr_rarefied <- map_dfr(combos, function(combo) {
  
  VLP_long %>%
    filter(NEXT_ID %in% lucky) %>%
    filter(Timepoint_new %in% combo) %>%
    mutate(combos = paste0(combo, collapse = ",")) %>%
    group_by(NEXT_ID, combos) %>%
    summarise(IV_size = n_distinct(vOTU), .groups = "drop") %>%
    right_join(smeta %>%
                 filter(NEXT_ID %in% lucky) %>%
                 group_by(NEXT_ID) %>%
                 summarise(n_tp = n())) %>%
    mutate(IV_size = replace_na(IV_size, 0))
  
  
})

results <- ppv_rarefied  %>%
  group_by(NEXT_ID) %>%
  summarise(PPV_min = round(min(PPV_size), 0),
            PPV_max = round(max(PPV_size), 0),
    PPV_size= round(mean(PPV_size), 0),
            .groups = "drop") %>%
  mutate(n_tp = 4) %>%
  left_join(ivr_rarefied %>%
              group_by(NEXT_ID) %>%
              summarise(IV_min = round(min(IV_size), 0),
                        IV_max = round(max(IV_size), 0),
                        ivr= round(mean(IV_size), 0),
                        .groups = "drop") %>%
              mutate(n_tp = 4)) %>%
  mutate(PPV_perc = PPV_size / ivr * 100)


all_kids_ppv <- results %>%
  bind_rows(PPV_VLP %>%
              filter(n_tp == 3))

testik_rarefied <- map_dfr(c("ivr", "PPV_size", "PPV_perc"), function(variable) {
  
  F1 <- as.formula(paste0(variable, " ~ n_tp"))
  
  saver <- t.test(F1, data = all_kids_ppv)
  
  data.frame(tested = variable,
             p_value = saver$p.value)
  
} ) # no more separation by N timepoints

new_summary <- all_kids_ppv %>%
  create_summary()

# tidy up:
rm(list = c("babah", "combos", "ivr_rarefied", "ppv_rarefied", "PPV_VLP", "results", "testik", "testik_rarefied"))

ivr_ppv <- all_kids_ppv %>%
  select(PPV_size, ivr, NEXT_ID) %>%
  pivot_longer(!NEXT_ID) %>%
  ggplot(aes(name, log10(value + 1), color = name, fill = name)) +
  ggrastr::rasterise(geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), size=0.01, alpha=1), dpi = 300)+
  geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.3, color = "black") +
  labs(x = "", y = "N vOTUs, log10") +
  scale_x_discrete(labels = c("All vOTUs\never\ndetected", "PPV")) +
  scale_fill_manual(values = met.brewer("Cassatt2")[c(2,9)]) +
  scale_color_manual(values = met.brewer("Cassatt2")[c(2,9)]) +
  theme_minimal() +
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.position = "none")
  
ggsave("05.PLOTS/06.DYNAMICS/IVR_PPV_log_scaled_absolute.pdf",
       ivr_ppv, "pdf", width = 5, height = 7, units = "cm", dpi = 300)

#############################################################
# 3.2 Analysis: PPV percentage & abundance per timepoint 
#############################################################

busya <- all_kids_ppv %>%
  left_join(smeta %>%
              select(NEXT_ID, Timepoint_new, vir_richness_cf) %>%
              pivot_wider(names_from = Timepoint_new, values_from = vir_richness_cf)) %>%
  select(NEXT_ID, ivr, M1, M3, M6, M12) %>%
  pivot_longer(!NEXT_ID)

summary(busya$value[busya$name == "M12"], na.rm = T)

bodel <- lm(value ~ name, data = busya)

TUKEY <- TukeyHSD(aov(bodel))

my_comparisons <- list( c("ivr", "M1"), 
                        c("ivr", "M3"), 
                        c("ivr", "M6"),
                        c("ivr", "M12"))

ivr_vs_richness <- busya %>%
  mutate(color_factor = ifelse(name == "ivr", "All vOTUs\never\ndetected", "Timepoint")) %>%
  mutate(name = factor(name, levels = c("ivr", "M1", "M3", "M6", "M12"), ordered = T)) %>%
  ggplot(aes(name, value, fill = color_factor, color = color_factor)) +
  ggrastr::rasterise(geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), size=0.01, alpha=1), dpi = 300)+
  geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.5, width = 0.3) +
  labs(x = "", y = "N vOTUs, log10") +
  scale_x_discrete(labels = c("All vOTUs ever\ndetected", "M1", "M3", "M6", "M12")) +
  theme_minimal() +
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.position = "none") +
  scale_fill_manual(values = c(met.brewer("Cassatt2")[2], met.brewer("Kandinsky")[2])) +
  scale_color_manual(values = c(met.brewer("Cassatt2")[2], met.brewer("Kandinsky")[2])) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size=3, p.adjust.method = "BH")

ggsave("05.PLOTS/06.DYNAMICS/IVR_vs_richness.pdf",
       ivr_vs_richness, "pdf", width = 11, height = 7, units = "cm", dpi = 300)

# timepoint-wise detected PPVs:
persistent_rich <- VLP_long %>%
  left_join(persistent_lookup, by = c("NEXT_ID", "vOTU")) %>%
  mutate(is_persistent = replace_na(is_persistent, FALSE)) %>%
  group_by(NEXT_ID, Timepoint_new) %>%
  summarise(
    richness = n(), 
    PPV_richness = sum(is_persistent),
    PPV_fraction = PPV_richness / richness * 100,
    .groups = "drop"
  )

summary(persistent_rich$PPV_fraction)

persistent_abundance <- VLP_long %>%
  left_join(persistent_lookup, by = c("NEXT_ID", "vOTU")) %>%
  mutate(is_persistent = replace_na(is_persistent, FALSE)) %>%
  group_by(NEXT_ID, Timepoint_new) %>%
  summarise(
    total_RPKM = sum(RPKM, na.rm = TRUE),
    RPKM_persistent = sum(ifelse(is_persistent == TRUE, RPKM, 0), na.rm = TRUE),
    PPV_abundance = RPKM_persistent / total_RPKM * 100
  ) %>%
  ungroup()

summary(persistent_abundance$PPV_abundance)

summary(lmer(PPV_abundance ~ Timepoint_new + (1|NEXT_ID), data = persistent_abundance)) #Timepoint_new.L  -17.387      2.165 401.923  -8.030 1.08e-14 ***

rich_ppv <- persistent_rich %>%
  ggplot(aes(Timepoint_new, PPV_fraction)) +
  ggrastr::rasterise(geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), color = "#A6B1E1", size=0.01, alpha=1), dpi = 300)+
  geom_boxplot(color = "#424874", fill = "#DCD6F7",outlier.shape = NA, alpha = 0.5, width = 0.3) +
  labs(x = "Timepoint", y = "PPV contribution to infant richness") +
  theme_minimal() +
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.position = "none") 

ggsave("05.PLOTS/06.DYNAMICS/Rich_PPV_VLP.pdf",
       rich_ppv, "pdf", width = 6, height = 7, units = "cm", dpi = 300)

abu_ppv <- persistent_abundance %>%
  ggplot(aes(Timepoint_new, PPV_abundance)) +
  ggrastr::rasterise(geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), color = "#A6B1E1", size=0.01, alpha=1), dpi = 300)+
  geom_boxplot(color = "#424874", fill = "#DCD6F7",outlier.shape = NA, alpha = 0.5, width = 0.3) +
  labs(x = "Timepoint", y = "PPV cumulative abundance") +
  theme_minimal() +
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.position = "none") 

ggsave("05.PLOTS/06.DYNAMICS/RAb_PPV_VLP.pdf",
       abu_ppv, "pdf", width = 6, height = 7, units = "cm", dpi = 300)
#############################################################
# 3.3 Analysis: PPV composition 
#############################################################
PPV_metadata <- persistent_lookup %>%
  group_by(vOTU) %>%
  summarise(prev = n(), .groups = "drop") %>%
  left_join(ETOF_vOTUr, by = c("vOTU" = "New_CID")) %>%
  separate(
    col    = tax_ictv_aai,
    into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  separate(
    col    = Host_taxonomy,
    into = paste0('Host_', c('Domain', 'Phylum',  'Class', 'Order', 'Family',  'Genus')),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  mutate(Host_Genus = ifelse(Host_Genus %in% c("g__Bacteroides", "g__Phocaeicola"), "Bacteroides_Phocaeicola", Host_Genus))

by_simple_host <- PPV_metadata %>%
  group_by(Host) %>%
  summarise(n = n()) # bacteriophages total: 2640
  
by_host_gen <- PPV_metadata %>%
  group_by(Host, genome) %>%
  summarise(n = n()) # bacteriophages total: 2640
  
PPV_RPKM <- PPV_metadata %>% 
  mutate(Order = ifelse(Order == "Unassigned", "Unclassified", Order),
         Host = ifelse(Host == "Unknown", "unknown", Host)) %>%
  group_by(genome, Host, Order) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(Host = factor(Host, levels = c("archaea", "bacteria", "invertebrates", "protists", "vertebrates", "unknown"), order = T))
  
ppv_met <- ggplot(PPV_RPKM, aes(x = Host, y = Order,
                              size = n, colour = genome)) +
  geom_point(alpha = 0.85) +
  geom_text(aes(label = n), colour = "black",
            size = 2, show.legend = FALSE) +
  scale_size_area(max_size = 6, name = "N PPV vOTUs", transform = "log10") +
  scale_colour_manual(values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                                 "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                                 "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2],
                                 "Unclassified" = "grey"), name = "Genome") +
  scale_x_discrete(position = "top") +
  labs(x = "Predicted host", y = "Taxonmic rank: order") +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 30, hjust = 0, vjust = -0.1),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.position  = "bottom",
        legend.box = "vertical",
        legend.box.spacing = unit(0, "cm"),
        legend.margin=unit(0, "cm")) +
  guides(color  = guide_legend(order = 1),
         size = guide_legend(order = 2))

ggsave("05.PLOTS/06.DYNAMICS/PPV_metadata_depicted.pdf",
       ppv_met, "pdf", width = 8, height = 9, units = "cm", dpi = 300)

phage_PPV <- PPV_metadata %>% 
  filter(Host == "bacteria") %>%
  group_by(Host_Genus) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(prop = round(n/2501*100, 1))

PPV_phages_RPKM <- persistent_lookup %>%
  left_join(PPV_metadata) %>% 
  mutate(Order = case_when(Order == "Unassigned" ~ "Unclassified",
                           Order == "Crassvirales" | Order == "Petitvirales" ~ "Crass or Petit",
                           .default = Order),
         Host = ifelse(Host == "Unknown", "unknown", Host)) %>%
  group_by(NEXT_ID, Order) %>%
  summarise(
    PPV_richness = sum(is_persistent),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = NEXT_ID,
              values_from = PPV_richness,
              values_fill = 0) %>%
  mutate(prevalence = rowSums(across(-c(Order)) > 0)) %>%
  relocate(prevalence)


PPV_euk_RPKM <- persistent_lookup %>%
  left_join(PPV_metadata) %>% 
  filter(Host %in% c("invertebrates", "vertebrates")) %>%
  group_by(NEXT_ID, Host, Family) %>%
  summarise(
    PPV_richness = sum(is_persistent),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = NEXT_ID,
              values_from = PPV_richness,
              values_fill = 0) %>%
  mutate(prevalence = rowSums(across(-c(Host, Family)) > 0)) %>%
  relocate(prevalence)

by_lfs <- PPV_metadata %>%
  filter(Host == "bacteria") %>%
  group_by(lifestyle) %>%
  summarise(n = n())

#############################################################
# 3.4 Analysis: PPV - M1-colonizers overlap
#############################################################
# identify persistent vOTUs from M1
persistent_lookup_M1 <- VLP_long %>%
  filter(Timepoint_new == "M1") %>%
  right_join(persistent_lookup, by=c("NEXT_ID", "vOTU")) %>%
  group_by(vOTU) %>%
  summarise(PPV = n_distinct(vOTU),
            M1_PPV = n_distinct(vOTU[!is.na(Timepoint_new)]))


rev <- persistent_lookup %>%
  left_join(VLP_long %>% 
              filter(Timepoint_new == "M1") %>% 
              mutate(detected_M1 = T) %>%
              select(NEXT_ID, vOTU, detected_M1),
            by = c("NEXT_ID", "vOTU") ) %>%
  mutate(detected_M1 = replace_na(detected_M1, F)) %>%
  group_by(NEXT_ID) %>%
  summarise(
    total_PPV_in_infant = n_distinct(vOTU),
    PPV_detected_at_M1_in_infant = sum(detected_M1),
    proportion_M1_among_PPV = PPV_detected_at_M1_in_infant / total_PPV_in_infant
  ) %>%
  summarise(
    mean_proportion = mean(proportion_M1_among_PPV),
    sd_proportion = sd(proportion_M1_among_PPV)
  )

####### are M1 enriched?
pair_level <- VLP_long %>%
  group_by(NEXT_ID, vOTU) %>%
  summarise(
    is_PPV     = n_distinct(Timepoint_new) >= 3,
    seen_at_M1 = any(Timepoint_new == "M1"),
    .groups = "drop"
  ) %>%
  mutate(
    PPV    = factor(is_PPV,     levels = c(TRUE, FALSE), labels = c("PPV", "TDV")),
    At_M1  = factor(seen_at_M1, levels = c(TRUE, FALSE), labels = c("Detected M1", "Not detected M1"))
  )

ft <- fisher.test(table(pair_level$PPV, pair_level$At_M1)) # p-value < 2.2e-16, odds ratio: 5.9

plot_data <- pair_level %>%
  count(PPV, At_M1) %>%
  group_by(PPV) %>%
  mutate(
    prop      = n / sum(n),
    grp_n     = sum(n)
  ) %>%
  ungroup() %>%
  mutate(bar_width = grp_n / nrow(pair_level))

M1_mosaic <- ggplot(plot_data, aes(x=PPV, y=prop, width=grp_n, fill=At_M1)) + 
  geom_bar(stat = "identity", position = "fill", colour = "black") + 
  facet_grid(~PPV, scales = "free_x", space = "free_x")  +
  geom_text(aes(label = scales::percent(prop)), position = position_stack(vjust = 0.5), size=2) +
  scale_fill_manual(values=c(met.brewer("Kandinsky")[1],met.brewer("Kandinsky")[2])) + 
  labs(x="", y="Proportion of (Infant × vOTU) pairs", fill="") + 
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text=element_text(size=8), 
        axis.text.y = element_blank(),
        axis.title=element_text(size=9),
        panel.spacing.x = unit(0, "npc"),
        legend.text = element_text(size=8, hjust=0.5),
        legend.title = element_text(size=9),
        legend.position = "bottom",
        legend.key.size=unit(0.7, "line"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0)
        )

ggsave('05.PLOTS/06.DYNAMICS/M1_vOTUs_enrichment_PPV.pdf',
       M1_mosaic, "pdf", width=7, height=8, units="cm", dpi = 300)

#############################################################
# 3.5 Analysis: PPV - mom-shared overlap
#############################################################
moms <- shared_to_mom %>%
  rownames_to_column(var = "vOTU") %>%
  pivot_longer(!vOTU, names_to = "Sequencing_ID") %>%
  mutate(mom_shared = ifelse(value !=0, T, F)) %>%
  right_join(smeta %>% select(Timepoint_new, NEXT_ID, Sequencing_ID)) %>%
  group_by(NEXT_ID, vOTU) %>%
  summarise(shared_mom = any(mom_shared), .groups = "drop")

pair_level2 <- pair_level %>%
  left_join(moms, by = c("NEXT_ID", "vOTU")) %>%
  mutate(shared_mom = replace_na(shared_mom, F)) %>%
  mutate(shared_mom = ifelse(shared_mom == T, "Yes", "No")) %>%
  mutate(shared_mom = factor(shared_mom, levels = c("Yes", "No"), ordered = T))

# 2x2: PPV ~ mom sharing
df_mom <- pair_level2 %>%
  count(PPV, shared_mom) %>%
  group_by(PPV) %>%
  mutate(prop = n / sum(n), grp_n = sum(n)) %>%
  ungroup() %>%
  mutate(bar_width = grp_n / nrow(pair_level))

ft2 <- fisher.test(table(pair_level2$shared_mom,
                  relevel(pair_level2$PPV, ref = "PPV")))

mom <- ggplot(df_mom, aes(x=PPV, y=prop, width=grp_n, fill=shared_mom)) + 
  geom_bar(stat = "identity", position = "fill", colour = "black") + 
  facet_grid(~PPV, scales = "free_x", space = "free_x")  +
  geom_text(aes(label = scales::percent(prop)), position = position_stack(vjust = 0.5), size=2) +
  scale_fill_manual(values=c(met.brewer("Kandinsky")[1],met.brewer("Kandinsky")[2])) + 
  labs(x="", y="Proportion of (Infant × vOTU) pairs", fill="Shared to mom") + 
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text=element_text(size=8), 
        axis.text.y = element_blank(),
        axis.title=element_text(size=9),
        panel.spacing.x = unit(0, "npc"),
        legend.text = element_text(size=8, hjust=0.5),
        legend.title = element_text(size=9),
        legend.position = "bottom",
        legend.key.size=unit(0.7, "line"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0))

ggsave('05.PLOTS/06.DYNAMICS/Mom_vOTUs_enrichment_PPV.pdf',
       mom, "pdf", width=7, height=8, units="cm", dpi = 300)

#############################################################
# 3.6 Analysis: PPV - lifestyle
#############################################################
bair_level <- VLP_long %>%
  left_join(ETOF_vOTUr, by = c("vOTU" = "New_CID")) %>%
  filter(lifestyle != "Unknown") %>%
  group_by(NEXT_ID, vOTU) %>%
  summarise(
    is_PPV     = n_distinct(Timepoint_new) >= 3,
    is_temperate = ifelse(lifestyle == "Temperate", T, F),
    .groups = "drop"
  ) %>%
  mutate(
    PPV    = factor(is_PPV,     levels = c(TRUE, FALSE), labels = c("PPV", "TDV")),
    temperate  = factor(is_temperate, levels = c(FALSE, TRUE), labels = c("Virulent", "Temperate"))
  )

ft <- fisher.test(table(bair_level$PPV, bair_level$temperate)) 

blot_data <- bair_level %>%
  count(PPV, temperate) %>%
  group_by(PPV) %>%
  mutate(
    prop      = n / sum(n),
    grp_n     = sum(n)
  ) %>%
  ungroup() %>%
  mutate(bar_width = grp_n / nrow(bair_level))

lfs <- ggplot(blot_data, aes(x=PPV, y=prop, width=grp_n, fill=temperate)) + 
  geom_bar(stat = "identity", position = "fill", colour = "black") + 
  facet_grid(~PPV, scales = "free_x", space = "free_x")  +
  geom_text(aes(label = scales::percent(prop)), position = position_stack(vjust = 0.5), size=2) +
  scale_fill_manual(values=c(met.brewer("Kandinsky")[1],met.brewer("Kandinsky")[2])) + 
  labs(x="", y="Proportion of (Infant × vOTU) pairs", fill="Lifestyle") + 
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text=element_text(size=8), 
        axis.text.y = element_blank(),
        axis.title=element_text(size=9),
        panel.spacing.x = unit(0, "npc"),
        legend.text = element_text(size=8, hjust=0.5),
        legend.title = element_text(size=9),
        legend.position = "bottom",
        legend.key.size=unit(0.7, "line"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0)
  )

ggsave('05.PLOTS/06.DYNAMICS/Virulent_enrichment_PPV.pdf',
       lfs, "pdf", width=7, height=8, units="cm", dpi = 300)
