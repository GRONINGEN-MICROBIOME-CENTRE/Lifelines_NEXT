setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# Here we work only with kids that have all 4 timepoints
# i.e., 84 kids (296 fecal samples) trying to study
# induction and lysogenization events
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
library(patchwork)
#############################################################
# 2. Load Input Data
#############################################################
long_smeta <- read.delim("06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt", sep='\t', header=T)

long_smeta <- long_smeta %>%
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) %>% # treating 2nd pregnancy as a separate family:
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)) %>%
  group_by(NEXT_ID, seq_type) %>%
  mutate(n_tp = n()) %>%
  ungroup() %>%
  filter(n_tp == 4)

# abundance tables etc:
# VLP
VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)
VLP <- VLP[,colnames(VLP) %in% long_smeta$Sequencing_ID]
VLP <- VLP[rowSums(VLP) > 0,]

# MGS
MGS <- read.table('06.CLEAN_DATA/02.FINAL/MGS_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)
MGS <- MGS[,colnames(MGS) %in% long_smeta$Sequencing_ID]
MGS <- MGS[rowSums(MGS) > 0,]

# MGS
HLV <- read.table('06.CLEAN_DATA/Intermediate/Holovirome_RPKM_1110samples_120997vOTUs.txt', sep='\t', header=T)
HLV <- HLV[,colnames(HLV) %in% long_smeta$Universal_ID]
HLV <- HLV[rowSums(HLV) > 0,]

# virus metadata
ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)
#############################################################
# 3.0 Analysis: selection of HLV-persistent temperate vOTUs
#############################################################
HLV_long <- HLV %>%
  rownames_to_column("vOTU") %>%
  pivot_longer(-vOTU, names_to = "Universal_ID", values_to = "RPKM") %>%
  filter(RPKM != 0) %>% # essential, otherwise some merging will run out of memory
  filter(vOTU %in% ETOF_vOTUr$New_CID[ETOF_vOTUr$lifestyle == "Temperate"]) %>%
  left_join(long_smeta %>% filter(seq_type == "VLP") %>% select(Sequencing_ID, NEXT_ID, Universal_ID, Timepoint_new), by = "Universal_ID")

# persistent temperate regardless of metavirome:
HLV_persistent_lookup <- HLV_long %>%
  group_by(NEXT_ID, vOTU) %>%
  summarise(n_tps_HLV = n_distinct(Timepoint_new), .groups = "drop") %>%
  filter(n_tps_HLV == 4) %>% # select 4-timepoint persisters
  mutate(is_persistent_HLV = TRUE) %>%
  filter(vOTU %in% ETOF_vOTUr$New_CID[ETOF_vOTUr$lifestyle == "Temperate"]) # select only temperate

length(unique(HLV_persistent_lookup$vOTU)) # 694

#############################################################
# 3.1 Analysis: VLP vs MGS detection of holovirome persisters
#############################################################
VLP_detection <- VLP %>%
  rownames_to_column("vOTU") %>%
  pivot_longer(-vOTU, names_to = "Sequencing_ID", values_to = "RPKM") %>%
  filter(RPKM != 0) %>% # essential, otherwise some merging will run out of memory
  left_join(long_smeta %>% filter(seq_type == "VLP") %>% select(Sequencing_ID, NEXT_ID, Universal_ID, Timepoint_new), by = "Sequencing_ID") %>%
  mutate(VLP_detected = T)

MGS_detection <- MGS %>%
  rownames_to_column("vOTU") %>%
  pivot_longer(-vOTU, names_to = "Sequencing_ID", values_to = "RPKM") %>%
  filter(RPKM != 0) %>% # essential, otherwise some merging will run out of memory
  left_join(long_smeta %>% filter(seq_type == "MGS") %>% select(Sequencing_ID, NEXT_ID, Universal_ID, Timepoint_new), by = "Sequencing_ID") %>%
  mutate(MGS_detected = T)

# table for all 4-time persisters:
across <- HLV_long %>%
  left_join(HLV_persistent_lookup %>% select(vOTU, NEXT_ID, is_persistent_HLV)) %>% # cleaning HLV_long to only contain persisters:
  filter(!is.na(is_persistent_HLV)) %>%
  select(-RPKM, -Sequencing_ID) %>%
  left_join(VLP_detection %>% select(-RPKM, -Sequencing_ID)) %>%
  left_join(MGS_detection %>% select(-RPKM, -Sequencing_ID)) %>%
  mutate(VLP_detected = replace_na(VLP_detected, F),
         MGS_detected = replace_na(MGS_detected, F))

# timepoint-wise classification of vOTU detection:
tp_classified <- across %>%
  mutate(tp_state = case_when(VLP_detected == 1 ~ "virion",   # VLP present -> virion
                              VLP_detected == 0 & MGS_detected == 1 ~ "prophage"))   # MGS only present -> prophage

# infant-wise classification of vOTUs:
votu_kid_labels <- tp_classified %>%
  group_by(vOTU, NEXT_ID) %>%
  summarise(n_induced = sum(tp_state == "virion"),
            n_dormant = sum(tp_state == "prophage"),
            n_tp = n(),
            .groups   = "drop") %>%
  mutate(category = case_when(
    n_induced > 0 & n_dormant == 0 ~ "always_induced",
    n_dormant > 0 & n_induced == 0 ~ "always_dormant",
    n_induced > 0 & n_dormant > 0  ~ "switcher"
    )
  )

# vOTU-wise classification of vOTUs:
global_labels <- votu_kid_labels %>%
  group_by(vOTU) %>%
  summarise(n_kids_dormant  = sum(category == "always_dormant"),
            n_kids_induced  = sum(category == "always_induced"),
            n_kids_switcher = sum(category == "switcher"),
            n_kids_total = n(),
            .groups = "drop") %>%
  mutate(global_category = case_when(n_kids_switcher > 0 ~ "switcher",
                                     n_kids_induced  > 0 & n_kids_dormant == 0 ~ "always_induced",
                                     n_kids_dormant  > 0 & n_kids_induced == 0 ~ "always_dormant",
                                     TRUE ~ "mixed"   # induced in some, dormant in others (no switcher)
    )
  )

table(global_labels$global_category) # 142 vOTU are PPV and always dormant

#############################################################
# 3.2 Analysis: visualization of detection patterns
#############################################################
all_patterns <- tp_classified %>%
  arrange(NEXT_ID, vOTU, Timepoint_new) %>%
  group_by(vOTU, NEXT_ID) %>%
  summarise(pattern = paste(tp_state, collapse = "_"), .groups = "drop") %>%
  count(pattern, sort = TRUE) %>%
  mutate(
    traj_group = case_when(
      pattern == "prophage_prophage_prophage_prophage" ~ "always_prophage",
      pattern == "prophage_prophage_prophage_virion" ~ "prophage_then_virion1",
      pattern == "prophage_prophage_virion_virion" ~ "prophage_then_virion2",
      pattern == "prophage_virion_virion_virion" ~ "prophage_then_virion3",
      pattern == "virion_virion_virion_virion" ~ "always_virion",
      pattern == "virion_virion_virion_prophage" ~ "virion_then_prophage3",
      pattern == "virion_virion_prophage_prophage" ~ "virion_then_prophage2",
      pattern == "virion_prophage_prophage_prophage" ~ "virion_then_prophage1",
      .default = pattern)) %>%
  mutate(category = case_when(traj_group == "always_prophage" ~ "Stable prophage",
                              grepl("prophage_then", traj_group) ~ "Induction",
                              grepl("virion_then", traj_group) ~ "Lysogenization",
                              traj_group == "always_virion" ~ "Stable replication",
                              .default = "Switcher"))
# ordering:
all_patterns_ord <- all_patterns %>%
  mutate(traj_group = factor(traj_group,
                             levels = c("always_prophage", "prophage_then_virion1", "prophage_then_virion2", "prophage_then_virion3",
                                        "always_virion", "virion_then_prophage3", "virion_then_prophage2", "virion_then_prophage1",
                                        all_patterns$traj_group[all_patterns$category == "Switcher"]), ordered = T)) 

pattern_order <- all_patterns_ord %>%
  arrange(desc(traj_group)) %>% # inverse order, otherwise does not look nice
  mutate(ymax = cumsum(n),
         ymin = lag(ymax, default = 0),
         ymid = (ymin + ymax) / 2)

tp_levels <- c("M1", "M3", "M6", "M12")

heatmap_data <- pattern_order %>%
  mutate(states = strsplit(pattern, "_")) %>%
  unnest_longer(states) %>%
  group_by(pattern) %>%
  mutate(timepoint = factor(tp_levels[row_number()], levels = tp_levels)) %>%
  ungroup() %>%
  rename(tp_state = states)

# category annotation positions
group_annot <- heatmap_data %>%
  group_by(category) %>%
  summarise(ymin = min(ymin), ymax = max(ymax), ymid = mean(ymid),
            .groups = "drop") %>%
  arrange(desc(ymin)) %>%
  mutate(category = factor(category, levels = unique(group_annot$category), ordered = T))


plot <- ggplot() +
  geom_tile(data = heatmap_data,
            aes(x = timepoint, y = ymid, height = n, fill = tp_state),
            color = "white", linewidth = 0.3) +
  geom_text(data = pattern_order,
            aes(x = 4.6, y = ymid, label = paste0("n=", n)),
            hjust = 0, size = 2, color = "grey30") +
  scale_fill_manual(values = c("virion" = "#B31312", "prophage" = "#D8D9CF"),
                    breaks = c("virion", "prophage"),
                    labels = c("Virion (VLP)", "Prophage (MGS only)")) +
  labs(x = "Timepoint", y = "Infant x temperate persistent vOTU", fill = "State") +
  ggnewscale::new_scale_fill() +
  geom_rect(data = group_annot,
            aes(xmin = -0.5, xmax = 0.52, ymin = ymin, ymax = ymax, fill = category),
            inherit.aes = FALSE, show.legend = F) +
  geom_text(data = group_annot,
            aes(x = 1, y = ymid, label = category),
            inherit.aes = FALSE, hjust = 0.5, vjust = -7, size = 2, angle = 90) +
  scale_fill_manual(values = met.brewer("OKeeffe1")[c(8,7,4,5,6)]) +
  geom_hline(data = group_annot %>% filter(ymin > 0),
             aes(yintercept = ymin),color = "white", linewidth = 1.5)  +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-0.5, 5), clip = "off") +
  theme_minimal() +
  theme(
    axis.text.y     = element_blank(),
    axis.text.x = element_text(size=8),
    axis.ticks.y    = element_blank(),
    panel.grid      = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 9),
    legend.title = element_text(size=9),
    legend.text = element_text(size=8)
  )

ggsave("05.PLOTS/06.DYNAMICS/HLV_temperate_PPV_VLP_MGS.pdf",
       plot, "pdf", width = 5.5, height = 10, units = "cm", dpi = 300)
