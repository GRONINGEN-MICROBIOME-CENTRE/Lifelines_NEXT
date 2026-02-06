setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

#############################################################
# MGS & VLP overlap and unique vOTUs bar graph
#############################################################

#############################################################
# 0. Used file source
#############################################################
## Extended_ETOF_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt -
## Extended filtered ETOF for 127,553 vOTUs identified in 2,220 
## MGS and VLP samples (above 3kbp); merged with VITAP and 
## geNomad assignemnts & lifestyle assignment from BACPHLIP

## HQ_vOTU_recoverable_by_source.txt - table with per-contig HQ vOTU recovery info by source,
## derived using the post-CheckV clustering information

## HQ_vOTU_recovery_stat.txt - summary stats of HQ vOTU recovery across sources
## dervied based on HQ_vOTU_recoverable_by_source.txt 
#############################################################
# 1. Functions
#############################################################

# ggplot2 theme used recurrently in a few plots:
theme_VLP_MGS_plots <- theme_minimal() +
  theme(legend.position = "right", 
        legend.text = element_text(size=8), 
        legend.title.position = "top",
        legend.title = element_text(size=10, hjust=0.5),
        axis.title.x = element_text(size=8, face="bold"), 
        legend.key.size = unit(0.8,"line"),
        legend.margin=margin(0,0,0,0),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
  )
#############################################################
# 2. Loading libraries
#############################################################
library(tidyverse)
library(dplyr)

library(ggplot2)
#############################################################
# 3. Load Input Data
#############################################################
ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

hq_by_recovery <- read.table('06.CLEAN_DATA/03.RESULTS_to_plot/HQ_vOTU_recoverable_by_source.txt', sep='\t', header = T)

#############################################################
# 2. Analysis: novel vOTU by recovery source
#############################################################
# vOTU members sharing:
overall <- hq_by_recovery %>%
  filter(!grepl('DB', vOTU_cluster_type)) %>%
  group_by(vOTU_cluster_type) %>%
  summarise(N_vOTUs = n(), .groups = "drop") %>%
  mutate(perc = N_vOTUs/nrow(hq_by_recovery) * 100) %>%
  rename(Part = vOTU_cluster_type) %>%
  mutate(Part = factor(Part, levels = c("NEXT_MGS", "NEXT_VLP+NEXT_MGS", "NEXT_VLP", ordered = T))) %>%
  arrange(Part) %>%
  mutate(xmin = c(0,
                  cumsum(N_vOTUs)[-n()])) %>%
  mutate(xmax = xmin + N_vOTUs,
         ymin = 0,
         ymax = 1)

recovery <- hq_by_recovery %>%
  filter(!grepl('DB', vOTU_cluster_type)) %>%
  group_by(recoverable_by) %>%
  summarise(N_vOTUs = n(), .groups = "drop") %>%
  mutate(hq_recovery_source = factor(recoverable_by,
                                     levels = c('MGS_QoL', 'MGS_QoL+VLP_QoL', 'VLP_QoL'),
                                     ordered = T)) %>%
  arrange(hq_recovery_source) %>%
  mutate(xmin = c(0,
                  cumsum(N_vOTUs)[-n()])) %>%
  mutate(xmax = xmin + N_vOTUs,
         ymin = 0,
         ymax = 1)
  
# lifestyle & genome annotation
# by genome type
by_genome <- hq_by_recovery %>%
  filter(!grepl('DB', vOTU_cluster_type)) %>%
  group_by(genome, recoverable_by) %>%
  summarise(N_vOTUs = n(), .groups = "drop") %>%
  rename(source = recoverable_by) %>%
  mutate(source = factor(source, levels = c("MGS_QoL", "MGS_QoL+VLP_QoL", "VLP_QoL"), ordered = TRUE)) %>%
  mutate(genome = factor(genome, levels = c('dsDNA', 'RNA', 'ssDNA', 'Unclassified'), ordered = TRUE)) %>%
  arrange(source, genome) %>%
  mutate(xmin = c(0,
                  cumsum(N_vOTUs)[-n()]) ) %>%
  mutate(xmax = xmin + N_vOTUs,
         ymin = 0,
         ymax = 1) 

# by the lifestyle type
by_lifestyle <- hq_by_recovery %>%
  left_join(ETOF_vOTUr %>% select(New_CID, lifestyle), by = c("Representative" = "New_CID")) %>%
  filter(!grepl('DB', vOTU_cluster_type)) %>%
  group_by(lifestyle, recoverable_by) %>%
  summarise(N_vOTUs = n(), .groups = "drop") %>%
  rename(source = recoverable_by) %>%
  mutate(source = factor(source, levels = c("MGS_QoL", "MGS_QoL+VLP_QoL", "VLP_QoL"), ordered = TRUE)) %>%
  mutate(lifestyle = factor(lifestyle, levels = c("Temperate", "Virulent", "Unknown"), ordered = TRUE)) %>%
  arrange(source, lifestyle) %>%
  mutate(xmin = c(0,
                  cumsum(N_vOTUs)[-n()]) ) %>%
  mutate(xmax = xmin + N_vOTUs,
         ymin = 0,
         ymax = 1) 
  

by_all <- ggplot() +
  geom_rect(data = overall, aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Part), color="black") +
  geom_rect(data = recovery, aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=hq_recovery_source), color="black") +
  geom_rect(data = by_genome, aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=genome), color="black") + 
  geom_rect(data = by_lifestyle, aes(xmin = xmin, xmax=xmax, ymin=ymin-3, ymax=ymax-3, fill=lifestyle), color="black") + 
  scale_fill_manual(labels = c("MGS unique", "Overlap", "VLP unique"),
                    values = c("NEXT_MGS"=MetBrewer::met.brewer("Kandinsky")[4],
                               "NEXT_VLP+NEXT_MGS"=MetBrewer::met.brewer("Kandinsky")[2],
                               "NEXT_VLP"=MetBrewer::met.brewer("Kandinsky")[1])) +
  labs(x = "Number of genomes", y = "", fill="Source") +
  # dependent recovery
  ggnewscale::new_scale_fill() +
  geom_rect(data = recovery, aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=hq_recovery_source), color="black") +
  scale_fill_manual(labels = c("MGS-dependent", "Dual source", "VLP-dependent"),
                    values = c("MGS_QoL"=MetBrewer::met.brewer("VanGogh1")[1],
                               "MGS_QoL+VLP_QoL"=MetBrewer::met.brewer("VanGogh1")[5],
                               "VLP_QoL"=MetBrewer::met.brewer("VanGogh1")[7])) +
  labs(fill="Recovery") +
  # Genome graph
  ggnewscale::new_scale_fill() +
  geom_rect(data = by_genome, aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=genome), color=NA) +
  scale_fill_manual(labels = c("dsDNA",  
                               "RNA", 
                               "ssDNA", "Unclassified"),
                    values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                               "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                               "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2],
                               "Unclassified"="darkgrey")) +
  geom_rect(data = recovery, aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2), fill = "transparent", color="black") +
  labs(fill="Genome") +
  # Lifestyle
  ggnewscale::new_scale_fill() +
  geom_rect(data = by_lifestyle, aes(xmin = xmin, xmax=xmax, ymin=ymin-3, ymax=ymax-3, fill=lifestyle), color=NA) +
  scale_fill_manual(labels = c("Temperate",  
                               "Virulent", 
                               "Unknown"),
                    values = c("Temperate"=MetBrewer::met.brewer("Thomas")[6],
                               "Virulent"=MetBrewer::met.brewer("Thomas")[5],
                               "Unknown"="darkgrey")) +
  geom_rect(data = recovery, aes(xmin = xmin, xmax=xmax, ymin=ymin-3, ymax=ymax-3), fill = "transparent", color="black") +
  labs(fill="Lifestyle   ") +
  scale_x_continuous(breaks = seq(0, 8500, by = 8500), limits = c(0, 8500), labels = abs) + 
  scale_y_continuous(limits = c(-3.2,1.5)) +
  theme_VLP_MGS_plots +
  ggpubr::geom_bracket(xmin = 0, xmax = 437, label = "Unique to MGS", y.position = 1.3, size = 0.55, label.size = 3, tip.length = 0.02) +
  ggpubr::geom_bracket(xmin = 4499, xmax = 8348, label = "Unique to VLP", y.position = 1.3, size=0.55, label.size = 3, tip.length = 0.02) +
  geom_text(data = overall, aes(x = (xmin + xmax)/2, y = (ymax + 0.1), label = paste(N_vOTUs)), size = 2.5, color = "black")

ggsave('05.PLOTS/05.VLP_MGS/Novel_recovery_genome_ls.png',
       by_all,  "png", width=20, height=12, units="cm", dpi = 300)
