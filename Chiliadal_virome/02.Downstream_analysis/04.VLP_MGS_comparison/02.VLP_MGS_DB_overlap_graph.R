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
ETOF_vOTUr <- read.table('06.CLEAN_DATA/Extended_ETOF_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

hq_by_recovery <- read.table('06.CLEAN_DATA/03.RESULTS_to_plot/HQ_vOTU_recoverable_by_source.txt', sep='\t', header = T)

recovery_summary <- read.table('06.CLEAN_DATA/03.RESULTS_to_plot/HQ_vOTU_recovery_stat.txt', sep='\t', header=T)

#############################################################
# 2. Analysis: novel vOTU by recovery source
#############################################################
# vOTU members sharing:
overall <- as.data.frame(table(ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality" & 
                                            !grepl('DB', ETOF_vOTUr$vOTU_cluster_type),]$vOTU_cluster_type))
colnames(overall) <- c('Part', 'N_vOTUs')

#overall$simple_new <- c('MGS_only', 'VLP_only', 'Overlap')
overall$Part <- factor(overall$Part, levels=c("NEXT_MGS", "NEXT_VLP+NEXT_MGS", "NEXT_VLP"), ordered = T)

overall <- overall[order(overall$Part),]

overall$xmin <- c(0,
                  cumsum(overall$N_vOTUs)[-nrow(overall)])
overall$xmax <- overall$xmin + overall$N_vOTUs

overall$ymin <- 0
overall$ymax <- 1


# vOTU HQ recovery:
recovery <- recovery_summary %>%
  filter(vOTU_cluster_type == "NEXT_VLP+NEXT_MGS") %>%
  transmute(
    hq_dep_VLP = HQ_VLP - HQ_MGS_VLP,
    hq_dep_MGS = HQ_MGS - HQ_MGS_VLP,
    dual_hq    = HQ_MGS_VLP
  ) %>%
  pivot_longer(
    everything(),
    names_to = "hq_recovery_source",
    values_to = "N_vOTUs"
  ) %>%
  bind_rows(
    recovery_summary %>%
      filter(vOTU_cluster_type %in% c("NEXT_VLP", "NEXT_MGS")) %>%
      transmute(
        hq_recovery_source = vOTU_cluster_type,
        N_vOTUs = N_HQ
      )
  )
recovery$hq_recovery_source <- factor(recovery$hq_recovery_source, levels=c('NEXT_MGS', 'hq_dep_MGS', 'dual_hq', 'hq_dep_VLP', 'NEXT_VLP'), ordered = T)
recovery <- recovery[order(recovery$hq_recovery_source),]

recovery$xmin <- c(0,
                  cumsum(recovery$N_vOTUs)[-nrow(recovery)])
recovery$xmax <- recovery$xmin + recovery$N_vOTUs

recovery$ymin <- 0
recovery$ymax <- 1

# lifestyle & genome annotation
interm <- ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality" & 
                       !grepl('DB', ETOF_vOTUr$vOTU_cluster_type), c("genome_simple", "POST_CHV_length", "lifestyle", "vOTU_cluster_type")]

# by genome type
by_genome <- as.data.frame(table(interm$genome_simple, interm$vOTU_cluster_type)) %>%
  rename(genome = Var1, source = Var2, N_vOTUs = Freq) %>%
  mutate(
    genome = as.character(genome),
    source = as.character(source),
    
    genome = ifelse(grepl("\\+", source), "NEXT_VLP+NEXT_MGS", genome),
    
    source = factor(source, levels = c("NEXT_MGS", "NEXT_VLP+NEXT_MGS", "NEXT_VLP"), ordered = TRUE),
    
    genome = factor(genome, levels =c('RNA', 'ssDNA', 'dsDNA', 'Unknown', 'NEXT_VLP+NEXT_MGS'), ordered = TRUE)
  ) %>%
  group_by(genome, source) %>%
  summarise(N_vOTUs = sum(N_vOTUs)) %>%
  arrange(source, genome) %>%
  mutate(
    genome = as.character(genome))

by_genome$xmin <- c(0,
                    cumsum(by_genome$N_vOTUs)[-nrow(by_genome)])

by_genome$xmax <- by_genome$xmin + by_genome$N_vOTUs

by_genome$ymin <- 0
by_genome$ymax <- 1

# by the lifestyle type
by_lifestyle <- as.data.frame(table(interm$lifestyle, interm$vOTU_cluster_type)) %>%
  rename(lifestyle = Var1, source = Var2, N_vOTUs = Freq) %>%
  mutate(
    lifestyle = as.character(lifestyle),
    source = as.character(source),
    lifestyle = ifelse(grepl("\\+", source), "NEXT_VLP+NEXT_MGS", lifestyle),
    #source = ifelse(grepl("\\+", source), "Overlap", source),
    source = factor(source, levels = c("NEXT_MGS", "NEXT_VLP+NEXT_MGS", "NEXT_VLP"), ordered = TRUE)
  ) %>%
  group_by(lifestyle, source) %>%
  summarise(N_vOTUs = sum(N_vOTUs)) %>%
  arrange(source)

by_lifestyle$xmin <- NA

by_lifestyle$xmin <- c(0,
                       cumsum(by_lifestyle$N_vOTUs)[-nrow(by_lifestyle)])

by_lifestyle$xmax <- by_lifestyle$xmin + by_lifestyle$N_vOTUs

by_lifestyle$ymin <- 0
by_lifestyle$ymax <- 1

# original coloring:

# ggplot() +
#   geom_rect(data = overall, aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Part), color="black") +
#   geom_rect(data = recovery, aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=hq_recovery_source), color="black") +
#   geom_rect(data = by_genome[by_genome$source=="NEXT_VLP+NEXT_MGS",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=genome), color="black") + 
#   geom_rect(data = by_lifestyle[by_lifestyle$source=="NEXT_VLP+NEXT_MGS",], aes(xmin = xmin, xmax=xmax, ymin=ymin-3, ymax=ymax-3, fill=lifestyle), color="black") + 
#   scale_fill_manual(labels = c("MGS unique", "Overlap", "VLP unique"),
#                     values = c("NEXT_MGS"=MetBrewer::met.brewer("VanGogh1")[1],
#                                "NEXT_VLP+NEXT_MGS"=MetBrewer::met.brewer("VanGogh1")[5],
#                                "NEXT_VLP"=MetBrewer::met.brewer("VanGogh1")[7])) +
#   labs(x = "Number of genomes", y = "", fill="Source") +
#   # dependent recovery
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = recovery[!recovery$hq_recovery_source %in% c('NEXT_MGS', 'NEXT_VLP'),], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=hq_recovery_source), color="black") +
#   scale_fill_manual(labels = c("MGS-dependent", "Dual source", "VLP-dependent"),
#                     values = c("hq_dep_MGS"=MetBrewer::met.brewer("VanGogh1")[2],
#                                "dual_hq"=MetBrewer::met.brewer("VanGogh1")[5],
#                                "hq_dep_VLP"=MetBrewer::met.brewer("VanGogh1")[6])) +
#   labs(fill="  Recovery") +
#   # Genome graph
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = by_genome[by_genome$source!="NEXT_VLP+NEXT_MGS",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=genome), color="black") +
#   scale_fill_manual(labels = c("dsDNA",  
#                                "RNA", 
#                                "ssDNA", "Unknown"),
#                     values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
#                                "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
#                                "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2],
#                                "Unknown"="darkgrey")) +
#   labs(fill="Genome") +
#   # Lifestyle
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = by_lifestyle[by_lifestyle$source!="NEXT_VLP+NEXT_MGS",], aes(xmin = xmin, xmax=xmax, ymin=ymin-3, ymax=ymax-3, fill=lifestyle), color="black") +
#   scale_fill_manual(labels = c("Temperate",  
#                                "Virulent"),
#                     values = c("Temperate"=MetBrewer::met.brewer("Thomas")[6],
#                                "Virulent"=MetBrewer::met.brewer("Thomas")[5])) +
#   labs(fill="Lifestyle") +
#   scale_x_continuous(breaks = seq(0, 9000, by = 9000), limits = c(0, 9000), labels = abs) + 
#   scale_y_continuous(limits = c(-3.2,1.5)) +
#   theme_VLP_MGS_plots +
#   ggpubr::geom_bracket(xmin = 0, xmax = 439, label = "Unique to MGS", y.position = 1.3, size = 0.55, label.size = 3, tip.length = 0.02) +
#   ggpubr::geom_bracket(xmin = 4535, xmax = 8580, label = "Unique to VLP", y.position = 1.3, size=0.55, label.size = 3, tip.length = 0.02) +
#   geom_text(data = overall, aes(x = (xmin + xmax)/2, y = (ymax + 0.1), label = paste(N_vOTUs)), size = 2.5, color = "black")
# 


exp_w_color <- ggplot() +
  geom_rect(data = overall, aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Part), color="black") +
  geom_rect(data = recovery, aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=hq_recovery_source), color="black") +
  geom_rect(data = by_genome[by_genome$source=="NEXT_VLP+NEXT_MGS",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=genome), color="black") + 
  geom_rect(data = by_lifestyle[by_lifestyle$source=="NEXT_VLP+NEXT_MGS",], aes(xmin = xmin, xmax=xmax, ymin=ymin-3, ymax=ymax-3, fill=lifestyle), color="black") + 
  scale_fill_manual(labels = c("MGS unique", "Overlap", "VLP unique"),
                    values = c("NEXT_MGS"=MetBrewer::met.brewer("VanGogh1")[1],
                               "NEXT_VLP+NEXT_MGS"=MetBrewer::met.brewer("VanGogh1")[5],
                               "NEXT_VLP"=MetBrewer::met.brewer("VanGogh1")[7])) +
  labs(x = "Number of genomes", y = "", fill="Source") +
  # dependent recovery
  ggnewscale::new_scale_fill() +
  geom_rect(data = recovery[!recovery$hq_recovery_source %in% c('NEXT_MGS', 'NEXT_VLP'),], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=hq_recovery_source), color="black") +
  scale_fill_manual(labels = c("MGS-dependent", "Dual source", "VLP-dependent"),
                    values = c("hq_dep_MGS"=MetBrewer::met.brewer("VanGogh1")[1],
                               "dual_hq"=MetBrewer::met.brewer("VanGogh1")[5],
                               "hq_dep_VLP"=MetBrewer::met.brewer("VanGogh1")[7])) +
  labs(fill="  Recovery") +
  # Genome graph
  ggnewscale::new_scale_fill() +
  geom_rect(data = by_genome[by_genome$source!="NEXT_VLP+NEXT_MGS",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=genome), color="black") +
  scale_fill_manual(labels = c("dsDNA",  
                               "RNA", 
                               "ssDNA", "Unknown"),
                    values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                               "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                               "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2],
                               "Unknown"="darkgrey")) +
  labs(fill="Genome") +
  # Lifestyle
  ggnewscale::new_scale_fill() +
  geom_rect(data = by_lifestyle[by_lifestyle$source!="NEXT_VLP+NEXT_MGS",], aes(xmin = xmin, xmax=xmax, ymin=ymin-3, ymax=ymax-3, fill=lifestyle), color="black") +
  scale_fill_manual(labels = c("Temperate",  
                               "Virulent"),
                    values = c("Temperate"=MetBrewer::met.brewer("Thomas")[6],
                               "Virulent"=MetBrewer::met.brewer("Thomas")[5])) +
  labs(fill="Lifestyle") +
  scale_x_continuous(breaks = seq(0, 9000, by = 9000), limits = c(0, 9000), labels = abs) + 
  scale_y_continuous(limits = c(-3.2,1.5)) +
  theme_VLP_MGS_plots +
  ggpubr::geom_bracket(xmin = 0, xmax = 439, label = "Unique to MGS", y.position = 1.3, size = 0.55, label.size = 3, tip.length = 0.02) +
  ggpubr::geom_bracket(xmin = 4535, xmax = 8580, label = "Unique to VLP", y.position = 1.3, size=0.55, label.size = 3, tip.length = 0.02) +
  geom_text(data = overall, aes(x = (xmin + xmax)/2, y = (ymax + 0.1), label = paste(N_vOTUs)), size = 2.5, color = "black")

ggsave('05.PLOTS/05.VLP_MGS/Novel_recovery_genome_ls.png',
       exp_w_color,  "png", width=20, height=12, units="cm", dpi = 300)
