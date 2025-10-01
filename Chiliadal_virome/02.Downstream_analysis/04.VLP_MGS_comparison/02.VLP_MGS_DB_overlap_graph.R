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
#############################################################
# 2. Analysis: plots merged DB by genome & source
#############################################################
overall <- as.data.frame(table(ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality",]$vOTU_cluster_type))
colnames(overall) <- c('Part', 'N_vOTUs')

overall$simple <- c('ignore', 'MGS_only', 'ignore',
                 'VLP_only', 'ignore', 'Overlap',
                 'Overlap')

overall$simple_new <- c('ignore', 'MGS_only', 'ignore',
                    'VLP_only', 'ignore', 'Overlap',
                    'ignore')

overall_wDB <- overall %>% 
  group_by(simple) %>% 
  summarise(N_vOTUs = sum(N_vOTUs))


overall_wDB <- overall_wDB[overall_wDB$simple!="ignore",]

overall_wDB$xmin <- c(0,
               cumsum(overall_wDB$N_vOTUs)[-nrow(overall_wDB)])

overall_wDB$xmax <- overall_wDB$xmin + overall_wDB$N_vOTUs

overall_wDB$ymin <- 0
overall_wDB$ymax <- 1

# for lifestyle & genome 
interm <- ETOF_vOTUr[ETOF_vOTUr$vOTU_cluster_type %in% c("NEXT_VLP", "NEXT_MGS", "NEXT_VLP+NEXT_MGS+DB_sum", "NEXT_VLP+NEXT_MGS") & 
                       ETOF_vOTUr$miuvig_quality=="High-quality", c("genome_simple", "POST_CHV_length", "lifestyle", "vOTU_cluster_type")]

# by genome type
by_genome <- as.data.frame(table(interm$genome_simple, interm$vOTU_cluster_type)) %>%
  rename(genome = Var1, source = Var2, N_vOTUs = Freq) %>%
  mutate(
    genome = as.character(genome),
    source = as.character(source),
    
    genome = ifelse(grepl("\\+", source), "Overlap", genome),
    source = ifelse(grepl("\\+", source), "Overlap", source),
    
    source = factor(source, levels = c("NEXT_MGS", "Overlap", "NEXT_VLP"), ordered = TRUE),
    
    genome = factor(genome, levels =c('RNA', 'ssDNA', 'dsDNA', 'Unknown', 'Overlap'), ordered = TRUE)
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

# as it ideally should look like:
by_genome_beauty <- by_genome
by_genome_beauty <- by_genome_beauty[-(2:4),]
by_genome_beauty[by_genome_beauty$source=="NEXT_MGS",]$genome <- "NEXT_MGS"
by_genome_beauty[by_genome_beauty$source=="NEXT_MGS",]$N_vOTUs <- 439
by_genome_beauty[by_genome_beauty$source=="NEXT_MGS",]$xmax <- 439

# by the lifestyle type

by_lifestyle <- as.data.frame(table(interm$lifestyle, interm$vOTU_cluster_type)) %>%
  rename(lifestyle = Var1, source = Var2, N_vOTUs = Freq) %>%
  mutate(
    lifestyle = as.character(lifestyle),
    source = as.character(source),
    lifestyle = ifelse(grepl("\\+", source), "Overlap", lifestyle),
    source = ifelse(grepl("\\+", source), "Overlap", source),
    source = factor(source, levels = c("NEXT_MGS", "Overlap", "NEXT_VLP"), ordered = TRUE)
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

######## merged plot

# RNA & ssDNA in MGS unique
# mergedplot <- ggplot() +
#   geom_rect(data = overall_wDB, aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=simple), color="black") +
#   geom_rect(data = by_genome[by_genome$source=="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=genome), color="black") + 
#   geom_rect(data = by_lifestyle[by_lifestyle$source=="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=lifestyle), color="black") + 
#   scale_fill_manual(labels = c("MGS unique", "Overlap", "VLP unique"),
#                     values = c("MGS_only"=MetBrewer::met.brewer("VanGogh2")[7],
#                                "Overlap"=MetBrewer::met.brewer("VanGogh2")[6],
#                                "VLP_only"=MetBrewer::met.brewer("VanGogh2")[5])) +
#   labs(x = "Number of genomes", y = "", fill="Source") +
#   # Genome graph
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = by_genome[by_genome$source!="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=genome), color="black") +
#   scale_fill_manual(labels = c("dsDNA",  
#                                "RNA", 
#                                "ssDNA", "Unknown"),
#                     values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
#                                "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
#                                "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2],
#                                "Unknown"=MetBrewer::met.brewer("VanGogh2")[9])) +
#   labs(x = "Number of genomes", y = "", fill="Genome") +
#   # Lifestyle graph
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = by_lifestyle[by_lifestyle$source!="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=lifestyle), color="black") +
#   scale_fill_manual(labels = c("Temperate",  
#                                "Virulent"),
#                     values = c("Temperate"=MetBrewer::met.brewer("VanGogh2")[8],
#                                "Virulent"=MetBrewer::met.brewer("VanGogh2")[3])) +
#   labs(x = "Number of genomes", y = "", fill="Lifestyle") +
#   
#   scale_x_continuous(breaks = seq(0, 14500, by = 14500), limits = c(0, 14500), labels = abs) + 
#   scale_y_continuous(limits = c(-2.2,1.5)) +
#   theme_VLP_MGS_plots +
#   ggpubr::geom_bracket(xmin = 0, xmax = 439, label = "Unique to MGS", y.position = 1.3, size = 0.55, label.size = 3, tip.length = 0.02) +
#   ggpubr::geom_bracket(xmin = 9985, xmax = 14030, label = "Unique to VLP", y.position = 1.3, size=0.55, label.size = 3, tip.length = 0.02) +
#   geom_text(data = overall_wDB, aes(x = (xmin + xmax)/2, y = (ymax + 0.1), label = paste(N_vOTUs)), size = 2.5, color = "black")
# 
# 
# ggsave('05.PLOTS/05.VLP_MGS/mergedplot.png',
#        mergedplot,  "png", width=20, height=6, units="cm", dpi = 300)

#### pretier in genome type for MGS:

beautymergedplot <- ggplot() +
  geom_rect(data = overall_wDB, aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=simple), color="black") +
  geom_rect(data = by_genome[by_genome$source=="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=genome), color="black") + 
  geom_rect(data = by_lifestyle[by_lifestyle$source=="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=lifestyle), color="black") + 
  scale_fill_manual(labels = c("MGS unique", "Overlap", "VLP unique"),
                    values = c("MGS_only"=MetBrewer::met.brewer("VanGogh2")[7],
                               "Overlap"=MetBrewer::met.brewer("VanGogh2")[6],
                               "VLP_only"=MetBrewer::met.brewer("VanGogh2")[5])) +
  labs(x = "Number of genomes", y = "", fill="Source") +
  # Genome graph
  ggnewscale::new_scale_fill() +
  geom_rect(data = by_genome_beauty[by_genome_beauty$source!="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=genome), color="black") +
  scale_fill_manual(labels = c("dsDNA",  
                               "RNA", 
                               "ssDNA", "Unknown"),
                    values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                               "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                               "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2],
                               "Unknown"=MetBrewer::met.brewer("VanGogh2")[9])) +
  labs(x = "Number of genomes", y = "", fill="Genome") +
  # beauty
  ggnewscale::new_scale_fill() +
  geom_rect(data = overall_wDB[overall_wDB$simple=="MGS_only",], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=simple), color="black") + 
  scale_fill_manual(labels=c("MGS unique"),values = c(
    "MGS_only"=MetBrewer::met.brewer("VanGogh2")[7])) +
  guides(fill = "none") +
  # Lifestyle graph
  ggnewscale::new_scale_fill() +
  geom_rect(data = by_lifestyle[by_lifestyle$source!="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=lifestyle), color="black") +
  scale_fill_manual(labels = c("Temperate",  
                               "Virulent"),
                    values = c("Temperate"=MetBrewer::met.brewer("VanGogh2")[8],
                               "Virulent"=MetBrewer::met.brewer("VanGogh2")[3])) +
  labs(x = "Number of genomes", y = "", fill="Lifestyle") +
  
  scale_x_continuous(breaks = seq(0, 14500, by = 14500), limits = c(0, 14500), labels = abs) + 
  scale_y_continuous(limits = c(-2.2,1.5)) +
  theme_VLP_MGS_plots +
  ggpubr::geom_bracket(xmin = 0, xmax = 439, label = "Unique to MGS", y.position = 1.3, size = 0.55, label.size = 3, tip.length = 0.02) +
  ggpubr::geom_bracket(xmin = 9985, xmax = 14030, label = "Unique to VLP", y.position = 1.3, size=0.55, label.size = 3, tip.length = 0.02) +
  geom_text(data = overall_wDB, aes(x = (xmin + xmax)/2, y = (ymax + 0.1), label = paste(N_vOTUs)), size = 2.5, color = "black")

ggsave('05.PLOTS/05.VLP_MGS/beautymergedplot.png',
       beautymergedplot,  "png", width=20, height=6, units="cm", dpi = 300)


#############################################################
# CONSIDERING ONLY NOVEL vOTU
#############################################################
overall_new <- overall[overall$simple_new!="ignore",]
overall_new$simple_new <- factor(overall_new$simple_new, 
                                 levels=c('MGS_only', 'Overlap', 'VLP_only'),
                                 ordered=T)

overall_new <- overall_new[order(overall_new$simple_new),]
overall_new$xmin <- NA
overall_new$xmin <- c(0,
                      cumsum(overall_new$N_vOTUs)[-3])

overall_new$xmax <- overall_new$xmin + overall_new$N_vOTUs

overall_new$ymin <- 0
overall_new$ymax <- 1

# by genome type
new_interm <- interm[interm$vOTU_cluster_type!="NEXT_VLP+NEXT_MGS+DB_sum",]

new_by_genome <- as.data.frame(table(new_interm$genome_simple, new_interm$vOTU_cluster_type)) %>%
  rename(genome = Var1, source = Var2, N_vOTUs = Freq) %>%
  mutate(
    genome = as.character(genome),
    source = as.character(source),
    
    genome = ifelse(grepl("\\+", source), "Overlap", genome),
    source = ifelse(grepl("\\+", source), "Overlap", source),
    
    source = factor(source, levels = c("NEXT_MGS", "Overlap", "NEXT_VLP"), ordered = TRUE),
    
    genome = factor(genome, levels =c('RNA', 'ssDNA', 'dsDNA', 'Unknown', 'Overlap'), ordered = TRUE)
  ) %>%
  group_by(genome, source) %>%
  summarise(N_vOTUs = sum(N_vOTUs)) %>%
  arrange(source, genome) %>%
  mutate(
    genome = as.character(genome))



new_by_genome$xmin <- c(0,
                    cumsum(new_by_genome$N_vOTUs)[-nrow(new_by_genome)])

new_by_genome$xmax <- new_by_genome$xmin + new_by_genome$N_vOTUs

new_by_genome$ymin <- 0
new_by_genome$ymax <- 1

# as it ideally should look like:
new_by_genome_beauty <- new_by_genome
new_by_genome_beauty <- new_by_genome_beauty[-(2:4),]
new_by_genome_beauty[new_by_genome_beauty$source=="NEXT_MGS",]$genome <- "NEXT_MGS"
new_by_genome_beauty[new_by_genome_beauty$source=="NEXT_MGS",]$N_vOTUs <- 439
new_by_genome_beauty[new_by_genome_beauty$source=="NEXT_MGS",]$xmax <- 439

# by the lifestyle type

new_by_lifestyle <- as.data.frame(table(new_interm$lifestyle, new_interm$vOTU_cluster_type)) %>%
  rename(lifestyle = Var1, source = Var2, N_vOTUs = Freq) %>%
  mutate(
    lifestyle = as.character(lifestyle),
    source = as.character(source),
    lifestyle = ifelse(grepl("\\+", source), "Overlap", lifestyle),
    source = ifelse(grepl("\\+", source), "Overlap", source),
    source = factor(source, levels = c("NEXT_MGS", "Overlap", "NEXT_VLP"), ordered = TRUE)
  ) %>%
  group_by(lifestyle, source) %>%
  summarise(N_vOTUs = sum(N_vOTUs)) %>%
  arrange(source)

new_by_lifestyle$xmin <- NA

new_by_lifestyle$xmin <- c(0,
                       cumsum(new_by_lifestyle$N_vOTUs)[-nrow(new_by_lifestyle)])

new_by_lifestyle$xmax <- new_by_lifestyle$xmin + new_by_lifestyle$N_vOTUs

new_by_lifestyle$ymin <- 0
new_by_lifestyle$ymax <- 1

######## merged plot

# RNA & ssDNA in MGS unique
# mergedplot_novel <- ggplot() +
#   geom_rect(data = overall_new, aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=simple), color="black") +
#   geom_rect(data = new_by_genome[new_by_genome$source=="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=genome), color="black") + 
#   geom_rect(data = new_by_lifestyle[new_by_lifestyle$source=="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=lifestyle), color="black") + 
#   scale_fill_manual(labels = c("MGS unique", "Overlap", "VLP unique"),
#                     values = c("MGS_only"=MetBrewer::met.brewer("VanGogh2")[7],
#                                "Overlap"=MetBrewer::met.brewer("VanGogh2")[6],
#                                "VLP_only"=MetBrewer::met.brewer("VanGogh2")[5])) +
#   labs(x = "Number of genomes", y = "", fill="Source") +
#   # Genome graph
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = new_by_genome[new_by_genome$source!="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=genome), color="black") +
#   scale_fill_manual(labels = c("dsDNA",  
#                                "RNA", 
#                                "ssDNA", "Unknown"),
#                     values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
#                                "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
#                                "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2],
#                                "Unknown"=MetBrewer::met.brewer("VanGogh2")[9])) +
#   labs(x = "Number of genomes", y = "", fill="Genome") +
#   # Lifestyle graph
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = new_by_lifestyle[new_by_lifestyle$source!="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=lifestyle), color="black") +
#   scale_fill_manual(labels = c("Temperate",  
#                                "Virulent"),
#                     values = c("Temperate"=MetBrewer::met.brewer("VanGogh2")[8],
#                                "Virulent"=MetBrewer::met.brewer("VanGogh2")[3])) +
#   labs(x = "Number of genomes", y = "", fill="Lifestyle") +
#   
#   scale_x_continuous(breaks = seq(0, 14500, by = 9000), limits = c(0, 9000), labels = abs) + 
#   scale_y_continuous(limits = c(-2.2,1.5)) +
#   theme_minimal() +
#   theme(legend.position = "right", 
#         legend.text = element_text(size=8), 
#         legend.title.position = "top",
#         legend.title = element_text(size=10, hjust=0.5),
#         axis.title.x = element_text(size=8, face="bold"), 
#         legend.key.size = unit(0.8,"line"),
#         legend.margin=margin(0,0,0,0),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.title.y=element_blank()
#   ) +
#   ggpubr::geom_bracket(xmin = 0, xmax = 439, label = "Unique to MGS", y.position = 1.3, size = 0.55, label.size = 3, tip.length = 0.02) +
#   ggpubr::geom_bracket(xmin = 4535, xmax = 8580, label = "Unique to VLP", y.position = 1.3, size=0.55, label.size = 3, tip.length = 0.02) +
#   geom_text(data = overall_new, aes(x = (xmin + xmax)/2, y = (ymax + 0.1), label = paste(N_vOTUs)), size = 2.5, color = "black")
# 
# 
# ggsave('05.PLOTS/05.VLP_MGS/mergedplot_novel_only.png',
#        mergedplot_novel,  "png", width=20, height=6, units="cm", dpi = 300)


# with beautified Genome part:
mergedplot_novel_beauty <- ggplot() +
  geom_rect(data = overall_new, aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=simple), color="black") +
  geom_rect(data = new_by_genome[new_by_genome$source=="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=genome), color="black") + 
  geom_rect(data = new_by_lifestyle[new_by_lifestyle$source=="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=lifestyle), color="black") + 
  scale_fill_manual(labels = c("MGS unique", "Overlap", "VLP unique"),
                    values = c("MGS_only"=MetBrewer::met.brewer("VanGogh2")[7],
                               "Overlap"=MetBrewer::met.brewer("VanGogh2")[6],
                               "VLP_only"=MetBrewer::met.brewer("VanGogh2")[5])) +
  labs(x = "Number of genomes", y = "", fill="Source") +
  # Genome graph
  ggnewscale::new_scale_fill() +
  geom_rect(data = new_by_genome[new_by_genome$source!="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=genome), color="black") +
  scale_fill_manual(labels = c("dsDNA",  
                               "RNA", 
                               "ssDNA", "Unknown"),
                    values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                               "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                               "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2],
                               "Unknown"=MetBrewer::met.brewer("VanGogh2")[9])) +
  labs(x = "Number of genomes", y = "", fill="Genome") +
  # beauty
  ggnewscale::new_scale_fill() +
  geom_rect(data = overall_new[overall_new$simple=="MGS_only",], aes(xmin = xmin, xmax=xmax, ymin=ymin-1, ymax=ymax-1, fill=simple), color="black") + 
  scale_fill_manual(labels=c("MGS unique"),values = c(
    "MGS_only"=MetBrewer::met.brewer("VanGogh2")[7])) +
  guides(fill = "none") +
  # Lifestyle graph
  ggnewscale::new_scale_fill() +
  geom_rect(data = new_by_lifestyle[new_by_lifestyle$source!="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin-2, ymax=ymax-2, fill=lifestyle), color="black") +
  scale_fill_manual(labels = c("Temperate",  
                               "Virulent"),
                    values = c("Temperate"=MetBrewer::met.brewer("VanGogh2")[8],
                               "Virulent"=MetBrewer::met.brewer("VanGogh2")[3])) +
  labs(x = "Number of genomes", y = "", fill="Lifestyle") +
  
  scale_x_continuous(breaks = seq(0, 14500, by = 9000), limits = c(0, 9000), labels = abs) + 
  scale_y_continuous(limits = c(-2.2,1.5)) +
  theme_minimal() +
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
  ) +
  ggpubr::geom_bracket(xmin = 0, xmax = 439, label = "Unique to MGS", y.position = 1.3, size = 0.55, label.size = 3, tip.length = 0.02) +
  ggpubr::geom_bracket(xmin = 4535, xmax = 8580, label = "Unique to VLP", y.position = 1.3, size=0.55, label.size = 3, tip.length = 0.02) +
  geom_text(data = overall_new, aes(x = (xmin + xmax)/2, y = (ymax + 0.1), label = paste(N_vOTUs)), size = 2.5, color = "black")

ggsave('05.PLOTS/05.VLP_MGS/mergedplot_novel_only_beautified.png',
       mergedplot_novel_beauty,  "png", width=20, height=6, units="cm", dpi = 300)

#############################################################
# Individual graph for overall overlap & unique for VLP
# and MGS comparison; redundant w merged one
# Includes DB-shared vOTUs
#############################################################
# p_simple_no_DB <- ggplot() +
#   geom_rect(data = overall_wDB, aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=simple), color="black") +
#   scale_fill_manual(labels = c("MGS unique", "Overlap", "VLP unique"),
#                     values = c("MGS_only"=MetBrewer::met.brewer("VanGogh2")[7],
#                                "Overlap"=MetBrewer::met.brewer("VanGogh2")[6],
#                                "VLP_only"=MetBrewer::met.brewer("VanGogh2")[5])) +
#   scale_x_continuous(breaks = seq(0, 14500, by = 14500), limits = c(0, 14500), labels = abs) +
#   scale_y_continuous(limits = c(-0.2,1.5)) +
#   theme_VLP_MGS_plots+
#   labs(x = "Number of genomes", y = "", fill="Source") + 
#   ggpubr::geom_bracket(xmin = 0, xmax = 439, label = "Unique to MGS", y.position = 1.1, size = 0.55, label.size = 3, tip.length = 0.02) +
#   ggpubr::geom_bracket(xmin = 9985, xmax = 14030, label = "Unique to VLP", y.position = 1.1, size=0.55, label.size = 3, tip.length = 0.02) +
#   geom_text(data = overall_wDB, aes(x = (xmin + xmax)/2, y = (ymin - 0.1), label = paste(N_vOTUs)), size = 2.5, color = "black")
# 
# ggsave('05.PLOTS/05.VLP_MGS/Simple_overlap_no_DB_no_partition.png',
#        p_simple_no_DB,  "png", width=20, height=6, units="cm", dpi = 300)

#############################################################
# Individual graph for overlap & unique for VLP and MGS 
# comparison; redundant w merged one; GENOME TYPE
# Includes DB-shared vOTUs
#############################################################
# p_genome_no_DB <- ggplot() +
#   geom_rect(data = by_genome[by_genome$source!="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=genome), color="black") +
#   scale_fill_manual(labels = c("dsDNA",  
#                                "RNA", 
#                                "ssDNA", "Unknown"),
#                     values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
#                                "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
#                                "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2],
#                                "Unknown"=MetBrewer::met.brewer("VanGogh2")[8])) +
#   labs(x = "Number of genomes", y = "", fill="Genome") + 
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = by_genome[by_genome$source=="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=genome), color="black") + 
#   scale_fill_manual(labels = c("Overlap"),
#                     values = c(
#                       "Overlap"=MetBrewer::met.brewer("VanGogh2")[6])) +
#   labs(x = "Number of genomes", y = "", fill="Source") +
#   scale_x_continuous(breaks = seq(0, 14500, by = 14500), limits = c(0, 14500), labels = abs) +
#   scale_y_continuous(limits = c(-0.2,1.5)) +
#   theme_VLP_MGS_plots +
#   ggpubr::geom_bracket(xmin = 0, xmax = 439, label = "Unique to MGS", y.position = 1.1, size = 0.55, label.size = 3, tip.length = 0.02) +
#   ggpubr::geom_bracket(xmin = 9985, xmax = 14030, label = "Unique to VLP", y.position = 1.1, size=0.55, label.size = 3, tip.length = 0.02) +
#   geom_text(data = overall_wDB, aes(x = (xmin + xmax)/2, y = (ymin - 0.1), label = paste(N_vOTUs)), size = 2.5, color = "black")
# 
# ggsave('05.PLOTS/05.VLP_MGS/Overlap_no_DB_genome_type.png',
#        p_genome_no_DB,  "png", width=20, height=6, units="cm", dpi = 300)

#############################################################
# Individual graph for overlap & unique for VLP and MGS 
# comparison; redundant w merged one; GENOME TYPE
# Includes DB-shared vOTUs
# IGNORES GENOME TYPE FOR MGS UNIQUE
#############################################################
# beautifiedp_genome_no_DB <- ggplot() +
#   geom_rect(data = by_genome_beauty[by_genome_beauty$source=="NEXT_VLP",], aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=genome), color="black") +
#   scale_fill_manual(labels = c("dsDNA",  
#                                "RNA", 
#                                "ssDNA", "Unknown"),
#                     values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
#                                "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
#                                "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2],
#                                "Unknown"=MetBrewer::met.brewer("VanGogh2")[8])) +
#   labs(x = "Number of genomes", y = "", fill="Genome") + 
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = by_genome_beauty[by_genome_beauty$source!="NEXT_VLP",], aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=genome), color="black") + 
#   scale_fill_manual(labels = c("MGS unique","Overlap"),
#                     values = c("NEXT_MGS"=MetBrewer::met.brewer("VanGogh2")[7],
#                                "Overlap"=MetBrewer::met.brewer("VanGogh2")[6])) +
#   labs(x = "Number of genomes", y = "", fill="Source") +
#   scale_x_continuous(breaks = seq(0, 14500, by = 14500), limits = c(0, 14500), labels = abs) +
#   scale_y_continuous(limits = c(-0.2,1.5)) +
#   theme_VLP_MGS_plots +
#   ggpubr::geom_bracket(xmin = 0, xmax = 439, label = "Unique to MGS", y.position = 1.1, size = 0.55, label.size = 3, tip.length = 0.02) +
#   ggpubr::geom_bracket(xmin = 9985, xmax = 14030, label = "Unique to VLP", y.position = 1.1, size=0.55, label.size = 3, tip.length = 0.02) +
#   geom_text(data = overall_wDB, aes(x = (xmin + xmax)/2, y = (ymin - 0.1), label = paste(N_vOTUs)), size = 2.5, color = "black")
# 
# 
# ggsave('05.PLOTS/05.VLP_MGS/Beautified_Overlap_no_DB_genome_type.png',
#        beautifiedp_genome_no_DB,  "png", width=20, height=6, units="cm", dpi = 300)

#############################################################
# Individual graph for overlap & unique for VLP and MGS 
# comparison; redundant w merged one; LIFESTYLE
# Includes DB-shared vOTUs
#############################################################
# p_ls_no_DB <- ggplot() +
#   geom_rect(data = by_lifestyle[by_lifestyle$source!="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=lifestyle), color="black") +
#   scale_fill_manual(labels = c("Temperate",  
#                                "Virulent"),
#                     values = c("Temperate"=MetBrewer::met.brewer("VanGogh2")[8],
#                                "Virulent"=MetBrewer::met.brewer("VanGogh2")[4])) +
#   labs(x = "Number of genomes", y = "", fill="Lifestyle") + 
#   ggnewscale::new_scale_fill() +
#   geom_rect(data = by_lifestyle[by_lifestyle$source=="Overlap",], aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=lifestyle), color="black") + 
#   scale_fill_manual(labels = c("Overlap"),
#                     values = c(
#                       "Overlap"=MetBrewer::met.brewer("VanGogh2")[6])) +
#   labs(x = "Number of genomes", y = "", fill="Source") +
#   scale_x_continuous(breaks = seq(0, 14500, by = 14500), limits = c(0, 14500), labels = abs) +
#   scale_y_continuous(limits = c(-0.2,1.5)) +
#   theme_VLP_MGS_plots +
#   ggpubr::geom_bracket(xmin = 0, xmax = 439, label = "Unique to MGS", y.position = 1.1, size = 0.55, label.size = 3, tip.length = 0.02) +
#   ggpubr::geom_bracket(xmin = 9985, xmax = 14030, label = "Unique to VLP", y.position = 1.1, size=0.55, label.size = 3, tip.length = 0.02) +
#   geom_text(data = overall_wDB, aes(x = (xmin + xmax)/2, y = (ymin - 0.1), label = paste(N_vOTUs)), size = 2.5, color = "black")
# 
# ggsave('05.PLOTS/05.VLP_MGS/Overlap_no_DB_lifestyle.png',
#        p_ls_no_DB,  "png", width=20, height=6, units="cm", dpi = 300)

