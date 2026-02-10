
setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")
#############################################################
# Generating supplementary figures, Supplementary Figure 1
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
library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)
library(MetBrewer)
library(ggridges)

vangogh1 <- met.brewer("VanGogh2")
#############################################################
# 2. Load Input Data
#############################################################
ETOF <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

votu_clustering <- read.table('06.CLEAN_DATA/NEXT_viral_clusters_MGS_VLP_long_format.txt', sep='\t', header=T)
votu_clustering <- votu_clustering %>%
  filter(Representative %in% ETOF$New_CID)

ETOF_all <- read.table('06.CLEAN_DATA/VLP_MGS_ETOF_full_rep.txt', sep='\t', header = T)
ETOF_only <- ETOF_all %>%
  filter(grepl('NEXT_', ETOF_all$New_CID)) %>%
  filter(New_CID %in% votu_clustering$Cluster_member) %>%
  left_join(votu_clustering, by = c("New_CID" = "Cluster_member")) %>%
  left_join(ETOF %>% select(New_CID, genome), by = c("Representative" = "New_CID")) %>%
  mutate( method = ifelse(grepl('NEXT_V', New_CID), 'VLP', 'MGS') )

esmeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)

#############################################################
# 3.1 Analysis: Supplementary figure 1a
#############################################################

length_by_genome <- ggplot(ETOF, aes(POST_CHV_length, fill = genome)) + 
  geom_histogram() + 
  facet_wrap(~ genome,nrow = 4, scales = "free_y" ) + 
  theme_bw() +
  labs(y = "N vOTU representatives", x = "Sequence length", fill = "Genome type") +
  theme(strip.background = element_rect(NA)) +
  scale_x_log10(breaks = c(3000, 30000, 300000), labels = c("3000", "30000", "300000")) +
  scale_fill_manual(values = c(vangogh1[1],
                               vangogh1[3],
                               vangogh1[5],
                               vangogh1[8]))

ggsave('05.PLOTS/03.GENERAL_STATS/Length_by_genome_vOTUrs.png',
       length_by_genome,  "png", width=10, height=10, units="cm", dpi = 300)

#############################################################
# 3.2 Analysis: Supplementary figure 1b
#############################################################

gen1b <- ETOF %>%
  group_by(genome) %>%
  summarise(sum=n()) %>%
  ggplot(aes(x = "by genome", y = sum, fill = genome)) +
  geom_bar(position = "stack",
           stat = "identity") +
  labs(y = "N vOTUs", fill = "Genome type") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        legend.key.size=unit(0.7, "line"),
        legend.key.spacing.y = unit(1, 'pt')) +
  guides(fill=guide_legend(ncol=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5)) +
  scale_fill_manual(values = met.brewer("Archambault"))
  
  

tax1b <- ETOF %>%
  separate(
    col    = tax_ictv_aai,
    into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  mutate(Class_aggr = ifelse( Class %in% c('Caudoviricetes', 'Malgrandaviricetes', 'Unclassified'), Class, "Other")) %>%
  group_by(Class_aggr) %>%
  summarise(sum=n()) %>%
  ggplot(aes(x = "by Class", y = sum, fill = Class_aggr)) +
  geom_bar(position = "stack",
           stat = "identity") +
  labs(fill = "Taxonomy (Class)") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        legend.key.size=unit(0.7, "line"),
        legend.key.spacing.y = unit(1, 'pt'))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5)) +
  scale_fill_manual(values = met.brewer("Austria"))
  

lfs1b <- ETOF %>%
  group_by(lifestyle) %>%
  summarise(sum=n()) %>%
  mutate(lifestyle = factor(lifestyle, 
                            levels = c('Temperate', 'Virulent', 'Unknown'), 
                            ordered = T)) %>%
  ggplot(aes(x = "by lifestyle", y = sum, fill = lifestyle)) +
  geom_bar(position = "stack",
           stat = "identity") +
  labs(fill = "Lifestyle") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        legend.key.size=unit(0.7, "line"),
        legend.key.spacing.y = unit(1, 'pt'))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5)) +
  scale_fill_manual(values = met.brewer("Cross"))
  
mvg1b <- ETOF %>%
  group_by(miuvig_quality) %>%
  summarise(sum=n()) %>%
  mutate(miuvig_quality = factor(miuvig_quality, 
                                 levels = c('High-quality', 'Genome-fragment'), 
                                 ordered = T)) %>%
  ggplot(aes(x = "by MIUViG\nquality", y = sum, fill = miuvig_quality)) +
  geom_bar(position = "stack",
           stat = "identity") +
  labs(fill = "MIUViG cateogory") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        legend.key.size=unit(0.7, "line"),
        legend.key.spacing.y = unit(1, 'pt'))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5)) +
  scale_fill_manual(values = met.brewer("Kandinsky"))

ckv1b <- ETOF %>%
  group_by(checkv_quality) %>%
  summarise(sum=n()) %>%
  mutate(checkv_quality = factor(checkv_quality, 
                                 levels = c('Complete', 'High-quality', 'Medium-quality', 'Low-quality', 'Not-determined'), 
                                 ordered = T)) %>%
  ggplot(aes(x = "by CheckV\nquality", y = sum, fill = checkv_quality)) +
  geom_bar(position = "stack",
           stat = "identity") +
  labs(fill = "CheckV cateogory") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        legend.key.size=unit(0.7, "line"),
        legend.key.spacing.y = unit(1, 'pt'))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5)) +
  scale_fill_manual(values = met.brewer("Johnson"))

FSUP1b <- (gen1b | tax1b | lfs1b | mvg1b | ckv1b) + plot_layout(guides = 'collect')

ggsave('05.PLOTS/03.GENERAL_STATS/ETOF_vOTU_stat.png',
       FSUP1b,  "png", width=14, height=13, units="cm", dpi = 300)

#############################################################
# 3.3 Analysis: Supplementary figure 1c
#############################################################
dat_text <- data.frame(
  label = c("Cohen's d = 1.1\np-value = 6.9e-122", 
            "Cohen's d = 1.7\np-value = 4.5e-292", 
            "Cohen's d = 0.1\np-value = 5.8e-3"),
  Type_reads   = c("Raw reads", "Human reads", "Clean reads"),
  x     = c(7.3e+07, 4e+07, 6e+07),
  y     = c(2.5, 2.5, 2.5)
)

reads1c <- esmeta %>%
  filter(raw_reads != max(raw_reads)) %>% # there is one outlier that collapses distribution
  select(raw_reads, human_reads, clean_reads, seq_type) %>%
  rename(`Raw reads` = raw_reads, `Human reads` = human_reads, `Clean reads` = clean_reads) %>%
  pivot_longer(cols = c(`Raw reads`, `Human reads`, `Clean reads`),
    names_to = "Type_reads",
    values_to = "N_reads") %>%
  ggplot(aes(x = N_reads, y = seq_type, fill = seq_type)) +
  geom_density_ridges() +
  facet_wrap(~Type_reads, nrow = 3, scales = "free", strip.position = "left") +
  scale_y_discrete(expand = expand_scale(mult = c(0.00, 1))) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(NA)) +
  labs(x = "Number of reads",
       fill="Sequencing") + 
  scale_fill_manual(values = c(met.brewer("Kandinsky")[1],
                               met.brewer("Kandinsky")[2])) +
  geom_text(inherit.aes = FALSE,
    data    = dat_text,
    mapping = aes(x = x, y = y, label = label), 
    size = 3
  )


contigs1c_text <- data.frame(
  label = c("Cohen's d = 0.5\np-value = 1.1e-46", 
            "Cohen's d = 0.7\np-value = 6.6e-86"),
  Type_contigs   = c("All contigs", "Contigs > 1kbp"),
  x     = c(4e+05, 5.0e+04),
  y     = c(3, 3)
)

contigs1c <- esmeta %>%
  filter(raw_reads != max(raw_reads)) %>%
  select(contigs_0_bp, contigs_1000_bp, seq_type) %>%
  rename(`All contigs` = contigs_0_bp, `Contigs > 1kbp` = contigs_1000_bp) %>%
  pivot_longer(
    cols = c(`All contigs`, `Contigs > 1kbp`),
    names_to = "Type_contigs",
    values_to = "N_contigs") %>%
  ggplot(aes(x = N_contigs, y = seq_type, fill = seq_type)) +
  geom_density_ridges() +
  facet_wrap(~Type_contigs, nrow = 2, scales = "free_x", strip.position = "left") +
  scale_y_discrete(expand = expand_scale(mult = c(0.00, 1.9))) +
  theme_bw() +  
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(NA)) +
  labs(x = "Number of contigs",
       fill="Sequencing") +
  scale_x_continuous(labels = scales::scientific) + 
  scale_fill_manual(values = c(met.brewer("Kandinsky")[1],
                               met.brewer("Kandinsky")[2])) +
  geom_text(inherit.aes = FALSE,
            data    = contigs1c_text,
            mapping = aes(x = x, y = y, label = label), 
            size = 3
  )


sc_en_text <- data.frame(
  label = c("Cohen's d = 2.7\np-value = 0"),
  calc_type   = c("Virus enrichment"),
  x     = 60,
  y     = 2.5
)

sc_en1c <- esmeta %>%
  select(sc_enrichment, seq_type) %>%
  mutate(calc_type = "Virus enrichment") %>%
  ggplot(aes(x = sc_enrichment, y = seq_type, fill = seq_type)) +
  geom_density_ridges() +
  facet_wrap(~calc_type, nrow = 1, scales = "free_x", strip.position = "left") +
  scale_y_discrete(expand = expand_scale(mult = c(0.00, 1))) +
  theme_bw() +  
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(NA)) +
  labs(x = "ViromeQC enrichment score",
       fill="Sequencing") + 
  scale_fill_manual(values = c(met.brewer("Kandinsky")[1],
                               met.brewer("Kandinsky")[2])) +
  geom_text(inherit.aes = FALSE,
            data    = sc_en_text,
            mapping = aes(x = x, y = y, label = label), 
            size = 3
  )

csc <- (contigs1c/sc_en1c) + plot_layout(heights = c(7,3))

FSUP1c <- (reads1c | csc) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

ggsave('05.PLOTS/03.GENERAL_STATS/Esmeta_VLP_MGS_stat.png',
       FSUP1c,  "png", width=16, height=13, units="cm", dpi = 300)


#############################################################
# 3.4 Analysis: Supplementary figure 1d (abandoned atm)
#############################################################
# 
# busya <- ETOF_only %>%
#   filter(miuvig_quality == "High-quality" & method %in% c('MGS', 'VLP')) %>%
#   group_by(method) %>%
#   summarise(n())
# 
# 
# ### also make a ridgeline? otherwise does not look good? split by genome? does it actually make sense to make it?? does it add much?
# ETOF_only %>%
#   filter(miuvig_quality == "High-quality" & method %in% c('MGS', 'VLP')) %>%
#   ggplot(aes(x=POST_CHV_length, fill = method)) +
#   geom_histogram()
# 
# boxplot(log10(ETOF_only$POST_CHV_length) ~ ETOF_only$method)
