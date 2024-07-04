## Code description
#######################################################################################################################################
## Script for the negative control sharing paper
## 
#######################################################################################################################################

## Load libraries
#######################################################################################################################################
.libPaths("/home1/p309176/R_libraries")
library(readr)
library(vegan)
library(dplyr)
library(stringr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
setwd("/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/plots")
#######################################################################################################################################

## Load the  data
#######################################################################################################################################
RPKM <- read.delim("../RPKM_counts_VLP_NCP.txt")
dim(RPKM)  # 306275   1301

extended_tof <- read.delim('../../VIR_DB/virus_contigs/MERGED_Extended_TOF_NCP')
row.names(extended_tof) <- extended_tof$New_CID
dim(extended_tof)  # 971583     41
extended_tof <- extended_tof %>%
  mutate(virus_group = ifelse(grepl("Duplodnaviria", taxonomy), "dsDNA", "Unclassified"),
         virus_group = ifelse(grepl("Inoviridae", taxonomy), "ssDNA", virus_group),
         virus_group = ifelse(grepl("Microviridae", taxonomy), "ssDNA", virus_group),
         virus_group = ifelse(grepl("Cressdnaviricota", taxonomy), "ssDNA", virus_group),
         virus_group = ifelse(grepl("Cossaviricota", taxonomy), "dsDNA", virus_group),
         virus_group = ifelse(grepl("Riboviria", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Nodaviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Tolivirales", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Retroviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Cystoviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Picornaviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Astroviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Caliciviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Varidnaviria", taxonomy), "dsDNA", virus_group),
         virus_group = ifelse(grepl("Anelloviridae", taxonomy), "ssDNA", virus_group),
         virus_group = ifelse(grepl("Portogloboviridae", taxonomy), "dsDNA", virus_group),
         virus_group = ifelse(grepl("Bicaudaviridae", taxonomy), "dsDNA", virus_group),
         host_group = ifelse(grepl("Caudoviricetes", taxonomy), "bacteria", "Unclassified"),
         host_group = ifelse(grepl("Herviviricetes", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Inoviridae", taxonomy), "bacteria", host_group),
         host_group = ifelse(grepl("Microviridae", taxonomy), "bacteria", host_group),
         host_group = ifelse(grepl("Genomoviridae", taxonomy), "plant", host_group),
         host_group = ifelse(grepl("Circoviridae", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Nanoviridae", taxonomy), "plant", host_group),
         host_group = ifelse(grepl("Geminiviridae", taxonomy), "plant", host_group),
         host_group = ifelse(grepl("Smacoviridae", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Cossaviricota", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Nodaviridae", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Tombusviridae", taxonomy), "plant", host_group),
         host_group = ifelse(grepl("Retroviridae", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Cystoviridae", taxonomy), "bacteria", host_group),
         host_group = ifelse(grepl("Picornaviridae", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Astroviridae", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Caliciviridae", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Phycodnaviridae", taxonomy), "plant", host_group),
         host_group = ifelse(grepl("Mimiviridae", taxonomy), "animal(giant_virus)", host_group),
         host_group = ifelse(grepl("Adenoviridae", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Corticoviridae", taxonomy), "bacteria", host_group),
         host_group = ifelse(grepl("Iridoviridae", taxonomy), "bacteria", host_group),
         host_group = ifelse(grepl("Lavidaviridae", taxonomy), "viral phages", host_group),
         host_group = ifelse(grepl("Adintoviridae", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Tectiviridae", taxonomy), "bacteria", host_group),
         host_group = ifelse(grepl("Sphaerolipoviridae", taxonomy), "archaea", host_group),
         host_group = ifelse(grepl("Autolykiviridae", taxonomy), "bacteria", host_group),
         host_group = ifelse(grepl("Marseilleviridae", taxonomy), "animal(giant_virus)", host_group),
         host_group = ifelse(grepl("Poxviridae", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Anelloviridae", taxonomy), "animal", host_group),
         host_group = ifelse(grepl("Portogloboviridae", taxonomy), "archaea", host_group),
         host_group = ifelse(grepl("Bicaudaviridae", taxonomy), "archaea", host_group)
  )


meta_all_with_qc_curated <- as.data.frame(read_tsv('../../metadata_with_qc_NCP.tsv'))
dim(meta_all_with_qc_curated)  # 1376   28

host_prediction <- read.csv('../../VIR_DB/host_prediction_w_neg_der95_NCP/results/MERGED_Host_prediction_to_genus_m90_v2.csv')
dim(host_prediction)  # 216679      5 

filtered_host_prediction <- host_prediction %>%
  group_by(Virus) %>%
  slice_max(order_by = Confidence.score, with_ties = FALSE) %>%
  ungroup()
filtered_host_prediction <- as.data.frame(filtered_host_prediction)
dim(filtered_host_prediction)  # 188122      5
#######################################################################################################################################

##############################################################################################################################################################################################################################################################################


# THIS BLOCK IS DEDICATED TO EXPLORING THE RESULTS, WHEN THE RPKM TABLE IS NOT FILTERED, HENCE ALL vOTUs INCLUDED


##############################################################################################################################################################################################################################################################################

## Preparing metadata
#######################################################################################################################################
meta_working <- meta_all_with_qc_curated[meta_all_with_qc_curated$Sample_name %in% colnames(RPKM), ]
#######################################################################################################################################

## Checking alpha diversity
#######################################################################################################################################
row.names(meta_working) <- meta_working$Sample_name
diversity_df <- as.data.frame(diversity(as.data.frame(t(RPKM)), index='shannon'))
meta_working <- merge(meta_working, diversity_df, by="row.names")
colnames(meta_working)[colnames(meta_working) == "diversity(as.data.frame(t(RPKM)), index = \"shannon\")"] <- "diversity"
meta_working$Row.names <- NULL
meta_working <- meta_working %>%
  mutate(ncvssample = ifelse(Type == "Neg_ctrl", "NC", "SAMPLE"),
         timepoint_type = ifelse(Type == "Neg_ctrl", "NC", NA),
         timepoint_type = ifelse(Type == "Infant" & Timepoint %in% c("M0", "M1", "M2", "M3", "M4"), "Infant (age < 5 months)", timepoint_type),
         timepoint_type = ifelse(Type == "Infant" & Timepoint %in% c("M6", "M12", "Y2-5"), "Infant (age > 5 months)", timepoint_type),
         timepoint_type = ifelse(Type == "Mother", "Mother", timepoint_type))


pdf('alpha_diversity_comparison.pdf', width = 12/2.54, height = 8/2.54)
ggplot(meta_working, aes(x=ncvssample, y=diversity)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(color = timepoint_type), size=0.1) +
  facet_grid(. ~ cohort) +
  scale_colour_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#003049")) +
  labs(y = "Shannon diversity") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))
dev.off()
#######################################################################################################################################

## Checking read number -> INCLUDED & DONE
#######################################################################################################################################

pdf('raw_read_comparison.pdf', width = 12/2.54, height = 8/2.54)
ggplot(meta_working, aes(x=ncvssample, y=in_reads_comb)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  scale_y_log10() +
  labs(y = "Raw reads combined") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))
dev.off()

pdf('clean_read_comparison.pdf', width = 12/2.54, height = 8/2.54)
ggplot(meta_working, aes(x=ncvssample, y=clean_reads_comb)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  scale_y_log10() +
  labs(y = "Clean reads combined") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))
dev.off()


#######################################################################################################################################

## Analysis for the vOTUs composition in NC vs samples: richness and ratio
#######################################################################################################################################

extended_tof <- extended_tof %>%
  mutate(checkv_quality_and_plasmid = ifelse(checkv_quality %in% c("Complete", "High-quality", "Medium-quality"), "CHM", NA),
         checkv_quality_and_plasmid = ifelse(checkv_quality %in% c("Not-determined", "Low-quality"), "LU", checkv_quality_and_plasmid),
         checkv_quality_and_plasmid = ifelse(plasmid == "Yes", "plasmid", checkv_quality_and_plasmid)
         )

RPKM_count <- RPKM
RPKM_count[RPKM_count > 0] <- 1

RPKM_count_ch1 <- merge(RPKM_count, extended_tof[colnames(extended_tof) %in% c("checkv_quality_and_plasmid")], by="row.names")
RPKM_count_ch1 <- RPKM_count_ch1[RPKM_count_ch1$checkv_quality_and_plasmid != "plasmid", ]
RPKM_count_ch1$Row.names <- NULL

summarized_df <- RPKM_count_ch1 %>%
  group_by(checkv_quality_and_plasmid) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df <- as.data.frame(t(summarized_df))
colnames(summarized_df) <- c("CHM", "LU")
summarized_df <- summarized_df[row.names(summarized_df) != "checkv_quality_and_plasmid", ]

summarized_df <- summarized_df %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df$CHM_LU_ratio <- summarized_df$CHM / (summarized_df$LU + summarized_df$CHM)


row.names(meta_working) <- meta_working$Sample_name
meta_working <- merge(meta_working, summarized_df, by="row.names")

meta_working$Row.names <- NULL

pdf('CHM_LU_ratio_comparison.pdf', width = 12/2.54, height = 8/2.54)
ggplot(meta_working[!(meta_working$Sample_name %in% c("bctrl4633v", "bctrl4654v", "bctrl4676v", "bctrl4699v")), ], aes(x=ncvssample, y=CHM_LU_ratio)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Richness ratio of the C/H/M quality to all detected vOTUs \n (excluding plasmids and orsay controls)") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))
dev.off()

pdf('CHM_number_comparison.pdf', width = 12/2.54, height = 8/2.54)
ggplot(meta_working[!(meta_working$Sample_name %in% c("bctrl4633v", "bctrl4654v", "bctrl4676v", "bctrl4699v")), ], aes(x=ncvssample, y=CHM)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Richness of the C/H/M quality vOTUs \n (excluding plasmids and orsay controls)") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))
dev.off()


#######################################################################################################################################


## Analysis for the vOTUs composition in NC vs samples: abundance and ratio
#######################################################################################################################################
RPKM_ch1 <- merge(RPKM, extended_tof[colnames(extended_tof) %in% c("checkv_quality_and_plasmid")], by="row.names")
RPKM_ch1 <- RPKM_ch1[RPKM_ch1$checkv_quality_and_plasmid != "plasmid", ]
RPKM_ch1$Row.names <- NULL

summarized_df1 <- RPKM_ch1 %>%
  group_by(checkv_quality_and_plasmid) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df1 <- as.data.frame(t(summarized_df1))
colnames(summarized_df1) <- c("CHM_abundance", "LU_abundance")
summarized_df1 <- summarized_df1[row.names(summarized_df1) != "checkv_quality_and_plasmid", ]

summarized_df1 <- summarized_df1 %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df1$CHM_LU_abundance_ratio <- summarized_df1$CHM_abundance / (summarized_df1$LU_abundance + summarized_df1$CHM_abundance)


row.names(meta_working) <- meta_working$Sample_name
meta_working <- merge(meta_working, summarized_df1, by="row.names")

meta_working$Row.names <- NULL

pdf('CHM_LU_abundance_ratio_comparison.pdf', width = 12/2.54, height = 8/2.54)
ggplot(meta_working[!(meta_working$Sample_name %in% c("bctrl4633v", "bctrl4654v", "bctrl4676v", "bctrl4699v")), ], aes(x=ncvssample, y=CHM_LU_abundance_ratio)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Abundance ratio of the C/H/M quality to all detected vOTUs \n (excluding plasmids and orsay controls)") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))
dev.off()

pdf('CHM_abundance_comparison.pdf', width = 12/2.54, height = 8/2.54)
ggplot(meta_working[!(meta_working$Sample_name %in% c("bctrl4633v", "bctrl4654v", "bctrl4676v", "bctrl4699v")), ], aes(x=ncvssample, y=CHM_abundance)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Abundance of the C/H/M quality vOTUs \n (excluding plasmids and orsay controls)") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))
dev.off()


#######################################################################################################################################


## Analysis for the virus groups richness
#######################################################################################################################################
CHM_viruses <- row.names(extended_tof[extended_tof$checkv_quality_and_plasmid == "CHM", ])
RPKM_count_ch2 <- RPKM_count[row.names(RPKM_count) %in% CHM_viruses, ]
RPKM_count_ch2 <- merge(RPKM_count_ch2, extended_tof[colnames(extended_tof) %in% c("virus_group")], by="row.names")
RPKM_count_ch2$Row.names <- NULL

summarized_df2 <- RPKM_count_ch2 %>%
  group_by(virus_group) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df2 <- as.data.frame(t(summarized_df2))
colnames(summarized_df2) <- c("RNA", "Unclassified", "dsDNA", "ssDNA")
summarized_df2 <- summarized_df2[row.names(summarized_df2) != "virus_group", ]

summarized_df2 <- summarized_df2 %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df2 <- summarized_df2[rowSums(summarized_df2) > 0, ]

fraction_per_sample <- function(row) {
  row / sum(row)
}

summarized_df2_frac <- as.data.frame(t(apply(summarized_df2, 1, fraction_per_sample)))

summarized_df2_frac$Sample_name <- row.names(summarized_df2_frac)
summarized_df2_frac <- merge(summarized_df2_frac, meta_working[, colnames(meta_working) %in% c("Sample_name", "ncvssample", "cohort")], all.x = T, by="Sample_name")
summarized_df2_frac$Sample_name <- NULL

summarized_df2_frac <- summarized_df2_frac %>%
  group_by(ncvssample, cohort) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

summarized_df2_frac <- as.data.frame(summarized_df2_frac)

## Redo for long format
summarized_df2_frac_melt <- melt(summarized_df2_frac, id.vars = c("ncvssample", "cohort"))
summarized_df2_frac_melt$variable <- as.factor(summarized_df2_frac_melt$variable)
summarized_df2_frac_melt$variable <- factor(summarized_df2_frac_melt$variable, levels = c("RNA", "ssDNA", "dsDNA", "Unclassified"))

## Mention in the plot, that only CHM viruses included here - no plasmids, no LU

pdf('mean_virus_group_richness_CHM.pdf', width = 12/2.54, height = 8/2.54)
ggplot(summarized_df2_frac_melt, aes(x = ncvssample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#a1c181", "#fcca46", "#fe7f2d", "#233d4d")) +
  facet_grid(. ~ cohort) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 6),  # Decrease y-axis font size
    axis.title = element_text(size = 8),  # Decrease axis titles font size
    legend.title = element_text(size = 8),  # Decrease legend title font size
    legend.text = element_text(size = 6),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "", y = "Mean fraction of the viral group richness \n (complete, high- and medium- quality vOTUs)", fill = "Viral group") #+
dev.off()

#######################################################################################################################################


## Analysis for the virus groups abundance
#######################################################################################################################################

RPKM_ch2 <- RPKM[row.names(RPKM) %in% CHM_viruses, ]
RPKM_ch2 <- merge(RPKM_ch2, extended_tof[colnames(extended_tof) %in% c("virus_group")], by="row.names")
RPKM_ch2$Row.names <- NULL

summarized_df2a <- RPKM_ch2 %>%
  group_by(virus_group) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df2a <- as.data.frame(t(summarized_df2a))
colnames(summarized_df2a) <- c("RNA", "Unclassified", "dsDNA", "ssDNA")
summarized_df2a <- summarized_df2a[row.names(summarized_df2a) != "virus_group", ]

summarized_df2a <- summarized_df2a %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df2a <- summarized_df2a[rowSums(summarized_df2a) > 0, ]

summarized_df2a_frac <- as.data.frame(t(apply(summarized_df2a, 1, fraction_per_sample)))

summarized_df2a_frac$Sample_name <- row.names(summarized_df2a_frac)
summarized_df2a_frac <- merge(summarized_df2a_frac, meta_working[, colnames(meta_working) %in% c("Sample_name", "ncvssample", "cohort")], all.x = T, by="Sample_name")
summarized_df2a_frac$Sample_name <- NULL

summarized_df2a_frac <- summarized_df2a_frac %>%
  group_by(ncvssample, cohort) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

summarized_df2a_frac <- as.data.frame(summarized_df2a_frac)

## Redo for long format
summarized_df2a_frac_melt <- melt(summarized_df2a_frac, id.vars = c("ncvssample", "cohort"))
summarized_df2a_frac_melt$variable <- as.factor(summarized_df2a_frac_melt$variable)
summarized_df2a_frac_melt$variable <- factor(summarized_df2a_frac_melt$variable, levels = c("RNA", "ssDNA", "dsDNA", "Unclassified"))

pdf('mean_virus_group_abundance_CHM.pdf', width = 12/2.54, height = 8/2.54)
ggplot(summarized_df2a_frac_melt, aes(x = ncvssample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#a1c181", "#fcca46", "#fe7f2d", "#233d4d")) +
  facet_grid(. ~ cohort) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 6),  # Decrease y-axis font size
    axis.title = element_text(size = 8),  # Decrease axis titles font size
    legend.title = element_text(size = 8),  # Decrease legend title font size
    legend.text = element_text(size = 6),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "", y = "Mean fraction of the viral group abundance \n (complete, high- and medium- quality vOTUs)", fill = "Viral group") #+
dev.off()

#######################################################################################################################################


## Analysis for the virus host richness
#######################################################################################################################################
RPKM_count_ch3 <- RPKM_count[row.names(RPKM_count) %in% CHM_viruses, ]
RPKM_count_ch3 <- merge(RPKM_count_ch3, extended_tof[colnames(extended_tof) %in% c("host_group")], by="row.names")
RPKM_count_ch3$Row.names <- NULL

summarized_df3 <- RPKM_count_ch3 %>%
  group_by(host_group) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df3 <- as.data.frame(t(summarized_df3))
colnames(summarized_df3) <- c("Unclassified", "animal", "archaea", "bacteria", "plant")
summarized_df3 <- summarized_df3[row.names(summarized_df3) != "host_group", ]

summarized_df3 <- summarized_df3 %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df3 <- summarized_df3[rowSums(summarized_df3) > 0, ]

summarized_df3_frac <- as.data.frame(t(apply(summarized_df3, 1, fraction_per_sample)))

summarized_df3_frac$Sample_name <- row.names(summarized_df3_frac)
summarized_df3_frac <- merge(summarized_df3_frac, meta_working[, colnames(meta_working) %in% c("Sample_name", "ncvssample", "cohort")], all.x = T, by="Sample_name")
summarized_df3_frac$Sample_name <- NULL

summarized_df3_frac <- summarized_df3_frac %>%
  group_by(ncvssample, cohort) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

summarized_df3_frac <- as.data.frame(summarized_df3_frac)

## Redo for long format
summarized_df3_frac_melt <- melt(summarized_df3_frac, id.vars = c("ncvssample", "cohort"))
summarized_df3_frac_melt$variable <- as.factor(summarized_df3_frac_melt$variable)
summarized_df3_frac_melt$variable <- factor(summarized_df3_frac_melt$variable, levels = c("bacteria", "archaea", "animal", "plant", "Unclassified"))

## Mention in the plot, that only CHM viruses included here - no plasmids, no LU

pdf('mean_host_group_richness_CHM.pdf', width = 12/2.54, height = 8/2.54)
ggplot(summarized_df3_frac_melt, aes(x = ncvssample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#a1c181", "#619b8a", "#fcca46", "#fe7f2d", "#233d4d")) +
  facet_grid(. ~ cohort) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 6),  # Decrease y-axis font size
    axis.title = element_text(size = 8),  # Decrease axis titles font size
    legend.title = element_text(size = 8),  # Decrease legend title font size
    legend.text = element_text(size = 6),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "", y = "Mean fraction of the host group richness \n (complete, high- and medium- quality vOTUs)", fill = "Host group") #+
dev.off()

#######################################################################################################################################


## Analysis for the virus host abundance
#######################################################################################################################################
RPKM_ch3 <- RPKM[row.names(RPKM) %in% CHM_viruses, ]
RPKM_ch3 <- merge(RPKM_ch3, extended_tof[colnames(extended_tof) %in% c("host_group")], by="row.names")
RPKM_ch3$Row.names <- NULL

summarized_df3a <- RPKM_ch3 %>%
  group_by(host_group) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df3a <- as.data.frame(t(summarized_df3a))
colnames(summarized_df3a) <- c("Unclassified", "animal", "archaea", "bacteria", "plant")
summarized_df3a <- summarized_df3a[row.names(summarized_df3a) != "host_group", ]

summarized_df3a <- summarized_df3a %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df3a <- summarized_df3a[rowSums(summarized_df3a) > 0, ]

summarized_df3a_frac <- as.data.frame(t(apply(summarized_df3a, 1, fraction_per_sample)))

summarized_df3a_frac$Sample_name <- row.names(summarized_df3a_frac)
summarized_df3a_frac <- merge(summarized_df3a_frac, meta_working[, colnames(meta_working) %in% c("Sample_name", "ncvssample", "cohort")], all.x = T, by="Sample_name")
summarized_df3a_frac$Sample_name <- NULL

summarized_df3a_frac <- summarized_df3a_frac %>%
  group_by(ncvssample, cohort) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

summarized_df3a_frac <- as.data.frame(summarized_df3a_frac)

## Redo for long format
summarized_df3a_frac_melt <- melt(summarized_df3a_frac, id.vars = c("ncvssample", "cohort"))
summarized_df3a_frac_melt$variable <- as.factor(summarized_df3a_frac_melt$variable)
summarized_df3a_frac_melt$variable <- factor(summarized_df3a_frac_melt$variable, levels = c("bacteria", "archaea", "animal", "plant", "Unclassified"))

pdf('mean_host_group_abundance_CHM.pdf', width = 12/2.54, height = 8/2.54)
ggplot(summarized_df3a_frac_melt, aes(x = ncvssample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#a1c181", "#619b8a", "#fcca46", "#fe7f2d", "#233d4d")) +
  facet_grid(. ~ cohort) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 6),  # Decrease y-axis font size
    axis.title = element_text(size = 8),  # Decrease axis titles font size
    legend.title = element_text(size = 8),  # Decrease legend title font size
    legend.text = element_text(size = 6),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "", y = "Mean fraction of the host group abundance \n (complete, high- and medium- quality vOTUs)", fill = "Host group") #+
dev.off()

#######################################################################################################################################



## Analysis for the phages host richness
#######################################################################################################################################
BU_viruses <- row.names(extended_tof[extended_tof$host_group %in% c("bacteria", "archaea", "Unclassified"), ])
RPKM_count_ch4 <- RPKM_count[row.names(RPKM_count) %in% CHM_viruses & row.names(RPKM_count) %in% BU_viruses, ]


row.names(filtered_host_prediction) <- filtered_host_prediction$Virus
filtered_host_prediction$Genus <- sub(".*;g__", "", filtered_host_prediction$Host.genus)
filtered_host_prediction$Genus[filtered_host_prediction$Genus == ""] <- "Unclassified"

RPKM_count_ch4 <- merge(RPKM_count_ch4, filtered_host_prediction[colnames(filtered_host_prediction) %in% c("Genus")], all.x = T, by="row.names")
RPKM_count_ch4$Genus[is.na(RPKM_count_ch4$Genus)] <- "Unclassified"
RPKM_count_ch4$Row.names <- NULL

summarized_df4 <- RPKM_count_ch4 %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df4 <- as.data.frame(t(summarized_df4))
colnames(summarized_df4) <- c(t(summarized_df4[1,]))
summarized_df4 <- summarized_df4[row.names(summarized_df4) != "Genus", ]

summarized_df4 <- summarized_df4 %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df4 <- summarized_df4[rowSums(summarized_df4) > 0, ]

summarized_df4_frac <- as.data.frame(t(apply(summarized_df4, 1, fraction_per_sample)))

summarized_df4_frac$Sample_name <- row.names(summarized_df4_frac)
summarized_df4_frac <- merge(summarized_df4_frac, meta_working[, colnames(meta_working) %in% c("Sample_name", "ncvssample", "cohort")], all.x = T, by="Sample_name")
summarized_df4_frac$Sample_name <- NULL

summarized_df4_frac <- summarized_df4_frac %>%
  group_by(ncvssample, cohort) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

summarized_df4_frac <- as.data.frame(summarized_df4_frac)

columns_other <- c()
for (col in names(summarized_df4_frac[3:ncol(summarized_df4_frac)])) {
  if (max(summarized_df4_frac[[col]], na.rm = TRUE) < 0.01) {
    columns_other <- c(columns_other, col)
  }
}

summarized_df4_frac <- summarized_df4_frac %>%
  mutate(Other = rowSums(select(., all_of(columns_other)), na.rm = TRUE))

summarized_df4_frac <- summarized_df4_frac %>%
  select(-all_of(columns_other))


## Redo for long format
summarized_df4_frac_melt <- melt(summarized_df4_frac, id.vars = c("ncvssample", "cohort"))
summarized_df4_frac_melt$variable <- as.factor(summarized_df4_frac_melt$variable)
summarized_df4_frac_melt$variable <- factor(summarized_df4_frac_melt$variable, levels = c("bacteria", "archaea", "animal", "plant", "Unclassified"))
summarized_df4_frac_melt$variable <- factor(summarized_df4_frac_melt$variable, 
                                            levels = levels(summarized_df4_frac_melt$variable)[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 30)])

pdf('mean_host_bacteria_richness_CHM.pdf', width = 12/2.54, height = 8/2.54)
ggplot(summarized_df4_frac_melt, aes(x = ncvssample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("cornflowerblue", "gold2", "turquoise3", "orchid", "steelblue", "burlywood", "purple3", "slategray4","yellow", "blue", "skyblue","darkgoldenrod1",
                               "darksalmon", "green", "pink", "tomato2", "rosybrown3", "snow2","plum2",
                               "palegreen3", "seagreen", "aquamarine", "yellow3", "chartreuse", "thistle", "orange2", 
                               "navajowhite", "firebrick", "purple", "violetred1", "red", "darkblue", "black")) +
  facet_grid(. ~ cohort) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 6),  # Decrease y-axis font size
    axis.title = element_text(size = 8),  # Decrease axis titles font size
    legend.title = element_text(size = 8),  # Decrease legend title font size
    legend.text = element_text(size = 6),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "", y = "Mean fraction of the host group abundance \n (complete, high- and medium- quality vOTUs)", fill = "Host group") #+
dev.off()
#######################################################################################################################################






































## Analysis for the vOTUs found in negative controls: Bray-curtis heat-map
#######################################################################################################################################
negative_controls <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" &
                                                meta_working$Sample_name %in% colnames(RPKM)]

RPKM_filtered_nc <- RPKM[, colnames(RPKM) %in% negative_controls]
RPKM_filtered_nc <- RPKM_filtered_nc[rowSums(RPKM_filtered_nc) > 0, colSums(RPKM_filtered_nc) > 0]

meta_working <- meta_working %>%
  mutate(nc_group = ifelse(grepl("bctrl1", Sample_name), "maqsood_buffer", NA),
         nc_group = ifelse(grepl("bctrl4", Sample_name), "maqsood_orsay", nc_group),
         nc_group = ifelse(grepl("ctl", Sample_name), "shah_buffer", nc_group),
         nc_group = ifelse(grepl("LN_7C08_VL_405", Sample_name), "garmaeva_buffer", nc_group),
         nc_group = ifelse(Sample_ID %in% c("LT1_D", "LT4_D", "LT5_D", "LT2_D", "LT3_D"), "liang_tube_DNA", nc_group),
         nc_group = ifelse(Sample_ID %in% c("LT1_R", "LT2_R", "LT3_R", "LT4_R", "LT5_R"), "liang_tube_RNA", nc_group),
         nc_group = ifelse(Sample_ID %in% c("LB1_D", "LB2_D", "LB3_D", "LB4_D", "LB5_D"), "liang_reagent_DNA", nc_group),
         nc_group = ifelse(Sample_ID %in% c("LB1_R", "LB2_R", "LB3_R", "LB4_R", "LB5_R"), "liang_reagent_RNA", nc_group),
         nc_group = ifelse(Sample_ID %in% c("LD1_D", "LD2_D", "LDH_D", "LDN1_D", "LDN2_D"), "liang_diaper_DNA", nc_group),
         nc_group = ifelse(Sample_ID %in% c("LD1_R", "LD2_R", "LDH_R", "LDN1_R", "LDN2_R"), "liang_diaper_RNA", nc_group),)  # NOT FINISHED!!!!!!!

# preapring phenotypes
meta_nc_nmds <- meta_working[meta_working$Sample_name %in% colnames(RPKM_filtered_nc), c("Sample_name", "cohort", "clean_reads_comb", "total_viruses_discovered")]
row.names(meta_nc_nmds) <- meta_nc_nmds$Sample_name
meta_nc_nmds <- merge(meta_nc_nmds, as.data.frame(t(as.data.frame(lapply(RPKM_filtered_nc, diversity, index = "shannon")))),
                      by="row.names")
row.names(meta_nc_nmds) <- meta_nc_nmds$Sample_name
meta_nc_nmds <- meta_nc_nmds[c("cohort", "clean_reads_comb", "total_viruses_discovered", "V1")]
names(meta_nc_nmds) <- c("cohort", "clean_reads_comb", "total_viruses_discovered", "shannon_div")

## Bray-curtis heat map analysis
nc_bray_dist_matrix <- as.matrix(vegdist(t(RPKM_filtered_nc), method="bray"))
nc_bray_dist_matrix_rev <- 1 - nc_bray_dist_matrix



annotation <- meta_nc_nmds %>%
  mutate(clean_reads_comb = NULL,
         total_viruses_discovered = NULL,
         shannon_div = NULL)

pdf('dist_mat_bray_curtis.pdf', width=12/2.54, height=10/2.54)
pheatmap(nc_bray_dist_matrix_rev, 
         annotation = annotation, 
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F, 
         color = colorRampPalette(c("white", "darkblue"))(100),
         fontsize = 5,
         main = "Negative controls Bray-Curtis Similarity")
dev.off()

#######################################################################################################################################










## Analysis for the vOTUs in uncontaminated vs contaminated samples: Bray-curtis box-plots
#######################################################################################################################################
viruses_in_nc <- row.names(RPKM_filtered_nc)

RPKM_viruses_found_in_nc <- RPKM[row.names(RPKM) %in% viruses_in_nc, ]
RPKM_viruses_found_in_nc <- RPKM_viruses_found_in_nc[rowSums(RPKM_viruses_found_in_nc) > 0, colSums(RPKM_viruses_found_in_nc) > 0]



RPKM_filtered_nc <- RPKM[, colnames(RPKM) %in% negative_controls]
RPKM_filtered_nc <- RPKM_filtered_nc[rowSums(RPKM_filtered_nc) > 0, colSums(RPKM_filtered_nc) > 0]


# preapring phenotypes
meta_nc_nmds <- meta_working[meta_working$Sample_name %in% colnames(RPKM_filtered_nc), c("Sample_name", "cohort", "clean_reads_comb", "total_viruses_discovered")]
row.names(meta_nc_nmds) <- meta_nc_nmds$Sample_name
meta_nc_nmds <- merge(meta_nc_nmds, as.data.frame(t(as.data.frame(lapply(RPKM_filtered_nc, diversity, index = "shannon")))),
                      by="row.names")
row.names(meta_nc_nmds) <- meta_nc_nmds$Sample_name
meta_nc_nmds <- meta_nc_nmds[c("cohort", "clean_reads_comb", "total_viruses_discovered", "V1")]
names(meta_nc_nmds) <- c("cohort", "clean_reads_comb", "total_viruses_discovered", "shannon_div")

## Bray-curtis heat map analysis
nc_bray_dist_matrix <- as.matrix(vegdist(t(RPKM_filtered_nc), method="bray"))
nc_bray_dist_matrix_rev <- 1 - nc_bray_dist_matrix

# heatmap.2(x=nc_bray_dist_matrix, dendrogram=F)

annotation <- meta_nc_nmds %>%
  mutate(clean_reads_comb = NULL,
         total_viruses_discovered = NULL,
         shannon_div = NULL)

pdf('dist_mat_bray_curtis.pdf', width=12/2.54, height=10/2.54)
pheatmap(nc_bray_dist_matrix_rev, 
         annotation = annotation, 
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F, 
         color = colorRampPalette(c("white", "darkblue"))(100),
         fontsize = 5,
         main = "Negative controls Bray-Curtis Similarity")
dev.off()



#######################################################################################################################################

























## Analysis for the virus taxonomy in the nc
#######################################################################################################################################

RPKM_filtered_nc_count <- RPKM_filtered_nc
RPKM_filtered_nc_count[RPKM_filtered_nc_count > 0] <- 1

RPKM_filtered_nc_count <- merge(RPKM_filtered_nc_count, extended_tof[colnames(extended_tof) %in% c("virus_group")], all.x = T, by="row.names")

CHM_viruses <- row.names(extended_tof[extended_tof$checkv_quality_and_plasmid %in% c("Complete", "High-quality", "Medium-quality"), ])
LU_viruses <- row.names(extended_tof[extended_tof$checkv_quality_and_plasmid %in% c("Not-determined", "Low-quality"), ])

# length(CHM_viruses)  # 910010
# length(LU_viruses)  # 38988 

RPKM_filtered_nc_count_CHM <- RPKM_filtered_nc_count[RPKM_filtered_nc_count$Row.names %in% CHM_viruses, ]
RPKM_filtered_nc_count_LU <- RPKM_filtered_nc_count[RPKM_filtered_nc_count$Row.names %in% LU_viruses, ]

RPKM_filtered_nc_count$Row.names <- NULL

summarized_by_virus_group_CHM <- RPKM_filtered_nc_count_CHM %>%
  group_by(virus_group) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

long_df_virus_group_CHM <- summarized_by_virus_group_CHM %>%
  pivot_longer(cols = -virus_group, 
               names_to = "sample", 
               values_to = "count")

total_counts_virus_group_CHM <- long_df_virus_group_CHM %>%
  group_by(sample) %>%
  summarise(total_count = sum(count))


pdf('summary_per_nc_by_virus_group_count_CHM.pdf', width = 12/2.54, height = 8/2.54)
ggplot(long_df_virus_group_CHM, aes(x = sample, y = count, fill = virus_group)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 4),  # Decrease y-axis font size
    axis.title = element_text(size = 5),  # Decrease axis titles font size
    legend.title = element_text(size = 5),  # Decrease legend title font size
    legend.text = element_text(size = 4),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "Sample", y = "Count", fill = "Viral group") #+
# geom_text(data = total_counts_host_CHM, aes(x = sample, y = total_counts_host_CHM, label = total_counts_host_CHM), 
#           angle = 45, size = 1, inherit.aes = FALSE)  # Position text above the bars on a single line
dev.off()


summarized_by_virus_group_LU <- RPKM_filtered_nc_count_LU %>%
  group_by(virus_group) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

dim(summarized_by_virus_group_LU)  # 180  44

long_df_virus_group_LU <- summarized_by_virus_group_LU %>%
  pivot_longer(cols = -virus_group, 
               names_to = "sample", 
               values_to = "count")

total_counts_virus_group_LU <- long_df_virus_group_LU %>%
  group_by(sample) %>%
  summarise(total_count = sum(count))


pdf('summary_per_nc_by_virus_group_count_LU.pdf', width = 12/2.54, height = 8/2.54)
ggplot(long_df_virus_group_LU, aes(x = sample, y = count, fill = virus_group)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 4),  # Decrease y-axis font size
    axis.title = element_text(size = 5),  # Decrease axis titles font size
    legend.title = element_text(size = 5),  # Decrease legend title font size
    legend.text = element_text(size = 4),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "Sample", y = "Count", fill = "Viral group") #+
# geom_text(data = total_counts_host_CHM, aes(x = sample, y = total_counts_host_CHM, label = total_counts_host_CHM), 
#           angle = 45, size = 1, inherit.aes = FALSE)  # Position text above the bars on a single line
dev.off()

#######################################################################################################################################

## Analysis for the host range (kindoms) in the nc
#######################################################################################################################################

RPKM_filtered_nc_count <- RPKM_filtered_nc
RPKM_filtered_nc_count[RPKM_filtered_nc_count > 0] <- 1

RPKM_filtered_nc_count <- merge(RPKM_filtered_nc_count, extended_tof[colnames(extended_tof) %in% c("host_group")], all.x = T, by="row.names")

CHM_viruses <- row.names(extended_tof[extended_tof$checkv_quality_and_plasmid %in% c("Complete", "High-quality", "Medium-quality"), ])
LU_viruses <- row.names(extended_tof[extended_tof$checkv_quality_and_plasmid %in% c("Not-determined", "Low-quality"), ])

# length(CHM_viruses)  # 910010
# length(LU_viruses)  # 38988 

RPKM_filtered_nc_count_CHM <- RPKM_filtered_nc_count[RPKM_filtered_nc_count$Row.names %in% CHM_viruses, ]
RPKM_filtered_nc_count_LU <- RPKM_filtered_nc_count[RPKM_filtered_nc_count$Row.names %in% LU_viruses, ]

RPKM_filtered_nc_count$Row.names <- NULL

summarized_by_host_group_CHM <- RPKM_filtered_nc_count_CHM %>%
  group_by(host_group) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

long_df_host_group_CHM <- summarized_by_host_group_CHM %>%
  pivot_longer(cols = -host_group, 
               names_to = "sample", 
               values_to = "count")

total_counts_host_group_CHM <- long_df_host_group_CHM %>%
  group_by(sample) %>%
  summarise(total_count = sum(count))


pdf('summary_per_nc_by_host_group_count_CHM.pdf', width = 12/2.54, height = 8/2.54)
ggplot(long_df_host_group_CHM, aes(x = sample, y = count, fill = host_group)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 4),  # Decrease y-axis font size
    axis.title = element_text(size = 5),  # Decrease axis titles font size
    legend.title = element_text(size = 5),  # Decrease legend title font size
    legend.text = element_text(size = 4),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "Sample", y = "Count", fill = "Host group") #+
# geom_text(data = total_counts_host_CHM, aes(x = sample, y = total_counts_host_CHM, label = total_counts_host_CHM), 
#           angle = 45, size = 1, inherit.aes = FALSE)  # Position text above the bars on a single line
dev.off()


summarized_by_host_group_LU <- RPKM_filtered_nc_count_LU %>%
  group_by(host_group) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

dim(summarized_by_host_group_LU)  # 180  44

long_df_host_group_LU <- summarized_by_host_group_LU %>%
  pivot_longer(cols = -host_group, 
               names_to = "sample", 
               values_to = "count")

total_counts_host_group_LU <- long_df_host_group_LU %>%
  group_by(sample) %>%
  summarise(total_count = sum(count))


pdf('summary_per_nc_by_host_group_count_LU.pdf', width = 12/2.54, height = 8/2.54)
ggplot(long_df_host_group_LU, aes(x = sample, y = count, fill = host_group)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 4),  # Decrease y-axis font size
    axis.title = element_text(size = 5),  # Decrease axis titles font size
    legend.title = element_text(size = 5),  # Decrease legend title font size
    legend.text = element_text(size = 4),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "Sample", y = "Count", fill = "Viral group") #+
# geom_text(data = total_counts_host_CHM, aes(x = sample, y = total_counts_host_CHM, label = total_counts_host_CHM), 
#           angle = 45, size = 1, inherit.aes = FALSE)  # Position text above the bars on a single line
dev.off()

#######################################################################################################################################






























## Analysis for the hosts for the viruses in the nc
#######################################################################################################################################

RPKM_filtered_nc_count <- RPKM_filtered_nc
RPKM_filtered_nc_count[RPKM_filtered_nc_count > 0] <- 1

filtered_host_prediction <- as.data.frame(filtered_host_prediction)
row.names(filtered_host_prediction) <- filtered_host_prediction$Virus

RPKM_filtered_nc_count <- merge(RPKM_filtered_nc_count, filtered_host_prediction[colnames(filtered_host_prediction) %in% c("Host.genus")], all.x = T, by="row.names")

CHM_viruses <- row.names(extended_tof[extended_tof$checkv_quality_and_plasmid %in% c("Complete", "High-quality", "Medium-quality"), ])
LU_viruses <- row.names(extended_tof[extended_tof$checkv_quality_and_plasmid %in% c("Not-determined", "Low-quality"), ])

# length(CHM_viruses)  # 910010
# length(LU_viruses)  # 38988 

RPKM_filtered_nc_count_CHM <- RPKM_filtered_nc_count[RPKM_filtered_nc_count$Row.names %in% CHM_viruses, ]
RPKM_filtered_nc_count_LU <- RPKM_filtered_nc_count[RPKM_filtered_nc_count$Row.names %in% LU_viruses, ]

RPKM_filtered_nc_count$Row.names <- NULL

summarized_by_host_CHM <- RPKM_filtered_nc_count_CHM %>%
  group_by(Host.genus) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

long_df_host_CHM <- summarized_by_host_CHM %>%
  pivot_longer(cols = -Host.genus, 
               names_to = "sample", 
               values_to = "count")
long_df_host_CHM$Host.genus <- sub(".*;g__", "", long_df_host_CHM$Host.genus)

total_counts_host_CHM <- long_df_host_CHM %>%
  group_by(sample) %>%
  summarise(total_count = sum(count))


pdf('summary_per_nc_by_host_count_CHM.pdf', width = 12/2.54, height = 8/2.54)
ggplot(long_df_host_CHM, aes(x = sample, y = count, fill = Host.genus)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 4),  # Decrease y-axis font size
    axis.title = element_text(size = 5),  # Decrease axis titles font size
    legend.title = element_text(size = 5),  # Decrease legend title font size
    legend.text = element_text(size = 4),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "Sample", y = "Count", fill = "Host genus") #+
  # geom_text(data = total_counts_host_CHM, aes(x = sample, y = total_counts_host_CHM, label = total_counts_host_CHM), 
  #           angle = 45, size = 1, inherit.aes = FALSE)  # Position text above the bars on a single line
dev.off()


summarized_by_host_LU <- RPKM_filtered_nc_count_LU %>%
  group_by(Host.genus) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

dim(summarized_by_host_LU)  # 180  44

long_df_host_LU <- summarized_by_host_LU %>%
  pivot_longer(cols = -Host.genus, 
               names_to = "sample", 
               values_to = "count")
long_df_host_LU$Host.genus <- sub(".*;g__", "", long_df_host_LU$Host.genus)

total_counts_host_LU <- long_df_host_LU %>%
  group_by(sample) %>%
  summarise(total_count = sum(count))


pdf('summary_per_nc_by_host_count_LU.pdf', width = 12/2.54, height = 8/2.54)
ggplot(long_df_host_LU, aes(x = sample, y = count, fill = Host.genus)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 4),  # Decrease y-axis font size
    axis.title = element_text(size = 5),  # Decrease axis titles font size
    legend.title = element_text(size = 1),  # Decrease legend title font size
    legend.text = element_text(size = 1),  # Decrease legend text font size
    legend.key.size = unit(0.1, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "Sample", y = "Count", fill = "Host genus") #+
# geom_text(data = total_counts_host_CHM, aes(x = sample, y = total_counts_host_CHM, label = total_counts_host_CHM), 
#           angle = 45, size = 1, inherit.aes = FALSE)  # Position text above the bars on a single line
dev.off()

#######################################################################################################################################





## Analysis for the vOTUs sharedness with the other cohorts
#######################################################################################################################################
negative_controls_garmaeva <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "garmaeva" & 
                                                            meta_working$Sample_name %in% colnames(RPKM)]

negative_controls_liang <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "liang" &
                                                         meta_working$Sample_name %in% colnames(RPKM)]

negative_controls_maqsood <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "maqsood" &
                                                           meta_working$Sample_name %in% colnames(RPKM)]

negative_controls_shah <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "shah" &
                                                        meta_working$Sample_name %in% colnames(RPKM)]

RPKM_filtered_w_dummy <- RPKM

RPKM_filtered_w_dummy$dummy_all <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls])

RPKM_filtered_w_dummy$dummy_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) == negative_controls_garmaeva]
RPKM_filtered_w_dummy$dummy_liang <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_liang])
RPKM_filtered_w_dummy$dummy_maqsood <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_maqsood])
RPKM_filtered_w_dummy$dummy_shah <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_shah])

RPKM_filtered_w_dummy[RPKM_filtered_w_dummy > 0] <- 1

samples <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl"]

temp_shared_w_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_garmaeva
temp_shared_w_maqsood <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_maqsood
temp_shared_w_liang <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_liang
temp_shared_w_shah <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_shah

temp_shared_w_garmaeva[temp_shared_w_garmaeva < 2] <- 0
temp_shared_w_maqsood[temp_shared_w_maqsood < 2] <- 0
temp_shared_w_liang[temp_shared_w_liang < 2] <- 0
temp_shared_w_shah[temp_shared_w_shah < 2] <- 0

table_for_plot_g <- as.data.frame(colSums(temp_shared_w_garmaeva) / 2)
table_for_plot_m <- as.data.frame(colSums(temp_shared_w_maqsood) / 2)
table_for_plot_l <- as.data.frame(colSums(temp_shared_w_liang) / 2)
table_for_plot_s <- as.data.frame(colSums(temp_shared_w_shah) / 2)

table_for_plot_g_m <- merge(table_for_plot_g, table_for_plot_m, by="row.names", all = T)
table_for_plot_l_s <- merge(table_for_plot_l, table_for_plot_s, by="row.names", all = T)
table_for_plot <- merge(table_for_plot_g_m, table_for_plot_l_s, by="Row.names", all = T)
row.names(table_for_plot) <- table_for_plot$Row.names
table_for_plot$Row.names <- NULL

RPKM_count <- RPKM
RPKM_count[RPKM_count > 0] <- 1

table_for_plot_add <- as.data.frame(colSums(RPKM_count))
table_for_plot <- merge(table_for_plot, table_for_plot_add, by='row.names', all.x = T)

colnames(table_for_plot) <- c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                              "Shared_w_nc_shah", "Viral_richness")

table_for_plot$Perc_shared_w_garmaeva <- table_for_plot$Shared_w_nc_garmaeva / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_maqsood <- table_for_plot$Shared_w_nc_maqsood / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_liang <- table_for_plot$Shared_w_nc_liang / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_shah <- table_for_plot$Shared_w_nc_shah / table_for_plot$Viral_richness *100

table_for_plot <- merge(table_for_plot, meta_working[c("Sample_name", "Type", "Timepoint", "cohort")], by="Sample_name", all.x = T)

table_for_plot_melt <- melt(table_for_plot, id = c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                                                   "Shared_w_nc_shah", "Viral_richness", "Type", "Timepoint", "cohort"))

table_for_plot_melt$variable <- factor(table_for_plot_melt$variable, levels = c("Perc_shared_w_garmaeva", "Perc_shared_w_liang", 
                                                                                "Perc_shared_w_maqsood", "Perc_shared_w_shah"))

## Plotting only sharedness for infants
pdf('shared_inf.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ ., scales = "free") +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

## Plotting only sharendess for mothers
pdf('shared_moms.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Mother", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ .) +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

variable.labs <- c("Garmaeva NCs", "Maqsood NCs", "Liang NCs", "Shah NCs")
names(variable.labs) <- c("Perc_shared_w_garmaeva", "Perc_shared_w_maqsood", "Perc_shared_w_liang", "Perc_shared_w_shah")

# Convert named vector to labeller
variable_labeller <- as_labeller(variable.labs)

## Plotting sharedness depending on the timepoints
pdf('shared_timepoints.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant" & table_for_plot_melt$cohort %in% c("garmaeva", "liang"), ], aes(x = Timepoint, y = value)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.13, size = 0.5) +
  ylim(0, 100) +
  facet_grid(rows = vars(variable), cols = vars(cohort), scales = "free", labeller = labeller(variable = variable_labeller)) +
  scale_x_discrete() + 
  labs(y = "%", x = "Timepoint") +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()
#######################################################################################################################################




##############################################################################################################################################################################################################################################################################


# THIS BLOCK IS DEDICATED TO EXPLORING THE RESULTS, WHEN THE RPKM TABLE IS FILTERED FOR PLASMIDS


##############################################################################################################################################################################################################################################################################


## File preparation
#######################################################################################################################################
virus_selection <- extended_tof$New_CID[extended_tof$plasmid == "No"]
length(virus_selection)  # 948998 (22,585 identified redundant "viruses" are actually plasmids)

RPKM_filtered <- RPKM[row.names(RPKM) %in% virus_selection, ]
RPKM_filtered <- RPKM_filtered[rowSums(RPKM_filtered) > 0, colSums(RPKM_filtered) > 0]
dim(RPKM_filtered)  # 298689   1300

meta_working <- meta_all_with_qc_curated[meta_all_with_qc_curated$Sample_name %in% colnames(RPKM_filtered), ]
#######################################################################################################################################

## Analysis for the vOTUs found in negative controls: Bray-curtis heat-map
#######################################################################################################################################
negative_controls <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" &
                                                meta_working$Sample_name %in% colnames(RPKM_filtered)]

RPKM_filtered_nc <- RPKM_filtered[, colnames(RPKM_filtered) %in% negative_controls]
RPKM_filtered_nc <- RPKM_filtered_nc[rowSums(RPKM_filtered_nc) > 0, colSums(RPKM_filtered_nc) > 0]


# preapring phenotypes
meta_nc_nmds <- meta_working[meta_working$Sample_name %in% colnames(RPKM_filtered_nc), c("Sample_name", "cohort", "clean_reads_comb", "total_viruses_discovered")]
row.names(meta_nc_nmds) <- meta_nc_nmds$Sample_name
meta_nc_nmds <- merge(meta_nc_nmds, as.data.frame(t(as.data.frame(lapply(RPKM_filtered_nc, diversity, index = "shannon")))),
                      by="row.names")
row.names(meta_nc_nmds) <- meta_nc_nmds$Sample_name
meta_nc_nmds <- meta_nc_nmds[c("cohort", "clean_reads_comb", "total_viruses_discovered", "V1")]
names(meta_nc_nmds) <- c("cohort", "clean_reads_comb", "total_viruses_discovered", "shannon_div")

## Bray-curtis heat map analysis
nc_bray_dist_matrix <- as.matrix(vegdist(t(RPKM_filtered_nc), method="bray"))
nc_bray_dist_matrix_rev <- 1 - nc_bray_dist_matrix

# heatmap.2(x=nc_bray_dist_matrix, dendrogram=F)

annotation <- meta_nc_nmds %>%
  mutate(clean_reads_comb = NULL,
         total_viruses_discovered = NULL,
         shannon_div = NULL)

pdf('dist_mat_bray_curtis_no_plasmids.pdf', width=12/2.54, height=10/2.54)
pheatmap(nc_bray_dist_matrix_rev, 
         annotation = annotation, 
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F, 
         color = colorRampPalette(c("white", "darkblue"))(100),
         fontsize = 5,
         main = "Negative controls Bray-Curtis Similarity")
dev.off()


#######################################################################################################################################

## Analysis for the vOTUs sharedness with the other cohorts
#######################################################################################################################################
negative_controls_garmaeva <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "garmaeva" & 
                                                         meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_liang <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "liang" &
                                                      meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_maqsood <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "maqsood" &
                                                        meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_shah <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "shah" &
                                                     meta_working$Sample_name %in% colnames(RPKM_filtered)]

RPKM_filtered_w_dummy <- RPKM_filtered

RPKM_filtered_w_dummy$dummy_all <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls])

RPKM_filtered_w_dummy$dummy_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) == negative_controls_garmaeva]
RPKM_filtered_w_dummy$dummy_liang <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_liang])
RPKM_filtered_w_dummy$dummy_maqsood <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_maqsood])
RPKM_filtered_w_dummy$dummy_shah <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_shah])

RPKM_filtered_w_dummy[RPKM_filtered_w_dummy > 0] <- 1

samples <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl"]

temp_shared_w_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_garmaeva
temp_shared_w_maqsood <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_maqsood
temp_shared_w_liang <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_liang
temp_shared_w_shah <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_shah

temp_shared_w_garmaeva[temp_shared_w_garmaeva < 2] <- 0
temp_shared_w_maqsood[temp_shared_w_maqsood < 2] <- 0
temp_shared_w_liang[temp_shared_w_liang < 2] <- 0
temp_shared_w_shah[temp_shared_w_shah < 2] <- 0

table_for_plot_g <- as.data.frame(colSums(temp_shared_w_garmaeva) / 2)
table_for_plot_m <- as.data.frame(colSums(temp_shared_w_maqsood) / 2)
table_for_plot_l <- as.data.frame(colSums(temp_shared_w_liang) / 2)
table_for_plot_s <- as.data.frame(colSums(temp_shared_w_shah) / 2)

table_for_plot_g_m <- merge(table_for_plot_g, table_for_plot_m, by="row.names", all = T)
table_for_plot_l_s <- merge(table_for_plot_l, table_for_plot_s, by="row.names", all = T)
table_for_plot <- merge(table_for_plot_g_m, table_for_plot_l_s, by="Row.names", all = T)
row.names(table_for_plot) <- table_for_plot$Row.names
table_for_plot$Row.names <- NULL

RPKM_count <- RPKM_filtered
RPKM_count[RPKM_count > 0] <- 1

table_for_plot_add <- as.data.frame(colSums(RPKM_count))
table_for_plot <- merge(table_for_plot, table_for_plot_add, by='row.names', all.x = T)

colnames(table_for_plot) <- c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                              "Shared_w_nc_shah", "Viral_richness")

table_for_plot$Perc_shared_w_garmaeva <- table_for_plot$Shared_w_nc_garmaeva / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_maqsood <- table_for_plot$Shared_w_nc_maqsood / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_liang <- table_for_plot$Shared_w_nc_liang / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_shah <- table_for_plot$Shared_w_nc_shah / table_for_plot$Viral_richness *100

table_for_plot <- merge(table_for_plot, meta_working[c("Sample_name", "Type", "Timepoint", "cohort")], by="Sample_name", all.x = T)

table_for_plot_melt <- melt(table_for_plot, id = c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                                                   "Shared_w_nc_shah", "Viral_richness", "Type", "Timepoint", "cohort"))

table_for_plot_melt$variable <- factor(table_for_plot_melt$variable, levels = c("Perc_shared_w_garmaeva", "Perc_shared_w_liang", 
                                                                                "Perc_shared_w_maqsood", "Perc_shared_w_shah"))

## Plotting only sharedness for infants
pdf('shared_inf_no_plasmids.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ ., scales = "free") +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

## Plotting only sharendess for mothers
pdf('shared_moms_no_plasmids.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Mother", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ .) +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

variable.labs <- c("Garmaeva NCs", "Maqsood NCs", "Liang NCs", "Shah NCs")
names(variable.labs) <- c("Perc_shared_w_garmaeva", "Perc_shared_w_maqsood", "Perc_shared_w_liang", "Perc_shared_w_shah")

# Convert named vector to labeller
variable_labeller <- as_labeller(variable.labs)

## Plotting sharedness depending on the timepoints
pdf('shared_timepoints_no_plasmids.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant" & table_for_plot_melt$cohort %in% c("garmaeva", "liang"), ], aes(x = Timepoint, y = value)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.13, size = 0.5) +
  ylim(0, 100) +
  facet_grid(rows = vars(variable), cols = vars(cohort), scales = "free", labeller = labeller(variable = variable_labeller)) +
  scale_x_discrete() + 
  labs(y = "%", x = "Timepoint") +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()





##############################################################################################################################################################################################################################################################################


# THIS BLOCK IS DEDICATED TO EXPLORING THE RESULTS, WHEN THE RPKM TABLE IS FILTERED FOR UNDETERMINED LOW QUALITY VIRUSES AND PLASMIDS


##############################################################################################################################################################################################################################################################################


#######################################################################################################################################
virus_selection <- extended_tof$New_CID[extended_tof$plasmid == "No" & extended_tof$checkv_quality %in% c("High-quality", "Medium-quality", "Complete")]
length(virus_selection)  # 38988 (932,595 identified redundant "viruses" are either viruses of low- undetermined quality or plasmids)

RPKM_filtered <- RPKM[row.names(RPKM) %in% virus_selection, ]
RPKM_filtered <- RPKM_filtered[rowSums(RPKM_filtered) > 0, colSums(RPKM_filtered) > 0]
dim(RPKM_filtered)  # 14519  1185

meta_working <- meta_all_with_qc_curated[meta_all_with_qc_curated$Sample_name %in% colnames(RPKM_filtered), ]
#######################################################################################################################################

## Analysis for the vOTUs found in negative controls: Bray-curtis heat-map
#######################################################################################################################################
negative_controls <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" &
                                                meta_working$Sample_name %in% colnames(RPKM_filtered)]

RPKM_filtered_nc <- RPKM_filtered[, colnames(RPKM_filtered) %in% negative_controls]
RPKM_filtered_nc <- RPKM_filtered_nc[rowSums(RPKM_filtered_nc) > 0, colSums(RPKM_filtered_nc) > 0]


# preapring phenotypes
meta_nc_nmds <- meta_working[meta_working$Sample_name %in% colnames(RPKM_filtered_nc), c("Sample_name", "cohort", "clean_reads_comb", "total_viruses_discovered")]
row.names(meta_nc_nmds) <- meta_nc_nmds$Sample_name
meta_nc_nmds <- merge(meta_nc_nmds, as.data.frame(t(as.data.frame(lapply(RPKM_filtered_nc, diversity, index = "shannon")))),
                      by="row.names")
row.names(meta_nc_nmds) <- meta_nc_nmds$Sample_name
meta_nc_nmds <- meta_nc_nmds[c("cohort", "clean_reads_comb", "total_viruses_discovered", "V1")]
names(meta_nc_nmds) <- c("cohort", "clean_reads_comb", "total_viruses_discovered", "shannon_div")

## Bray-curtis heat map analysis
nc_bray_dist_matrix <- as.matrix(vegdist(t(RPKM_filtered_nc), method="bray"))
nc_bray_dist_matrix_rev <- 1 - nc_bray_dist_matrix

# heatmap.2(x=nc_bray_dist_matrix, dendrogram=F)

annotation <- meta_nc_nmds %>%
  mutate(clean_reads_comb = NULL,
         total_viruses_discovered = NULL,
         shannon_div = NULL)

pdf('dist_mat_bray_curtis_no_plasmids_no_lowundetermined_quality.pdf', width=12/2.54, height=10/2.54)
pheatmap(nc_bray_dist_matrix_rev, 
         annotation = annotation, 
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F, 
         color = colorRampPalette(c("white", "darkblue"))(100),
         fontsize = 5,
         main = "Negative controls Bray-Curtis Similarity")
dev.off()

#######################################################################################################################################

## Analysis for the vOTUs sharedness with the other cohorts
#######################################################################################################################################
negative_controls_garmaeva <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "garmaeva" & 
                                                         meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_liang <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "liang" &
                                                      meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_maqsood <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "maqsood" &
                                                        meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_shah <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "shah" &
                                                     meta_working$Sample_name %in% colnames(RPKM_filtered)]

RPKM_filtered_w_dummy <- RPKM_filtered

RPKM_filtered_w_dummy$dummy_all <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls])

RPKM_filtered_w_dummy$dummy_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) == negative_controls_garmaeva]
RPKM_filtered_w_dummy$dummy_liang <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_liang])
RPKM_filtered_w_dummy$dummy_maqsood <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_maqsood])
RPKM_filtered_w_dummy$dummy_shah <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_shah])

RPKM_filtered_w_dummy[RPKM_filtered_w_dummy > 0] <- 1

samples <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl"]

temp_shared_w_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_garmaeva
temp_shared_w_maqsood <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_maqsood
temp_shared_w_liang <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_liang
temp_shared_w_shah <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_shah

temp_shared_w_garmaeva[temp_shared_w_garmaeva < 2] <- 0
temp_shared_w_maqsood[temp_shared_w_maqsood < 2] <- 0
temp_shared_w_liang[temp_shared_w_liang < 2] <- 0
temp_shared_w_shah[temp_shared_w_shah < 2] <- 0

table_for_plot_g <- as.data.frame(colSums(temp_shared_w_garmaeva) / 2)
table_for_plot_m <- as.data.frame(colSums(temp_shared_w_maqsood) / 2)
table_for_plot_l <- as.data.frame(colSums(temp_shared_w_liang) / 2)
table_for_plot_s <- as.data.frame(colSums(temp_shared_w_shah) / 2)

table_for_plot_g_m <- merge(table_for_plot_g, table_for_plot_m, by="row.names", all = T)
table_for_plot_l_s <- merge(table_for_plot_l, table_for_plot_s, by="row.names", all = T)
table_for_plot <- merge(table_for_plot_g_m, table_for_plot_l_s, by="Row.names", all = T)
row.names(table_for_plot) <- table_for_plot$Row.names
table_for_plot$Row.names <- NULL

RPKM_count <- RPKM_filtered
RPKM_count[RPKM_count > 0] <- 1

table_for_plot_add <- as.data.frame(colSums(RPKM_count))
table_for_plot <- merge(table_for_plot, table_for_plot_add, by='row.names', all.x = T)

colnames(table_for_plot) <- c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                              "Shared_w_nc_shah", "Viral_richness")

table_for_plot$Perc_shared_w_garmaeva <- table_for_plot$Shared_w_nc_garmaeva / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_maqsood <- table_for_plot$Shared_w_nc_maqsood / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_liang <- table_for_plot$Shared_w_nc_liang / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_shah <- table_for_plot$Shared_w_nc_shah / table_for_plot$Viral_richness *100

table_for_plot <- merge(table_for_plot, meta_working[c("Sample_name", "Type", "Timepoint", "cohort")], by="Sample_name", all.x = T)

table_for_plot_melt <- melt(table_for_plot, id = c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                                                   "Shared_w_nc_shah", "Viral_richness", "Type", "Timepoint", "cohort"))

table_for_plot_melt$variable <- factor(table_for_plot_melt$variable, levels = c("Perc_shared_w_garmaeva", "Perc_shared_w_liang", 
                                                                                "Perc_shared_w_maqsood", "Perc_shared_w_shah"))

## Plotting only sharedness for infants
pdf('shared_inf_no_plasmids_no_lowundetermined_quality.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ ., scales = "free") +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

## Plotting only sharendess for mothers
pdf('shared_moms_no_plasmids_no_lowundetermined_quality.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Mother", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ .) +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

variable.labs <- c("Garmaeva NCs", "Maqsood NCs", "Liang NCs", "Shah NCs")
names(variable.labs) <- c("Perc_shared_w_garmaeva", "Perc_shared_w_maqsood", "Perc_shared_w_liang", "Perc_shared_w_shah")

# Convert named vector to labeller
variable_labeller <- as_labeller(variable.labs)

## Plotting sharedness depending on the timepoints
pdf('shared_timepoints_no_plasmids_no_lowundetermined_quality.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant" & table_for_plot_melt$cohort %in% c("garmaeva", "liang"), ], aes(x = Timepoint, y = value)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.13, size = 0.5) +
  ylim(0, 100) +
  facet_grid(rows = vars(variable), cols = vars(cohort), scales = "free", labeller = labeller(variable = variable_labeller)) +
  scale_x_discrete() + 
  labs(y = "%", x = "Timepoint") +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()





## Analysis for the vOTUs sharedness with the other cohorts: relative abundance
#######################################################################################################################################
# Get RA table

col_sums <- colSums(RPKM_filtered)
RA_filtered <- sweep(RPKM_filtered, 2, col_sums, FUN = "/")




RA_filtered_garmaeva_nc <- RA_filtered[row.names(RA_filtered) %in% ]















## Analysis for the vOTUs sharedness with different samples
#######################################################################################################################################

test <- RPKM_filtered_w_dummy[RPKM_filtered_w_dummy$dummy_all > 0, ]
test1 <- test[c("dummy_garmaeva", "dummy_liang", "dummy_maqsood", "dummy_shah")]
test1$shared_neg_ctrl <- rowSums(test1)
test1$shared_all_samples <- rowSums(test[colnames(test) %in% 
                                           meta_working$Sample_name[meta_working$Type != "Neg_ctrl"] ])
test1$shared_infants <- rowSums(test[colnames(test) %in% 
                                       meta_working$Sample_name[meta_working$Type == "Infant"] ])
test1$shared_mothers <- rowSums(test[colnames(test) %in% 
                                       meta_working$Sample_name[meta_working$Type == "Mother"] ])


test1$shared_infants_l <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type == "Infant" & meta_working$cohort == "liang"] ])

test1$shared_infants_s <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type == "Infant" & meta_working$cohort == "shah"] ])

test1$shared_infants_m <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type == "Infant" & meta_working$cohort == "maqsood"] ])

test1$shared_infants_g <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type == "Infant" & meta_working$cohort == "garmaeva"] ])

test1$shared_mother_m <- rowSums(test[colnames(test) %in% 
                                        meta_working$Sample_name[meta_working$Type == "Mother" & meta_working$cohort == "maqsood"] ])

test1$shared_mothers_g <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type == "Mother" & meta_working$cohort == "garmaeva"] ])


test1$shared_samples_l <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "liang"] ])

test1$shared_samples_s <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "shah"] ])

test1$shared_samples_m <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "maqsood"] ])

test1$shared_samples_g <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "garmaeva"] ])



## To do
#######################################################################################################################################
# Tune line 113
# Tune the low-unknown quality virus host prediction plot
# color the dots in the violin plot by timepoint
#######################################################################################################################################


# Redo all the analysis for the following






# New summary things -> are going directly after the data upload

extended_tof <- read.delim("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\data\\MERGED_Extended_TOF_NCP")





# REMOVED PARTS
# Create the plot and save it as a PDF
# pdf('summary_per_nc_absolute.pdf', width = 12/2.54, height = 8/2.54)
# ggplot(long_df, aes(x = sample, y = count, fill = checkv_quality_and_plasmid)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("#ef476f", "#f78c6b", "#06d6a0", "#ffd166", "#118ab2", "#073b4c")) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 4),  # Rotate x-ticks and decrease font size
#     axis.text.y = element_text(size = 4),  # Decrease y-axis font size
#     axis.title = element_text(size = 5),  # Decrease axis titles font size
#     legend.title = element_text(size = 5),  # Decrease legend title font size
#     legend.text = element_text(size = 4),  # Decrease legend text font size
#     legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
#   ) +
#   labs(x = "Sample", y = "Count", fill = "Virus Type") +
#   geom_text(data = total_counts, aes(x = sample, y = total_count, label = total_count), 
#             angle = 45, size = 1, inherit.aes = FALSE)  # Position text above the bars on a single line
# dev.off()
# 
# long_df$count <- long_df$count + 1
# 
# pdf('summary_per_nc_log.pdf', width = 12/2.54, height = 8/2.54)
# ggplot(long_df, aes(x = sample, y = count, fill = checkv_quality_and_plasmid)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c("#ef476f", "#f78c6b", "#06d6a0", "#ffd166", "#118ab2", "#073b4c")) +
#   scale_y_log10() +  # Set y-axis to log scale
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 4),  # Rotate x-ticks and decrease font size
#     axis.text.y = element_text(size = 4),  # Decrease y-axis font size
#     axis.title = element_text(size = 5),  # Decrease axis titles font size
#     legend.title = element_text(size = 5),  # Decrease legend title font size
#     legend.text = element_text(size = 4),  # Decrease legend text font size
#     legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
#   ) +
#   labs(x = "Sample", y = "Count", fill = "Virus Type") +
#   geom_text(data = total_counts, aes(x = sample, y = total_count, label = total_count), 
#             angle = 45, size = 1, inherit.aes = FALSE)  # Position text above the bars on a single line
# dev.off()
