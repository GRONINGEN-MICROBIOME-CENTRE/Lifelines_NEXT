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
library(ggtext)
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
         host_group = ifelse(grepl("Caudoviricetes", taxonomy), "Bacteria", "Unclassified"),
         host_group = ifelse(grepl("Herviviricetes", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Inoviridae", taxonomy), "Bacteria", host_group),
         host_group = ifelse(grepl("Microviridae", taxonomy), "Bacteria", host_group),
         host_group = ifelse(grepl("Genomoviridae", taxonomy), "Plant", host_group),
         host_group = ifelse(grepl("Circoviridae", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Nanoviridae", taxonomy), "Plant", host_group),
         host_group = ifelse(grepl("Geminiviridae", taxonomy), "Plant", host_group),
         host_group = ifelse(grepl("Smacoviridae", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Cossaviricota", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Nodaviridae", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Tombusviridae", taxonomy), "Plant", host_group),
         host_group = ifelse(grepl("Retroviridae", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Cystoviridae", taxonomy), "Bacteria", host_group),
         host_group = ifelse(grepl("Picornaviridae", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Astroviridae", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Caliciviridae", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Phycodnaviridae", taxonomy), "Plant", host_group),
         host_group = ifelse(grepl("Mimiviridae", taxonomy), "Animal(giant_virus)", host_group),
         host_group = ifelse(grepl("Adenoviridae", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Corticoviridae", taxonomy), "Bacteria", host_group),
         host_group = ifelse(grepl("Iridoviridae", taxonomy), "Bacteria", host_group),
         host_group = ifelse(grepl("Lavidaviridae", taxonomy), "viral phages", host_group),
         host_group = ifelse(grepl("Adintoviridae", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Tectiviridae", taxonomy), "Bacteria", host_group),
         host_group = ifelse(grepl("Sphaerolipoviridae", taxonomy), "Archaea", host_group),
         host_group = ifelse(grepl("Autolykiviridae", taxonomy), "Bacteria", host_group),
         host_group = ifelse(grepl("Marseilleviridae", taxonomy), "Animal(giant_virus)", host_group),
         host_group = ifelse(grepl("Poxviridae", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Anelloviridae", taxonomy), "Animal", host_group),
         host_group = ifelse(grepl("Portogloboviridae", taxonomy), "Archaea", host_group),
         host_group = ifelse(grepl("Bicaudaviridae", taxonomy), "Archaea", host_group)
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
  mutate(ncvssample = ifelse(Type == "Neg_ctrl", "NCs", "SAMPLES"),  # Changed for SAMPLES and NCs here; original - SAMPLE & NC
         timepoint_type = ifelse(Type == "Neg_ctrl", "NC", NA),
         timepoint_type = ifelse(Type == "Infant" & Timepoint %in% c("M0", "M1", "M2", "M3", "M4"), "Infant (age < 5 months)", timepoint_type),
         timepoint_type = ifelse(Type == "Infant" & Timepoint %in% c("M6", "M12", "Y2-5"), "Infant (age > 5 months)", timepoint_type),
         timepoint_type = ifelse(Type == "Mother", "Mother", timepoint_type))


pdf('alpha_diversity_comparison.pdf', width = 14/2.54, height = 8/2.54)
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

pdf('raw_read_comparison.pdf', width = 14/2.54, height = 8/2.54)
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

pdf('clean_read_comparison.pdf', width = 14/2.54, height = 8/2.54)
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

## Analysis for the vOTUs composition in NC vs samples: richness and ratio detected (CHM & palsmid)
#######################################################################################################################################

extended_tof <- extended_tof %>%
  mutate(checkv_quality_and_plasmid = ifelse(checkv_quality %in% c("Complete", "High-quality", "Medium-quality"), "CHM", NA),
         checkv_quality_and_plasmid = ifelse(checkv_quality %in% c("Not-determined", "Low-quality"), "LU", checkv_quality_and_plasmid),
         checkv_quality_and_plasmid = ifelse(plasmid == "Yes", "plasmid", checkv_quality_and_plasmid)
         )

RPKM_count <- RPKM
RPKM_count[RPKM_count > 0] <- 1

## Small code for the richness check in pilot - beginning
#############################################################################
RPKM_count_ch0 <- RPKM_count[, colnames(RPKM_count) %in% meta_working$Sample_name[meta_working$cohort == "garmaeva"]]
RPKM_count_ch0 <- RPKM_count_ch0[rowSums(RPKM_count_ch0) > 0, colSums(RPKM_count_ch0) > 0]
RPKM_count_ch0 <- merge(RPKM_count_ch0, extended_tof[colnames(extended_tof) %in% c("checkv_quality_and_plasmid")], by="row.names")
RPKM_count_ch0$Row.names <- NULL


summarized_df0 <- RPKM_count_ch0 %>%
  group_by(checkv_quality_and_plasmid) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

total_row <- summarized_df0 %>%
  summarise(across(-checkv_quality_and_plasmid, sum))
total_row <- cbind(checkv_quality_and_plasmid = "total", total_row)
summarized_df0 <- bind_rows(summarized_df0, total_row)

summarized_df0 <- as.data.frame(t(summarized_df0))
colnames(summarized_df0) <- c("CHM", "LU", "plasmid", "total")
summarized_df0 <- summarized_df0[row.names(summarized_df0) != "checkv_quality_and_plasmid", ]
row.names(meta_working) <- meta_working$Sample_name
summarized_df0 <- merge(summarized_df0, meta_working[c("Type", "Timepoint")], all.x=T, by="row.names")
summarized_df0a <- summarized_df0 %>%
  mutate(Timepoint = ifelse(Type == "Mother", "Mother", Timepoint))

summarized_df0a$Timepoint <- as.factor(summarized_df0a$Timepoint)
summarized_df0a$Timepoint <- factor(summarized_df0a$Timepoint, levels = c( "M1", "M2", "M3", "M6", "M12", "Mother", NA))

summarized_df0a[c("CHM", "LU", "plasmid", "total")] <- summarized_df0a[c("CHM", "LU", "plasmid", "total")] %>%
  mutate(across(where(is.character), ~ as.numeric(.)))


pdf('Total_number_pilot.pdf', width = 12/2.54, height = 8/2.54)
ggplot(summarized_df0a, aes(x=Timepoint, y=total)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter() +
  labs(y = "Richness of the total vOTUs in pilot") +
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


pdf('CHM_number_pilot.pdf', width = 12/2.54, height = 8/2.54)
ggplot(summarized_df0a, aes(x=Timepoint, y=CHM)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter() +
  labs(y = "Richness of the C/H/M quality vOTUs in pilot") +
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

pdf('Plasmid_number_pilot.pdf', width = 12/2.54, height = 8/2.54)
ggplot(summarized_df0a, aes(x=Timepoint, y=plasmid)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter() +
  labs(y = "Richness of the plasmid vOTUs in pilot") +
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

#############################################################################
## Small code for the richness check in pilot - end

RPKM_count_ch1 <- merge(RPKM_count, extended_tof[colnames(extended_tof) %in% c("checkv_quality_and_plasmid")], by="row.names")
# RPKM_count_ch1 <- RPKM_count_ch1[RPKM_count_ch1$checkv_quality_and_plasmid != "plasmid", ]
RPKM_count_ch1$Row.names <- NULL

summarized_df <- RPKM_count_ch1 %>%
  group_by(checkv_quality_and_plasmid) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df <- as.data.frame(t(summarized_df))
colnames(summarized_df) <- c("CHM", "LU", "plasmid")
summarized_df <- summarized_df[row.names(summarized_df) != "checkv_quality_and_plasmid", ]

summarized_df <- summarized_df %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df$CHM_LU_ratio <- summarized_df$CHM / (summarized_df$LU + summarized_df$CHM)
summarized_df$plamid_vOTU_ratio <- summarized_df$plasmid / (summarized_df$plasmid + summarized_df$LU + summarized_df$CHM)

row.names(meta_working) <- meta_working$Sample_name
meta_working <- merge(meta_working, summarized_df, by="row.names")

meta_working$Row.names <- NULL

pdf('CHM_LU_ratio_comparison.pdf', width = 14/2.54, height = 8/2.54)
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

pdf('CHM_number_comparison.pdf', width = 14/2.54, height = 8/2.54)
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



pdf('plasmid_vOTU_ratio_comparison.pdf', width = 14/2.54, height = 8/2.54)
ggplot(meta_working, aes(x=ncvssample, y=plamid_vOTU_ratio)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Richness ratio of the plasmids to all detected vOTUs \n (including plasmids)") +
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


pdf('plasmid_number_comparison.pdf', width = 14/2.54, height = 8/2.54)
ggplot(meta_working, aes(x=ncvssample, y=plasmid)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Richness of the plasmids") +
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


## Analysis for the vOTUs composition in NC vs samples: abundance and ratio (CHM & palsmid)
#######################################################################################################################################
RPKM_ch1 <- merge(RPKM, extended_tof[colnames(extended_tof) %in% c("checkv_quality_and_plasmid")], by="row.names")
# RPKM_ch1 <- RPKM_ch1[RPKM_ch1$checkv_quality_and_plasmid != "plasmid", ]
RPKM_ch1$Row.names <- NULL

summarized_df1 <- RPKM_ch1 %>%
  group_by(checkv_quality_and_plasmid) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df1 <- as.data.frame(t(summarized_df1))
colnames(summarized_df1) <- c("CHM_abundance", "LU_abundance", "plasmid_abundance")  # Check the order of the columns!!!
summarized_df1 <- summarized_df1[row.names(summarized_df1) != "checkv_quality_and_plasmid", ]

summarized_df1 <- summarized_df1 %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df1$CHM_LU_abundance_ratio <- summarized_df1$CHM_abundance / (summarized_df1$LU_abundance + summarized_df1$CHM_abundance)
summarized_df1$plasmid_vOTU_abundance_ratio <- summarized_df1$plasmid_abundance / (summarized_df1$plasmid_abundance + summarized_df1$LU_abundance + summarized_df1$CHM_abundance)

row.names(meta_working) <- meta_working$Sample_name
meta_working <- merge(meta_working, summarized_df1, by="row.names")

meta_working$Row.names <- NULL

pdf('CHM_LU_abundance_ratio_comparison.pdf', width = 14/2.54, height = 8/2.54)
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

pdf('CHM_abundance_comparison.pdf', width = 14/2.54, height = 8/2.54)
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


pdf('plasmid_vOTU_abundance_ratio_comparison.pdf', width = 14/2.54, height = 8/2.54)
ggplot(meta_working, aes(x=ncvssample, y=plasmid_vOTU_abundance_ratio)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Abundance ratio of the plasmids to all detected vOTUs \n (including plasmids)") +
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

pdf('plasmid_abundance_comparison.pdf', width = 14/2.54, height = 8/2.54)
ggplot(meta_working, aes(x=ncvssample, y=plasmid_abundance)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Abundance of the plasmids") +
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


## Analysis for the vOTUs composition in NC vs samples: richness and ratio discovered (CHM & palsmid)
#######################################################################################################################################

get_study_and_sample_of_origin <- function(cid) {
  parts <- strsplit(cid, "_")[[1]]
  num_parts <- length(parts)
  sample_of_origin <- paste(parts[1:(num_parts - 6)], collapse = "_")
  split_parts <- strsplit(sample_of_origin, "_", fixed = TRUE)[[1]]
  study_of_origin <- split_parts[1]
  sample_of_origin <- paste(split_parts[-1], collapse = "_")
  
  return(c(study_of_origin, sample_of_origin))
}

# Applying the function to create the new columns
split_results <- t(sapply(extended_tof$New_CID, get_study_and_sample_of_origin))

# Adding the new columns to the dataframe
extended_tof$study_of_origin <- split_results[, 1]
extended_tof$sample_of_origin <- split_results[, 2]

summarized_discovered_vOTUs <- extended_tof %>% group_by(sample_of_origin, study_of_origin, checkv_quality_and_plasmid) %>% summarise(count=n())

summarized_discovered_vOTUs <- summarized_discovered_vOTUs %>%
  pivot_wider(names_from = checkv_quality_and_plasmid, values_from = count, values_fill = list(count = 0))

summarized_discovered_vOTUs <- as.data.frame(summarized_discovered_vOTUs)

colnames(summarized_discovered_vOTUs) <- c("Sample_name", "cohort", "CHM_discovered", "LU_discovered", "plasmid_discovered")

summarized_discovered_vOTUs$CHM_LU_richness_discovered_ratio <- summarized_discovered_vOTUs$CHM_discovered / (summarized_discovered_vOTUs$LU_discovered + summarized_discovered_vOTUs$CHM_discovered)
summarized_discovered_vOTUs$plasmid_richness_discovered_ratio <- summarized_discovered_vOTUs$plasmid_discovered / (summarized_discovered_vOTUs$plasmid_discovered + summarized_discovered_vOTUs$LU_discovered + summarized_discovered_vOTUs$CHM_discovered)
row.names(summarized_discovered_vOTUs) <- summarized_discovered_vOTUs$Sample_name

row.names(meta_working) <- meta_working$Sample_name
meta_working <- merge(meta_working, summarized_discovered_vOTUs[, !colnames(summarized_discovered_vOTUs) %in% c("cohort", "Sample_name")], all.x = T, by="row.names")
meta_working$Row.names <- NULL

pdf('CHM_LU_richness_discovered_ratio_comparison.pdf', width = 14/2.54, height = 8/2.54)
ggplot(meta_working[!(meta_working$Sample_name %in% c("bctrl4633v", "bctrl4654v", "bctrl4676v", "bctrl4699v")), ], aes(x=ncvssample, y=CHM_LU_richness_discovered_ratio)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Richness ratio of the discovered C/H/M quality vOTUs to all discovered vOTUs \n (excluding plasmids and orsay controls)") +
  theme_bw() +
  theme(plot.title = element_text(size=8),
        strip.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 4))
dev.off()

pdf('CHM_richness_discovered_comparison.pdf', width = 14/2.54, height = 8/2.54)
ggplot(meta_working[!(meta_working$Sample_name %in% c("bctrl4633v", "bctrl4654v", "bctrl4676v", "bctrl4699v")), ], aes(x=ncvssample, y=CHM_discovered)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Richness of the discovered C/H/M quality vOTUs \n (excluding plasmids and orsay controls)") +
  theme_bw() +
  theme(plot.title = element_text(size=8),
        strip.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 4))
dev.off()


pdf('plasmid_vOTU_richness_discovered_ratio_comparison.pdf', width = 14/2.54, height = 8/2.54)
ggplot(meta_working, aes(x=ncvssample, y=plasmid_richness_discovered_ratio)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Richness ratio of the discovered plasmids to all discovered vOTUs \n (including plasmids)") +
  theme_bw() +
  theme(plot.title = element_text(size=8),
        strip.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 4))
dev.off()

pdf('plasmid_richness_discovered_comparison.pdf', width = 14/2.54, height = 8/2.54)
ggplot(meta_working, aes(x=ncvssample, y=plasmid_discovered)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Richness of the discovered plasmidss") +
  theme_bw() +
  theme(plot.title = element_text(size=8),
        strip.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 4))
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

pdf('mean_virus_group_richness_CHM.pdf', width = 14/2.54, height = 8/2.54)
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

pdf('mean_virus_group_abundance_CHM.pdf', width = 14/2.54, height = 8/2.54)
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
colnames(summarized_df3) <- c("Unclassified", "Animal", "Archaea", "Bacteria", "Plant")
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
summarized_df3_frac_melt$variable <- factor(summarized_df3_frac_melt$variable, levels = c("Bacteria", "Archaea", "Animal", "Plant", "Unclassified"))

## Mention in the plot, that only CHM viruses included here - no plasmids, no LU

pdf('mean_host_group_richness_CHM.pdf', width = 14/2.54, height = 8/2.54)
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
colnames(summarized_df3a) <- c("Unclassified", "Animal", "Archaea", "Bacteria", "Plant")
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
summarized_df3a_frac_melt$variable <- factor(summarized_df3a_frac_melt$variable, levels = c("Bacteria", "Archaea", "Animal", "Plant", "Unclassified"))

pdf('mean_host_group_abundance_CHM.pdf', width = 14/2.54, height = 8/2.54)
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
BU_viruses <- row.names(extended_tof[extended_tof$host_group %in% c("Bacteria", "Archaea", "Unclassified"), ])
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
summarized_df4_frac_melt$variable <- factor(summarized_df4_frac_melt$variable, levels = c("Bacteria", "Archaea", "Animal", "Plant", "Unclassified"))
summarized_df4_frac_melt$variable <- factor(summarized_df4_frac_melt$variable, 
                                            levels = levels(summarized_df4_frac_melt$variable)[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 30)])

pdf('mean_host_bacteria_richness_CHM.pdf', width = 15/2.54, height = 8/2.54)
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

## Bray-Curtis boxplots
#######################################################################################################################################

meta_sharedness <- meta_working %>%
  mutate(Subject_ID = ifelse(grepl("kid", Sample_name), Sample_name, Subject_ID)) %>%
  filter(!is.na(Subject_ID)) %>%
  group_by(Subject_ID) %>%
  slice_sample(n = 1) %>%
  ungroup()

meta_sharedness_wnc <- rbind(meta_sharedness, meta_working[meta_working$Type == "Neg_ctrl", ])

RPKM_one_tppi <- RPKM[, colnames(RPKM) %in% meta_sharedness_wnc$Sample_name]  # one timepoint per person (randomly picked)
RPKM_one_tppi <- RPKM_one_tppi[rowSums(RPKM_one_tppi) > 0, colSums(RPKM_one_tppi) > 0]

bray_dist_matrix_one_tppi_wnc <- as.matrix(vegdist(t(RPKM_one_tppi), method="bray"))
bray_dist_matrix_one_tppi_wnc_rev <- 1 - bray_dist_matrix_one_tppi_wnc

upper_tri <- upper.tri(bray_dist_matrix_one_tppi_wnc_rev)

# Get row and column indices for upper triangle
row_indices <- row(bray_dist_matrix_one_tppi_wnc_rev)[upper_tri]
col_indices <- col(bray_dist_matrix_one_tppi_wnc_rev)[upper_tri]

# Extract distances corresponding to upper triangle indices
distances <- bray_dist_matrix_one_tppi_wnc_rev[upper_tri]

# Create a dataframe
df_distances <- data.frame(
  Sample1 = rownames(bray_dist_matrix_one_tppi_wnc_rev)[row_indices],
  Sample2 = colnames(bray_dist_matrix_one_tppi_wnc_rev)[col_indices],
  Distance = distances
)

# Ensure Sample1 < Sample2 in each row (order does not matter)
df_distances <- df_distances[order(pmin(df_distances$Sample1, df_distances$Sample2)),
]  # Sorting by Sample1, Sample2

meta_sharedness_wnc <- as.data.frame(meta_sharedness_wnc)
meta_sharedness_wnc$Sample1 <- meta_sharedness_wnc$Sample_name
meta_sharedness_wnc$Sample2 <- meta_sharedness_wnc$Sample_name

df_distances <- merge(df_distances, meta_sharedness_wnc[c("Sample1", "ncvssample", "cohort")], by="Sample1", all.x=T)
df_distances <- merge(df_distances, meta_sharedness_wnc[c("Sample2", "ncvssample", "cohort")], by="Sample2", all.x=T)
colnames(df_distances) <- c("Sample1", "Sample2", "Distance", "Sample1_ncvssample",  "Sample1_cohort", "Sample2_ncvssample",  "Sample2_cohort")

df_distances_nc <- df_distances[!(df_distances$Sample1_ncvssample == "SAMPLES" & df_distances$Sample2_ncvssample == "SAMPLES"), ]
df_distances_nc <- df_distances_nc %>%
  mutate(category = ifelse(Sample1_ncvssample == "NCs" & Sample2_ncvssample == "SAMPLES", "Samples", NA),
         category = ifelse(Sample2_ncvssample == "NCs" & Sample1_ncvssample == "SAMPLES", "Samples", category),
         category = ifelse(Sample2_ncvssample == "NCs" & Sample1_ncvssample == "NCs", "NCs", category),
         type_cohort = ifelse(Sample1_cohort == Sample2_cohort, "same", "different"),
         cohort_nc = ifelse(type_cohort == "same" & category == "NCs",  Sample1_cohort, NA),
         cohort_nc = ifelse(category == "Samples" & Sample1_ncvssample == "NCs",  Sample1_cohort, cohort_nc),
         cohort_nc = ifelse(category == "Samples" & Sample2_ncvssample == "NCs",  Sample2_cohort, cohort_nc))

switched_df <- df_distances_nc %>%
  filter(type_cohort == "different" & category == "NCs") %>%
  mutate(Sample3 = Sample1, 
         Sample1 = Sample2,
         Sample2 = Sample3,
         Sample3 = NULL,
         Sample3_cohort = Sample1_cohort,
         Sample1_cohort = Sample2_cohort,
         Sample2_cohort = Sample3_cohort,
         Sample3_cohort = NULL)

df_distances_nc <- rbind(df_distances_nc, switched_df)

df_distances_nc <- df_distances_nc %>%
  mutate(cohort_nc = ifelse(category == "NCs" & type_cohort == "different",  Sample1_cohort, cohort_nc))


pdf('bray-curtis_similarity_boxplot.pdf', width=12/2.54, height=10.5/2.54)
ggplot(df_distances_nc, aes(x = type_cohort, y = Distance)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outliers = FALSE) +
  ylim(0, 1) +
  facet_grid(cohort_nc ~ category, scales = "free", labeller = labeller(
    cohort_nc = c(
      "garmaeva" = "NCs Garmaeva *et al.* <br>",
      "liang" = "NCs Liang *et al.* <br>",
      "maqsood" = "NCs Maqsood *et al.* <br>",
      "shah" = "NCs Shah *et al.* <br>"
    ),
    category = c(
      "NCs" = "NCs and NCs",
      "Samples" = "NCs and Samples"
    )
  )) +
  labs(y = "Bray-Curtis similarity", x = "Type of cohort \n (samples from the same or different cohort compared)") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    strip.background = element_rect(fill = "transparent"),
  )
dev.off()

#######################################################################################################################################


## Jackrat boxplots
#######################################################################################################################################

RPKM_one_tppi_count <- RPKM_one_tppi
RPKM_one_tppi_count[RPKM_one_tppi_count > 0] <- 1

jaccard_dist_matrix_one_tppi_wnc <- as.matrix(vegdist(t(RPKM_one_tppi), method="jaccard"))
jaccard_dist_matrix_one_tppi_wnc_rev <- 1 - jaccard_dist_matrix_one_tppi_wnc

upper_tri <- upper.tri(jaccard_dist_matrix_one_tppi_wnc_rev)

# Get row and column indices for upper triangle
row_indices <- row(jaccard_dist_matrix_one_tppi_wnc_rev)[upper_tri]
col_indices <- col(jaccard_dist_matrix_one_tppi_wnc_rev)[upper_tri]

# Extract distances corresponding to upper triangle indices
distances <- jaccard_dist_matrix_one_tppi_wnc_rev[upper_tri]

# Create a dataframe
df_distances <- data.frame(
  Sample1 = rownames(jaccard_dist_matrix_one_tppi_wnc_rev)[row_indices],
  Sample2 = colnames(jaccard_dist_matrix_one_tppi_wnc_rev)[col_indices],
  Distance = distances
)

# Ensure Sample1 < Sample2 in each row (order does not matter)
df_distances <- df_distances[order(pmin(df_distances$Sample1, df_distances$Sample2)),
]  # Sorting by Sample1, Sample2

meta_sharedness_wnc <- as.data.frame(meta_sharedness_wnc)
meta_sharedness_wnc$Sample1 <- meta_sharedness_wnc$Sample_name
meta_sharedness_wnc$Sample2 <- meta_sharedness_wnc$Sample_name

df_distances <- merge(df_distances, meta_sharedness_wnc[c("Sample1", "ncvssample", "cohort")], by="Sample1", all.x=T)
df_distances <- merge(df_distances, meta_sharedness_wnc[c("Sample2", "ncvssample", "cohort")], by="Sample2", all.x=T)
colnames(df_distances) <- c("Sample1", "Sample2", "Distance", "Sample1_ncvssample",  "Sample1_cohort", "Sample2_ncvssample",  "Sample2_cohort")

df_distances_nc <- df_distances[!(df_distances$Sample1_ncvssample == "SAMPLES" & df_distances$Sample2_ncvssample == "SAMPLES"), ]
df_distances_nc <- df_distances_nc %>%
  mutate(category = ifelse(Sample1_ncvssample == "NCs" & Sample2_ncvssample == "SAMPLES", "Samples", NA),
         category = ifelse(Sample2_ncvssample == "NCs" & Sample1_ncvssample == "SAMPLES", "Samples", category),
         category = ifelse(Sample2_ncvssample == "NCs" & Sample1_ncvssample == "NCs", "NCs", category),
         type_cohort = ifelse(Sample1_cohort == Sample2_cohort, "same", "different"),
         cohort_nc = ifelse(type_cohort == "same" & category == "NCs",  Sample1_cohort, NA),
         cohort_nc = ifelse(category == "Samples" & Sample1_ncvssample == "NCs",  Sample1_cohort, cohort_nc),
         cohort_nc = ifelse(category == "Samples" & Sample2_ncvssample == "NCs",  Sample2_cohort, cohort_nc))

switched_df <- df_distances_nc %>%
  filter(type_cohort == "different" & category == "NCs") %>%
  mutate(Sample3 = Sample1, 
         Sample1 = Sample2,
         Sample2 = Sample3,
         Sample3 = NULL,
         Sample3_cohort = Sample1_cohort,
         Sample1_cohort = Sample2_cohort,
         Sample2_cohort = Sample3_cohort,
         Sample3_cohort = NULL)

df_distances_nc <- rbind(df_distances_nc, switched_df)

df_distances_nc <- df_distances_nc %>%
  mutate(cohort_nc = ifelse(category == "NCs" & type_cohort == "different",  Sample1_cohort, cohort_nc))


pdf('jaccart_similarity_boxplot.pdf', width=12/2.54, height=10.5/2.54)
ggplot(df_distances_nc, aes(x = type_cohort, y = Distance)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outliers = FALSE) +
  ylim(0, 1) +
  facet_grid(cohort_nc ~ category, scales = "free", labeller = labeller(
    cohort_nc = c(
      "garmaeva" = "NCs Garmaeva *et al.* <br>",
      "liang" = "NCs Liang *et al.* <br>",
      "maqsood" = "NCs Maqsood *et al.* <br>",
      "shah" = "NCs Shah *et al.* <br>"
    ),
    category = c(
      "NCs" = "NCs and NCs",
      "Samples" = "NCs and Samples"
    )
  )) +
  labs(y = "Bray-Curtis similarity", x = "Type of cohort \n (samples from the same or different cohort compared)") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    strip.background = element_rect(fill = "transparent"),
  )
dev.off()

#######################################################################################################################################

## Analysis for the vOTUs sharedness with the other cohorts: add the sharedness with the samples: wo dummy 
#######################################################################################################################################
negative_controls <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" &  meta_working$Sample_name %in% colnames(RPKM)]

RPKM_filtered <- RPKM_count[, colnames(RPKM_count) %in% meta_sharedness_wnc$Sample_name]
RPKM_filtered$dummy_all <- rowSums(RPKM_filtered[, colnames(RPKM_filtered) %in% negative_controls])
RPKM_filtered_vOTUsshared <- RPKM_filtered[RPKM_filtered$dummy_all > 0, ]
samples_0_shared_to_any_NC <- colnames(RPKM_filtered_vOTUsshared)[colSums(RPKM_filtered_vOTUsshared) == 0]
RPKM_filtered$dummy_all <- NULL

samples_not_shared_garmaeva <- sample(samples_0_shared_to_any_NC[samples_0_shared_to_any_NC %in% meta_sharedness_wnc$Sample_name[meta_sharedness_wnc$cohort == "garmaeva"]], 1)
samples_not_shared_liang <- sample(samples_0_shared_to_any_NC[samples_0_shared_to_any_NC %in% meta_sharedness_wnc$Sample_name[meta_sharedness_wnc$cohort == "liang"]], 26)
samples_not_shared_maqsood <- sample(samples_0_shared_to_any_NC[samples_0_shared_to_any_NC %in% meta_sharedness_wnc$Sample_name[meta_sharedness_wnc$cohort == "maqsood"]], 8)
samples_not_shared_shah <- sample(samples_0_shared_to_any_NC[samples_0_shared_to_any_NC %in% meta_sharedness_wnc$Sample_name[meta_sharedness_wnc$cohort == "shah"]], 8)

nc_samples0share <- c(negative_controls, samples_not_shared_garmaeva, samples_not_shared_liang, samples_not_shared_maqsood, samples_not_shared_shah)
non_nc_samples0share <- setdiff(names(RPKM_filtered), nc_samples0share)

# Initialize an empty dataframe for the results
result <- data.frame(matrix(ncol = length(nc_samples0share), nrow = length(non_nc_samples0share)))
names(result) <- nc_samples0share
row.names(result) <- non_nc_samples0share

# Calculate the percent similarity
for (non_dummy_col in non_nc_samples0share) {
  for (dummy_col in nc_samples0share) {
    # Number of viruses present in the non-dummy column
    present_in_non_dummy <- sum(RPKM_filtered[[non_dummy_col]] == 1)
    
    if (present_in_non_dummy > 0) {
      # Number of viruses present in both non-dummy and dummy column
      present_in_both <- sum(RPKM_filtered[[non_dummy_col]] == 1 & RPKM_filtered[[dummy_col]] == 1)
      
      # Calculate the percent similarity
      similarity <- present_in_both / present_in_non_dummy * 100
    } else {
      similarity <- NA  # or 0 if you prefer
    }
    
    # Store the result
    result[non_dummy_col, dummy_col] <- similarity
  }
}

result$Sample_name <- row.names(result)
row.names(result) <- NULL
table_for_plot <- merge(result, meta_working[c("Sample_name", "Type", "cohort")], by="Sample_name", all.x = T)

table_for_plot_melt <- melt(table_for_plot, id = c("Sample_name", "Type", "cohort"))
table_for_plot_melt$Sample_name <- NULL
table_for_plot_melt$Sample_name <- table_for_plot_melt$variable
table_for_plot_melt <- merge(table_for_plot_melt, meta_working[c("Sample_name", "ncvssample", "cohort")], by="Sample_name", all.x = T, suffixes = c("", "_origin"))

table_for_plot_melt$cohort_origin <- factor(table_for_plot_melt$cohort_origin, levels = c("garmaeva", "liang",
                                                                                          "maqsood", "shah"))

## Plotting only sharedness for infants
pdf('shared_infv2.pdf', width=12/2.54, height=10.5/2.54)
ggplot(subset(table_for_plot_melt, Type == "Infant"), aes(x = ncvssample, y = value)) +
  #  geom_violin(trim = FALSE, fill = "#C8ACD6", alpha=0.3) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outliers = FALSE) +
  ylim(0, 100) +
  facet_grid(cohort ~ cohort_origin, scales = "free", labeller = labeller(
    cohort_origin = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"
    ),
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>",
      "shah" = "Samples Shah *et al.* <br>"
    )
  )) +
  labs(y = "% shared vOTUs") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    plot.title = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6), 
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "transparent")
  )
dev.off()



pdf('shared_momv2.pdf', width=12/2.54, height=7/2.54)
ggplot(subset(table_for_plot_melt, Type == "Mother"), aes(x = shared_w, y = value)) +
  #  geom_violin(trim = FALSE, fill = "#C8ACD6", alpha=0.3) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outliers = FALSE) +
  ylim(0, 100) +
  facet_grid(cohort ~ cohort_origin, scales = "free", labeller = labeller(
    cohort_origin = c(
      "garmaeva" = "NC Garmaeva *et al.* <br>",
      "liang" = "NC Liang *et al.* <br>",
      "maqsood" = "NC Maqsood *et al.* <br>",
      "shah" = "NC Shah *et al.* <br>"
    ),
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>",
      "shah" = "Samples Shah *et al.* <br>"
    )
  )) +
  labs(y = "% shared vOTUs") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    plot.title = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6), 
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "transparent")
  )
dev.off()

#######################################################################################################################################
