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
library(VennDiagram)
library(lmerTest)
library(patchwork)
library(MuMIn)
#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
setwd("/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/plots")
#######################################################################################################################################

## Load the  data
#######################################################################################################################################
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
         host_group = ifelse(grepl("Caudoviricetes", taxonomy), "Prokaryotes", "Unclassified"),
         host_group = ifelse(grepl("Herviviricetes", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Inoviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Microviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Genomoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Circoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Nanoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Geminiviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Smacoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Cossaviricota", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Nodaviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Tombusviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Retroviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Cystoviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Picornaviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Astroviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Caliciviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Phycodnaviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Mimiviridae", taxonomy), "Eukaryotes(giant_virus)", host_group),
         host_group = ifelse(grepl("Adenoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Corticoviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Iridoviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Lavidaviridae", taxonomy), "viral phages", host_group),
         host_group = ifelse(grepl("Adintoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Tectiviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Sphaerolipoviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Autolykiviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Marseilleviridae", taxonomy), "Eukaryotes(giant_virus)", host_group),
         host_group = ifelse(grepl("Poxviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Anelloviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Portogloboviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Bicaudaviridae", taxonomy), "Prokaryotes", host_group),
         checkv_quality_and_plasmid = ifelse(checkv_quality %in% c("Complete", "High-quality", "Medium-quality"), "CHM", NA),
         checkv_quality_and_plasmid = ifelse(checkv_quality %in% c("Not-determined", "Low-quality"), "LU", checkv_quality_and_plasmid),
         checkv_quality_and_plasmid = ifelse(plasmid == "Yes", "plasmid", checkv_quality_and_plasmid)
  )

RPKM_contigs_keep <- extended_tof$New_CID[extended_tof$POST_CHV_length >= 1000 & (extended_tof$viral_genes >= extended_tof$host_genes) & extended_tof$plasmid=="No"]

RPKM_initial <- read.delim("../RPKM_counts_VLP_NCP.txt")
dim(RPKM_initial)  # 306275   1301

RPKM <- RPKM_initial[row.names(RPKM_initial) %in% RPKM_contigs_keep, ]
dim(RPKM)  # 193970   1301
RPKM <- RPKM[rowSums(RPKM) > 0, colSums(RPKM) > 0]
dim(RPKM)  # 193970   1291

RPKM_count <- RPKM
RPKM_count[RPKM_count > 0] <- 1

meta_all_with_qc_curated <- as.data.frame(read_tsv('../../metadata_with_qc_NCPv2.tsv'))
dim(meta_all_with_qc_curated)  # 1376   28
meta_all_with_qc_curated <- meta_all_with_qc_curated %>%
  mutate(Subject_ID = ifelse(grepl("kid", Sample_name), Sample_name, Subject_ID),
         nc_subject_group = ifelse(grepl("bctrl1", Sample_name), "NC_maqsood_buffer", "SAMPLE"),
         nc_subject_group = ifelse(grepl("bctrl4", Sample_name), "NC_maqsood_orsay", nc_subject_group),
         nc_subject_group = ifelse(grepl("ctl", Sample_name), "NC_shah_buffer", nc_subject_group),
         nc_subject_group = ifelse(grepl("LN_7C08_VL_405", Sample_name), "NC_garmaeva_buffer", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LT1_D", "LT4_D", "LT5_D", "LT2_D", "LT3_D"), "NC_liang_tube", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LT1_R", "LT2_R", "LT3_R", "LT4_R", "LT5_R"), "NC_liang_tube", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LB1_D", "LB2_D", "LB3_D", "LB4_D", "LB5_D", "LNCB3_D", "LNCB2_D", "LNCB1_D"), "NC_liang_reagent", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LB1_R", "LB2_R", "LB3_R", "LB4_R", "LB5_R", "LNCB1_R", "LNCB3_R", "LNCB2_R"), "NC_liang_reagent", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LD1_D", "LD2_D", "LDH_D", "LDN1_D", "LDN2_D"), "NC_liang_diaper", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LD1_R", "LD2_R", "LDH_R", "LDN1_R", "LDN2_R"), "NC_liang_diaper", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LMDNC_D", "LMDNC_R"), "NC_liang_buffer", nc_subject_group),
         rna_dna = ifelse(grepl("_D", Sample_ID), "DNA", NA),
         rna_dna = ifelse(grepl("_R", Sample_ID), "RNA", rna_dna),
         ncvssample = ifelse(Type == "Neg_ctrl", "NCs", "SAMPLES"),  # Changed for SAMPLES and NCs here; original - SAMPLE & NC
         nc_subject_group = ifelse(ncvssample == "SAMPLES", Subject_ID, nc_subject_group),
         timepoint_type = ifelse(Type == "Neg_ctrl", "NC", NA),
         timepoint_type = ifelse(Type == "Infant" & Timepoint %in% c("M0", "M1", "M2", "M3", "M4"), "Infant (age < 5 months)", timepoint_type),
         timepoint_type = ifelse(Type == "Infant" & Timepoint %in% c("M6", "M12", "Y2-5"), "Infant (age > 5 months)", timepoint_type),
         timepoint_type = ifelse(Type == "Mother", "Mother", timepoint_type)
  )

meta_all_with_qc_curated$nc_subject_group = as.factor(meta_all_with_qc_curated$nc_subject_group)
meta_all_with_qc_curated$cohort = as.factor(meta_all_with_qc_curated$cohort)

meta_working <- meta_all_with_qc_curated[meta_all_with_qc_curated$Sample_name %in% colnames(RPKM), ]

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

## Adding alpha diversity and richness to meta_working
#######################################################################################################################################
row.names(meta_working) <- meta_working$Sample_name
diversity_df <- as.data.frame(diversity(as.data.frame(t(RPKM)), index='shannon'))
meta_working <- merge(meta_working, diversity_df, by="row.names")
colnames(meta_working)[colnames(meta_working) == "diversity(as.data.frame(t(RPKM)), index = \"shannon\")"] <- "diversity"
meta_working$Row.names <- NULL
richness_table <- as.data.frame(colSums(RPKM_count))
row.names(meta_working) <- meta_working$Sample_name
meta_working <- merge(meta_working, richness_table, by="row.names")
colnames(meta_working)[colnames(meta_working) == "colSums(RPKM_count)"] <- "richness"
meta_working$Row.names <- NULL
#######################################################################################################################################

## Analysis for the vOTUs composition in NC vs samples: richness and CHM ratio detected
#######################################################################################################################################
RPKM_count_ch1 <- merge(RPKM_count, extended_tof[colnames(extended_tof) %in% c("checkv_quality_and_plasmid")], by="row.names")
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
#######################################################################################################################################

## Analysis for the vOTUs composition in NC vs samples: richness and CHM ratio discovered
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

summarized_discovered_vOTUs <- extended_tof %>%
  filter(POST_CHV_length >= 1000 & plasmid == "No") %>%
  group_by(sample_of_origin, study_of_origin, checkv_quality_and_plasmid) %>%
  summarise(count = n())

summarized_discovered_vOTUs <- summarized_discovered_vOTUs %>%
  pivot_wider(names_from = checkv_quality_and_plasmid, values_from = count, values_fill = list(count = 0))

summarized_discovered_vOTUs <- as.data.frame(summarized_discovered_vOTUs)
colnames(summarized_discovered_vOTUs) <- c("Sample_name", "cohort", "CHM_discovered", "LU_discovered")

summarized_discovered_vOTUs$CHM_LU_richness_discovered_ratio <- summarized_discovered_vOTUs$CHM_discovered / (summarized_discovered_vOTUs$LU_discovered + summarized_discovered_vOTUs$CHM_discovered)
row.names(summarized_discovered_vOTUs) <- summarized_discovered_vOTUs$Sample_name

row.names(meta_working) <- meta_working$Sample_name
meta_working <- merge(meta_working, summarized_discovered_vOTUs[, !colnames(summarized_discovered_vOTUs) %in% c("cohort", "Sample_name")], all.x = T, by="row.names")
meta_working$Row.names <- NULL
#######################################################################################################################################

## Analysis for the virus groups richness + STATS
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
summarized_df2_frac <- merge(summarized_df2_frac, meta_working[, colnames(meta_working) %in% c("Sample_name", "ncvssample", "cohort", "nc_subject_group")], all.x = T, by="Sample_name")
summarized_df2_frac$Sample_name <- NULL

summarized_df2_frac_stats <- melt(summarized_df2_frac, id.vars = c("ncvssample", "cohort", "nc_subject_group"))

summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "liang" & summarized_df2_frac_stats$variable == "RNA", ]))
summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "liang" & summarized_df2_frac_stats$variable == "dsDNA", ]))
summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "liang" & summarized_df2_frac_stats$variable == "ssDNA", ]))

summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "maqsood" & summarized_df2_frac_stats$variable == "RNA", ]))
summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "maqsood" & summarized_df2_frac_stats$variable == "dsDNA", ]))
summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "maqsood" & summarized_df2_frac_stats$variable == "ssDNA", ]))

summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "shah" & summarized_df2_frac_stats$variable == "RNA", ]))
summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "shah" & summarized_df2_frac_stats$variable == "dsDNA", ]))
summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "shah" & summarized_df2_frac_stats$variable == "ssDNA", ]))

p.adjust(c(0.421, 0.0229, 0.0628, 0.6409, 0.895524, 5.75e-13, 0.245, 0.955), method="BH")

summarized_df2_frac <- summarized_df2_frac %>%
  group_by(ncvssample, cohort) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

summarized_df2_frac <- as.data.frame(summarized_df2_frac)

## Redo for long format
summarized_df2_frac_melt <- melt(summarized_df2_frac, id.vars = c("ncvssample", "cohort"))
summarized_df2_frac_melt$variable <- as.factor(summarized_df2_frac_melt$variable)
summarized_df2_frac_melt$variable <- factor(summarized_df2_frac_melt$variable, levels = c("RNA", "ssDNA", "dsDNA", "Unclassified"))
#######################################################################################################################################

## Analysis for the virus host richness -> prokaryotes vs eukaryotes + stats
#######################################################################################################################################
RPKM_count_ch3 <- RPKM_count[row.names(RPKM_count) %in% CHM_viruses, ]
RPKM_count_ch3 <- merge(RPKM_count_ch3, extended_tof[colnames(extended_tof) %in% c("host_group")], by="row.names")
RPKM_count_ch3$Row.names <- NULL

summarized_df3 <- RPKM_count_ch3 %>%
  group_by(host_group) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df3 <- as.data.frame(t(summarized_df3))
colnames(summarized_df3) <- c("Eukaryotes", "Prokaryotes", "Unclassified")
summarized_df3 <- summarized_df3[row.names(summarized_df3) != "host_group", ]

summarized_df3 <- summarized_df3 %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df3 <- summarized_df3[rowSums(summarized_df3) > 0, ]

summarized_df3_frac <- as.data.frame(t(apply(summarized_df3, 1, fraction_per_sample)))

summarized_df3_frac$Sample_name <- row.names(summarized_df3_frac)
summarized_df3_frac <- merge(summarized_df3_frac, meta_working[, colnames(meta_working) %in% c("Sample_name", "ncvssample", "cohort", "nc_subject_group")], all.x = T, by="Sample_name")
summarized_df3_frac$Sample_name <- NULL

summarized_df3_frac_stats <- melt(summarized_df3_frac, id.vars = c("ncvssample", "cohort", "nc_subject_group"))

summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "liang" & summarized_df3_frac_stats$variable == "Eukaryotes", ]))
summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "liang" & summarized_df3_frac_stats$variable == "Prokaryotes", ]))

summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "maqsood" & summarized_df3_frac_stats$variable == "Eukaryotes", ]))
summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "maqsood" & summarized_df3_frac_stats$variable == "Prokaryotes", ]))

summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "shah" & summarized_df3_frac_stats$variable == "Eukaryotes", ]))
summary(lmer(value ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "shah" & summarized_df3_frac_stats$variable == "Prokaryotes", ]))
p.adjust(c(0.235, 0.195, 0.355, 0.159, 0.3878, 5.05e-05), method="BH")

# Answering if prokaryotes dominating eukaryotes

summary(lmer(Prokaryotes ~ Eukaryotes + (1|nc_subject_group), REML = F, data = summarized_df3_frac[summarized_df3_frac_stats$cohort == "liang", ]))
summary(lmer(Prokaryotes ~ Eukaryotes + (1|nc_subject_group), REML = F, data = summarized_df3_frac[summarized_df3_frac_stats$cohort == "maqsood", ]))
summary(lmer(Prokaryotes ~ Eukaryotes + (1|nc_subject_group), REML = F, data = summarized_df3_frac[summarized_df3_frac_stats$cohort == "shah", ]))
p.adjust(c(2e-16, 2.16e-11, 2e-16), method="BH")

summarized_df3_frac <- summarized_df3_frac %>%
  group_by(ncvssample, cohort) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

summarized_df3_frac <- as.data.frame(summarized_df3_frac)

## Redo for long format
summarized_df3_frac_melt <- melt(summarized_df3_frac, id.vars = c("ncvssample", "cohort"))
summarized_df3_frac_melt$variable <- as.factor(summarized_df3_frac_melt$variable)
summarized_df3_frac_melt$variable <- factor(summarized_df3_frac_melt$variable, levels = c("Prokaryotes", "Eukaryotes", "Unclassified"))
#######################################################################################################################################

## Analysis for the phages host richness
#######################################################################################################################################
BU_viruses <- row.names(extended_tof[extended_tof$host_group %in% c("Prokaryotes", "Unclassified"), ])
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

summarized_df4_frac_melt$variable <- factor(summarized_df4_frac_melt$variable, 
                                            levels = levels(summarized_df4_frac_melt$variable)[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 31)])

#######################################################################################################################################

## Patching the Figure 1
#######################################################################################################################################
figure_1A <- ggplot(meta_working, aes(x=ncvssample, y=clean_reads_comb)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  scale_y_log10() +
  labs(y = "Number of clean reads", tag="a") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.position = "none",
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=6))  # delete the background for facets (as in %shared)

figure_1B <- ggplot(meta_working[!(meta_working$Sample_name %in% c("bctrl4633v", "bctrl4654v", "bctrl4676v", "bctrl4699v")), ], aes(x=ncvssample, y=CHM_LU_richness_discovered_ratio)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=0.75, shape = 21, stroke = 0.1, color = "white") +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=6),
    legend.position = "bottom") +
  labs(y = "Proportion of at least medium quality sequences", fill="Timepoint", tag="b") # guides = collect (under the first two plots)

figure_1C <- ggplot(summarized_df2_frac_melt, aes(x = ncvssample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#a1c181", "#fcca46", "#fe7f2d", "#233d4d")) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  labs(tag="c") +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=6),
    legend.position = "bottom"
  ) +
  labs(x = "", y = "Mean fraction of the viral group richness \n (at least medium quality)", fill = "Viral group")

figure_1D <- ggplot(summarized_df3_frac_melt, aes(x = ncvssample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#619b8a", "#f4a261", "#233d4d")) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  labs(tag="d")+
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=6),
    legend.position = "bottom"
  ) +
  labs(x = "", y = "Mean fraction of the host group richness \n (at least medium quality vOTUs)", fill = "Host group")

figure_1E <- ggplot(summarized_df4_frac_melt, aes(x = ncvssample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("cornflowerblue", "gold2", "turquoise3", "orchid", "steelblue", "burlywood", "purple3", "slategray4","yellow", "blue", "skyblue","darkgoldenrod1",
                               "darksalmon", "green", "pink", "tomato2", "rosybrown3", "snow2","plum2", "aquamarine3",
                               "palegreen3", "seagreen", "aquamarine", "yellow3", "chartreuse", "thistle", "orange2", 
                               "navajowhite", "firebrick", "purple", "violetred1", "red", "darkblue", "black")) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  labs(tag="e") +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6, face = "italic"),
    legend.key.size = unit(0.3, "cm"),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=6),
    legend.position = "bottom"
  ) +
  labs(x = "", y = "Mean fraction of the host genus richness \n (at least medium quality vOTUs)", fill = "Host genus")

# Combine the plots using patchwork
combined_plot <- (figure_1A + figure_1B) / (figure_1C + figure_1D) / figure_1E

# Save the combined plot as a PDF
ggsave("combined_figure.png", combined_plot, width = 22/2.54, height = 30/2.54)

#######################################################################################################################################

## Stats calculation for part 1
#######################################################################################################################################

summary(lmer(clean_reads_comb ~ cohort + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated))
summary(lmer(contigs_1000 ~ cohort + clean_reads_comb + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated))
summary(lmer(total_viruses_discovered ~ cohort + clean_reads_comb + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated))
summary(lmer(total_viruses_discovered ~ cohort + clean_reads_comb + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$ncvssample == "NCs", ]))
summary(lmer(CHM_LU_richness_discovered_ratio ~ cohort + clean_reads_comb + (1|nc_subject_group), REML = F, data = meta_working))

summary(lmer(clean_reads_comb ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "liang", ]))
summary(lmer(clean_reads_comb ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "maqsood", ]))
summary(lmer(clean_reads_comb ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "shah", ]))
p.adjust(c(0.684, 0.0383, 0.725), method="BH")

summary(lmer(contigs_1000 ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "liang", ]))
summary(lmer(contigs_1000 ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "maqsood", ]))
summary(lmer(contigs_1000 ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "shah", ]))
p.adjust(c(0.0366, 0.802, 0.360), method="BH")

summary(lmer(total_viruses_discovered ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "liang", ]))
summary(lmer(total_viruses_discovered ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "maqsood", ]))
summary(lmer(total_viruses_discovered ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "shah", ]))
p.adjust(c(0.000301, 0.645, 0.368), method="BH")

summary(lmer(CHM_LU_richness_discovered_ratio ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "liang", ]))
summary(lmer(CHM_LU_richness_discovered_ratio ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "maqsood", ]))
summary(lmer(CHM_LU_richness_discovered_ratio ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "shah", ]))
p.adjust(c(0.0742, 0.020469, 0.8654), method="BH")

neg_ctrl_data <- meta_all_with_qc_curated %>%
  filter(Type == "Neg_ctrl") %>%
  pull(total_viruses_discovered)

mean_value <- mean(neg_ctrl_data, na.rm = TRUE)

standard_error <- sd(neg_ctrl_data, na.rm = TRUE) / sqrt(length(neg_ctrl_data))

alpha <- 0.05
t_value <- qt(1 - alpha/2, df = length(neg_ctrl_data) - 1)
confidence_interval <- t_value * standard_error
lower_bound <- mean_value - confidence_interval
upper_bound <- mean_value + confidence_interval

cat("Mean:", mean_value, "\n")
cat("95% Confidence Interval: [", lower_bound, ", ", upper_bound, "]\n")

summary(lmer(diversity ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "liang", ]))
summary(lmer(diversity ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "maqsood", ]))
summary(lmer(diversity ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "shah", ]))
p.adjust(c(0.00923, 0.0696, 0.00723), method="BH")

summary(lmer(richness ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "liang", ]))
summary(lmer(richness ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "maqsood", ]))
summary(lmer(richness ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "shah", ]))
p.adjust(c(0.00244, 0.540, 0.00723), method="BH")


#######################################################################################################################################

## Bray-Curtis boxplots: all samples; log scale
#######################################################################################################################################

bray_dist_matrix_full <- as.matrix(vegdist(t(RPKM), method="bray"))
bray_dist_matrix_full_rev <- 1 - bray_dist_matrix_full

#write.table(bray_dist_matrix_full_rev, "bray_dist_matrix_full_rev.tsv", sep='\t', row.names=T, col.names=T, quote=F)

upper_tri <- upper.tri(bray_dist_matrix_full_rev)

# Get row and column indices for upper triangle
row_indices <- row(bray_dist_matrix_full_rev)[upper_tri]
col_indices <- col(bray_dist_matrix_full_rev)[upper_tri]

# Extract distances corresponding to upper triangle indices
distances <- bray_dist_matrix_full_rev[upper_tri]

# Create a dataframe
df_distances <- data.frame(
  Sample1 = rownames(bray_dist_matrix_full_rev)[row_indices],
  Sample2 = colnames(bray_dist_matrix_full_rev)[col_indices],
  Distance = distances
)

# Ensure Sample1 < Sample2 in each row (order does not matter)
df_distances <- df_distances[order(pmin(df_distances$Sample1, df_distances$Sample2)),
]  # Sorting by Sample1, Sample2

meta_working <- as.data.frame(meta_working)
meta_working$Sample1 <- meta_working$Sample_name
meta_working$Sample2 <- meta_working$Sample_name

df_distances <- merge(df_distances, meta_working[c("Sample1", "ncvssample", "cohort")], by="Sample1", all.x=T)
df_distances <- merge(df_distances, meta_working[c("Sample2", "ncvssample", "cohort")], by="Sample2", all.x=T)
colnames(df_distances) <- c("Sample1", "Sample2", "Distance", "Sample1_ncvssample",  "Sample1_cohort", "Sample2_ncvssample",  "Sample2_cohort")

df_distances_nc <- df_distances[!(df_distances$Sample1_ncvssample == "SAMPLES" & df_distances$Sample2_ncvssample == "SAMPLES"), ]
df_distances_nc <- df_distances_nc %>%
  mutate(Sample1_cohort = as.character(Sample1_cohort),
         Sample2_cohort = as.character(Sample2_cohort)) %>%
  mutate(category = ifelse(Sample1_ncvssample == "NCs" & Sample2_ncvssample == "SAMPLES", "Samples", NA),
         category = ifelse(Sample2_ncvssample == "NCs" & Sample1_ncvssample == "SAMPLES", "Samples", category),
         category = ifelse(Sample2_ncvssample == "NCs" & Sample1_ncvssample == "NCs", "NCs", category),
         type_cohort = ifelse(Sample1_cohort == Sample2_cohort, "same", "different"),
         cohort_nc = ifelse(type_cohort == "same" & category == "NCs", Sample1_cohort, NA),
         cohort_nc = ifelse(category == "Samples" & Sample1_ncvssample == "NCs", Sample1_cohort, cohort_nc),
         cohort_nc = ifelse(category == "Samples" & Sample2_ncvssample == "NCs", Sample2_cohort, cohort_nc))

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

# Permutation analysis


# Initialize a list to store the results for each combination
results <- list()

# Define the unique values in cohort_nc
cohort_values <- c("garmaeva", "liang", "maqsood", "shah")

# Loop over each cohort_nc value
for (cohort in cohort_values) {
  
  # Filter the dataframe for the current cohort
  df_filtered <- df_distances_nc[df_distances_nc$cohort_nc == cohort, ]
  
  # Apply the two conditions for Sample1_ncvssample
  conditions <- list(
    "Condition_1" = df_filtered[df_filtered$Sample1_ncvssample == df_filtered$Sample2_ncvssample, ],
    "Condition_2" = df_filtered[df_filtered$Sample1_ncvssample != df_filtered$Sample2_ncvssample, ]
  )
  
  # Loop over each condition
  for (condition_name in names(conditions)) {
    
    # Skip the combination cohort == "garmaeva" and Condition_1
    if (cohort == "garmaeva" && condition_name == "Condition_1") {
      next
    }
    
    # Get the filtered data for the current condition
    df_condition <- conditions[[condition_name]]
    
    # Perform the initial Wilcoxon test to get the baseline p-value
    baseline_result <- wilcox.test(Distance ~ type_cohort, data = df_condition)
    baseline_pvalue <- baseline_result$p.value
    
    # Set up the permutation test
    set.seed(123)  # Set seed for reproducibility
    n_permutations <- 10000
    ppermute <- numeric(n_permutations)
    
    # Perform permutations
    for (i in 1:n_permutations) {
      # Shuffle the "Distance" column
      shuffled_distances <- sample(df_condition$Distance)
      
      # Create a new dataframe with the shuffled "Distance" column
      df_shuffled <- df_condition
      df_shuffled$Distance <- shuffled_distances
      
      # Perform Wilcoxon test on the shuffled data
      result <- wilcox.test(Distance ~ type_cohort, data = df_shuffled)
      
      # Store the p-value in the ppermute vector
      ppermute[i] <- result$p.value
    }
    
    # Calculate the final p-value
    final_pvalue <- sum(ppermute <= baseline_pvalue) / n_permutations
    
    # Save the results in the list
    results[[paste(cohort, condition_name, sep = "_")]] <- list(
      baseline_pvalue = baseline_pvalue,
      final_pvalue = final_pvalue
    )
  }
}

# Output the results
results

#CHECK FOR 0 VALUES FIRST!!!
min_nonzero <- min(df_distances_nc$Distance[df_distances_nc$Distance > 0])
df_distances_nc_for_plot <- df_distances_nc
df_distances_nc_for_plot$Distance <- df_distances_nc_for_plot$Distance + (min_nonzero / 2)

df_distances_nc_for_plot$type_cohort <- factor(df_distances_nc_for_plot$type_cohort, levels = c("same", "different"))
df_distances_nc_for_plot$cohort_nc <- factor(df_distances_nc_for_plot$cohort_nc, levels = c("garmaeva", "liang", "maqsood", "shah"))

#######################################################################################################################################

## Analysis for the vOTUs sharedness with the other cohorts: correlation between sharedness with own vs different cohort
#######################################################################################################################################
RPKM_filtered_w_dummy <- RPKM_count

RPKM_filtered_w_dummy$dummy_non_garmaeva_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_liang, negative_controls_maqsood, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_liang_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_maqsood, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_maqsood_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_liang, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_shah_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_liang, negative_controls_maqsood)])

RPKM_filtered_w_dummy$dummy_garmaeva_NC <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) == negative_controls_garmaeva]
RPKM_filtered_w_dummy$dummy_liang_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_liang])
RPKM_filtered_w_dummy$dummy_maqsood_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_maqsood])
RPKM_filtered_w_dummy$dummy_shah_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_shah])

RPKM_filtered_w_dummy[RPKM_filtered_w_dummy > 0] <- 1
RPKM_filtered_w_dummy <- RPKM_filtered_w_dummy[, !colnames(RPKM_filtered_w_dummy) %in% c(negative_controls)]

dummy_cols <- grep("^dummy", names(RPKM_filtered_w_dummy), value = TRUE)
non_dummy_cols <- setdiff(names(RPKM_filtered_w_dummy), dummy_cols)

# Initialize an empty dataframe for the results
result <- data.frame(matrix(ncol = length(dummy_cols), nrow = length(non_dummy_cols)))
names(result) <- dummy_cols
row.names(result) <- non_dummy_cols

# Calculate the percent similarity
for (non_dummy_col in non_dummy_cols) {
  for (dummy_col in dummy_cols) {
    # Number of viruses present in the non-dummy column
    present_in_non_dummy <- sum(RPKM_filtered_w_dummy[[non_dummy_col]] == 1)
    
    if (present_in_non_dummy > 0) {
      # Number of viruses present in both non-dummy and dummy column
      present_in_both <- sum(RPKM_filtered_w_dummy[[non_dummy_col]] == 1 & RPKM_filtered_w_dummy[[dummy_col]] == 1)
      
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
table_for_plot_cor <- merge(result, meta_working[c("Sample_name", "Type", "cohort", "Subject_ID", "Timepoint", "nc_subject_group")], by="Sample_name", all.x = T)
table_for_plot_cor <- table_for_plot_cor %>%
  mutate(same_cohort_NC = ifelse(cohort == "maqsood", dummy_maqsood_NC, NA),
         same_cohort_NC = ifelse(cohort == "liang", dummy_liang_NC, same_cohort_NC),
         same_cohort_NC = ifelse(cohort == "shah", dummy_shah_NC, same_cohort_NC),
         same_cohort_NC = ifelse(cohort == "garmaeva", dummy_garmaeva_NC, same_cohort_NC),
         different_cohort_NC = ifelse(cohort == "maqsood", dummy_non_maqsood_NC, NA),
         different_cohort_NC = ifelse(cohort == "liang", dummy_non_liang_NC, different_cohort_NC),
         different_cohort_NC = ifelse(cohort == "shah", dummy_non_shah_NC, different_cohort_NC),
         different_cohort_NC = ifelse(cohort == "garmaeva", dummy_non_garmaeva_NC, different_cohort_NC),
         Timepoint_numeric = as.integer(gsub("M", "", Timepoint)),
         Timepoint_numeric = ifelse(grepl("Y2-5", Timepoint), 24, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "Mtrim3", 0, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M0", 3, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M1", 4, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M2", 5, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M3", 6, Timepoint_numeric))

RPKM_filtered_w_dummy <- RPKM

RPKM_filtered_w_dummy$dummy_non_garmaeva_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_liang, negative_controls_maqsood, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_liang_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_maqsood, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_maqsood_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_liang, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_shah_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_liang, negative_controls_maqsood)])

RPKM_filtered_w_dummy$dummy_garmaeva_NC <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) == negative_controls_garmaeva]
RPKM_filtered_w_dummy$dummy_liang_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_liang])
RPKM_filtered_w_dummy$dummy_maqsood_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_maqsood])
RPKM_filtered_w_dummy$dummy_shah_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_shah])

RPKM_filtered_w_dummy <- RPKM_filtered_w_dummy[, !colnames(RPKM_filtered_w_dummy) %in% c(negative_controls)]

dummy_cols <- grep("^dummy", names(RPKM_filtered_w_dummy), value = TRUE)
non_dummy_cols <- setdiff(names(RPKM_filtered_w_dummy), dummy_cols)

# Initialize an empty dataframe for the results
result <- data.frame(matrix(ncol = length(dummy_cols), nrow = length(non_dummy_cols)))
names(result) <- dummy_cols
row.names(result) <- non_dummy_cols

# Calculate the percent similarity
for (non_dummy_col in non_dummy_cols) {
  for (dummy_col in dummy_cols) {
    # Number of viruses present in the non-dummy column
    present_in_non_dummy <- sum(RPKM_filtered_w_dummy[[non_dummy_col]])
    
    if (present_in_non_dummy > 0) {
      # Number of viruses present in both non-dummy and dummy column
      present_in_both <- sum(RPKM_filtered_w_dummy[[non_dummy_col]][RPKM_filtered_w_dummy[[non_dummy_col]] > 0 & RPKM_filtered_w_dummy[[dummy_col]] > 0])
      
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
table_for_plot_cor_a <- merge(result, meta_working[c("Sample_name", "Type", "cohort", "Subject_ID", "Timepoint", "nc_subject_group")], by="Sample_name", all.x = T)
table_for_plot_cor_a <- table_for_plot_cor_a %>%
  mutate(same_cohort_NC = ifelse(cohort == "maqsood", dummy_maqsood_NC, NA),
         same_cohort_NC = ifelse(cohort == "liang", dummy_liang_NC, same_cohort_NC),
         same_cohort_NC = ifelse(cohort == "shah", dummy_shah_NC, same_cohort_NC),
         same_cohort_NC = ifelse(cohort == "garmaeva", dummy_garmaeva_NC, same_cohort_NC),
         different_cohort_NC = ifelse(cohort == "maqsood", dummy_non_maqsood_NC, NA),
         different_cohort_NC = ifelse(cohort == "liang", dummy_non_liang_NC, different_cohort_NC),
         different_cohort_NC = ifelse(cohort == "shah", dummy_non_shah_NC, different_cohort_NC),
         different_cohort_NC = ifelse(cohort == "garmaeva", dummy_non_garmaeva_NC, different_cohort_NC),
         Timepoint_numeric = as.integer(gsub("M", "", Timepoint)),
         Timepoint_numeric = ifelse(grepl("Y2-5", Timepoint), 24, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "Mtrim3", 0, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M0", 3, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M1", 4, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M2", 5, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M3", 6, Timepoint_numeric))


table_for_plot_cor_combine <- merge(table_for_plot_cor, table_for_plot_cor_a[c("Sample_name", "same_cohort_NC", "different_cohort_NC")], by="Sample_name", 
                                    suffixes = c("_presence","_abundance"))

calculate_mean_ci <- function(data_vector, confidence_level = 0.95) {
  # Remove NA values
  data_vector <- na.omit(data_vector)
  
  # Calculate the mean
  mean_value <- mean(data_vector)
  
  # Calculate the standard error of the mean
  standard_error <- sd(data_vector) / sqrt(length(data_vector))
  
  # Calculate the t-value for the given confidence level
  alpha <- 1 - confidence_level
  t_value <- qt(1 - alpha/2, df = length(data_vector) - 1)
  
  # Calculate the confidence interval
  confidence_interval <- t_value * standard_error
  lower_bound <- mean_value - confidence_interval
  upper_bound <- mean_value + confidence_interval
  
  # Return a list with the mean and confidence interval
  return(list(mean = mean_value, lower_bound = lower_bound, upper_bound = upper_bound))
}

# Example usage
result <- calculate_mean_ci(table_for_plot_cor_combine$same_cohort_NC_abundance)
cat("Mean:", result$mean, "\n")
cat("95% Confidence Interval: [", result$lower_bound, ", ", result$upper_bound, "]\n")

result <- calculate_mean_ci(table_for_plot_cor_combine$same_cohort_NC_presence)
cat("Mean:", result$mean, "\n")
cat("95% Confidence Interval: [", result$lower_bound, ", ", result$upper_bound, "]\n")

result <- calculate_mean_ci(table_for_plot_cor_combine$different_cohort_NC_presence)
cat("Mean:", result$mean, "\n")
cat("95% Confidence Interval: [", result$lower_bound, ", ", result$upper_bound, "]\n")

summary(lmer(same_cohort_NC_abundance ~ cohort + (1|nc_subject_group), REML = F, data = table_for_plot_cor_combine))

model_perc_shared_predic <- lmer(different_cohort_NC_presence ~ same_cohort_NC_presence + Type + cohort + (1 | Subject_ID), data=table_for_plot_cor_combine)
summary(model_perc_shared_predic)
r_squared <- r.squaredGLMM(model_perc_shared_predic)
print(r_squared)

table_for_plot_cor_combine_melt <- melt(table_for_plot_cor_combine[c("Sample_name", "cohort", "same_cohort_NC_presence", "different_cohort_NC_presence", 
                                                                     "same_cohort_NC_abundance", "different_cohort_NC_abundance")], 
                                        id.vars = c("Sample_name", "cohort"))

table_for_plot_cor_combine_melt <- table_for_plot_cor_combine_melt %>%
  mutate(pres_abun = ifelse(grepl("presence", variable), "presence", "abundance"),
         same_diff = ifelse(grepl("same", variable), "same", "different"))

table_for_plot_cor_combine_melt$same_diff <- factor(table_for_plot_cor_combine_melt$same_diff, levels = c("same", "different"))
table_for_plot_cor_combine_melt$pres_abun <- factor(table_for_plot_cor_combine_melt$pres_abun, levels = c("presence", "abundance"))

#######################################################################################################################################

## Analysis for the vOTUs sharedness along the timepoints: with dummy, only Liang and Garmaeva cohorts, no cross-cohort comparison
#######################################################################################################################################
filter_RPKM_table <- function(RPKM_data, samples) {
  filtered_RPKM <- RPKM_data[, colnames(RPKM_data) %in% samples]
  filtered_RPKM <- filtered_RPKM[rowSums(filtered_RPKM) > 0, colSums(filtered_RPKM) > 0]
  return(filtered_RPKM)
}


calculate_similarity <- function(RPKM_data, nc_samples, non_nc_samples) {
  result <- data.frame(matrix(ncol = length(nc_samples), nrow = length(non_nc_samples)))
  names(result) <- nc_samples
  row.names(result) <- non_nc_samples
  
  for (non_dummy_col in non_nc_samples) {
    for (dummy_col in nc_samples) {
      present_in_non_dummy <- sum(RPKM_data[[non_dummy_col]] == 1)
      
      if (present_in_non_dummy > 0) {
        present_in_both <- sum(RPKM_data[[non_dummy_col]] == 1 & RPKM_data[[dummy_col]] == 1)
        similarity <- present_in_both / present_in_non_dummy * 100
      } else {
        similarity <- NA
      }
      
      result[non_dummy_col, dummy_col] <- similarity
    }
  }
  
  result$Sample_name <- row.names(result)
  row.names(result) <- NULL
  return(result)
}

calculate_similarity_a <- function(RPKM_data, nc_samples, non_nc_samples) {
  result <- data.frame(matrix(ncol = length(nc_samples), nrow = length(non_nc_samples)))
  names(result) <- nc_samples
  row.names(result) <- non_nc_samples
  
  for (non_dummy_col in non_nc_samples) {
    for (dummy_col in nc_samples) {
      abundance_in_non_dummy <- sum(RPKM_data[[non_dummy_col]])
      
      if (abundance_in_non_dummy > 0) {
        abundace_present_in_both <- sum(RPKM_data[[non_dummy_col]][RPKM_data[[non_dummy_col]] > 0 & RPKM_data[[dummy_col]] > 0])  
        similarity <- abundace_present_in_both / abundance_in_non_dummy * 100
      } else {
        similarity <- NA
      }
      
      result[non_dummy_col, dummy_col] <- similarity
    }
  }
  
  result$Sample_name <- row.names(result)
  row.names(result) <- NULL
  return(result)
}

negative_controls <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$Type == "Neg_ctrl" & !is.na(meta_all_with_qc_curated$Type)]
samples_g <- c(meta_working$Sample_name[meta_working$Type %in% c("Infant", "Neg_ctrl") & meta_working$cohort == "garmaeva"])

RPKM_filtered_g <- filter_RPKM_table(RPKM_count, samples_g)

nc_samples_g <- colnames(RPKM_filtered_g)[colnames(RPKM_filtered_g) %in% negative_controls]
non_nc_samples_g <- setdiff(names(RPKM_filtered_g), nc_samples_g)

result_g <- calculate_similarity(RPKM_filtered_g, nc_samples_g, non_nc_samples_g)
table_for_plot_g <- merge(result_g, meta_working[c("Sample_name", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb")], by="Sample_name", all.x = TRUE)
table_for_plot_g$Timepoint <- factor(table_for_plot_g$Timepoint, levels = c("M1", "M2", "M3", "M6", "M12"))
table_for_plot_g$Dataset <- "RPKM_count"
colnames(table_for_plot_g) <- c("Sample_name", "value", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb", "Dataset")
table_for_plot_g$Subject_ID <- as.factor(table_for_plot_g$Subject_ID)
table_for_plot_g$Timepoint_numeric <- as.integer(gsub("M", "", table_for_plot_g$Timepoint))

summary(lmer(value ~ Timepoint_numeric + (1|Subject_ID), REML = F, data = table_for_plot_g))

RPKM_filtered_ga <- filter_RPKM_table(RPKM, samples_g)

result_ga <- calculate_similarity_a(RPKM_filtered_ga, nc_samples_g, non_nc_samples_g)
table_for_plot_ga <- merge(result_ga, meta_working[c("Sample_name", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb")], by="Sample_name", all.x = TRUE)
table_for_plot_ga$Timepoint <- factor(table_for_plot_ga$Timepoint, levels = c("M1", "M2", "M3", "M6", "M12"))
table_for_plot_ga$Dataset <- "RPKM"
colnames(table_for_plot_ga) <- c("Sample_name", "value", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb", "Dataset")
table_for_plot_ga$Subject_ID <- as.factor(table_for_plot_ga$Subject_ID)
table_for_plot_ga$Timepoint_numeric <- as.integer(gsub("M", "", table_for_plot_ga$Timepoint))

summary(lmer(value ~ Timepoint_numeric + (1|Subject_ID), REML = F, data = table_for_plot_ga))

samples_l <- c(meta_working$Sample_name[(meta_working$Subcohort == "Discovery" & !is.na(meta_working$Subcohort) 
                                         & meta_working$cohort == "liang" & meta_working$Type == "Infant") | 
                                          (meta_working$cohort == "liang" & meta_working$Type == "Neg_ctrl")])

RPKM_filtered_l <- filter_RPKM_table(RPKM_count, samples_l)

RPKM_filtered_l$liang_nc_dummy <- rowSums(RPKM_filtered_l[, colnames(RPKM_filtered_l) %in% negative_controls], na.rm = TRUE)
RPKM_filtered_l$liang_nc_dummy[RPKM_filtered_l$liang_nc_dummy > 0] <- 1
RPKM_filtered_l <- RPKM_filtered_l[, !colnames(RPKM_filtered_l) %in% negative_controls]

nc_samples_l <- "liang_nc_dummy"
non_nc_samples_l <- setdiff(names(RPKM_filtered_l), nc_samples_l)

result_l <- calculate_similarity(RPKM_filtered_l, nc_samples_l, non_nc_samples_l)

table_for_plot_l <- merge(result_l, meta_working[c("Sample_name", "Sample_ID", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb")], by="Sample_name", all.x = T)
table_for_plot_l <- table_for_plot_l %>%
  mutate(nucleic_acid = ifelse(grepl("_D", Sample_ID), "DNA", "RNA"))
table_for_plot_l <- table_for_plot_l[table_for_plot_l$nucleic_acid == "DNA", ]
table_for_plot_l$nucleic_acid <- NULL
table_for_plot_l$Sample_ID <- NULL
table_for_plot_l$Dataset <- "RPKM_count"
table_for_plot_l$Timepoint <- factor(table_for_plot_l$Timepoint, levels = c("M0", "M1", "M4"))
colnames(table_for_plot_l) <- c("Sample_name", "value", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb", "Dataset")
table_for_plot_l$Subject_ID <- as.factor(table_for_plot_l$Subject_ID)
table_for_plot_l$Timepoint_numeric <- as.integer(gsub("M", "", table_for_plot_l$Timepoint))

summary(lmer(value ~ Timepoint_numeric + (1|Subject_ID), REML = F, data = table_for_plot_l))

RPKM_filtered_la <- filter_RPKM_table(RPKM, samples_l)

RPKM_filtered_la$liang_nc_dummy <- rowSums(RPKM_filtered_la[, colnames(RPKM_filtered_la) %in% negative_controls], na.rm = TRUE)
RPKM_filtered_la$liang_nc_dummy[RPKM_filtered_la$liang_nc_dummy > 0] <- 1
RPKM_filtered_la <- RPKM_filtered_la[, !colnames(RPKM_filtered_la) %in% negative_controls]

result_la <- calculate_similarity_a(RPKM_filtered_la, nc_samples_l, non_nc_samples_l)
table_for_plot_la <- merge(result_la, meta_working[c("Sample_name", "Sample_ID", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb")], by="Sample_name", all.x = T)
table_for_plot_la <- table_for_plot_la %>%
  mutate(nucleic_acid = ifelse(grepl("_D", Sample_ID), "DNA", "RNA"))
table_for_plot_la <- table_for_plot_la[table_for_plot_la$nucleic_acid == "DNA", ]
table_for_plot_la$nucleic_acid <- NULL
table_for_plot_la$Sample_ID <- NULL
table_for_plot_la$Dataset <- "RPKM"
table_for_plot_la$Timepoint <- factor(table_for_plot_la$Timepoint, levels = c("M0", "M1", "M4"))
colnames(table_for_plot_la) <- c("Sample_name", "value", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb", "Dataset")
table_for_plot_la$Subject_ID <- as.factor(table_for_plot_la$Subject_ID)
table_for_plot_la$Timepoint_numeric <- as.integer(gsub("M", "", table_for_plot_la$Timepoint))

summary(lmer(value ~ Timepoint_numeric + (1|Subject_ID), REML = F, data = table_for_plot_la))

combined_table_for_plot <- rbind(table_for_plot_g, table_for_plot_ga, table_for_plot_l, table_for_plot_la)
combined_table_for_plot$Timepoint <- factor(combined_table_for_plot$Timepoint, levels = c("M0", "M1", "M2", "M3", "M4", "M6", "M12"))
combined_table_for_plot$Dataset <- factor(combined_table_for_plot$Dataset, levels = c("RPKM_count", "RPKM"))

#CHECK FOR 0 VALUES FIRST!!!
min_nonzero <- min(combined_table_for_plot$value[combined_table_for_plot$value > 0])
combined_table_for_plot_log <- combined_table_for_plot
combined_table_for_plot_log$value <- combined_table_for_plot_log$value + (min_nonzero / 2)

#######################################################################################################################################

## Patching the Figure 2
#######################################################################################################################################
figure_2A <- ggplot(df_distances_nc_for_plot, aes(x = type_cohort, y = Distance)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outliers = FALSE) +
  scale_y_log10() +  # Apply log scale to the y-axis
  #  coord_cartesian(ylim = c(new_min_y, 1)) +  # Set limits in the original scale
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
  labs(y = "log10(1 - (Bray-Curtis dissimilarity))", x = "Type of cohort \n (samples from the same or different cohort compared)", tag="a") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    strip.background = element_rect(fill = "transparent")
  )


figure_2B <- ggplot(table_for_plot_cor_combine_melt, aes(x = same_diff, y = value)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outliers = FALSE) +
  ylim(0, 100) +
  facet_grid(cohort ~ pres_abun, scales = "free", labeller = labeller(
    pres_abun = c(
      "presence" = "% shared cotigs (presence)",
      "abundance" = "% shared cotigs (abundance)"
    ),
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>",
      "shah" = "Samples Shah *et al.* <br>"
    )
  )) +
  labs(tag="b") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    plot.title = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8), 
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "transparent")
  )

figure_2С <- ggplot(combined_table_for_plot_log, aes(x = Timepoint, y = value)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outliers = FALSE) +
  scale_y_log10() +
  facet_grid(Dataset ~ cohort, scales = "free", labeller = labeller(
    Dataset = c(
      "RPKM" = "Abundance",
      "RPKM_count" = "Presence"
    ),
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>"
    )
  )) +
  labs(y = "log10(% shared vOTUs)", tag="c") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "transparent")
  )

figure_2D <- ggplot(table_for_plot_cor_combine, aes(x=different_cohort_NC_presence, y=same_cohort_NC_presence)) +
  geom_point(size = 1, aes(color=Type), alpha=0.5) +  # Adjusted point size and added transparency
  geom_smooth(method=lm, color="black", fill="lightgray", se=TRUE, linewidth=0.5) +  # Use linewidth instead of size
  scale_color_manual(values=c("#5B84B1FF", "#FC766AFF")) +
  labs(x = "% shared vOTUs with NCs \n from different studies", y = "% shared vOTUs with NCs \n from same study") +
  labs(tag="d") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom"
  )


# Combine the plots using patchwork
combined_plot2 <- (figure_2A + figure_2B) / (figure_2С + figure_2D)

# Save the combined plot as a PDF
ggsave("combined_figure2.png", combined_plot2, width = 28/2.54, height = 24/2.54)

#######################################################################################################################################

## Checking viral clusters -> supplementary
#######################################################################################################################################
negative_controls_all <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$Type == "Neg_ctrl"]
viral_clusters <- as.data.frame(read_tsv('/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/initial_dereplication/NCP_viral_clusters.tsv',
                                         col_names = FALSE))

colnames(viral_clusters) <- c("Representative", "vOTUs_included")



RPKM_ncs <- RPKM[, colnames(RPKM) %in% negative_controls_all]
RPKM_ncs <- RPKM_ncs[rowSums(RPKM_ncs) > 0, colSums(RPKM_ncs) > 0]
dim(RPKM_ncs)  # 7376   43

vc_detected_in_nc <- row.names(RPKM_ncs)

contains_negative_control <- function(text, negative_controls_all) {
  any(sapply(negative_controls_all, grepl, text))
}

filtered_viral_clusters <- viral_clusters %>%
  filter(sapply(vOTUs_included, contains_negative_control, negative_controls_all))

vc_discovered_in_nc <- filtered_viral_clusters$Representative

venn.plot <- venn.diagram(
  x = list(Vector1 = vc_discovered_in_nc, Vector2 = vc_detected_in_nc),
  category.names = c("VC containing \n discovered vOTUs \n in NCs", "VC containing \n detected vOTUs \n in NCs"),
  filename = NULL,
  output = TRUE,
  height = 2000,    # Adjusting height to make the picture smaller
  width = 2000,     # Adjusting width to make the picture smaller
  resolution = 300, # Higher resolution for better quality
  lwd = 2,
  col = c("black", "black"),
  fill = c("lightblue", "lightcoral"),
  alpha = 0.5,
  cex = 0.8,        # Decrease font size for printed numbers
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.8,    # Decrease font size for category names
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.col = "black", # Set category label color to black
  cat.pos = c(-20, 20),  # Adjusted position of category labels
  cat.dist = c(0.05, 0.05),  # Adjusted distance of category labels from the circles
  cat.just = list(c(0.5, 0), c(0.5, 0)), # Adjusted justification of category labels
  ext.pos = 30,
  ext.dist = -0.1,
  ext.length = 0.85,
  ext.line.lwd = 2,
  ext.line.lty = "dashed"
)

# Save the Venn diagram to a PDF file
pdf(file = "venn_vc_in_nc_detected_vs_discovered.pdf")
grid.draw(venn.plot)
dev.off()

#######################################################################################################################################


