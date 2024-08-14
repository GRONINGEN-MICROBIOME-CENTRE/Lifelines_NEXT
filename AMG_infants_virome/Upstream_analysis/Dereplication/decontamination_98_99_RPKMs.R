## Code description
#######################################################################################################################################
## Script for the negative control sharing paper: decontaminating strategies
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
library(ggrepel)
library(ggtext)
library(VennDiagram)
library(lmerTest)
#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
setwd("/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R")
#######################################################################################################################################

## Load the  data
#######################################################################################################################################
extended_tof <- read.delim('../VIR_DB/virus_contigs/MERGED_Extended_TOF_NCP')
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

RPKM_initial_98_cov75 <- read.delim("RPKM_counts_VLP_wo_nc98_der95_bc75_NCP.txt")
RPKM_initial_98_cov75 <- RPKM_initial_98_cov75[row.names(RPKM_initial_98_cov75) %in% RPKM_contigs_keep, ]
RPKM_initial_98_cov75 <- RPKM_initial_98_cov75[rowSums(RPKM_initial_98_cov75) > 0, colSums(RPKM_initial_98_cov75) > 0]
dim(RPKM_initial_98_cov75)  # 192922   1289

RPKM_initial_98_cov98 <- read.delim("RPKM_counts_VLP_wo_nc98_der95_bc98_NCP.txt")
RPKM_initial_98_cov98 <- RPKM_initial_98_cov98[row.names(RPKM_initial_98_cov98) %in% RPKM_contigs_keep, ]
RPKM_initial_98_cov98 <- RPKM_initial_98_cov98[rowSums(RPKM_initial_98_cov98) > 0, colSums(RPKM_initial_98_cov98) > 0]
dim(RPKM_initial_98_cov98)  # 150827   1244

RPKM_initial_99_cov75 <- read.delim("RPKM_counts_VLP_wo_nc99_der95_bc75_NCP.txt")
RPKM_initial_99_cov75 <- RPKM_initial_99_cov75[row.names(RPKM_initial_99_cov75) %in% RPKM_contigs_keep, ]
RPKM_initial_99_cov75 <- RPKM_initial_99_cov75[rowSums(RPKM_initial_99_cov75) > 0, colSums(RPKM_initial_99_cov75) > 0]
dim(RPKM_initial_99_cov75)  # 193118   1290

RPKM_initial_99_cov99 <- read.delim("RPKM_counts_VLP_wo_nc99_der95_bc99_NCP.txt")
RPKM_initial_99_cov99 <- RPKM_initial_99_cov99[row.names(RPKM_initial_99_cov99) %in% RPKM_contigs_keep, ]
RPKM_initial_99_cov99 <- RPKM_initial_99_cov99[rowSums(RPKM_initial_99_cov99) > 0, colSums(RPKM_initial_99_cov99) > 0]
dim(RPKM_initial_99_cov99)  # 137764   1235

meta_all_with_qc_curated <- as.data.frame(read_tsv('../metadata_with_qc_NCPv2.tsv'))
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
#######################################################################################################################################

## Do final cleaning for the RPKM tables
#######################################################################################################################################
## 98%
negative_controls <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$Type == "Neg_ctrl"]

negative_controls_garmaeva <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$Type == "Neg_ctrl" & meta_all_with_qc_curated$cohort == "garmaeva"]
negative_controls_liang <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$Type == "Neg_ctrl" & meta_all_with_qc_curated$cohort == "liang"]
negative_controls_maqsood <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$Type == "Neg_ctrl" & meta_all_with_qc_curated$cohort == "maqsood"]
negative_controls_shah <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$Type == "Neg_ctrl" & meta_all_with_qc_curated$cohort == "shah"]

samples_garmaeva <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$cohort == "garmaeva"]
samples_liang <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$cohort == "liang"]
samples_maqsood <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$cohort == "maqsood"]
samples_shah <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$cohort == "shah"]

RPKM_NCs_and_samplesshared_98_cov98 <- RPKM_initial_98_cov98

RPKM_NCs_and_samplesshared_98_cov98$dummy_garmaeva_NC <- RPKM_NCs_and_samplesshared_98_cov98[, colnames(RPKM_NCs_and_samplesshared_98_cov98) == negative_controls_garmaeva]
RPKM_NCs_and_samplesshared_98_cov98$dummy_liang_NC <- rowSums(RPKM_NCs_and_samplesshared_98_cov98[, colnames(RPKM_NCs_and_samplesshared_98_cov98) %in% negative_controls_liang])
RPKM_NCs_and_samplesshared_98_cov98$dummy_maqsood_NC <- rowSums(RPKM_NCs_and_samplesshared_98_cov98[, colnames(RPKM_NCs_and_samplesshared_98_cov98) %in% negative_controls_maqsood])
RPKM_NCs_and_samplesshared_98_cov98$dummy_shah_NC <- rowSums(RPKM_NCs_and_samplesshared_98_cov98[, colnames(RPKM_NCs_and_samplesshared_98_cov98) %in% negative_controls_shah])

# Function to perform the analysis with a check for sample existence
analyze_pairs <- function(samples_vector, dummy_column, df) {
  results <- list()
  
  # Check if the dummy column exists
  if (!(dummy_column %in% colnames(df))) {
    stop(paste("Dummy column", dummy_column, "not found in the dataframe"))
  }
  
  for (sample in samples_vector) {
    # Check if the sample column exists in the dataframe
    if (sample %in% colnames(df)) {
      selected_rows <- rownames(df)[
        df[[sample]] > 0 & 
          df[[dummy_column]] > 0
      ]
      results[[sample]] <- selected_rows
    } else {
      message(paste("Sample", sample, "not found in the dataframe. Skipping..."))
    }
  }
  
  return(results)
}

ANI98_cov98_results <- list()

# Perform analysis for each pair
ANI98_cov98_results$garmaeva <- analyze_pairs(samples_garmaeva, "dummy_garmaeva_NC", RPKM_NCs_and_samplesshared_98_cov98)
ANI98_cov98_results$liang <- analyze_pairs(samples_liang, "dummy_liang_NC", RPKM_NCs_and_samplesshared_98_cov98)
ANI98_cov98_results$maqsood <- analyze_pairs(samples_maqsood, "dummy_maqsood_NC", RPKM_NCs_and_samplesshared_98_cov98)
ANI98_cov98_results$shah <- analyze_pairs(samples_shah, "dummy_shah_NC", RPKM_NCs_and_samplesshared_98_cov98)

RPKM_fin_98_cov75 <- RPKM_initial_98_cov75

for (pair in names(ANI98_cov98_results)) {
  for (sample in names(ANI98_cov98_results[[pair]])) {
    rows_to_zero <- ANI98_cov98_results[[pair]][[sample]]
    
    # Check if the sample (column) and row exist in RPKM_fin_98_cov75
    if (sample %in% colnames(RPKM_fin_98_cov75)) {
      for (row in rows_to_zero) {
        if (row %in% rownames(RPKM_fin_98_cov75)) {
          RPKM_fin_98_cov75[row, sample] <- 0
        }
      }
    }
  }
}

RPKM_fin_98_cov75 <- RPKM_fin_98_cov75[rowSums(RPKM_fin_98_cov75) > 0, colSums(RPKM_fin_98_cov75) > 0]
dim(RPKM_fin_98_cov75)
RPKM_fin_98_cov75 <- RPKM_fin_98_cov75[, !colnames(RPKM_fin_98_cov75) %in% negative_controls]
RPKM_fin_98_cov75 <- RPKM_fin_98_cov75[rowSums(RPKM_fin_98_cov75) > 0, colSums(RPKM_fin_98_cov75) > 0]
dim(RPKM_fin_98_cov75)

## 99%
RPKM_NCs_and_samplesshared_99_cov99 <- RPKM_initial_99_cov99

RPKM_NCs_and_samplesshared_99_cov99$dummy_garmaeva_NC <- RPKM_NCs_and_samplesshared_99_cov99[, colnames(RPKM_NCs_and_samplesshared_99_cov99) == negative_controls_garmaeva]
RPKM_NCs_and_samplesshared_99_cov99$dummy_liang_NC <- rowSums(RPKM_NCs_and_samplesshared_99_cov99[, colnames(RPKM_NCs_and_samplesshared_99_cov99) %in% negative_controls_liang])
RPKM_NCs_and_samplesshared_99_cov99$dummy_maqsood_NC <- rowSums(RPKM_NCs_and_samplesshared_99_cov99[, colnames(RPKM_NCs_and_samplesshared_99_cov99) %in% negative_controls_maqsood])
RPKM_NCs_and_samplesshared_99_cov99$dummy_shah_NC <- rowSums(RPKM_NCs_and_samplesshared_99_cov99[, colnames(RPKM_NCs_and_samplesshared_99_cov99) %in% negative_controls_shah])

ANI99_cov99_results <- list()

# Perform analysis for each pair
ANI99_cov99_results$garmaeva <- analyze_pairs(samples_garmaeva, "dummy_garmaeva_NC", RPKM_NCs_and_samplesshared_99_cov99)
ANI99_cov99_results$liang <- analyze_pairs(samples_liang, "dummy_liang_NC", RPKM_NCs_and_samplesshared_99_cov99)
ANI99_cov99_results$maqsood <- analyze_pairs(samples_maqsood, "dummy_maqsood_NC", RPKM_NCs_and_samplesshared_99_cov99)
ANI99_cov99_results$shah <- analyze_pairs(samples_shah, "dummy_shah_NC", RPKM_NCs_and_samplesshared_99_cov99)

RPKM_fin_99_cov75 <- RPKM_initial_99_cov75

for (pair in names(ANI99_cov99_results)) {
  for (sample in names(ANI99_cov99_results[[pair]])) {
    rows_to_zero <- ANI99_cov99_results[[pair]][[sample]]
    
    # Check if the sample (column) and row exist in RPKM_fin_99_cov75
    if (sample %in% colnames(RPKM_fin_99_cov75)) {
      for (row in rows_to_zero) {
        if (row %in% rownames(RPKM_fin_99_cov75)) {
          RPKM_fin_99_cov75[row, sample] <- 0
        }
      }
    }
  }
}

RPKM_fin_99_cov75 <- RPKM_fin_99_cov75[rowSums(RPKM_fin_99_cov75) > 0, colSums(RPKM_fin_99_cov75) > 0]
dim(RPKM_fin_99_cov75)
RPKM_fin_99_cov75 <- RPKM_fin_99_cov75[, !colnames(RPKM_fin_99_cov75) %in% negative_controls]
RPKM_fin_99_cov75 <- RPKM_fin_99_cov75[rowSums(RPKM_fin_99_cov75) > 0, colSums(RPKM_fin_99_cov75) > 0]
dim(RPKM_fin_99_cov75)
#######################################################################################################################################

## Saving the tables 
#######################################################################################################################################
write.table(RPKM_fin_98_cov75, "RPKM_counts_VLP_NCP_fin_98_cov75.txt", sep='\t')
write.table(RPKM_fin_99_cov75, "RPKM_counts_VLP_NCP_fin_99_cov75.txt", sep='\t')
#######################################################################################################################################