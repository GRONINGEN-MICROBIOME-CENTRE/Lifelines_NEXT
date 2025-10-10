### LIFELINES NEXT ASSOCIATION PHENOTYPES & TAXONOMIC SPECIES INFANTS ####
### AUTHOR:  TRISHLA SINHA, ALEXANDER KURULSHIKOV
### ORIGINAL SCRIPT: 3rd AUGUST, 2023
### LAST UPDATE: 23rd September, 2025 

# Load packages 
library(tidyverse)
library(wesanderson)
library(reshape2)
library(foreach)
library(data.table)
library(ggplot2)
library(dplyr)
library(lmerTest)
library(pheatmap)
library(mmrm)
library (mgcv) 

# Load all functions from script : 
source("Functions_associations_phenotypes_LLNEXT.R")

# Load metadata 
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"
metadata_infants$Timepoint_categorical=factor(metadata_infants$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))

# Loading taxanomic data 
metaphlan_clr = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")


#################### Processing cross-sectional phenotypes ###############################

# Loading cross-sectional phenotypes 
cross_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_2023_11_15.txt")
cross_selection <-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_selection_AFTER_correlation_&_correction_10_07_2024.txt")


# Selecting data relevant for infant 
infant_cross_selection<-cross_selection[cross_selection$infant_microbiome_selection==1,]
column_names <- infant_cross_selection[[1]]
infant_cross_phenotypes <- cross_phenotypes %>% select(all_of(column_names))
infant_metadata_cross_phenotypes<-left_join(metadata_infants, infant_cross_phenotypes)
row.names(infant_metadata_cross_phenotypes)<-infant_metadata_cross_phenotypes$NG_ID

# For this specific data remove duplicates for categorical timepoints 
infant_metadata_cross_phenotypes_1 <- infant_metadata_cross_phenotypes %>%
  distinct(SAMPLE_ID, .keep_all = TRUE) 

# Remove phenotypes with too many NA's 
infant_metadata_cross_phenotypes_2 <-drop_phenotypes(infant_metadata_cross_phenotypes_1, 200, 5)
names (infant_metadata_cross_phenotypes_2)[2]<-"NEXT_ID"
infant_metadata_cross_phenotypes_2$NEXT_ID=as.factor(infant_metadata_cross_phenotypes_2$NEXT_ID)


# Matching ID's 
infant_metadata_cross_phenotypes_2<- infant_metadata_cross_phenotypes_2[rownames(infant_metadata_cross_phenotypes_2) %in% rownames(metaphlan_clr), ]
metaphlan_clr <- metaphlan_clr [rownames(metaphlan_clr ) %in% rownames(infant_metadata_cross_phenotypes_2), ]

# Select and normalize phenotypes 
metadata_columns <- c("NG_ID", "NEXT_ID","Modified_NEXT_ID_without_preg_number", "days_from_first_collection","FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                      "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                      "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                      "exact_age_days_at_collection", "exact_age_months_at_collection",
                      "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                      "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                      "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                      "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint")


phenos <- infant_metadata_cross_phenotypes_2[, !(colnames(infant_metadata_cross_phenotypes_2) %in% metadata_columns)] 
phenos<-run_invrank_dataFrame(phenos)
infant_metadata_cross_phenotypes_2$BATCH_NUMBER=as.factor(infant_metadata_cross_phenotypes_2$BATCH_NUMBER)

# Select and normalize covariates (Base model with only correction for technical variables)
covariates_base <- infant_metadata_cross_phenotypes_2[,c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER")]
covariates_normalized_base <- run_invrank_dataFrame(covariates_base)
results_base <- run_mmrm_analysis(metaphlan_clr, phenos, 
                                  time = infant_metadata_cross_phenotypes_2$exact_age_months_at_collection, 
                                  time_cat = infant_metadata_cross_phenotypes_2$Timepoint_categorical, 
                                  NEXT_ID = infant_metadata_cross_phenotypes_2$NEXT_ID, 
                                  covariates = covariates_normalized_base)
results_cross_trait_base <- as.data.frame(results_base$trait)

# Correction for mode of delivery
covariates_delivery <- infant_metadata_cross_phenotypes_2[,c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary")]
covariates_normalized_delivery <- run_invrank_dataFrame(covariates_delivery)
results_delivery <- run_mmrm_analysis(metaphlan_clr, phenos, 
                                      time = infant_metadata_cross_phenotypes_2$exact_age_months_at_collection, 
                                      time_cat = infant_metadata_cross_phenotypes_2$Timepoint_categorical, 
                                      NEXT_ID = infant_metadata_cross_phenotypes_2$NEXT_ID, 
                                      covariates = covariates_normalized_delivery)
results_cross_trait_delivery <- as.data.frame(results_delivery$trait)


# Correction for mode of delivery and ever never breastfed
covariates_delivery_breastfed <- infant_metadata_cross_phenotypes_2[,c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")]
covariates_normalized_delivery_breastfed <- run_invrank_dataFrame(covariates_delivery_breastfed)
results_delivery_breastfed <- run_mmrm_analysis(metaphlan_clr, phenos, 
                                                time = infant_metadata_cross_phenotypes_2$exact_age_months_at_collection, 
                                                time_cat = infant_metadata_cross_phenotypes_2$Timepoint_categorical, 
                                                NEXT_ID = infant_metadata_cross_phenotypes_2$NEXT_ID, 
                                                covariates = covariates_normalized_delivery_breastfed)
results_cross_trait_delivery_breastfed <- as.data.frame(results_delivery_breastfed$trait)


rename_except_cross <- function(df, suffix) {
  cols_excluded <- c("outcome", "trait","levels")
  colnames(df) <- ifelse(colnames(df) %in% cols_excluded, colnames(df), paste0(colnames(df), suffix))
  return(df)
}

results_cross_trait_base <- rename_except_cross (results_cross_trait_base, "_cor_base")
results_cross_trait_delivery  <- rename_except_cross (results_cross_trait_delivery, "_cor_delivery")
results_cross_trait_delivery_breastfed <- rename_except_cross (results_cross_trait_delivery_breastfed , "_cor_delivery_breastfed")


combined_results_cross <- reduce(
  list(results_cross_trait_base, results_cross_trait_delivery, results_cross_trait_delivery_breastfed),
  full_join,
  by = c("outcome", "trait", "levels")
)
combined_results_cross$trait_group <- sapply(strsplit(as.character(combined_results_cross$trait), "_"), function(x) paste(x[1:2], collapse = "_"))
combined_results_cross$type_phenotype<-"cross_sectional_phenotype"

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/taxa/")

combined_results_cross$FDR_base<-p.adjust(combined_results_cross$p_cor_base, method = "fdr")
combined_results_cross$FDR_delivery<-p.adjust(combined_results_cross$p_cor_delivery, method = "fdr")
combined_results_cross$FDR_delivery_feeding<-p.adjust(combined_results_cross$p_cor_delivery_breastfed, method = "fdr")

# Supplementary table S27
write.table(combined_results_cross, "taxa_SGB_cross_phenotypes_results_23_09.txt", sep = "\t", row.names = F)

#################### Processing dynamic phenotypes ###############################
# Load dynamic phenotypes 
dynamic_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_longitudinal_2023_09_29.txt")
dynamic_selection<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_longitudinal_selection_AFTER_correlation_&_correction_10_07_2024.txt")

# Merging files and selecting infant relevant data
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"
metadata_infants$BATCH_NUMBER<-as.factor(metadata_infants$BATCH_NUMBER)
metadata_infants$Timepoint_categorical=factor(metadata_infants$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))


# Selecting data relevant for infant 
infant_dynamic_selection<-dynamic_selection[dynamic_selection$infant_microbiome_selection==1,]
column_names <- infant_dynamic_selection[[1]]
infant_dynamic_phenotypes <- dynamic_phenotypes %>% select(all_of(column_names))
infant_metadata_dynamic_phenotypes<-left_join(metadata_infants, infant_dynamic_phenotypes)
row.names(infant_metadata_dynamic_phenotypes)<-infant_metadata_dynamic_phenotypes$NG_ID

# For this specific data remove duplicates for categorical timepoints 
infant_metadata_dynamic_phenotypes_1 <- infant_metadata_dynamic_phenotypes %>%
  distinct(SAMPLE_ID, .keep_all = TRUE) 

# Remove phenotypes with too many NA's 
infant_metadata_dynamic_phenotypes_2 <-drop_phenotypes(infant_metadata_dynamic_phenotypes_1, 200, 5)

delivery_feeding_mode<-infant_metadata_cross_phenotypes_1[,c( "SAMPLE_ID", "birth_deliverybirthcard_mode_binary","infant_ffq_ever_never_breastfed") ]
infant_metadata_dynamic_phenotypes_2<-left_join(infant_metadata_dynamic_phenotypes_2, delivery_feeding_mode)

names (infant_metadata_dynamic_phenotypes_2)[2]<-"NEXT_ID"
infant_metadata_dynamic_phenotypes_2$NEXT_ID=as.factor(infant_metadata_dynamic_phenotypes_2$NEXT_ID)
row.names(infant_metadata_dynamic_phenotypes_2)<-infant_metadata_dynamic_phenotypes_2$NG_ID
infant_metadata_dynamic_phenotypes_2<- infant_metadata_dynamic_phenotypes_2[rownames(infant_metadata_dynamic_phenotypes_2) %in% rownames(metaphlan_clr), ]

# Taxa analysis 
metadata_columns <- c("NG_ID", "NEXT_ID","Modified_NEXT_ID_without_preg_number", "days_from_first_collection","FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                      "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                      "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                      "exact_age_days_at_collection", "exact_age_months_at_collection",
                      "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                      "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                      "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                      "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint", "birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")

result_dynamic_TAXA_base <- gam_function(metaphlan_clr, infant_metadata_dynamic_phenotypes_2,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER"))
result_dynamic_TAXA_cor_delivery <- gam_function(metaphlan_clr, infant_metadata_dynamic_phenotypes_2,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary"))
result_dynamic_TAXA_cor_delivery_feeding<-gam_function(metaphlan_clr, infant_metadata_dynamic_phenotypes_2,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary","infant_ffq_ever_never_breastfed" ))



rename_except_dynamic <- function(df, suffix) {
  cols_excluded <- c("outcome", "trait", "effect.level")
  colnames(df) <- ifelse(colnames(df) %in% cols_excluded, colnames(df), paste0(colnames(df), suffix))
  return(df)
}

result_dynamic_TAXA_base_re <- rename_except_dynamic(result_dynamic_TAXA_base, "_cor_base")
result_dynamic_TAXA_cor_delivery_re  <- rename_except_dynamic(result_dynamic_TAXA_cor_delivery, "_cor_delivery")
result_dynamic_TAXA_cor_delivery_feeding_re <- rename_except_dynamic(result_dynamic_TAXA_cor_delivery_feeding , "_cor_delivery_breastfed")



combined_results_dynamic <- reduce(
  list(result_dynamic_TAXA_base_re , result_dynamic_TAXA_cor_delivery_re, result_dynamic_TAXA_cor_delivery_feeding_re),
  full_join,
  by = c("outcome", "trait", "effect.level")
)

combined_results_dynamic$trait_group <- sapply(strsplit(as.character(combined_results_dynamic$trait), "_"), function(x) paste(x[1:2], collapse = "_"))
combined_results_dynamic$type_phenotype<-"dynamic_phenotype"

combined_results_dynamic$FDR_base<-p.adjust(combined_results_dynamic$p_cor_base, method = "fdr")
combined_results_dynamic$FDR_delivery<-p.adjust(combined_results_dynamic$p_cor_delivery, method = "fdr")
combined_results_dynamic$FDR_delivery_feeding<-p.adjust(combined_results_dynamic$p_cor_delivery_breastfed, method = "fdr")

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/taxa/")
# Supplementary table S28
write.table(combined_results_dynamic, "taxa_SGB_dynamic_phenotypes_results_infants_24_09_2025.txt", sep = "\t", row.names = F)
