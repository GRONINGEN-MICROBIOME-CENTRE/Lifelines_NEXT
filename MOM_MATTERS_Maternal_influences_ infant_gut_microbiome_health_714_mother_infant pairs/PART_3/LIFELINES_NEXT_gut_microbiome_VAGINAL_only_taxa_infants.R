### LIFELINES NEXT GUT MICROBIOME ANALYSIS TAXA INFANTS VAGINAL DELIVERY ONLY ####
### AUTHOR:  TRISHLA SINHA
### ORIGINAL SCRIPT: 3rd AUGUST, 2023
### LAST UPDATE: 26th September, 2025

rm (list = ls())

#load packages 
library(tidyverse)
library(mmrm)
library(wesanderson)
library(reshape2)

source("Functions_associations_phenotypes_LLNEXT.R")


# Load metadata and phenotypes
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
cross_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_2023_11_15.txt")
cross_selection <-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_selection_AFTER_correlation_&_correction_10_07_2024.txt")


# Merging files and selecting infant relevant data
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"

# Selecting data relevant for vaginal delivery 
infant_cross_selection <- cross_selection %>% 
  filter(only_run_vaginal_delivery == 1) %>% 
  na.omit()
column_names <- infant_cross_selection[[1]]
infant_cross_phenotypes <- cross_phenotypes %>% select(all_of(column_names))
infant_metadata_cross_phenotypes<-left_join(metadata_infants, infant_cross_phenotypes)
row.names(infant_metadata_cross_phenotypes)<-infant_metadata_cross_phenotypes$NG_ID


# Select only vaginal delivery 
infant_metadata_cross_phenotypes <- infant_metadata_cross_phenotypes %>% 
  filter(!is.na(birth_deliverybirthcard_mode_binary) & birth_deliverybirthcard_mode_binary == "VG")

infant_metadata_cross_phenotypes$birth_deliverybirthcard_mode_binary=NULL

# For this specific data remove duplicates for categorical timepoints 
infant_metadata_cross_phenotypes_1 <- infant_metadata_cross_phenotypes %>%
  distinct(SAMPLE_ID, .keep_all = TRUE) 

# Remove phenotypes with too many NA's 
infant_metadata_cross_phenotypes_2 <-drop_phenotypes(infant_metadata_cross_phenotypes_1, 200, 5)
names (infant_metadata_cross_phenotypes_2)[2]<-"NEXT_ID"
infant_metadata_cross_phenotypes_2$NEXT_ID=as.factor(infant_metadata_cross_phenotypes_2$NEXT_ID)

#Taxa
taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")

# Analysis

# Select outcome 
taxa_1 <- taxa[match(rownames(infant_metadata_cross_phenotypes_2), rownames(taxa)),]


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
results_base <- run_mmrm_analysis(taxa_1, phenos, 
                                  time = infant_metadata_cross_phenotypes_2$exact_age_months_at_collection, 
                                  time_cat = infant_metadata_cross_phenotypes_2$Timepoint_categorical, 
                                  NEXT_ID = infant_metadata_cross_phenotypes_2$NEXT_ID, 
                                  covariates = covariates_normalized_base)
results_cross_trait_base <- as.data.frame(results_base$trait)
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/taxa/")
write.table(results_cross_trait_base, "SGB_associations_VAGINAL_ONLY_27_09_2025.txt", sep = "\t", row.names = F)


