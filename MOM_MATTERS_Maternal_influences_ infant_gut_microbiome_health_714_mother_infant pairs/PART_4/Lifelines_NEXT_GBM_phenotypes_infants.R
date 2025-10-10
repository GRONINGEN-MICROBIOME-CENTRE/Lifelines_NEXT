### LIFELINES NEXT ASSOCIATION PHENOTYPES & GUT-BRAIN MODULES ####
### AUTHOR:  TRISHLA SINHA, ALEXANDER KURULSHIKOV
### ORIGINAL SCRIPT: 9th June, 2025
### LAST UPDATE: 1st October, 2025


# Load packages 
library(tidyverse)
library(mmrm)
library(wesanderson)
library(dplyr)
library(tidyverse)
library (mgcv) 
library(reshape2)
library(foreach)


# Load all functions from script : Functions_associations_phenotypes_LLNEXT 

# Load metadata and cross sectional phenotypes
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
cross_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_2023_11_15.txt")
cross_selection <-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_selection_AFTER_correlation_&_correction_10_07_2024.txt")


# Merging files and selecting infant relevant data
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"
metadata_infants$Timepoint_categorical=factor(metadata_infants$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))

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

# Analysis

# Select outcome and match with phenotypic data 

GBM<-read.delim("~/Desktop/LLNEXT/Analysis/GBM/NEXT_GBM_merged_clean.tsv")
row.names(GBM)<-GBM$ID
GBM$ID=NULL
GBM<- na.omit(GBM)

# Setting a prevalence cut-off 
infant_NEXT_ID<-infant_metadata_cross_phenotypes_2 %>%
  select(NEXT_ID)
infant_GBM_all<-merge(infant_NEXT_ID,GBM, by="row.names" )
row.names(infant_GBM_all)<-infant_GBM_all$Row.names
infant_GBM_all$Row.names=NULL
unique_counts <- sapply(infant_GBM_all, function(x) length(unique(infant_GBM_all$NEXT_ID[x >0.001]))) 
infant_GBM_all_filt <- infant_GBM_all[, unique_counts >= 0.3*length(unique(infant_GBM_all$NEXT_ID)) ] # Setting a 30% cut-off on prevalence 
infant_GBM_all_filt$NEXT_ID=NULL


## Loading microbiome data 
metaphlan = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
metaphlan_clr = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
metaphlan_infants = metaphlan[match(as.character(metadata_infants$NG_ID),rownames(metaphlan)),grep("t__",colnames(metaphlan))]
colnames(metaphlan_infants) = sub(".*s__","",colnames(metaphlan_infants))
metaphlan_infants = metaphlan_infants[,colnames(metaphlan_clr)] # Filtering to only the columns in the filtered file 


metaphlan_infants<- metaphlan_infants[rownames(metaphlan_infants) %in% rownames(infant_GBM_all_filt), ]
infant_GBM_all_filt<- infant_GBM_all_filt[rownames(infant_GBM_all_filt) %in% rownames(metaphlan_infants), ]
infant_metadata_cross_phenotypes_2<- infant_metadata_cross_phenotypes_2[rownames(infant_metadata_cross_phenotypes_2) %in% rownames(infant_GBM_all_filt), ]

# Data transformations 
GBM_alr <-do_clr_externalWeighting(infant_GBM_all_filt,metaphlan_infants)
GBM_alr<-as.data.frame(GBM_alr)
GBM_alr_null <-nullify_zeros(GBM_alr,infant_GBM_all_filt)


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
results_base <- run_mmrm_analysis(GBM_alr_null, phenos, 
                                  time = infant_metadata_cross_phenotypes_2$exact_age_months_at_collection, 
                                  time_cat = infant_metadata_cross_phenotypes_2$Timepoint_categorical, 
                                  NEXT_ID = infant_metadata_cross_phenotypes_2$NEXT_ID, 
                                  covariates = covariates_normalized_base)
results_cross_trait_base <- as.data.frame(results_base$trait)


# Correction for mode of delivery
covariates_delivery <- infant_metadata_cross_phenotypes_2[,c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary")]
covariates_normalized_delivery <- run_invrank_dataFrame(covariates_delivery)
results_delivery <- run_mmrm_analysis(GBM_alr_null, phenos, 
                                      time = infant_metadata_cross_phenotypes_2$exact_age_months_at_collection, 
                                      time_cat = infant_metadata_cross_phenotypes_2$Timepoint_categorical, 
                                      NEXT_ID = infant_metadata_cross_phenotypes_2$NEXT_ID, 
                                      covariates = covariates_normalized_delivery)
results_cross_trait_delivery <- as.data.frame(results_delivery$trait)


# Correction for mode of delivery and ever never breastfed
covariates_delivery_breastfed <- infant_metadata_cross_phenotypes_2[,c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")]
covariates_normalized_delivery_breastfed <- run_invrank_dataFrame(covariates_delivery_breastfed)
results_delivery_breastfed <- run_mmrm_analysis(GBM_alr_null, phenos, 
                                                time = infant_metadata_cross_phenotypes_2$exact_age_months_at_collection, 
                                                time_cat = infant_metadata_cross_phenotypes_2$Timepoint_categorical, 
                                                NEXT_ID = infant_metadata_cross_phenotypes_2$NEXT_ID, 
                                                covariates = covariates_normalized_delivery_breastfed)
results_cross_trait_delivery_breastfed <- as.data.frame(results_delivery_breastfed$trait)


rename_except <- function(df, suffix) {
  cols_excluded <- c("outcome", "trait","levels")
  colnames(df) <- ifelse(colnames(df) %in% cols_excluded, colnames(df), paste0(colnames(df), suffix))
  return(df)
}

results_cross_trait_base_re <- rename_except(results_cross_trait_base, "_cor_base")
results_cross_trait_delivery_re  <- rename_except(results_cross_trait_delivery, "_cor_delivery")
results_cross_trait_delivery_breastfed_re <- rename_except(results_cross_trait_delivery_breastfed , "_cor_delivery_breastfed")


combined_results_cross <- reduce(
  list(results_cross_trait_base_re, results_cross_trait_delivery_re, results_cross_trait_delivery_breastfed_re),
  full_join,
  by = c("outcome", "trait", "levels")
)
combined_results_cross$trait_group <- sapply(strsplit(as.character(combined_results_cross$trait), "_"), function(x) paste(x[1:2], collapse = "_"))
combined_results_cross$type_phenotype<-"cross_sectional_phenotype"


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/GBM")
description<-read.delim("MGBM_pathway_descriptions.txt")
names (description)[1]<-"outcome"
combined_results_cross<-left_join(combined_results_cross, description)
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/GBM/")
combined_results_cross$FDR_base<-p.adjust(combined_results_cross$p_cor_base, method = "fdr")
combined_results_cross$FDR_delivery<-p.adjust(combined_results_cross$p_cor_delivery, method = "fdr")
combined_results_cross$FDR_delivery_feeding<-p.adjust(combined_results_cross$p_cor_delivery_breastfed, method = "fdr")
# Supplementary table S31
write.table(combined_results_cross, "GBM_cross_phenotypes_results_all_01_10_2025.txt", sep = "\t", row.names = F)


# Load metadata and dynamic phenotypes
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/metadata_basic_phenotypes.txt")
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
names (infant_metadata_dynamic_phenotypes_2)[2]<-"NEXT_ID"
infant_metadata_dynamic_phenotypes_2$NEXT_ID=as.factor(infant_metadata_dynamic_phenotypes_2$NEXT_ID)
infant_metadata_dynamic_phenotypes_2<- infant_metadata_dynamic_phenotypes_2[rownames(infant_metadata_dynamic_phenotypes_2) %in% rownames(infant_GBM_all_filt), ]

# Select outcome and covariates 
covariate_columns=c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER")

# GBM analysis 
metadata_columns <- c("NG_ID", "NEXT_ID","Modified_NEXT_ID_without_preg_number", "days_from_first_collection","FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                      "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                      "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                      "exact_age_days_at_collection", "exact_age_months_at_collection",
                      "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                      "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                      "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                      "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint", "birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")

result_dynamic_GBM_base <- gam_function(GBM_alr_null, infant_metadata_dynamic_phenotypes_2,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER"))
result_dynamic_GBM_cor_delivery <- gam_function(GBM_alr_null, infant_metadata_dynamic_phenotypes_2,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary"))
result_dynamic_GBM_cor_delivery_feeding<-gam_function(GBM_alr_null, infant_metadata_dynamic_phenotypes_2,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary","infant_ffq_ever_never_breastfed" ))



rename_except <- function(df, suffix) {
  cols_excluded <- c("outcome", "trait", "effect.level")
  colnames(df) <- ifelse(colnames(df) %in% cols_excluded, colnames(df), paste0(colnames(df), suffix))
  return(df)
}

result_dynamic_GBM_base_re <- rename_except(result_dynamic_GBM_base, "_cor_base")
result_dynamic_GBM_cor_delivery_re <- rename_except(result_dynamic_GBM_cor_delivery, "_cor_delivery")
result_dynamic_GBM_cor_delivery_feeding_re <- rename_except(result_dynamic_GBM_cor_delivery_feeding, "_cor_delivery_breastfed")


combined_results_dynamic <- reduce(
  list(result_dynamic_GBM_base_re, result_dynamic_GBM_cor_delivery_re, result_dynamic_GBM_cor_delivery_feeding_re),
  full_join,
  by = c("outcome", "trait", "effect.level")
)

combined_results_dynamic$trait_group <- sapply(strsplit(as.character(combined_results_dynamic$trait), "_"), function(x) paste(x[1:2], collapse = "_"))
combined_results_dynamic$type_phenotype<-"dynamic_phenotype"


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/GBM")
description<-read.delim("MGBM_pathway_descriptions.txt")
names (description)[1]<-"outcome"
combined_results_dynamic<-left_join(combined_results_dynamic, description)
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/GBM/")

combined_results_dynamic$FDR_base<-p.adjust(combined_results_dynamic$p_cor_base, method = "fdr")
combined_results_dynamic$FDR_delivery<-p.adjust(combined_results_dynamic$p_cor_delivery, method = "fdr")
combined_results_dynamic$FDR_delivery_feeding<-p.adjust(combined_results_dynamic$p_cor_delivery_breastfed, method = "fdr")
# Supplementary table S32
write.table(combined_results_dynamic, "GBM_dynamic_phenotypes_results_all_01_10_2025.txt", sep = "\t", row.names = F)



#### Plots of significant results ######


all<-merge(infant_metadata_cross_phenotypes_2,GBM ,by="row.names" )
all$Row.names=NULL
row.names(all)<-all$NG_ID

# Parity 
ggplot(all, aes(x = mother_birthcard_parity, y = MGB036)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Scatter Plot of mother_birthcard_parity vs. MGB036",
    x = "Mother Birthcard Parity",
    y = "MGB036"
  ) +
  theme_minimal()

# Feeding 
clean_data <- all %>%
  filter(!is.na(infant_birthcard_feeding_mode_after_birth)) %>%
  mutate(infant_birthcard_feeding_mode_after_birth = factor(
    infant_birthcard_feeding_mode_after_birth,
    levels = c("BF", "MF", "FF")  
  ))

ggplot(clean_data, aes(x = infant_birthcard_feeding_mode_after_birth, y = MGB043, fill = infant_birthcard_feeding_mode_after_birth)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  
  geom_jitter(aes(color = infant_birthcard_feeding_mode_after_birth), width = 0.2, alpha = 0.6, size = 2) +
  labs(
    title = "Acetate synthesis I",
    x = "Feeding Mode",
    y = "MGB043 Abundance"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

head (clean_data)


### DYNAMIC phenotypes 

all<-merge(infant_metadata_dynamic_phenotypes_2,GBM_alr_null ,by="row.names" )
all$Row.names=NULL
row.names(all)<-all$NG_ID

# Feeding complex 
clean_data <- all %>%
  filter(!is.na(infant_ffq_feeding_mode_complex)) %>%
  mutate(infant_ffq_feeding_mode_complex = factor(
    infant_ffq_feeding_mode_complex,
    levels = c("BF", "MF", "FF")  
  ))

ggplot(clean_data, aes(x = infant_ffq_feeding_mode_complex, y = MGB031, fill = infant_ffq_feeding_mode_complex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Hide default outliers
  geom_jitter(aes(color = infant_ffq_feeding_mode_complex), width = 0.2, alpha = 0.6, size = 2) +
  labs(
    title = "17-beta-Estradiol degradation",
    x = "Feeding Mode",
    y = "MGB031 Abundance"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# Feeding simple 
clean_data <- all %>%
  filter(!is.na(infant_ffq_feeding_mode_complex)) %>%
  mutate(infant_ffq_feeding_mode_complex = factor(
    infant_ffq_feeding_mode_complex,
    levels = c("BF", "MF", "FF")  
  ))

ggplot(clean_data, aes(x = infant_ffq_feeding_mode_complex, y = MGB031, fill = infant_ffq_feeding_mode_complex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Hide default outliers
  geom_jitter(aes(color = infant_ffq_feeding_mode_complex), width = 0.2, alpha = 0.6, size = 2) +
  labs(
    title = "17-beta-Estradiol degradation",
    x = "Feeding Mode",
    y = "MGB031 Abundance"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggplot(clean_data, aes(
  x = infant_ffq_feeding_mode_complex,
  y = MGB043,
  fill = infant_ffq_feeding_mode_complex
)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Hide default outliers
  geom_jitter(aes(color = infant_ffq_feeding_mode_complex),
              width = 0.2, alpha = 0.6, size = 2) +
  labs(
    title = "",
    x = "Feeding Mode",
    y = "Acetate synthesis I"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  facet_wrap(~ Timepoint_categorical)  # Separate plot per timepoint

clean_data <- all %>%
  filter(!is.na(infant_ffq_feeding_mode_simple)) %>%
  mutate(infant_ffq_feeding_mode_simple = factor(
    infant_ffq_feeding_mode_simple,
    levels = c("excl_BF", "excl_FF")  
  ))




acetate_synthesis<-ggplot(clean_data, aes(
  x = infant_ffq_feeding_mode_simple,
  y = MGB043,
  fill = infant_ffq_feeding_mode_simple
)) +
  
  geom_jitter(aes(color = infant_ffq_feeding_mode_simple),
              width = 0.2, alpha = 0.6, size = 2) +
  labs(
    title = "",
    x = "Feeding Mode",
    y = "Acetate synthesis I"
  ) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  facet_wrap(~ Timepoint_categorical, nrow = 1)

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_3_new")
ggsave("acetate_synthesis_by_feeding_mode.pdf", plot = acetate_synthesis,
       width = 4, height = 4, units = "in")

