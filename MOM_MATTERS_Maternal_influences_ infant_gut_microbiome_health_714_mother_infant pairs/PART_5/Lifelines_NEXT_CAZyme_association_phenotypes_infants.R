### LIFELINES NEXT ASSOCIATION PHENOTYPES & CAZYMES ####
### AUTHOR:  TRISHLA SINHA, ALEXANDER KURULSHIKOV
### ORIGINAL SCRIPT: 9th June, 2025
### LAST UPDATE: 28th September, 2025 

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
source("Functions_associations_phenotypes_LLNEXT_new.R")

# Function to prepare cazyme data 
Prepare_cazymes = function(File = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes/NEXT_cayman_infant.txt", Prevalence_min = 0.4, Transform = "log" ){
  caz = read_tsv(File)
  caz %>% filter(! feature %in% c("total_reads", "filtered_reads")) -> caz
  #Get unique counts + ambigous counts normalized
  caz %>% select(ID, feature, combined_rpkm ) %>% spread(key = feature, value = combined_rpkm) -> caz
  #Change ID
  caz$ID = caz$ID 
  #Get summary
  summary_as_df <- function(x) {
    result <- summary(x)
    data.frame(stat = names(result), value = as.vector(result))
  }
  summary_df <- caz %>%
    select(-ID) %>%
    purrr::map(summary_as_df) %>%
    bind_rows(.id = "Cazyme") %>% spread(stat, value) %>% as_tibble()
  summary_df = summary_df %>% mutate( `NA's` = ifelse( is.na(`NA's`), 0,  `NA's` ), N = dim(caz)[1], Present= N - `NA's`)
  summary_df = summary_df %>% mutate(Present_perc = Present/N ) %>% mutate(Keep = ifelse(Present_perc>=Prevalence_min, T, F ) )
  #Make NA into 0
  caz %>% mutate_all(~ replace(., is.na(.), 0)) -> caz
  #Remove Cazymes with high NA prop
  caz %>% select(c("ID", summary_df$Cazyme[summary_df$Keep==T] )) -> caz_anal
  if (Transform == "log"){
    Log_transform(caz_anal)  ->caz_anal_tranf
  } else if (Transform == "clr"){
    CLR_transformation(caz_anal)  -> caz_anal_tranf
  } else { caz_anal_tranf = NULL}
  
  return(list("caz" = caz, "caz_filtered" = caz_anal, "caz_tranf"=caz_anal_tranf  ,"Summary_caz" = summary_df ) )
}


Log_transform = function(df){
  colnames(df)[1] = "ID"
  PS = Pseudocount(df)
  df %>% select(-ID) -> df2
  log10(df2 + PS) %>% as_tibble() %>% mutate(ID = df$ID, .before=1) %>% return()
}
Pseudocount = function(df_tara){
  df_tara %>% select(-ID) %>% as.matrix() -> PS
  min(PS[PS!=0])/2 %>% return()
}
CLR_transformation = function( df_tara ){
  df_2 = select(df_tara, -ID)
  df_2 = df_2 + Pseudocount(df_tara)
  matrix_data = as.matrix(df_2)
  
  clr_matrix <- apply(matrix_data, 1, function(row){
    exp(mean(log(row)) ) -> geometric_mean
    log(row / geometric_mean)
  } ) %>% t()
  clr_matrix %>% as_tibble() %>% mutate(ID =df_tara$ID  , .before=1) %>% return()
  
}
Filter_prevalence = function(df_tara, Prevalence = 0.4 ){
  df_tara %>% select(-ID) %>% apply(2, function(x){ mean(x!=0) }  ) -> Prevalences
  Prevalences[Prevalences>=Prevalence] %>% names() -> KEEP
  df_tara %>% select(c("ID",KEEP) ) %>% return()
}


# Load metadata 
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"
metadata_infants$Timepoint_categorical=factor(metadata_infants$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))


# Loading & processing cazyme data 

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes")
caz_tables = Prepare_cazymes(File = "NEXT_cayman_infant.txt", Prevalence_min = 0.3, Transform = "NULL" )
caz_filtered = as.data.frame (caz_tables[["caz_filtered"]])
row.names(caz_filtered)<-caz_filtered$ID
caz_filtered$ID<-NULL

# Loading microbiome data 
metaphlan = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
metaphlan_clr = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
metaphlan_infants = metaphlan[match(as.character(metadata_infants$NG_ID),rownames(metaphlan)),grep("t__",colnames(metaphlan))]
colnames(metaphlan_infants) = sub(".*s__","",colnames(metaphlan_infants))
metaphlan_infants = metaphlan_infants[,colnames(metaphlan_clr)] # Filtering to only the columns in the filtered file 
metaphlan_infants <- metaphlan_infants[rownames(metaphlan_infants) %in% rownames(caz_filtered), ]

# Data transformations 
cazyme_alr = do_clr_externalWeighting(caz_filtered,metaphlan_infants)
cazyme_alr<-as.data.frame(cazyme_alr)
cazyme_alr_null <-nullify_zeros(cazyme_alr,caz_filtered)


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
infant_metadata_cross_phenotypes_2<- infant_metadata_cross_phenotypes_2[rownames(infant_metadata_cross_phenotypes_2) %in% rownames(cazyme_alr), ]
cazyme_alr_null<- cazyme_alr_null[rownames(cazyme_alr_null) %in% rownames(infant_metadata_cross_phenotypes_2), ]


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
results_base <- run_mmrm_analysis(cazyme_alr_null, phenos, 
                                  time = infant_metadata_cross_phenotypes_2$exact_age_months_at_collection, 
                                  time_cat = infant_metadata_cross_phenotypes_2$Timepoint_categorical, 
                                  NEXT_ID = infant_metadata_cross_phenotypes_2$NEXT_ID, 
                                  covariates = covariates_normalized_base)
results_cross_trait_base <- as.data.frame(results_base$trait)

# Correction for mode of delivery
covariates_delivery <- infant_metadata_cross_phenotypes_2[,c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary")]
covariates_normalized_delivery <- run_invrank_dataFrame(covariates_delivery)
results_delivery <- run_mmrm_analysis(cazyme_alr_null, phenos, 
                                      time = infant_metadata_cross_phenotypes_2$exact_age_months_at_collection, 
                                      time_cat = infant_metadata_cross_phenotypes_2$Timepoint_categorical, 
                                      NEXT_ID = infant_metadata_cross_phenotypes_2$NEXT_ID, 
                                      covariates = covariates_normalized_delivery)
results_cross_trait_delivery <- as.data.frame(results_delivery$trait)


# Correction for mode of delivery and ever never breastfed
covariates_delivery_breastfed <- infant_metadata_cross_phenotypes_2[,c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")]
covariates_normalized_delivery_breastfed <- run_invrank_dataFrame(covariates_delivery_breastfed)
results_delivery_breastfed <- run_mmrm_analysis(cazyme_alr_null, phenos, 
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


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/cazymes")

combined_results_cross$FDR_base<-p.adjust(combined_results_cross$p_cor_base, method = "fdr")
combined_results_cross$FDR_delivery<-p.adjust(combined_results_cross$p_cor_delivery, method = "fdr")
combined_results_cross$FDR_delivery_feeding<-p.adjust(combined_results_cross$p_cor_delivery_breastfed, method = "fdr")

# Supplementary Table S44
write.table(combined_results_cross, "cazyme_cross_phenotypes_results_all_infant_29_09_2025.txt", sep = "\t", row.names = F)


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
infant_metadata_dynamic_phenotypes_2<- infant_metadata_dynamic_phenotypes_2[rownames(infant_metadata_dynamic_phenotypes_2) %in% rownames(cazyme_alr_null), ]

# Cazyme analysis 
metadata_columns <- c("NG_ID", "NEXT_ID","Modified_NEXT_ID_without_preg_number", "days_from_first_collection","FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                      "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                      "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                      "exact_age_days_at_collection", "exact_age_months_at_collection",
                      "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                      "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                      "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                      "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint", "birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")

result_dynamic_CAZ_base <- gam_function(cazyme_alr_null, infant_metadata_dynamic_phenotypes_2,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER"))
result_dynamic_CAZ_cor_delivery <- gam_function(cazyme_alr_null, infant_metadata_dynamic_phenotypes_2,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary"))
result_dynamic_CAZ_cor_delivery_feeding<-gam_function(cazyme_alr_null, infant_metadata_dynamic_phenotypes_2,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "birth_deliverybirthcard_mode_binary","infant_ffq_ever_never_breastfed" ))



rename_except_dynamic <- function(df, suffix) {
  cols_excluded <- c("outcome", "trait", "effect.level")
  colnames(df) <- ifelse(colnames(df) %in% cols_excluded, colnames(df), paste0(colnames(df), suffix))
  return(df)
}

result_dynamic_CAZ_base_re <- rename_except_dynamic(result_dynamic_CAZ_base, "_cor_base")
result_dynamic_CAZ_cor_delivery_re  <- rename_except_dynamic(result_dynamic_CAZ_cor_delivery, "_cor_delivery")
result_dynamic_CAZ_cor_delivery_feeding_re <- rename_except_dynamic(result_dynamic_CAZ_cor_delivery_feeding , "_cor_delivery_breastfed")



combined_results_dynamic <- reduce(
  list(result_dynamic_CAZ_base_re , result_dynamic_CAZ_cor_delivery_re, result_dynamic_CAZ_cor_delivery_feeding_re),
  full_join,
  by = c("outcome", "trait", "effect.level")
)

combined_results_dynamic$trait_group <- sapply(strsplit(as.character(combined_results_dynamic$trait), "_"), function(x) paste(x[1:2], collapse = "_"))
combined_results_dynamic$type_phenotype<-"dynamic_phenotype"


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/cazymes/")
combined_results_dynamic$FDR_base<-p.adjust(combined_results_dynamic$p_cor_base, method = "fdr")
combined_results_dynamic$FDR_delivery<-p.adjust(combined_results_dynamic$p_cor_delivery, method = "fdr")
combined_results_dynamic$FDR_delivery_feeding<-p.adjust(combined_results_dynamic$p_cor_delivery_breastfed, method = "fdr")

# Supplementary Table S45
write.table(combined_results_dynamic, "cazyme_dynamic_phenotypes_results_all_30_09_2025.txt", sep = "\t", row.names = F)


## Associations with HMOs in breast milk in exclusively breastfed infants  ##
hmo <-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/hmo/HMOs_matched_infant_MGS_breastfeeding_only_early_timepoints.txt")
row.names(hmo)<-hmo$NG_ID
hmo %>% filter(NG_ID %in% row.names(cazyme_alr_null) ) -> hmo
hmo_trans<-run_invrank_dataFrame(hmo)
cazyme_alr_null<- cazyme_alr_null[rownames(cazyme_alr_null) %in% rownames(hmo_trans), ]

cazymes_RO <- c(
  "CBM32", "CBM35", "CBM58", "CBM61", "CBM68", "CBM71", "GH1", "GH2", "GH4", 
  "GH5", "GH5_11", "GH5_13", "GH5_16", "GH5_20", "GH5_24", "GH5_30", "GH5_32", "GH5_33", 
  "GH5_39", "GH5_42", "GH5_49", "GH5_50", "GH5_51", "GH5_54", "GH5_56", "GH13", "GH13_2", 
  "GH13_20", "GH13_31", "GH13_34", "GH13_35", "GH13_36", "GH13_40", "GH15", "GH16", "GH16_2", 
  "GH16_3", "GH16_6", "GH16_7", "GH16_8", "GH16_9", "GH16_10", "GH16_22", "GH16_23", "GH18", 
  "GH20", "GH27", "GH29", "GH30", "GH30_5", "GH31", "GH33", "GH35", "GH36", 
  "GH39", "GH42", "GH43", "GH43_3", "GH43_8", "GH43_13", "GH43_15", "GH43_24", "GH43_25", 
  "GH43_28", "GH43_30", "GH43_31", "GH43_32", "GH43_34", "GH43_37", "GH49", "GH53", "GH57", 
  "GH66", "GH70", "GH87", "GH95", "GH97", "GH98", "GH110", "GH112", "GH117", 
  "GH136", "GH147", "GH151", "GH159", "GH165", "GH173"
)

# Only choosing cazymes responsible for resistant oligosaccha
cazyme_alr_null <- cazyme_alr_null[, colnames(cazyme_alr_null) %in% cazymes_RO]

# Select outcome and covariates 
covariate_columns=c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER")

# Cazyme analysis 
metadata_columns <- c("NG_ID", "NEXT_ID","Modified_NEXT_ID_without_preg_number", "days_from_first_collection","FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                      "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                      "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                      "exact_age_days_at_collection", "exact_age_months_at_collection",
                      "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                      "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                      "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                      "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint")

result_dynamic_CAZ_HMO_base <- gam_function(cazyme_alr_null, hmo_trans,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER"))
result_dynamic_CAZ_HMO_base$trait_group <- sapply(strsplit(as.character(result_dynamic_CAZ_HMO_base$trait), "_"), function(x) paste(x[1:2], collapse = "_"))
result_dynamic_CAZ_HMO_base<-na.omit(result_dynamic_CAZ_HMO_base)
result_dynamic_CAZ_HMO_base$FDR<-p.adjust(result_dynamic_CAZ_HMO_base$p, method = "fdr")


write.table(result_dynamic_CAZ_HMO_base, "NEXT_cayman_cazymes_hmo_GAM_01_10_2025.txt", sep = "\t", row.names = F)

# Associations of cazymes with dietary patterns (both using factor analysis and Latent Profile Analysis) at M6, M9 & M12 
dietary_patterns_factor<-read.delim("~/Desktop/LLNEXT/Analysis/dietary_clusters/infant_factor_analysis_dietary_clusters_micro.txt")
row.names(dietary_patterns_factor)<-dietary_patterns_factor$NG_ID
cazyme_alr_null<- cazyme_alr_null[rownames(cazyme_alr_null) %in% rownames(dietary_patterns_factor), ]




caz_linear_model_cor_feeding_kcal_dietary_patterns_factor_M6<- linear_model_taxa_cor_feeding_kcal_per_timepoint (dietary_patterns_factor,"NG_ID", cazyme_alr_null, 
                                                                                                                 c("infant_ffq_f1_veg_potatoes_m6",
                                                                                                                      "infant_ffq_f2_fruits_m6",
                                                                                                                      "infant_ffq_f3_tea_cheese_sweet_spreads_m6",
                                                                                                                      "infant_ffq_f4_fats_oils_nut_butters_m6",
                                                                                                                      "infant_ffq_f5_pasta_juice_meat_m6"), tp="M6")



caz_linear_model_cor_feeding_kcal_dietary_patterns_factor_M9<- linear_model_taxa_cor_feeding_kcal_per_timepoint (dietary_patterns_factor, "NG_ID", cazyme_alr_null, 
                                                                                                                 c("infant_ffq_f1_sweet_desserts_m9",
                                                                                                                            "infant_ffq_f2_wholegrain_bread_fats_oils_spreads_m9",
                                                                                                                            "infant_ffq_f3_cooked_foods_m9",
                                                                                                                            "infant_ffq_f4_starchy_foods_legumes_m9",
                                                                                                                            "infant_ffq_f5_biscuits_sweet_drinks_m9"), "M9")



caz_linear_model_cor_feeding_kcal_dietary_patterns_factor_M12<- linear_model_taxa_cor_feeding_kcal_per_timepoint (dietary_patterns_factor, "NG_ID", cazyme_alr_null, 
                                                                                                                  c("infant_ffq_f1_starchy_foods_legumes_m12",
                                                                                                                            "infant_ffq_f2_sweet_drinks_biscuits_pudding_m12",
                                                                                                                            "infant_ffq_f3_juice_dairy_drinks_m12",
                                                                                                                            "infant_ffq_f4_sweet_desserts_m12",
                                                                                                                            "infant_ffq_f5_protein_fats_oils_m12"), "M12")



# Merging all dataframes of factor analysis together for FDR correction 
all_data_frames_factor <- list(
  caz_linear_model_cor_feeding_kcal_dietary_patterns_factor_M6,
  caz_linear_model_cor_feeding_kcal_dietary_patterns_factor_M9,
  caz_linear_model_cor_feeding_kcal_dietary_patterns_factor_M12
)

merged_results_dietary_clusters_factor <- bind_rows(all_data_frames_factor )
merged_results_dietary_clusters_factor $FDR<-p.adjust (merged_results_dietary_clusters_factor $P, method = "fdr")


dietary_patterns_LPA<-read.delim("~/Desktop/LLNEXT/Analysis/dietary_clusters/infant_LPA_dietary_clusters_new.txt")

row.names(dietary_patterns_LPA)<-dietary_patterns_LPA$NG_ID
cazyme_alr_null<- cazyme_alr_null[rownames(cazyme_alr_null) %in% rownames(dietary_patterns_LPA), ]


caz_linear_model_cor_feeding_kcal_dietary_patterns_LPA_M6<- linear_model_taxa_cor_feeding_kcal_per_timepoint (dietary_patterns_LPA, "NG_ID", cazyme_alr_null, c("cluster_kday", "cluster_gday"), "M6")

caz_linear_model_cor_feeding_kcal_dietary_patterns_LPA_M9<- linear_model_taxa_cor_feeding_kcal_per_timepoint (dietary_patterns_LPA, "NG_ID", cazyme_alr_null, c("cluster_kday", "cluster_gday"), "M9")

caz_linear_model_cor_feeding_kcal_dietary_patterns_LPA_M12<- linear_model_taxa_cor_feeding_kcal_per_timepoint (dietary_patterns_LPA, "NG_ID", cazyme_alr_null, c("cluster_kday", "cluster_gday"), "M12")


# Merging all dataframes of LPA analysis together for FDR correction 
all_data_frames_LPA <- list(
  caz_linear_model_cor_feeding_kcal_dietary_patterns_LPA_M6,
  caz_linear_model_cor_feeding_kcal_dietary_patterns_LPA_M9,
  caz_linear_model_cor_feeding_kcal_dietary_patterns_LPA_M12
)
merged_results_dietary_clusters_LPA <- bind_rows(all_data_frames_LPA )
merged_results_dietary_clusters_LPA$FDR<-p.adjust (merged_results_dietary_clusters_LPA$P, method = "fdr")







#### Plotting results of interest ######


# Dynamic phenotypes 
all<-merge(infant_metadata_dynamic_phenotypes_2,cazyme_alr_null ,by="row.names" )
all$Row.names=NULL
row.names(all)<-all$NG_ID
all$Timepoint_categorical=factor(all$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))

infant_colors<-c("#b5dd88","#41c0b4", "#4397bb", "#eca4c9", "#cb4563","#a42097", "#390962")


ggplot(all, aes(x = infant_growth_weight_kg, y = GH89, color = Timepoint_categorical)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +  
  scale_color_manual(values = infant_colors) +
  labs(
    title = "Infant Growth vs. GH89 by Timepoint",
    x = "Infant Growth (kg)",
    y = "GH89",
    color = "Timepoint"
  ) +
  theme_minimal() +
  guides(
    color = guide_legend(override.aes = list(fill = NA))  
  )


ggplot(all, aes(x = infant_growth_weight_kg, y = PL2, color = Timepoint_categorical)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +  
  scale_color_manual(values = infant_colors) +
  labs(
    title = "Infant Growth vs. PL2 by Timepoint",
    x = "Infant Growth (kg)",
    y = "PL2",
    color = "Timepoint"
  ) +
  theme_minimal() +
  guides(
    color = guide_legend(override.aes = list(fill = NA))  
  )


ggplot(all, aes(x = infant_growth_weight_kg, y = PL17_1, color = Timepoint_categorical)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +  
  scale_color_manual(values = infant_colors) +
  labs(
    title = "Infant Growth vs. PL17_1 by Timepoint",
    x = "Infant Growth (kg)",
    y = "PL17_1",
    color = "Timepoint"
  ) +
  theme_minimal() +
  guides(
    color = guide_legend(override.aes = list(fill = NA))  
  )


results_df <- result_dynamic_CAZ_HMO_base %>%
  mutate(trait_full = ifelse(is.na(effect.level), trait,
                             paste0(trait, "_", effect.level)))


results_df <- results_df %>%
  mutate(label = ifelse(p < 0.05, "*", ""))

# Step 3: Create matrix of Estimate values for heatmap
heatmap_matrix <- results_df %>%
  select(outcome, trait_full, Estimate) %>%
  pivot_wider(names_from = trait_full, values_from = Estimate) %>%
  column_to_rownames("outcome") %>%
  as.matrix()

# Step 4: Create matrix of labels (with significance stars)
label_matrix <- results_df %>%
  select(outcome, trait_full, label) %>%
  pivot_wider(names_from = trait_full, values_from = label) %>%
  column_to_rownames("outcome") %>%
  as.matrix()


pheatmap(heatmap_matrix,
         display_numbers = label_matrix,   # Shows numbers + stars
         number_color = "black",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "CBM32-HMO Associations: Estimates with Significance",
         border_color = NA)


################ HEATMAPS TO DIEATARY PATTERNS ##################


mat_factor <- merged_results_dietary_clusters_factor %>%
  mutate(Feature_Time = paste0(Feature, "_", Timepoint)) %>%
  select(Outcome, Feature_Time, Estimate) %>%
  pivot_wider(names_from = Feature_Time, values_from = Estimate) %>%
  column_to_rownames("Outcome") %>%
  as.matrix()


mat_factor <- mat_factor[rowSums(is.na(mat_factor)) < ncol(mat_factor), ]
mat_factor <- mat_factor[, colSums(is.na(mat_factor)) < nrow(mat_factor)]


pheatmap(mat_factor,
         main = "CAZyme ~ Dietary Patterns (Factor Analysis)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6,
         fontsize_col = 7,
         border_color = NA,
         na_col = "grey90")

# ---- LPA HEATMAP ----

mat_LPA <- merged_results_dietary_clusters_LPA %>%
  mutate(Feature_Time = paste0(Feature, "_", Timepoint)) %>%
  select(Outcome, Feature_Time, Estimate) %>%
  pivot_wider(names_from = Feature_Time, values_from = Estimate) %>%
  column_to_rownames("Outcome") %>%
  as.matrix()

# Optional: remove all-NA rows/columns
mat_LPA <- mat_LPA[rowSums(is.na(mat_LPA)) < ncol(mat_LPA), ]
mat_LPA <- mat_LPA[, colSums(is.na(mat_LPA)) < nrow(mat_LPA)]

# Plot heatmap for LPA
pheatmap(mat_LPA,
         main = "CAZyme ~ Dietary Patterns (LPA)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6,
         fontsize_col = 7,
         border_color = NA,
         na_col = "grey90")
