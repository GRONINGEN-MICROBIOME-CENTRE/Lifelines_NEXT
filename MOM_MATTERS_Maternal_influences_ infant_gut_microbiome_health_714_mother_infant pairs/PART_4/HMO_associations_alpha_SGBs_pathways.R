################## Associations HMOs with bacterial taxa and metabolic pathways ####### 
# Author: Trishla Sinha 
# Last update: 7th October, 2025

# Load all functions from script : 
source("Functions_associations_phenotypes_LLNEXT_new.R")

## Associations with HMOs in breast milk in exclusively breastfed infants  ##
hmo_all <-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/hmo/HMOs_matched_infant_MGS_breastfeeding_only_early_timepoints.txt")
row.names(hmo_all)<-hmo_all$NG_ID

hmo_colnames <- names(hmo_all)[grepl("^mother_milk_HMO", names(hmo_all))]
# Tested HMO's 
# "mother_milk_HMO_milk_group"     "mother_milk_HMO_Le"             "mother_milk_HMO_Se"            
# "mother_milk_HMO_2FL_ugml"       "mother_milk_HMO_3FL_ugml"       "mother_milk_HMO_LDFT_ugml"     
# "mother_milk_HMO_3GL_ugml"       "mother_milk_HMO_6GL_ugml"       "mother_milk_HMO_A_tetra_ugml"  
# "mother_milk_HMO_3SL_ugml"       "mother_milk_HMO_LNT_ugml"       "mother_milk_HMO_LNnT_ugml"     
# "mother_milk_HMO_6SL_ugml"       "mother_milk_HMO_3F3SL_ugml"     "mother_milk_HMO_LNFP_V_ugml"   
# "mother_milk_HMO_LNnFP_V_ugml"   "mother_milk_HMO_LNFP_I_ugml"    "mother_milk_HMO_LNFP_III_ugml" 
# "mother_milk_HMO_LNFP_II_ugml"   "mother_milk_HMO_LNnDFH_ugml"    "mother_milk_HMO_LSTb_ugml"     
# "mother_milk_HMO_LNDH_I_ugml"    "mother_milk_HMO_LSTc_ugml"      "mother_milk_HMO_LNH_ugml"      
# "mother_milk_HMO_MFLNH_III_ugml" "mother_milk_HMO_DSLNT_ugml"     "mother_milk_HMO_DFLNHa_ugml"   
# "mother_milk_HMO_Total_ugml"     "mother_milk_HMO_Fuc_ugml"       "mother_milk_HMO_Neut_ugml"     
# "mother_milk_HMO_Sia_ugml"    

# Inverse transformation 
hmo_all[hmo_colnames] <- run_invrank_dataFrame(hmo_all[hmo_colnames])


# With infant alpha diversity 
alpha<-hmo_all %>% select(shannon)

# With Infant SGB's 
taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
taxa <- taxa[match(rownames(hmo_all), rownames(taxa)),]

# With pathways 
pathways<-read.delim("~/Desktop/LLNEXT/Analysis/pathways/metacyc_infants_ALR_zeroTreated_20250722.tsv")
pathways <- pathways[match(rownames(hmo_all), rownames(pathways)),]


# Cazyme analysis 
metadata_columns <- c("NG_ID", "NEXT_ID","Modified_NEXT_ID_without_preg_number", "days_from_first_collection","FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                      "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                      "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                      "exact_age_days_at_collection", "exact_age_months_at_collection",
                      "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                      "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                      "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                      "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint")

result_HMO_alpha <- gam_function(alpha, hmo_all,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER"))
result_HMO_taxa <- gam_function(taxa, hmo_all,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER"))
result_HMO_pathway <- gam_function(pathways, hmo_all,metadata_columns, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER"))


result_HMO_alpha$FDR<-p.adjust (result_HMO_alpha$p, method = "fdr")
result_HMO_taxa$FDR<-p.adjust (result_HMO_taxa$p, method = "fdr")
result_HMO_pathway$FDR<-p.adjust (result_HMO_pathway$p, method = "fdr")
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/hmo")

write.table(result_HMO_alpha, "result_HMO_alpha_oct_2025.txt", sep = "\t", row.names = F)
write.table(result_HMO_taxa, "result_HMO_taxa_SGB_oct_2025.txt", sep = "\t", row.names = F)
write.table(result_HMO_pathway, "result_HMO_pathway_oct_2025.txt", sep = "\t", row.names = F)

