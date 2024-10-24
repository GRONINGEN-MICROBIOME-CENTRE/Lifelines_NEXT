### LIFELINES NEXT GUT MICROBIOME ANALYSIS ALPHA DIVERSITY INFANTS DYNAMIC PHENOTYPES ####
### AUTHOR:  TRISHLA SINHA
### ORIGINAL SCRIPT: 13TH JUL]Y, 2023
### LAST UPDATE: 10th July, 2024



rm (list = ls())

#load packages 
library(tidyverse)
library (mgcv) 
library(reshape2)
library(foreach)

##### FUNCTIONS ##########

# To remove phenotypes with too many NA's and too little variance 
drop_phenotypes <- function(phenotype_table, non_NA_values, minimum_variance) {
  cleaned_phenotypes <- phenotype_table #create a variable with phenotypic data
  discarded_phenotypes <- NULL
  for (i in 1:ncol(cleaned_phenotypes)) { 
    # check if the number of non-NA values is below the threshold
    if (length(which(!is.na(cleaned_phenotypes[,i]))) < non_NA_values) { 
      print(paste(colnames(cleaned_phenotypes)[i], "has less than n non-NA values"))
      discarded_phenotypes <- c(discarded_phenotypes, colnames(cleaned_phenotypes[i]))
      # check if the variance is below the threshold
    } else if (as.vector(sort(table(cleaned_phenotypes[,i]), decreasing = T) [1]) > length(which(!is.na(cleaned_phenotypes[,i]))) * (1-minimum_variance*0.01)) {
      print(paste(colnames(cleaned_phenotypes)[i], "has less than n% variance"))
      discarded_phenotypes <- c(discarded_phenotypes, colnames(cleaned_phenotypes[i]))
      # the phenotypes that do not fullfill the previous statements should be included
    } else {
      print(paste(colnames(cleaned_phenotypes)[i], ": accepted"))
    }
  } 
  cleaned_phenotypes[, discarded_phenotypes] <- NULL #drop non-selected phenotypes
  return(cleaned_phenotypes)
}

# To inverse rank tansform numeric data 
run_invrank_dataFrame = function(data){
  invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))} 
  data_types = unlist(lapply(data,class))
  newdata = data
  for(i in 1:ncol(data)){
    if(data_types[i]=="numeric"|data_types[i]=="integer") {
      newdata[,i] = invrank(data[,i])
    }
  }
  newdata
}

gam_function <- function(taxa, phenotypes) {
  metadata_columns <- c("NG_ID", "NEXT_ID","Modified_NEXT_ID_without_preg_number", "days_from_first_collection","FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                        "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                        "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                        "exact_age_days_at_collection", "exact_age_months_at_collection",
                        "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                        "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                        "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                        "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint")
  
  
  phenotypes2correlate <- phenotypes[, !(colnames(phenotypes) %in% metadata_columns)] 
  phenotypes2correlate[sapply(phenotypes2correlate, is.character)] <- lapply(phenotypes2correlate[sapply(phenotypes2correlate, is.character)],  #convert character columns to factors
                                                                             as.factor)
  
  covariates <-phenotypes[, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "infant_ffq_energy_kcalday_item_sum_recal", "infant_ffq_feeding_mode_complex")] 
  
  phenotypes2correlate<-run_invrank_dataFrame(phenotypes2correlate) 
  covariates<- run_invrank_dataFrame( covariates)
  
  
  result = foreach(i = 1:ncol(taxa), .combine = rbind) %:% 
    foreach(j = 1:ncol(phenotypes2correlate), .combine = rbind) %do% {
      taxa_name = colnames(taxa)[i]
      phenotype_name = colnames(phenotypes2correlate)[j]
      print(paste("Processing taxa:", taxa_name, "- Phenotype:", phenotype_name))
      data.fit <- data.frame(bac = taxa[,i],
                          trait = phenotypes2correlate[,j],
                          Time = phenotypes$exact_age_months_at_collection,
                          RD = covariates$clean_reads_FQ_1,
                          Batch =covariates$BATCH_NUMBER,
                          DNAcon = covariates$dna_conc,
                          kcal=covariates$infant_ffq_energy_kcalday_item_sum_recal,
                          feed=covariates$infant_ffq_feeding_mode_complex,
                          NEXT_ID = phenotypes$NEXT_ID)
      data.fit$NEXT_ID <- factor(data.fit$NEXT_ID)
      data.fit <- data.fit[complete.cases(data.fit),]
      
      mod_gam1 <- tryCatch({
        bam(bac ~ trait + RD + DNAcon + Batch + feed+ kcal + s(Time) + s(NEXT_ID, bs = 're'), data = data.fit, discrete = T)
      }, warning = function(w) {
        return(NULL)  
      }, error = function(e) {
        return(NULL)  
      })
      
      if (is.null(mod_gam1)) {
        converged = FALSE
        beta = statistics[,1]
        SE = statistics[,2]
        P = statistics[,4]
        N = nrow(data.fit)
      } else {
        summary.results <- summary(mod_gam1)
        statistics = summary.results$p.table[2:(which(rownames(summary.results$p.table) == 'RD') - 1),, drop = FALSE]
        converged = TRUE
        beta = statistics[,1]
        SE = statistics[,2]
        P = statistics[,4]
        N = nrow(data.fit)
      }
      
      output = data.frame(
        species = colnames(taxa)[i],
        variable = colnames(phenotypes2correlate)[j],
        effect.level = if(!is.null(statistics)) sub("^trait", "", rownames(statistics)) else NA,
        beta = beta,
        SE = SE,
        P = P,
        N = N,
        converged = converged
      )
      
      output
    }
  result
}

# Load metadata and phenotypes
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_03_08_2023.txt")
dynamic_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_longitudinal_2023_09_29.txt")
dynamic_selection<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_longitudinal_selection_AFTER_correlation_2023_11_13.txt")

# Merging files and selecting infant relevant data
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"
metadata_infants$BATCH_NUMBER<-as.factor(metadata_infants$BATCH_NUMBER)

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

# Select only ffq's & metadta 
columns_with_infant_ffq <- grep("infant_ffq", names(infant_metadata_dynamic_phenotypes_2), value = TRUE)


metadata_columns <- c("NG_ID", "NEXT_ID", "FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                      "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                      "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                      "exact_age_days_at_collection", "exact_age_months_at_collection",
                      "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                      "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                      "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                      "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint", 
                      "infant_ffq_energy_kcalday_item_sum_recal", "infant_ffq_feeding_mode_complex")


all_relevant_columns <- union(metadata_columns, columns_with_infant_ffq )


food_phenotypes <- infant_metadata_dynamic_phenotypes_2[, all_relevant_columns]
food_phenotypes$infant_ffq_vitD_quant_amount=NULL
food_phenotypes$infant_ffq_vitK_quant_amount=NULL
food_phenotypes$infant_ffq_cmpltry_food_cdc=NULL

#Taxa
taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")

# Analysis

# Select outcome 
taxa_1 <- taxa[match(rownames(infant_metadata_dynamic_phenotypes_2), rownames(taxa)),]
alpha<-infant_metadata_dynamic_phenotypes_2 %>% select(shannon)

# Select merged phenotype & metadata file 
phenotypes<-food_phenotypes


result_taxa <- gam_function(taxa_1, phenotypes)
result_taxa$FDR<-p.adjust (result_taxa$P, method = "fdr")
result_taxa$trait_group <- sapply(strsplit(as.character(result_taxa$variable), "_"), function(x) paste(x[1:2], collapse = "_"))


result_alpha <- gam_function(alpha, phenotypes)
result_alpha$FDR<-p.adjust (result_alpha$P, method = "fdr")
result_alpha$trait_group <- sapply(strsplit(as.character(result_alpha$variable), "_"), function(x) paste(x[1:2], collapse = "_"))


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results")

write.table(result_alpha, "./alpha_diversity/alpha_diversity_GAM_infants_dynamic_phenotypes_cor_feeding.txt", sep = "\t", row.names = F)
write.table(result_taxa, "./taxa/taxa_GAM_infants_dynamic_phenotypes_cor_feeding.txt", sep = "\t", row.names = F)


