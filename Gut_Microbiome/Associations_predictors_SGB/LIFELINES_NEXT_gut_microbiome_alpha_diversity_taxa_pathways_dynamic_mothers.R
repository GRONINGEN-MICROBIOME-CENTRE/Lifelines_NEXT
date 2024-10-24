### LIFELINES NEXT GUT MICROBIOME ANALYSIS ALPHA DIVERSITY & TAXA MOTHERS DYNAMIC PHENOTYPES ####
### AUTHOR:  Trishla Sinha & Alexander Kurilshikov
### ORIGINAL SCRIPT: 12th September, 2023
### LAST UPDATE: 26th of June, 2024


rm (list = ls())

#load packages 
library(tidyverse)
library (mgcv) 
library(reshape2)
library(foreach)
library(wesanderson)

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
  phenotypes2correlate[sapply(phenotypes2correlate, is.character)] <- lapply(phenotypes2correlate[sapply(phenotypes2correlate, is.character)], as.factor)
  
  covariates <- phenotypes[, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER")] 
  
  phenotypes2correlate <- run_invrank_dataFrame(phenotypes2correlate) 
  covariates <- run_invrank_dataFrame(covariates)
  
  result = foreach(i = 1:ncol(taxa), .combine = rbind) %:% 
    foreach(j = 1:ncol(phenotypes2correlate), .combine = rbind) %do% {
      print(paste(i, j))
      data.fit <- data.frame(bac = taxa[,i],
                             trait = phenotypes2correlate[,j],
                             Time = phenotypes$exact_age_months_at_collection,
                             RD = covariates$clean_reads_FQ_1,
                             Batch = covariates$BATCH_NUMBER,
                             DNAcon = covariates$dna_conc,
                             NEXT_ID = phenotypes$NEXT_ID)
      data.fit$NEXT_ID <- factor(data.fit$NEXT_ID)
      data.fit <- data.fit[complete.cases(data.fit),]
      
      mod_gam1 <- tryCatch({
        bam(bac ~ trait + RD + DNAcon + Batch + s(Time) + s(NEXT_ID, bs = 're'), data = data.fit, discrete = TRUE)
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
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
dynamic_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_longitudinal_2023_09_29.txt")
dynamic_selection<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_longitudinal_selection_AFTER_correlation_&_correction_17_04_2024.txt")

# Merging files and selecting mother relevant data
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_mothers<-metadata[metadata$Type=="mother", ]

metadata_mothers$BATCH_NUMBER<-as.factor(metadata_mothers$BATCH_NUMBER)

# Selecting data relevant for mother 
mother_dynamic_selection<-dynamic_selection[dynamic_selection$mother_microbiome_selection==1,]
column_names <- mother_dynamic_selection[[1]]
mother_dynamic_phenotypes <- dynamic_phenotypes %>% select(all_of(column_names))
mother_dynamic_phenotypes<-mother_dynamic_phenotypes[, c(9:51)]

mother_dynamic_phenotypes <- mother_dynamic_phenotypes %>%
  distinct(SAMPLE_ID, .keep_all = TRUE)

mother_metadata_dynamic_phenotypes<-left_join(metadata_mothers, mother_dynamic_phenotypes)
row.names(mother_metadata_dynamic_phenotypes)<-mother_metadata_dynamic_phenotypes$NG_ID

# For this specific data remove duplicates for categorical timepoints 
mother_metadata_dynamic_phenotypes_1 <- mother_metadata_dynamic_phenotypes %>%
  distinct(SAMPLE_ID, .keep_all = TRUE) 

# Remove phenotypes with too many NA's 
mother_metadata_dynamic_phenotypes_2 <-drop_phenotypes(mother_metadata_dynamic_phenotypes_1, 200, 5)
mother_metadata_dynamic_phenotypes_2$NEXT_ID=as.factor(mother_metadata_dynamic_phenotypes_2$NEXT_ID)

#Taxa
taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/NEXT_metaphlan_4_CLR_transformed_fil_30_percent_SGB_mothers_03_08_2023.txt")

# Pathways 
pathways<-read.delim("~/Desktop/LLNEXT/Analysis/pathways/NEXT_humann_pathways_CLR_transformed_fil_30_percent_26_06_2024.txt")

# Analysis

# Select outcome 
taxa_1 <- taxa[match(rownames(mother_metadata_dynamic_phenotypes_2), rownames(taxa)),]
alpha<-mother_metadata_dynamic_phenotypes_2 %>% select(shannon)
pathways_1 <- pathways[match(rownames(mother_metadata_dynamic_phenotypes_2), rownames(pathways)),]

# Select merged phenotype & metadata file 
phenotypes<-mother_metadata_dynamic_phenotypes_2


# Associations with alpha diversity
result_alpha <- gam_function(alpha, phenotypes)
result_alpha<-result_alpha[result_alpha$converged=="TRUE",]
result_alpha$FDR<-p.adjust (result_alpha$P, method = "fdr")
result_alpha$trait_group <- sapply(strsplit(as.character(result_alpha$variable), "_"), function(x) paste(x[1:2], collapse = "_"))

# Associations with taxa (at SGB level)
result_taxa <- gam_function(taxa_1, phenotypes)
result_taxa<-result_taxa[result_taxa$converged=="TRUE",]
result_taxa$FDR<-p.adjust (result_taxa$P, method = "fdr")
result_taxa$trait_group <- sapply(strsplit(as.character(result_taxa$variable), "_"), function(x) paste(x[1:2], collapse = "_"))

# Associations with pathways
result_pathways <- gam_function(pathways_1, phenotypes)
result_pathways<-result_pathways[result_pathways$converged=="TRUE",]
result_pathways$FDR<-p.adjust (result_pathways$P, method = "fdr")
result_pathways$trait_group <- sapply(strsplit(as.character(result_pathways$variable), "_"), function(x) paste(x[1:2], collapse = "_"))


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results")

write.table(result_alpha, "./alpha_diversity/mothers/alpha_diversity_GAM_mothers_dynamic_phenotypes.txt", sep = "\t", row.names = F)
write.table(result_taxa, "./taxa/mothers/taxa_GAM_mothers_dynamic_phenotypes.txt", sep = "\t", row.names = F)
write.table(result_pathways, "./pathways/mothers/pathways_GAM_mothers_dynamic_phenotypes.txt", sep = "\t", row.names = F)

mother_metadata_dynamic_phenotypes_2$Timepoint_categorical<-factor(mother_metadata_dynamic_phenotypes_2$Timepoint_categorical, levels = c("P12", "P28", "B", "M3"))

