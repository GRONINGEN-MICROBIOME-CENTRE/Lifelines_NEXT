################### PERMANOVA EACH TIMEPOINT INFANTS AND MOTHERS ##########################
# Author: T.Sinha 
# Last updated: 26th of August, 2024

library(tidyverse)
library(vegan)
library(foreach)

rm (list = ls())

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


##############################
# Mother 
##############################
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/phenotypes")
metadata_mothers<-read.delim("mother_phenotypes_all_microbiome_fil.txt")
metadata_mothers <- subset(metadata_mothers, Timepoint_categorical != "M1" & Timepoint_categorical != "M2")
row.names(metadata_mothers)<-metadata_mothers$NG_ID

taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
mother_taxa<-taxa[row.names(taxa)%in% rownames(metadata_mothers),] 
mother_taxa<-mother_taxa[match(row.names(metadata_mothers),row.names(mother_taxa)),]

# Filtering data 
mother_NEXT_ID<-metadata_mothers %>%
  select(next_id_mother)
mother_taxa_all<-merge(mother_NEXT_ID,mother_taxa, by="row.names" )
row.names(mother_taxa_all)<-mother_taxa_all$Row.names
mother_taxa_all$Row.names=NULL
unique_counts <- sapply(mother_taxa_all, function(x) length(unique(mother_taxa_all$next_id_mother[x >0.001]))) # Here I am counting the unique elements in the NEXT_ID column where the corresponding value in each column (i.e., x) is greater than the given cut-off. 
mother_taxa_all_filt <- mother_taxa_all[, unique_counts >= 0.3*length(unique(mother_taxa_all$next_id_mother)) ] # Setting a 20% cut-off on prevalence (3rd august version) 
mother_taxa_all_filt$next_id_mother=NULL

mother_taxa_SGB<-mother_taxa_all_filt[,grep("t__",colnames(mother_taxa_all_filt))]
my_pseudocount_normal=min(mother_taxa_SGB[mother_taxa_SGB!=0])/2# 
distance=vegdist(mother_taxa_SGB, method = "aitchison", pseudocount=my_pseudocount_normal) 


Call_adonis = function(phenotypes,Phenotype_list, Distance, perm=10000, cores=8){
  Distance = as.matrix(Distance)
  adonis_results<- data.frame(matrix(ncol =5, nrow= ncol(phenotypes)))  
  colnames(adonis_results) <- c("Phenotype","Df", "F", "R2", "p-value")
  Phenotype_list = Phenotype_list[!is.na(Phenotype_list)]
  adon<-foreach(i=1:length(Phenotype_list),.combine=rbind)%do%{
    #This is the phentype vector
    Pheno = Phenotype_list[i]
    print(Pheno)
    as_vector(phenotypes[Pheno]) %>% as.vector() -> Values
    #Filter NAs
    r<-row.names(phenotypes)[!is.na(Values)]
    phenotypes[r,] -> phenos2
    #Match with Distance
    distmat_cleaned <- Distance[r,r]
    #Run adonis
    FR = formula( paste0("distmat_cleaned ~ clean_reads_FQ_1 + dna_conc + BATCH_NUMBER +", Pheno) )
    ad1<-adonis2(FR , phenos2, permutations=perm,parallel=cores, na.action = na.fail)
                 
    adonis_results[i,] <- c(Pheno,ad1$Df[4], ad1$F[4], ad1$R2[4],ad1$'Pr(>F)'[4])
  }
  
  adonis_results %>% drop_na() %>% return()
}
ResultsAdonis = tibble()

not_to_test<- c("NG_ID", "NEXT_ID", "FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2","next_id_mother","Modified_NEXT_ID_without_preg_number" ,"days_from_first_collection"  ,
                                   "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                                   "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                                   "exact_age_days_at_collection", "exact_age_months_at_collection",
                                   "Timepoint_categorical", "SAMPLE_ID",
                                   "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                                   "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                                   "shannon", "BATCH_NUMBER", "next_id_partner", "sibling_number", "timepoint")
              


for (Timepoint in unique(metadata_mothers$Timepoint_categorical)){
  print(Timepoint)
  phenotypes<-metadata_mothers %>% filter(Timepoint_categorical == Timepoint )
  #Put function stuff for phenotype selection. Phenolist
  PP=drop_phenotypes(phenotypes, non_NA_values = 100, minimum_variance = 5)
  PP %>% select(-one_of(not_to_test)) %>% colnames() -> PhenosTest
  #input selected phenotypes, as a vector, to Phenotype_list
  Call_adonis(phenotypes, Phenotype_list=PhenosTest ,perm=10000, Distance=distance) %>% mutate(Timepoint = Timepoint) %>% rbind(ResultsAdonis, .) -> ResultsAdonis
}


ResultsAdonis$Timepoint <- factor(ResultsAdonis$Timepoint, levels = c("P12", "P28", "B", "M3"))
  
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint")
write.table(ResultsAdonis, "result_adonis_per_timepoint_mothers_10000_perm.txt", sep = "\t", row.names = F)

##############################
# Infant 
##############################
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/phenotypes")
metadata_infants<-read.delim("infant_phenotypes_all_microbiome_fil.txt")
row.names(metadata_infants)<-metadata_infants$NG_ID

taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
infant_taxa<-taxa[row.names(taxa)%in% rownames(metadata_infants),] 
infant_taxa<-infant_taxa[match(row.names(metadata_infants),row.names(infant_taxa)),]

# Filtering data 
infant_NEXT_ID<-metadata_infants %>%
  select(next_id_infant)
infant_taxa_all<-merge(infant_NEXT_ID,infant_taxa, by="row.names" )
row.names(infant_taxa_all)<-infant_taxa_all$Row.names
infant_taxa_all$Row.names=NULL
unique_counts <- sapply(infant_taxa_all, function(x) length(unique(infant_taxa_all$next_id_infant[x >0.1]))) # Here I am counting the unique elements in the NEXT_ID column where the corresponding value in each column (i.e., x) is greater than the given cut-off. 
infant_taxa_all_filt <- infant_taxa_all[, unique_counts >= 0.1*length(unique(infant_taxa_all$next_id_infant)) ] # Setting a 20% cut-off on prevalence (3rd august version) 
infant_taxa_all_filt$next_id_infant=NULL

infant_taxa_SGB<-infant_taxa_all_filt[,grep("t__",colnames(infant_taxa_all_filt))]
my_pseudocount_normal=min(infant_taxa_SGB[infant_taxa_SGB!=0])/2# 
distance=vegdist(infant_taxa_SGB, method = "aitchison", pseudocount=my_pseudocount_normal) 


Call_adonis = function(phenotypes,Phenotype_list, Distance, perm=10000, cores=8){
  Distance = as.matrix(Distance)
  adonis_results<- data.frame(matrix(ncol =5, nrow= ncol(phenotypes)))  
  colnames(adonis_results) <- c("Phenotype","Df", "F", "R2", "p-value")
  Phenotype_list = Phenotype_list[!is.na(Phenotype_list)]
  adon<-foreach(i=1:length(Phenotype_list),.combine=rbind)%do%{
    #This is the phentype vector
    Pheno = Phenotype_list[i]
    print(Pheno)
    as_vector(phenotypes[Pheno]) %>% as.vector() -> Values
    #Filter NAs
    r<-row.names(phenotypes)[!is.na(Values)]
    phenotypes[r,] -> phenos2
    #Match with Distance
    distmat_cleaned <- Distance[r,r]
    #Run adonis
    FR = formula( paste0("distmat_cleaned ~ clean_reads_FQ_1 + dna_conc + BATCH_NUMBER +", Pheno) )
    ad1<-adonis2(FR , phenos2, permutations=perm,parallel=cores, na.action = na.fail)
    
    adonis_results[i,] <- c(Pheno,ad1$Df[4], ad1$F[4], ad1$R2[4],ad1$'Pr(>F)'[4])
  }
  
  adonis_results %>% drop_na() %>% return()
}
ResultsAdonis = tibble()

not_to_test<- c("NG_ID", "NEXT_ID", "FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2","next_id_infant", "next_id_mother","Modified_NEXT_ID_without_preg_number" ,"days_from_first_collection"  ,
                "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                "exact_age_days_at_collection", "exact_age_months_at_collection",
                "Timepoint_categorical", "SAMPLE_ID",
                "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                "shannon", "BATCH_NUMBER", "next_id_partner", "sibling_number", "timepoint")



for (Timepoint in unique(metadata_infants$Timepoint_categorical)){
  print(Timepoint)
  phenotypes<-metadata_infants %>% filter(Timepoint_categorical == Timepoint )
  PP=drop_phenotypes(phenotypes, non_NA_values = 50, minimum_variance = 5)
  PP %>% select(-one_of(not_to_test)) %>% colnames() -> PhenosTest
  Call_adonis(phenotypes, Phenotype_list=PhenosTest ,perm=10000, Distance=distance) %>% mutate(Timepoint = Timepoint) %>% rbind(ResultsAdonis, .) -> ResultsAdonis
}


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint")
write.table(ResultsAdonis, "result_adonis_per_timepoint_infants_10000_perm.txt", sep = "\t", row.names = F)



