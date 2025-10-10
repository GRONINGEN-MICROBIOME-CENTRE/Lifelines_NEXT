################### PERMANOVA TCAM INFANTS ##########################
# Author: T.Sinha 
# Last updated: 28th September, 2025

library(tidyverse)
library(vegan)
library(foreach)

# Load functions

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


Call_adonis = function(phenotypes, Phenotype_list, Distance, perm=10000, cores=8, covariates="base") {
  Distance = as.matrix(Distance)
  adonis_results <- data.frame(matrix(ncol = 6, nrow = length(Phenotype_list)))
  colnames(adonis_results) <- c("Phenotype", "Df", "F", "R2", "p-value", "N")
  Phenotype_list = Phenotype_list[!is.na(Phenotype_list)]
  
  # Base formula 
  base_formula <- "clean_reads_FQ_1 + dna_conc + BATCH_NUMBER"
  
  # Add covariates if required 
  if (covariates == "cor_delivery") {
    base_formula <- paste0(base_formula, " + birth_deliverybirthcard_mode_binary")
  } else if (covariates == "cor_delivery_feeding") {
    base_formula <- paste0(base_formula, " + birth_deliverybirthcard_mode_binary + infant_ffq_ever_never_breastfed")
  }
  
  adon <- foreach(i = 1:length(Phenotype_list), .combine = rbind) %do% {
    Pheno = Phenotype_list[i]
    print(Pheno)
    Values <- as.vector(as_vector(phenotypes[[Pheno]]))
    
    # Filter NAs
    r <- row.names(phenotypes)[!is.na(Values)]
    phenos2 <- phenotypes[r, ]
    
    # Match with Distance
    distmat_cleaned <- Distance[r, r]
    
    # Construct formula
    FR <- as.formula(paste0("distmat_cleaned ~ ", base_formula, " + ", Pheno))
    
    # Run adonis2
    ad1 <- adonis2(FR, phenos2, permutations = perm, parallel = cores, na.action = na.fail, by="margin")
    
    # Find the row corresponding to the phenotype term
    pheno_row <- which(rownames(ad1) == Pheno)
    residual_row <- which(rownames(ad1) == "Residual")
    
    if (length(pheno_row) == 0) {
      warning(paste("Phenotype", Pheno, "not found in adonis2 result."))
      adonis_results[i, ] <- c(Pheno, NA, NA, NA, NA)
    } else {
      adonis_results[i, ] <- c(Pheno, ad1$Df[pheno_row], ad1$F[pheno_row], ad1$R2[pheno_row], ad1$'Pr(>F)'[pheno_row], ad1$Df[residual_row])
    }
  }
  
  return(adonis_results %>% drop_na())
}


##############################
# Infant 
##############################
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/phenotypes")
metadata_infants<-read.delim("infant_phenotypes_all_microbiome_fil.txt")
metadata_infants <- metadata_infants %>%
  filter(!is.na(birth_deliverybirthcard_mode_binary) & !is.na(infant_ffq_ever_never_breastfed))
metadata_infants<-metadata_infants[,c(1:131)]
metadata_infants <- metadata_infants %>%
  distinct(next_id_infant, .keep_all = TRUE)
row.names(metadata_infants)<-metadata_infants$next_id_infant
metadata_infants<-drop_phenotypes(metadata_infants, non_NA_values = 50, minimum_variance = 5)


TCAM<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint/TCAM_infant_M1_M3_M6_M12.txt", sep = ",")
row.names(TCAM)<-TCAM$X
TCAM$X=NULL

infant_TCAM<-TCAM[row.names(TCAM)%in% rownames(metadata_infants),] 
metadata_infants<-metadata_infants[row.names(metadata_infants)%in% rownames(TCAM),] 
metadata_infants <- metadata_infants[match(rownames(infant_TCAM), rownames(metadata_infants)), ]


distance <- vegdist(infant_TCAM, method = "euclidean" )

not_to_test<- c("NG_ID", "NEXT_ID", "FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2","next_id_infant", "next_id_mother","Modified_NEXT_ID_without_preg_number" ,"days_from_first_collection"  ,
                "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                "exact_age_days_at_collection", "exact_age_months_at_collection","exact_age_years_at_collection",
                "Timepoint_categorical", "SAMPLE_ID",
                "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                "shannon", "BATCH_NUMBER", "next_id_partner", "sibling_number", "timepoint")

# Running base model 

PhenosTest <- metadata_infants %>% select(-one_of(not_to_test)) %>% colnames()

# Run base model without correction 
ResultsAdonis_base <- Call_adonis(metadata_infants, 
                             Phenotype_list = PhenosTest,
                             perm = 10000,
                             Distance = distance,
                             covariates = "base")

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint")
write.table(ResultsAdonis_base, "result_adonis_TCAM_INFANTS_BY_MARGIN_10000_perm_no_cor_28_09_2025.txt", sep = "\t", row.names = F)

# Running model corrected for delivery mode 
ResultsAdonis_cor_delivery <- Call_adonis(metadata_infants, 
                             Phenotype_list = PhenosTest,
                             perm = 10000,
                             Distance = distance,
                             covariates = "cor_delivery")
write.table(ResultsAdonis_cor_delivery, "result_adonis_TCAM_INFANTS_BY_MARGIN_10000_perm_cor_delivery_28_09_2025.txt", sep = "\t", row.names = F)


# Running model corrected for delivery mode & feeding 
ResultsAdonis_cor_delivery_feeding<- Call_adonis(metadata_infants, 
                             Phenotype_list = PhenosTest,
                             perm = 10000,
                             Distance = distance,
                             covariates = "cor_delivery_feeding")
write.table(ResultsAdonis_cor_delivery_feeding, "result_adonis_TCAM_INFANTS_BY_MARGIN_10000_perm_cor_delivery_feeding_28_09_2025.txt", sep = "\t", row.names = F)







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
  PP=drop_phenotypes(phenotypes, non_NA_values = 100, minimum_variance = 5)
  PP %>% select(-one_of(not_to_test)) %>% colnames() -> PhenosTest
  #input selected phenotypes, as a vector, to Phenotype_list
  Call_adonis(phenotypes, Phenotype_list=PhenosTest ,perm=10000, Distance=distance) %>% mutate(Timepoint = Timepoint) %>% rbind(ResultsAdonis, .) -> ResultsAdonis
}


ResultsAdonis$Timepoint <- factor(ResultsAdonis$Timepoint, levels = c("P12", "P28", "B", "M3"))

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint")
write.table(ResultsAdonis, "result_adonis_per_timepoint_mothers_10000_perm.txt", sep = "\t", row.names = F)



