################################################################################
##### TCAM analysis
### Author(s): Nataliia Kuzub, Asier Fernández-Pato, Trishla Sinha
### Last updated: 30th August, 2024
################################################################################

################################################################################
# Loading libraries
################################################################################
library(vegan)
library(dplyr)
library(ggplot2)
library(readr)
library(ggExtra)

# Set working directory 
setwd("Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint/")

################################################################################
# Loading functions
################################################################################

do_preprocessing_infants <- function(TCAM, cross_phenotypes, cross_selection, timepoints){
  # TCAM table preprocesing: output of the TCAM.ipynb
  rownames(TCAM) <- TCAM$X
  TCAM <- subset(TCAM, select = -c(X))
  
  # Fitering out the samples excluded for this analysis
  cross_phenotypes <- cross_phenotypes[cross_phenotypes$next_id_infant %in% row.names(TCAM), ]
  cross_phenotypes$mother_birthcard_parity <- ifelse(is.na(cross_phenotypes$mother_birthcard_parity), NA, 
                                                     ifelse(cross_phenotypes$mother_birthcard_parity >= 2, 2, cross_phenotypes$mother_birthcard_parity))
  row.names(cross_phenotypes) <- cross_phenotypes$next_id_infant
  
  # Selecting data relevant for infant 
  infant_cross_selection <- cross_selection[cross_selection$TCAM_infants == 1 & !is.na(cross_selection$TCAM_infants),]
  column_names <- infant_cross_selection[[1]]
  infant_cross_phenotypes <- cross_phenotypes %>% select(all_of(column_names))
  
  # infant_cross_phenotypes_filtered <-drop_phenotypes(infant_cross_phenotypes, non_NA_values, minimum_variance)
  
  # Preparing phenotypes for associations
  to_exclude <- c("FAMILY", "next_id_mother", "next_id_infant", "next_id_partner", "infant_relations", "sibling_number", "twin_pair")
  infant_cross_phenotypes_filtered <- infant_cross_phenotypes[,!(colnames(infant_cross_phenotypes) %in% to_exclude)]
  
  
  return(list(TCAM, infant_cross_phenotypes_filtered))
}



do_preprocessing_mothers <- function(TCAM, cross_phenotypes, cross_selection, timepoints){
  # TCAM table preprocesing: output of the TCAM.ipynb
  rownames(TCAM) <- TCAM$X
  TCAM <- subset(TCAM, select = -c(X))
  
  # Fitering out the samples excluded for this analysis
  cross_phenotypes <- cross_phenotypes[cross_phenotypes$next_id_mother %in% row.names(TCAM), ]
  cross_phenotypes$mother_birthcard_parity <- ifelse(is.na(cross_phenotypes$mother_birthcard_parity), NA, 
                                                     ifelse(cross_phenotypes$mother_birthcard_parity >= 2, 2, cross_phenotypes$mother_birthcard_parity))
  row.names(cross_phenotypes) <- cross_phenotypes$next_id_mother
  
  # Selecting data relevant for infant 
  mother_cross_selection <- cross_selection[cross_selection$TCAM_mothers == 1 & !is.na(cross_selection$TCAM_mothers),]
  column_names <- mother_cross_selection[[1]]
  mother_cross_phenotypes <- cross_phenotypes %>% select(all_of(column_names))
  
  
  # Preparing phenotypes for associations
  to_exclude <- c("FAMILY", "next_id_mother", "next_id_infant", "next_id_partner", "infant_relations", "sibling_number", "twin_pair")
  mother_cross_phenotypes_filtered <- mother_cross_phenotypes[,!(colnames(mother_cross_phenotypes) %in% to_exclude)]
  
  
  return(list(TCAM, mother_cross_phenotypes_filtered))
}


do_permanova_simple <- function(processed_data, permutation_num){
  # Generate the Aitchison distance matrices
  Distance <- vegdist(processed_data[[1]], method = "euclidean" )
  Distance <- data.frame(as.matrix((Distance)))
  
  # Run PERMANOVA to estimate the effect of phenotypes on overall bacterial composition
  adonis_results_raw <- data.frame(Df = numeric(ncol(processed_data[[2]])),
                                   F = numeric(ncol(processed_data[[2]])),
                                   R2 = numeric(ncol(processed_data[[2]])),
                                   p_value = numeric(ncol(processed_data[[2]])))
  
  rownames(adonis_results_raw) <- colnames(processed_data[[2]])
  
  for (i in 1:ncol(processed_data[[2]])) {
    # Subset distance matrix and phenotypes, removing rows with NAs
    na_rows <- which(is.na(processed_data[[2]][, i]))
    distmat_cleaned <- Distance[!(rownames(Distance) %in% row.names(processed_data[[2]])[na_rows]),
                                !(colnames(Distance) %in% row.names(processed_data[[2]])[na_rows])]
    phenos2 <- processed_data[[2]][row.names(processed_data[[2]]) %in% row.names(distmat_cleaned), ]
    phenos2 <- phenos2[ match(rownames(distmat_cleaned), rownames(phenos2) ), ]
    # Run adonis2
    ad1 <- adonis2(distmat_cleaned ~ phenos2[[i]],
                   permutations = permutation_num,
                   parallel = 8,
                   na.action = na.fail,
                   by = "terms")
    # Save results
    adonis_results_raw[i, ] <- c(ad1$Df[1], ad1$F[1], ad1$R2[1], ad1$"Pr(>F)"[1])
  }
  
  return(adonis_results_raw)
}


################################################################################
# Loading raw data
################################################################################

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/phenotypes")
cross_selection_corrected <- read.delim("masterfile_cross_sectional_selection_AFTER_correlation_&_correction_10_07_2024.txt")
cross_selection_corrected[cross_selection_corrected$variable_name == "birth_delivery_mode_complex", "TCAM"] <- 0
#cross_selection_corrected[cross_selection_corrected$variable_name == "birth_delivery_mode_simple", "TCAM"] <- 0

cross_phenotypes <- read.delim("masterfile_cross_sectional_2023_11_15.txt")
cross_selection_mothers <- read.delim("masterfile_cross_sectional_selection_AFTER_correlation_&_correction_10_07_2024.txt")


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint")

TCAM_infant_M1_M3_M6_M12 <- read.delim("TCAM_infant_M1_M3_M6_M12.txt", sep = ",")
TCAM_mother_P12_B_M3 <- read.delim("TCAM_mother_P12_B_M3.txt", sep = ",")


################################################################################
# Preprocessing raw data
################################################################################from the selection after the correlation


processed_infant_M1_M3_M6_M12 <- do_preprocessing_infants(TCAM_infant_M1_M3_M6_M12, cross_phenotypes, cross_selection_corrected, c('M1', 'M3', 'M6', 'M12'))

processed_mother_P12_B_M3 <- do_preprocessing_mothers(TCAM_mother_P12_B_M3, cross_phenotypes, cross_selection_mothers, c('P12', 'B', 'M3'))

################################################################################
# PERMANOVA analysis
################################################################################

permanova_simple_infant_M1_M3_M6_M12 <- do_permanova_simple(processed_infant_M1_M3_M6_M12, 10000)
permanova_simple_infant_M1_M3_M6_M12$Phenotype <-row.names(permanova_simple_infant_M1_M3_M6_M12)
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint")
write.table(permanova_simple_infant_M1_M3_M6_M12, "result_adonis_tcam_infants_10000_perm.txt", sep = "\t", row.names = F)


permanova_simple_mother_P12_B_M3 <- do_permanova_simple(processed_mother_P12_B_M3, 10000)
permanova_simple_mother_P12_B_M3$Phenotype<-row.names(permanova_simple_mother_P12_B_M3)
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint")
write.table(permanova_simple_mother_P12_B_M3, "result_adonis_tcam_mothers_10000_perm.txt", sep = "\t", row.names = F)

