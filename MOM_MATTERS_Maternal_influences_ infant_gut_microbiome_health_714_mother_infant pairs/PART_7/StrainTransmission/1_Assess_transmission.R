.libPaths(c("/groups/umcg-llnext/tmp01/umcg-sandreusanchez/Strain_transmission/Packages", .libPaths()))

library(tidyverse)
library(patchwork)
library(ggpubr)
library(cutpointr)
library(lme4)



args <- commandArgs(trailingOnly = TRUE)
file = args[1]
print(paste0("Initializing: ", file))


file_new = str_remove(basename(file), "\\.txt") 
SGB_name = str_remove(file_new, "\\_DistMat") 

print(SGB_name)


print("Reading file and metadata")
strain = read.delim(file)
metadata<-read.delim("metadata/LLNEXT_metadata_03_01_2024.txt")
cross_phenos<-read.delim("metadata/masterfile_cross_sectional_2023_11_15.txt")

prioritized_phenos <- cross_phenos %>%
  select(next_id_infant,
        mother_birthcardself_gestational_age_weeks, 
         birth_deliverybirthcard_mode_binary,
         birth_deliverybirthcard_place_delivery_simple,
         birth_delivery_mode_simple,
         birth_delivery_mode_complex,
         mother_birthcard_parity,
         family_pets_any, 
         infant_birthcard_feeding_mode_after_birth) %>% as_tibble()
dynamic_phenos<-read.delim("metadata/masterfile_longitudinal_2023_09_29.txt")

dynamic_prioritized_phenos<-dynamic_phenos %>%
  select(next_id_infant,
         SAMPLE_ID,
         infant_ffq_feeding_mode_simple,
         infant_ffq_feeding_mode_complex) %>% as_tibble()

names (dynamic_prioritized_phenos)[2]<-"sample_id_infant"


thresholds_mireia = read_tsv("Thresholds_vallesColomere_Jan21.tsv")


Get_strain_sharing = function(sp3, SGB_name ,Plot = F){
  print("Running strain sharing function")	  
  # Normalize distances 
  sp3$nGD <- sp3$dist / (max(sp3$dist))
  
  nGD_training <- rbind(sp3 %>% filter(same_individual == "same_individual") %>% filter(time_day_diff <= 180) %>%
                          group_by(NEXT_ID_short_1) %>% arrange(NEXT_ID_short_1, time_day_diff) %>% slice_head(n = 1) %>% ungroup(),
                        sp3 %>% filter(same_individual != "same_individual",  related == "unrelated" ) %>%
                          group_by(NEXT_ID_short_1, NEXT_ID_short_2) %>% slice_head(n = 1) %>% ungroup())
  
  
  res_youden <- cutpointr(data = nGD_training, x = nGD, class = same_individual, method = maximize_metric, metric = youden, na.rm=T) # Youdens cut-off 
  quantile_5pc <- nGD_training %>% filter(related == "unrelated") %>% pull(nGD) %>% quantile(0.05) # 5th quantile of unrelated individual distribution
  quantile_5pc<-unname (quantile_5pc) # was a named number 
  number_same_individuals <- sum(nGD_training$same_individual == "same_individual", na.rm = TRUE) # Number of within sample pairs 
  quantile_3pc <- nGD_training %>% filter(related == "unrelated") %>% pull(nGD) %>% quantile(0.03)# 3rd percentile of the unrelated individual distribution
  quantile_3pc<-unname(quantile_3pc) # was a named number 
  
  # Set the ideal threshold 
  
  # Condition 1: If res_youden$optimal_cutpoint is lower than quantile_5pc
  if (res_youden$optimal_cutpoint < quantile_5pc) {
    final_threshold <- res_youden$optimal_cutpoint
  } else if (number_same_individuals < 50) {
    # Condition 2: If number_same_individuals is less than 50
    final_threshold <- quantile_3pc
  } else if (number_same_individuals >= 50 & quantile_5pc < res_youden$optimal_cutpoint) {
    # Condition 3: If number_same_individuals is >= 50 and quantile_5pc is less than res_youden$optimal_cutpoint
    final_threshold <- quantile_5pc
  } else {
    # A default value or an error handling could be placed here if none of the conditions are met.
    final_threshold <- NA # or some other default or error value
  }
  
  print (paste0("Threshold for SGB: ", final_threshold) )
  
  thresholds_mireia %>% filter(SGB==SGB_name) -> t_m
  if (dim(t_m)[1] != 1 ){ threshold_valles = NA 
  } else { quantile(nGD_training$nGD , t_m$used_nGD_score_percentile) -> threshold_valles }
  
  sp4<-sp3 %>% mutate(Strain_sharing = ifelse(nGD <= final_threshold , "yes", "no" ) ) 
  
  
  output_table <- tibble(file_name = SGB_name, final_threshold = final_threshold, number_same_individuals=number_same_individuals, Threshold_valles = threshold_valles)
  print(output_table)
  if (Plot == T){
    pdf(paste0("plots_thresholds/", SGB_name, ".pdf"))
    plot<-ggplot(data = nGD_training) +
      geom_density(aes(x = nGD, fill = same_individual), alpha = 0.5) +
      geom_vline(aes(xintercept = res_youden$optimal_cutpoint, color = "youden")) +
      geom_vline(aes(xintercept = quantile_5pc, color = "quantile"), linetype = "dotted", lwd = 1) +
      theme_minimal() + xlab("StrainPhlAn nGD") + ylab("") +
      ggtitle(paste("Distribution for", SGB_name)) +
      theme(legend.title = element_blank(), legend.position = "bottom") +
      scale_color_manual(name = "Statistics", values = c(youden = "blue", quantile = "red"))
      dev.off() 
  }
  
  return(list(output_table, sp4))
}

print("Processing distance matrix file")

strain[upper.tri(strain)] <- NA # Setting upper triangle of distance matrix to NA, helpful later to not have duplicates
strain$Sample_ID_1=row.names(strain)
  
# Change format of data frame to long 
sp1=strain %>% pivot_longer(!Sample_ID_1, names_to = "Sample_ID_2", values_to = "dist")
# Remove NA's (as we have the same pair distances twice A to B and B to A) (this was the upper half of the diagonal)
sp2=sp1 %>% filter(Sample_ID_1!=Sample_ID_2)
sp3=sp2 %>% drop_na()
  
# Make sure both SAMPLE1 and SAMPLE2 are 12 characters (currently long names )
sp3$Sample_ID_1 <- substr(sp3$Sample_ID_1, 1, 12)
sp3$Sample_ID_2 <- substr(sp3$Sample_ID_2, 1, 12)
  
print("Merging with metadata")
# Merging with the metdata file 
sp3 <- left_join(sp3 %>% select(everything()),
                   metadata %>% select(Sample_ID_1 = NG_ID,
                                 NEXT_ID_long_1 = NEXT_ID,
                                 NEXT_ID_short_1=Modified_NEXT_ID_without_preg_number,
                                 FAMILY_1 = FAMILY,
                                 days_from_first_collection_1 = days_from_first_collection,
                                 Timepoint_categorical_1 = Timepoint_categorical,
                                 exact_age_months_at_collection_1 = exact_age_months_at_collection,
                                 Type_1=Type))
  
sp3 <- left_join(sp3 %>% select(everything()),
                   metadata %>% select(Sample_ID_2 = NG_ID,
                                       NEXT_ID_2 = NEXT_ID,
                                       NEXT_ID_short_2=Modified_NEXT_ID_without_preg_number,
                                       FAMILY_2 = FAMILY,
                                       days_from_first_collection_2 = days_from_first_collection,
                                       Timepoint_categorical_2 = Timepoint_categorical,
                                       exact_age_months_at_collection_2 = exact_age_months_at_collection,
                                       Type_2=Type))

# Computing time difference between sample (important for longitudinal samples, because some samples have many years between them ie. mothers of first and second pregnancy)
sp3$time_day_diff <- abs(sp3$days_from_first_collection_1 - sp3$days_from_first_collection_2)

  
sp3= sp3 %>% filter(Sample_ID_1!=Sample_ID_2)
  
# Annotating pairs of samples. Are they related? Are they from the same individual?
sp3$same_individual <- ifelse( sp3$NEXT_ID_short_1 ==  sp3$NEXT_ID_short_2, "same_individual", "different_individual")
sp3$related <- ifelse( sp3$ FAMILY_1 ==  sp3$ FAMILY_2, "related", "unrelated")
 
################################################################
#1. Find tresholds in training data and apply to all the dataset
################################################################
  
Distance_output = paste0("Distance_tables/", SGB_name, ".rds")
Threshold_file = paste0("Distance_tables/", SGB_name, "_threshold.tsv")
  
#Find trehsold: Distances are normalized, youden is used
Results_f = Get_strain_sharing(sp3, SGB_name, Plot=T)
sp4 = Results_f[[2]]
output_table =  Results_f[[1]]
print("Saving output")
write_rds(sp4, Distance_output)
write_tsv(output_table, Threshold_file)
