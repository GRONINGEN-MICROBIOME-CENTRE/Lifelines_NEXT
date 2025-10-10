############################ STRAIN TRANSMISSION ###########
# Authors: Sergio Andreu Sanchez, Trishla Sinha, Alex Kurilshikov
# Last update: 26st of Nov, 2024 
# This script to strain transmission analysis from StrainPhlAn data

library(tidyverse)
library(patchwork)
library(ggpubr)
library(cutpointr)
library(lme4)
library(reshape2)


# 1.1 Strain transmission evaluation --------------------------------------
# adopted from https://github.com/GRONINGEN-MICROBIOME-CENTRE/Gut_Microbiome/strain_analysis/StrainTransmission/1_Assess_transmission.R

## Functions 

Get_strain_sharing = function(sp3, Plot = T){
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
  
  ## This is the comparison with estimations from https://www.nature.com/articles/s41586-022-05620-1
  ## Skipped in the demo code. 
  #thresholds_mireia %>% filter(SGB==SGB_name) -> t_m
  #if (dim(t_m)[1] != 1 ){ threshold_valles = NA 
  #} else { quantile(nGD_training$nGD , t_m$used_nGD_score_percentile) -> threshold_valles }
  
  sp4<-sp3 %>% mutate(Strain_sharing = ifelse(nGD <= final_threshold , "yes", "no" ) ) 
  
  
  output_table <- tibble(file_name = "random_SGB", final_threshold = final_threshold, number_same_individuals=number_same_individuals)#, Threshold_valles = threshold_valles)
  print(output_table)
  if (Plot == T){
    #pdf(paste0("plots_thresholds/", SGB_name, ".pdf"))
    plot<-ggplot(data = nGD_training) +
      geom_density(aes(x = nGD, fill = same_individual), alpha = 0.5) +
      geom_vline(aes(xintercept = res_youden$optimal_cutpoint, color = "youden")) +
      geom_vline(aes(xintercept = quantile_5pc, color = "quantile"), linetype = "dotted", lwd = 1) +
      theme_minimal() + xlab("StrainPhlAn nGD") + ylab("") +
      ggtitle(paste("Distribution for", "random_SGB")) +
      theme(legend.title = element_blank(), legend.position = "bottom") +
      scale_color_manual(name = "Statistics", values = c(youden = "blue", quantile = "red"))
    #dev.off() 
  }
  print(plot)
  return(list(output_table, sp4))
}


## Preparing input table  

metadata = read.table("simulated_metadata.txt")
tree = read.tree("RAxMLtree_example.tre")
dist1 = as.matrix(cophenetic.phylo(tree))
diag(dist1) = NA
dist1[upper.tri(dist1)] = NA
dist_long = melt(dist1)
dist_long = dist_long[!is.na(dist_long$value),]

colnames(dist_long)[3] = "dist"

colnames(dist_long)[1:2] = c("NEXT_ID_short_1","NEXT_ID_short_2")
dist_long$Type_1 = metadata[dist_long$NEXT_ID_short_1,"Type"]
dist_long$Type_2 = metadata[dist_long$NEXT_ID_short_2,"Type"]
dist_long$Family_1 = metadata[dist_long$NEXT_ID_short_1,"Family"]
dist_long$Family_2 = metadata[dist_long$NEXT_ID_short_2,"Family"]
dist_long$Timepoint_1 = metadata[dist_long$NEXT_ID_short_1,"Timepoint"]
dist_long$Timepoint_2 = metadata[dist_long$NEXT_ID_short_2,"Timepoint"]

dist_long$Time_1 = sapply(dist_long$Timepoint_1, \(x) switch(x,W2 = 14,M1 = 30,M3 = 90,M6 = 185, M12 = 360))
dist_long$Time_2 = sapply(dist_long$Timepoint_2, \(x) switch(x,W2 = 14,M1 = 30,M3 = 90,M6 = 185, M12 = 360))
dist_long$time_day_diff = abs(dist_long$Time_1 - dist_long$Time_2)
dist_long$same_individual = ifelse(
  (dist_long$Family_1 == dist_long$Family_2) &
    (dist_long$Type_1 == dist_long$Type_2) , "same_individual","different_individual")

dist_long$related = ifelse(
  dist_long$NEXT_ID_short_1 == dist_long$NEXT_ID_short_2 , "related","unrelated")

## Running function

result_strain = Get_strain_sharing(dist_long,Plot = T)

# 1.2 Strain transmission association to phenotypes ---------------------------
# adopted from https://github.com/GRONINGEN-MICROBIOME-CENTRE/Lifelines_NEXT/Gut_Microbiome/strain_analysis/StrainTransmission/2_MotherInfantTranmission_Association.R


