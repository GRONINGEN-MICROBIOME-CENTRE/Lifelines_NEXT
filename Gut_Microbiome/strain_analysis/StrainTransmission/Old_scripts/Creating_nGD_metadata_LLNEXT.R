####### Creating a file with normalized phylogenetic distances along with relevant metaddata ######## 
### Author: Trishla Sinha
### Made: 3rd January, 2023
### Last updated: 21st January, 2023
### Input data: Distance matrices of each SGB (only those present in atleast 100 samples in LLNEXT)
### Output data: 
#1)Table for each each SGG with pair-wise distances (of each infant/mother at a particular timepoint with another infant/mother at a particular timepoint) with annotation whether they share the same strain or not
#2) Output table with the name of SGB with its final chose threshold with the number of same-individual pairs the SGB was found in 

#rm(list=ls())
vars_to_remove <- setdiff(ls(), "strain")
rm(list = vars_to_remove)

library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(ggpubr)
library(cutpointr)
library(lme4)


# Cleaning all the linkage files and deriving family ID's 

setwd("~/Desktop/LLNEXT/Analysis/strainphlan_NEXT/distance_matrices/")
folder = "~/Desktop/LLNEXT/Analysis/strainphlan_NEXT/distance_matrices/" 
file_list = list.files(path=folder,pattern="*.txt") # Names of all the species you have distance matrices for 
file_new<-str_remove(file_list, "\\.txt") 
file_short<-str_remove(file_new, "\\_DistMat") 
#strain=lapply(file_list, function(x)read.delim(x))
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_03_01_2024.txt")
cross_phenos<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_2023_11_15.txt")
prioritized_phenos <- cross_phenos %>%
  select(next_id_infant,
        mother_birthcardself_gestational_age_weeks, 
         birth_deliverybirthcard_mode_binary,
         birth_deliverybirthcard_place_delivery_simple,
         birth_delivery_mode_simple,
         birth_delivery_mode_complex,
         mother_birthcard_parity,
         family_pets_any, 
         infant_birthcard_feeding_mode_after_birth)
dynamic_phenos<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_longitudinal_2023_09_29.txt")

dynamic_prioritized_phenos<-dynamic_phenos %>%
  select(next_id_infant,
         SAMPLE_ID,
         infant_ffq_feeding_mode_simple,
         infant_ffq_feeding_mode_complex)

names (dynamic_prioritized_phenos)[2]<-"sample_id_infant"

output_table <- data.frame(file_name = character(), final_threshold = numeric(), number_same_individuals=numeric (), stringsAsFactors = FALSE) # To carry the final thresholds for each SGB
#i=113 Largest strain, know it works 

Result_all_SGB = tibble()
for(i in 328:length(strain)){  
  sp=strain[[i]]
  sp[upper.tri(sp)] <- NA # Setting upper triangle of distance matrix to NA, helpful later to not have duplicates
  sp$Sample_ID_1=row.names(sp)
  
  # Change format of data frame to long 
  sp1=sp %>%
    pivot_longer(!Sample_ID_1, names_to = "Sample_ID_2", values_to = "dist")

  
  # Remove NA's (as we have the same pair distances twice A to B and B to A) (this was the upper half of the diagonal)
  sp2=sp1 %>%
    filter(Sample_ID_1!=Sample_ID_2)
  
  sp3=sp2 %>%
    drop_na()
  
  # Make sure both SAMPLE1 and SAMPLE2 are 12 characters (currently long names )
  sp3$Sample_ID_1 <- substr(sp3$Sample_ID_1, 1, 12)
  sp3$Sample_ID_2 <- substr(sp3$Sample_ID_2, 1, 12)
  
  # Merging with the metdata file 
  
  sp3 <- left_join(sp3 %>% select(everything()),
                   metadata %>% select(Sample_ID_1 = NG_ID,
                                 NEXT_ID_long_1 = NEXT_ID,
                                 NEXT_ID_short_1=Modified_NEXT_ID_without_preg_number,
                                 FAMILY_1 = FAMILY,
                                 days_from_first_collection_1 = days_from_first_collection,
                                 Timepoint_categorical_1 = Timepoint_categorical,
                                 Type_1=Type))
  
  sp3 <- left_join(sp3 %>% select(everything()),
                   metadata %>% select(Sample_ID_2 = NG_ID,
                                       NEXT_ID_2 = NEXT_ID,
                                       NEXT_ID_short_2=Modified_NEXT_ID_without_preg_number,
                                       FAMILY_2 = FAMILY,
                                       days_from_first_collection_2 = days_from_first_collection,
                                       Timepoint_categorical_2 = Timepoint_categorical,
                                       Type_2=Type))

  # Computing time difference between sample (important for longitudinal samples, because some samples have many years between them ie. mothers of first and second pregnancy)
  sp3$time_day_diff <- abs(sp3$days_from_first_collection_1 - sp3$days_from_first_collection_2)
  
  
  sp3=sp3 %>%
    filter(Sample_ID_1!=Sample_ID_2)
  
  # Annotating pairs of samples. Are they related? Are they from the same individual?
  sp3$same_individual <- ifelse( sp3$NEXT_ID_short_1 ==  sp3$NEXT_ID_short_2, "same_individual", "different_individual")
  sp3$related <- ifelse( sp3$ FAMILY_1 ==  sp3$ FAMILY_2, "related", "unrelated")
 
  # Normalize distances 
 sp3$nGD <- sp3$dist / (max(sp3$dist))
 
 nGD_training <- rbind(sp3 %>% filter(same_individual == "same_individual") %>%
                         filter(time_day_diff <= 180) %>%
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

 print (final_threshold)
 
sp4<-sp3 %>% mutate(Strain_sharing = ifelse(nGD <= final_threshold , "yes", "no" ) ) 

#file_output_location <- paste0("~/Desktop/LLNEXT/Analysis/strainphlan_NEXT/distance_matrices/strain_sharing_all_distances/", file_short[i], ".txt")
#write.table(sp4, file = file_output_location, row.names = FALSE, sep = "\t", quote = FALSE)

output_table <- rbind(output_table, data.frame(file_name = file_short[i], final_threshold = final_threshold, number_same_individuals=number_same_individuals))
 
 pdf(paste0(folder, "plots_thresholds/", file_short[i], ".pdf"))
 
 
 plot<-ggplot(data = nGD_training) +
   geom_density(aes(x = nGD, fill = same_individual), alpha = 0.5) +
   geom_vline(aes(xintercept = res_youden$optimal_cutpoint, color = "youden")) +
   geom_vline(aes(xintercept = quantile_5pc, color = "quantile"), linetype = "dotted", lwd = 1) +
   theme_minimal() + xlab("StrainPhlAn nGD") + ylab("") +
   ggtitle(paste("Distribution for", file_short[i])) +
   theme(legend.title = element_blank(), legend.position = "bottom") +
   scale_color_manual(name = "Statistics", values = c(youden = "blue", quantile = "red"))
 print (plot)
 # Try here with histograms 
 dev.off() 
 
 sp5<-sp4 %>% filter(related == "related")
 
 # First working with mother-infant pairs to make conclusions about strain transmission over time 
 mother_infant_pairs <- sp5 %>% 
   filter((Type_1 == "mother" & Type_2 == "infant") | (Type_1 == "infant" & Type_2 == "mother"))
 #Just continue if there are enough mother_infant pairs
 if (dim(mother_infant_pairs)[1] <= 2) { next }
 
 
 
 # Define whether mother's timepoint is pre or post-birth
 mother_infant_pairs <- mother_infant_pairs %>%
   mutate(mother_timeframe = case_when(
     (Type_1 == "mother" & (Timepoint_categorical_1 %in% c("P12", "P28", "B"))) | 
       (Type_2 == "mother" & (Timepoint_categorical_2 %in% c("P12", "P28", "B"))) ~ "pre_birth",
     (Timepoint_categorical_1 %in% c("M1", "M2", "M3")) | 
       (Timepoint_categorical_2 %in% c("M1", "M2", "M3")) ~ "post_birth",
     TRUE ~ NA_character_  
   ))
 
 # Define exact timepoint categorical of infant 
 mother_infant_pairs$infant_timepoint_categorical <- ifelse(
   mother_infant_pairs$Type_1 == "infant", 
   mother_infant_pairs$Timepoint_categorical_1, 
   ifelse(
     mother_infant_pairs$Type_2 == "infant", 
     mother_infant_pairs$Timepoint_categorical_2, 
     NA  
   )
 )
 
 # Define exact timepoint categorical of mother 
 mother_infant_pairs$mother_timepoint_categorical <- ifelse(
   mother_infant_pairs$Type_1 == "mother", 
   mother_infant_pairs$Timepoint_categorical_1, 
   ifelse(
     mother_infant_pairs$Type_2 == "mother", 
     mother_infant_pairs$Timepoint_categorical_2, 
     NA  
   )
 )
 
 # Define early or late timepoint  of infant 
 mother_infant_pairs$infant_early_late_timepoint <- ifelse(
   mother_infant_pairs$infant_timepoint_categorical %in% c("W2", "M1", "M2", "M3"),
   "early",
   ifelse(
     mother_infant_pairs$infant_timepoint_categorical %in% c("M6", "M9", "M12"),
     "late",
     NA  
   )
 )
 
 timepoint_order <- c("W2", "M1", "M2", "M3", "M6", "M9", "M12")
 
 mother_infant_pairs$infant_timepoint_categorical <- factor(
   mother_infant_pairs$infant_timepoint_categorical,
   levels = timepoint_order 
 )
 
 timepoint_order_mothers<-c("P12", "P28", "B", "M1", "M2", "M3")
 mother_infant_pairs$mother_timepoint_categorical <- factor(
   mother_infant_pairs$mother_timepoint_categorical,
   levels = timepoint_order_mothers 
 )
 
 # For strain sharing with mother and different infant timepoints take into account only the birth timepoint and if not available the P28 or the M3 timepoint
 # Check this 
 subset_mother_infant_pairs <- mother_infant_pairs %>%
   filter(Type_1 == "mother" | Type_2 == "mother") %>%
   group_by(FAMILY_1) %>%
   mutate(
     selected_timepoint = ifelse("B" %in% mother_timepoint_categorical, "B",
                                 ifelse("P28" %in% mother_timepoint_categorical, "P28",
                                        ifelse("M3" %in% mother_timepoint_categorical, "M3", NA_character_)))
   ) %>%
   ungroup()%>%
   filter(!is.na(selected_timepoint)) %>% filter( as.character(selected_timepoint) == mother_timepoint_categorical)
 
 # Merging with phenotypes (first pulling out the infant ID)
 
subset_mother_infant_pairs$next_id_infant <- ifelse(
subset_mother_infant_pairs$Type_1 == "infant",
subset_mother_infant_pairs$NEXT_ID_short_1,
   ifelse(
subset_mother_infant_pairs$Type_2 == "infant",
subset_mother_infant_pairs$NEXT_ID_short_2,
     NA
   )
 )
 

subset_mother_infant_pairs<-left_join(subset_mother_infant_pairs, prioritized_phenos)
subset_mother_infant_pairs$sample_id_infant<-paste0(subset_mother_infant_pairs$next_id_infant,"_", subset_mother_infant_pairs$infant_timepoint_categorical)
subset_mother_infant_pairs<-left_join(subset_mother_infant_pairs, dynamic_prioritized_phenos)
subset_mother_infant_pairs$infant_birthcard_feeding_mode_after_birth <- factor(
  subset_mother_infant_pairs$infant_birthcard_feeding_mode_after_birth,
  levels = c("BF", "MF", "FF")
)
subset_mother_infant_pairs$infant_ffq_feeding_mode_complex <- factor(
  subset_mother_infant_pairs$infant_ffq_feeding_mode_complex,
  levels = c("BF", "MF", "FF")
)
# Visualize strain sharing over time 
pdf(paste0(folder, "plots_strain_sharing_timepoint/", file_short[i], ".pdf"))
dynamic_title_timepoint <- paste("Percentage of", file_short[i], "Sharing Between Mothers Over Time")
data_summary_timepoint <- subset_mother_infant_pairs %>%
  group_by(infant_timepoint_categorical) %>%
  summarise(count = n(),
            yes_count = sum(Strain_sharing == "yes"),
            percentage = (yes_count / count) * 100) %>%
  ungroup()
strain_sharing_mother_infant_timepoint <- ggplot(data_summary_timepoint, aes(x = infant_timepoint_categorical, y = percentage)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = paste("n=", count, sep=""), y = percentage), vjust = -0.5, position = position_dodge(0.9)) +
  labs(title = dynamic_title_timepoint,
       x = "Timepoint",
       y = "% Strain Sharing between mother & infant") +
  theme_minimal()
print (strain_sharing_mother_infant_timepoint)
dev.off() 

# Visualize strain sharing with mode delivery and strain sharing 
pdf(paste0(folder, "plots_strain_sharing_mode_delivery/", file_short[i], ".pdf"))
dynamic_title_mode_delivery <- paste("Percentage of", file_short[i], "Sharing Between Mothers by Timepoint and Mode of Delivery")

data_summary_mode_delivery <- subset_mother_infant_pairs %>%
  filter(!is.na(birth_delivery_mode_simple)) %>%
  group_by(infant_timepoint_categorical, birth_delivery_mode_simple) %>%
  summarise(count = n(),
            yes_count = sum(Strain_sharing == "yes"),
            percentage = (yes_count / count) * 100) %>%
  ungroup()

strain_sharing_mother_infant_plot_mode_delivery<-ggplot(data_summary_mode_delivery, aes(x = infant_timepoint_categorical, y = percentage, fill = birth_delivery_mode_simple)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = paste("n=", count, sep=""), y = percentage), vjust = -0.5, position = position_dodge(0.9)) +
  labs(title = dynamic_title_mode_delivery,
       x = "Timepoint",
       y = "% Strain Sharing between mother & infant") +
  theme_minimal()
print ( strain_sharing_mother_infant_plot_mode_delivery)

dev.off() 

# Visualize strain sharing with feeding after birth and strain sharing 
pdf(paste0(folder, "plots_strain_sharing_feeding_after_birth/", file_short[i], ".pdf"))
dynamic_title_feeding_mode <- paste("Percentage of", file_short[i], "Sharing Between Mothers by Timepoint and Feeding Mode After Birth")
data_summary_feeding_mode <- subset_mother_infant_pairs %>%
  filter(!is.na(infant_birthcard_feeding_mode_after_birth)) %>%
  group_by(infant_timepoint_categorical, infant_birthcard_feeding_mode_after_birth) %>%
  summarise(count = n(),
            yes_count = sum(Strain_sharing == "yes"),
            percentage = (yes_count / count) * 100) %>%
  ungroup()

strain_sharing_mother_infant_plot_feeding_mode <- ggplot(data_summary_feeding_mode, aes(x = infant_timepoint_categorical, y = percentage, fill = infant_birthcard_feeding_mode_after_birth)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = paste("n=", count, sep=""), y = percentage), vjust = -0.5, position = position_dodge(0.9)) +
  labs(title = dynamic_title_feeding_mode,
       x = "Timepoint",
       y = "% Strain Sharing between mother & infant") +
  theme_minimal()
print(strain_sharing_mother_infant_plot_feeding_mode)
dev.off() 


pdf(paste0(folder, "plots_strain_sharing_place_delivery/", file_short[i], ".pdf"))

data_summary_place_delivery <- subset_mother_infant_pairs %>%
  filter(birth_delivery_mode_simple == "VG", !is.na(birth_deliverybirthcard_place_delivery_simple)) %>%
  group_by(infant_timepoint_categorical, birth_deliverybirthcard_place_delivery_simple) %>%
  summarise(count = n(),
            yes_count = sum(Strain_sharing == "yes"),
            percentage = (yes_count / count) * 100) %>%
  ungroup()

 
 dynamic_title_place_delivery <- paste("Percentage of", file_short[i], "Sharing Between Mothers by Timepoint and Place of Delivery")
 
 strain_sharing_mother_infant_plot_place_delivery<- ggplot(data_summary_place_delivery, aes(x = infant_timepoint_categorical, y = percentage, fill = birth_deliverybirthcard_place_delivery_simple)) +
   geom_bar(stat = "identity", position = position_dodge()) +
   geom_text(aes(label = paste("n=", count, sep=""), y = percentage), vjust = -0.5, position = position_dodge(0.9)) +
   labs(title = dynamic_title_place_delivery,
        x = "Timepoint",
        y = "% Strain Sharing between mother & infant") +
   theme_minimal()

 
 print ( strain_sharing_mother_infant_plot_place_delivery)
 
 dev.off() 
 
 #Now for Feeding mode complex
 
 pdf(paste0(folder, "plots_strain_sharing_infant_ffq_feeding_mode_complex/", file_short[i], ".pdf"))
 
 data_summary_feeding_mode_complex <- subset_mother_infant_pairs %>%
   filter( !is.na(infant_ffq_feeding_mode_complex)) %>%
   group_by(infant_timepoint_categorical, infant_ffq_feeding_mode_complex) %>%
   summarise(count = n(),
             yes_count = sum(Strain_sharing == "yes"),
             percentage = (yes_count / count) * 100) %>%
   ungroup()
 

 dynamic_title_feeding_mode_complex <- paste("Percentage of", file_short[i], "Sharing Between Mothers by Timepoint and Complex Feeding Mode")
 
 strain_sharing_mother_infant_plot_feeding_mode_complex <- ggplot(data_summary_feeding_mode_complex, aes(x = infant_timepoint_categorical, y = percentage, fill = infant_ffq_feeding_mode_complex)) +
   geom_bar(stat = "identity", position = position_dodge()) +
   geom_text(aes(label = paste("n=", count, sep=""), y = percentage), vjust = -0.5, position = position_dodge(0.9)) +
   labs(title = dynamic_title_feeding_mode_complex,
        x = "Timepoint",
        y = "% Strain Sharing between mother & infant") +
   theme_minimal()
 

 print(strain_sharing_mother_infant_plot_feeding_mode_complex)

 dev.off()
 
 #associatiate phenotype with strain share yes/no. Each row is a mother-infant distance. Noa ccounting for siblings
 #dataframe for association: subset_mother_infant_pairs
 subset_mother_infant_pairs %>% mutate(Strain_sharing = as.factor(Strain_sharing)) -> subset_mother_infant_pairs
 Results_association_strainSharing = tibble()
 #1. loop through phenotypes if at least 10 mother-infant pairs with at least 3 timepoints.
 #If there is only one timepoint per mother, then if there are >= 3 entries per baby, it means it has at least three available timepoints
 subset_mother_infant_pairs %>% group_by(next_id_infant) %>% summarise(N = n()) %>% filter(N>=3) -> Babies_with_3
 N_babies_with_3 = dim(Babies_with_3)[1] 
 if (N_babies_with_3 > 10){
   Pheno_list = c("infant_timepoint_categorical",
   "mother_birthcardself_gestational_age_weeks","birth_deliverybirthcard_mode_binary",          
   "birth_deliverybirthcard_place_delivery_simple", "birth_delivery_mode_simple",                   
   "birth_delivery_mode_complex", "mother_birthcard_parity",                      
   "family_pets_any","infant_birthcard_feeding_mode_after_birth", "infant_ffq_feeding_mode_simple",               
   "infant_ffq_feeding_mode_complex")     
   for ( Pheno in Pheno_list){
     subset_mother_infant_pairs_it = subset_mother_infant_pairs
     subset_mother_infant_pairs_it %>% filter(!is.na(!!sym(Pheno))) -> subset_mother_infant_pairs_it
     print(paste0("Iterating through pheno: ", Pheno))
     #2. build the model
     if (Pheno != "infant_timepoint_categorical"){
      Model0 =  as.formula(paste0( "Strain_sharing ~ (1|next_id_infant) + infant_timepoint_categorical" ))
      Model = as.formula(paste0( "Strain_sharing ~ (1|next_id_infant) + infant_timepoint_categorical + ", Pheno ))
     } else { 
       Model0 = as.formula(paste0( "Strain_sharing ~ (1|next_id_infant)"))
       Model = as.formula(paste0( "Strain_sharing ~ (1|next_id_infant) +", Pheno)) 
      }
     #glmer(formula = Model, family = binomial(link = "logit"), data = subset_mother_infant_pairs_it) -> Result_model
     Result_model <- tryCatch({
       glmer(formula = Model, family = binomial(link = "logit"), data = subset_mother_infant_pairs_it)
     }, error = function(e) {
       print(paste0("Model fitting failed for phenotype: ", Pheno, " with error: ", e$message))
       return(NA)
     }) 
       
     if (class(Result_model) == "logical") {
       next
     }
     
     as.data.frame(summary(Result_model)$coefficients) -> Result_coefficients
     rownames(Result_coefficients)[grepl(Pheno, rownames(Result_coefficients))] -> Rows_to_take
     P_anova = NA
     if (length(Rows_to_take) > 2){ 
       Result_model0 <- tryCatch({
         glmer(formula = Model0, family = binomial(link = "logit"), data = subset_mother_infant_pairs_it)
       }, error = function(e) {
         print(paste0("Model fitting failed for phenotype: ", Pheno, " with error: ", e$message))
         return(NA)
       }) 
       
       if (class(Result_model0) == "logical") {
         next
       }
       
       c
       P_anova$`Pr(>Chisq)`[2] -> P_anova
     }
     Result_coefficients[Rows_to_take,] %>% rownames_to_column("Phenotype") %>% as_tibble() %>% mutate(P_ANOVA = P_anova, Overall_name = Pheno) %>% rbind(Results_association_strainSharing, . ) -> Results_association_strainSharing
   }
   Results_association_strainSharing %>% mutate(SGB = file_short[i], .before=1) -> Results_association_strainSharing
 }
 
 Results_association_strainSharing %>% rbind(Result_all_SGB, . ) -> Result_all_SGB
}

link_SGB_species<-read.delim("../link_SGB_species_names.txt")
Result_all_SGB_species<-left_join(Result_all_SGB, link_SGB_species)
write.table(Result_all_SGB_species, file = "~/Desktop/LLNEXT/Analysis/strainphlan_NEXT/Result_all_SGB_species.txt", row.names = FALSE, sep = "\t", quote = F)

write.table(output_table, file = "~/Desktop/LLNEXT/Analysis/strainphlan_NEXT/thresholds_output_NEXT_TS.txt", row.names = FALSE, sep = "\t", quote = F)

