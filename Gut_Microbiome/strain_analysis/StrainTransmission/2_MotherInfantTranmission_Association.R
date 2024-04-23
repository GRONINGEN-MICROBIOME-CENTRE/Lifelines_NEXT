library(tidyverse)
library(patchwork)
library(ggpubr)
library(lme4)

#1. Get RDS file from command line
args <- commandArgs(trailingOnly = TRUE)
file = args[1]
sp4 = read_rds(file)
#2. Get SGB name
SGB_name = str_remove(basename(file), "\\.rds") 

#Functions
source("Common_function.R") #Mother_infant_formatting
Associate_phenos = function(subset_mother_infant_pairs, SGB_name){
  #associatiate phenotype with strain share yes/no. Each row is a mother-infant distance. Noa ccounting for siblings
  #dataframe for association: subset_mother_infant_pairs
  subset_mother_infant_pairs %>% mutate(Strain_sharing = as.factor(Strain_sharing)) -> subset_mother_infant_pairs
  Results_association_strainSharing = tibble()
  #1. loop through phenotypes if at least 10 mother-infant pairs with at least 3 timepoints.
  #If there is only one timepoint per mother, then if there are >= 3 entries per baby, it means it has at least three available timepoints
  #subset_mother_infant_pairs %>% group_by(next_id_infant) %>% summarise(N = n()) %>% filter(N>=3) -> Babies_with_3
  #N_babies_with_3 = dim(Babies_with_3)[1] 
  dim(subset_mother_infant_pairs)[1] -> N_samples
  if (N_samples < 20){ return( tibble(Phenotype = NA, Estimate =NA, `Std. Error` =NA, `z value`=NA, `Pr(>|z|)`=NA, P_ANOVA = NA, Overall_name=NA, Model_converge=NA ) ) }
   Pheno_list = c("mother_birthcard_duration_water_broken_min" , "birth_birthcard_total_duration_delivery_min", "mother_birthcardself_gestational_age_weeks","birth_deliverybirthcard_mode_binary", "birth_deliverybirthcard_place_delivery_simple", "birth_delivery_mode_simple",  "birth_delivery_mode_complex", "mother_birthcard_parity",                      
   "family_pets_any","infant_birthcard_feeding_mode_after_birth", "infant_ffq_feeding_mode_simple",  "infant_ffq_feeding_mode_complex")     
   for ( Pheno in Pheno_list){
     subset_mother_infant_pairs_it = subset_mother_infant_pairs
     if (Pheno %in%  c("mother_birthcard_duration_water_broken_min" , "birth_birthcard_total_duration_delivery_min", "birth_deliverybirthcard_place_delivery_simple") ){ subset_mother_infant_pairs_it %>% filter(birth_deliverybirthcard_mode_binary == "VG") -> subset_mother_infant_pairs_it }
     subset_mother_infant_pairs_it %>% filter(!is.na(!!sym(Pheno))) -> subset_mother_infant_pairs_it
     
     #Check number of transmitted / not transmitted
     subset_mother_infant_pairs_it %>% group_by(Strain_sharing) %>% summarise(N = n()) %>% filter(Strain_sharing == "yes") -> NSharing
     if (NSharing$N < 10){ print(paste0("Not enough strain sharing ", Pheno) ); next } 
     #Check counts per level if categorical
     CLASS = class(as_vector( subset_mother_infant_pairs_it[,Pheno] ) )
     if ( CLASS != "numeric"){ 
         subset_mother_infant_pairs_it %>% group_by(!!sym(Pheno),Strain_sharing) %>% summarise(N = n()) -> D ; D %>% filter(N < 5) -> Problem
         if ( dim(Problem)[1] != 0  ){ print(paste0("Not enough samples per group in ",Pheno )); next }
     }
     #	

     subset_mother_infant_pairs_it[Pheno]

     print(paste0("Iterating through pheno: ", Pheno))
     #2. build the model
     if (! Pheno %in% c("infant_timepoint_categorical", "Age_days_infant") ){
      Model0 =  as.formula(paste0( "Strain_sharing ~ (1|next_id_infant) + infant_timepoint_categorical + mother_timepoint_categorical" ))
      Model = as.formula(paste0( "Strain_sharing ~ (1|next_id_infant) + infant_timepoint_categorical + mother_timepoint_categorical +  ", Pheno ))
      if  (CLASS == "numeric"){  subset_mother_infant_pairs_it[, Pheno] = scale(subset_mother_infant_pairs_it[, Pheno])  }
     } else { 
       Model0 = as.formula(paste0( "Strain_sharing ~ mother_timepoint_categorical + (1|next_id_infant)"))
       Model = as.formula(paste0( "Strain_sharing ~ mother_timepoint_categorical + (1|next_id_infant) +", Pheno)) 
     }
     #glmer(formula = Model, family = binomial(link = "logit"), data = subset_mother_infant_pairs_it) -> Result_model
     Result_model <- tryCatch({
        Fit_glmer(Model, subset_mother_infant_pairs_it)
        #glmer(formula = Model, family = binomial(link = "logit"), data = subset_mother_infant_pairs_it)
     }, error = function(e) {
       print(paste0("Model fitting failed for phenotype: ", Pheno, " with error: ", e$message))
       return(NA)
     }) 
       
     if (class(Result_model) == "logical") {
       next
     }

     if (Result_model[[2]] == T){ next }
     Result_model[[1]] -> Result_model

     #Check if convergence
     Convergence = summary(Result_model)$optinfo$conv$lme4
     if (length(Convergence) == 0 ){ Convergence = "yes"} else { Convergence = "no" }
     
     as.data.frame(summary(Result_model)$coefficients) -> Result_coefficients
     rownames(Result_coefficients)[grepl(Pheno, rownames(Result_coefficients))] -> Rows_to_take
     P_anova = NA
     if (length(Rows_to_take) > 2){ 
       Result_model0 <- tryCatch({
         glmer(formula = Model0, family = binomial(link = "logit"),data = subset_mother_infant_pairs_it)
       }, error = function(e) {
         print(paste0("Model fitting failed for phenotype: ", Pheno, " with error: ", e$message))
         return(NA)
       }) 
       
       if (class(Result_model0) == "logical") {
         next
       }
       anova(Result_model,Result_model0) -> P_anova
       P_anova$`Pr(>Chisq)`[2] -> P_anova
     }

     print("Preparing plot")
     plot_name = paste0("plots_associations/sharing_", SGB_name,"_vs_",Pheno , ".pdf")
     if ( CLASS != "numeric"){ 
     Plot_stuff(subset_mother_infant_pairs,Pheno, plot_name, Plot_type="Bar", Do_facet = F )
     } else {
	Plot_stuff(subset_mother_infant_pairs,Pheno, plot_name, Plot_type="Box", Do_facet = F )	
	}

     Result_coefficients[Rows_to_take,] %>% rownames_to_column("Phenotype") %>% as_tibble() %>% mutate(P_ANOVA = P_anova, Overall_name = Pheno, Model_converge = Convergence) %>% rbind(Results_association_strainSharing, . ) -> Results_association_strainSharing
   }
   Results_association_strainSharing %>% mutate(SGB = SGB_name, .before=1) -> Results_association_strainSharing
  return(Results_association_strainSharing)
}
Associated_phenos_age = function(subset_mother_infant_pairs, SGB_name){
  print("Running association with age")
  subset_mother_infant_pairs %>% group_by(next_id_infant) %>% summarise(N = n()) %>% filter(N>=2) -> Babies_with_3
  N_babies_with_3 = dim(Babies_with_3)[1]
  if (N_babies_with_3 < 10){ print("Not enough longitudinal samples with three timepoints"); return( tibble(SGB=SGB_name,Phenotype = "infant_timepoint_categorical", Estimate =NA, `Std. Error` =NA, `z value`=NA, `Pr(>|z|)`=NA, P_ANOVA = NA, Overall_name="infant_timepoint_categorical", Model_converge=NA ) ) } 
  subset_mother_infant_pairs %>% group_by(infant_timepoint_categorical) %>% summarise(N = n() ) -> Nub
  Nub %>% filter(N >= 5) -> Nub
  if (dim(Nub)[1] < 3){ print("Not enough longitudinal samples"); return( tibble(SGB=SGB_name, Phenotype = "infant_timepoint_categorical", Estimate =NA, `Std. Error` =NA, `z value`=NA, `Pr(>|z|)`= NA, P_ANOVA = NA, Overall_name="infant_timepoint_categorical", Model_converge=NA ) ) }
  subset_mother_infant_pairs %>% filter( infant_timepoint_categorical %in% Nub$infant_timepoint_categorical) -> subset_mother_infant_pairs

  Model = "as.factor(Strain_sharing) ~  Age_days_infant  + mother_timepoint_categorical + (1|next_id_infant)"
  Res = Fit_glmer(Model, subset_mother_infant_pairs)
  if (Res[[2]] == T){
	return( tibble(SGB=SGB_name, Phenotype = "infant_timepoint_categorical", Estimate =NA, `Std. Error` =NA, `z value`=NA, `Pr(>|z|)`= NA, P_ANOVA = NA, Overall_name="infant_timepoint_categorical", Model_converge="no" ) )
  }
  Res[[1]] -> Result_model   
	
  #glmer( as.factor(Strain_sharing) ~  Age_days_infant  + mother_timepoint_categorical + (1|next_id_infant) , subset_mother_infant_pairs ,  family = binomial(link = "logit") ) %>% summary() -> Result_model
  
  Convergence = summary(Result_model)$optinfo$conv$lme4
  if (length(Convergence) == 0 ){ Convergence = "yes"} else { Convergence = "no" }
  as.data.frame(summary(Result_model)$coefficients) -> Result_coefficients

  plot_name = paste0("plots_associations/sharing_", SGB_name, "_vs_Age_days_infant.pdf")
  Plot_stuff(subset_mother_infant_pairs,  "infant_timepoint_categorical" , plot_name, Plot_type="Bar" )

  Result_coefficients["Age_days_infant",] %>% rownames_to_column("Phenotype") %>% as_tibble() %>% mutate(P_ANOVA = NA, Overall_name = "Age_days_infant", Model_converge = Convergence) %>% mutate(SGB = SGB_name, .before=1) -> Res
  return(Res)
}


print("Formatting")
#4. Format
Mother_infant_formatting(sp4) ->  subset_mother_infant_pairs
print("Running associations")
#5. Run associations
Results_age = Associated_phenos_age(subset_mother_infant_pairs, SGB_name)
Result_all_phenos =  Associate_phenos(subset_mother_infant_pairs, SGB_name)
print("Saving")
#6. Save output
Output_v = paste0("Results/", SGB_name, "_summaryLMERTEST_COV_infanttime+mothertime+nextid.tsv")

print(Results_age)
print(Result_all_phenos)
write_tsv(rbind(Results_age,Result_all_phenos), Output_v)









