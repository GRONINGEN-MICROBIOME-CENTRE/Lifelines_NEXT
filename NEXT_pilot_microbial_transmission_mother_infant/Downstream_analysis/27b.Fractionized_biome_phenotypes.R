setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore which external factors may explain the 
# variation in the stability of the infant microbiome and
# virome
#############################################################

##############################
# Functions
##############################

mixed_models_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list, consider_time) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
    # print (c(Bug, Prevalence))
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      
      if (consider_time=='time_as_covariate') {
        Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + infant_mode_delivery + Age_months + (1|Individual_ID)"), collapse="" )) 
      } else { # else is mainly for associating entities with time alone
        Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + (1|Individual_ID)"), collapse="" )) 
      }
      
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      
      if (consider_time=='time_as_covariate') {
        Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + infant_mode_delivery + Age_months + ",pheno2, "+ (1|Individual_ID)"), collapse="" ))
      } else { # else is mainly for associating entities with time alone
        Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + ",pheno2, "+ (1|Individual_ID)"), collapse="" ))
      }
      
      lmer(Model2, df_pheno, REML = F) -> resultmodel2
      M = "Mixed"
      as.data.frame(anova(resultmodel2, base_model))['resultmodel2','Pr(>Chisq)']->p_simp
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
  
}

linear_model_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      
      Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + ",pheno2), collapse="" ))
      lm(Model2, df_pheno) -> resultmodel2
      
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = Summ_simple$`Pr(>|t|)`, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
}
##############################
# Loading libraries
##############################
library(vegan)
library(dplyr)
library(tibble)
library(lme4)
library(RLRsim)
library(lmerTest)

library(reshape2)
library(ggplot2)
library(ggforce)

##############################
# Input data
##############################

MGS_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata <- MGS_metadata[MGS_metadata$Type=='Infant',]
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels = c("M1", "M2", "M3",
                                                                    "M6", "M12"), ordered = T)
row.names(MGS_metadata) <- MGS_metadata$Short_sample_ID


PB_structure <- read.table('03a.RESULTS/Fractions_of_personal_infant_virome_N_and_abundance.txt', sep='\t', header=T)
MGS_metadata <- MGS_metadata[MGS_metadata$Short_sample_ID %in% PB_structure$Sample,]
PB_structure$N_species <- MGS_metadata$viral_richness[match(MGS_metadata$Short_sample_ID, PB_structure$Sample)]

if (sum( rowSums(PB_structure[ , grep( 'N_', colnames(PB_structure[,-length(PB_structure)]) ) ] )!=PB_structure$N_species)==0) {
  PB_structure$perc_PPB <- PB_structure$N_PPB/PB_structure$N_species
  PB_structure$perc_TDB <- PB_structure$N_TDB/PB_structure$N_species
  PB_structure$perc_Singleton <- PB_structure$N_Singletons/PB_structure$N_species
}

PB_structure_perc <- PB_structure[,grep('perc', colnames(PB_structure))]
# CLR-transformation
my_pseudocount_normal=min(PB_structure_perc[PB_structure_perc!=0])/2
PB_structure_perc_CLR<-decostand(PB_structure_perc, "clr", pseudocount=my_pseudocount_normal)
PB_structure_perc_CLR$Short_sample_ID <- row.names(PB_structure_perc_CLR)

PB_structure <- merge(PB_structure, MGS_metadata, by='row.names')
row.names(PB_structure) <- PB_structure$Row.names
PB_structure$Row.names <- NULL

PB_structure_perc <- PB_structure[,grep('ab', colnames(PB_structure))]
##############################
# ANALYSIS
##############################
PB_structure_melt <- melt( PB_structure_perc_CLR, 'Short_sample_ID' )
PB_structure_melt$Individual_ID <- MGS_metadata$Individual_ID[match(MGS_metadata$Short_sample_ID,PB_structure_melt$Short_sample_ID)]
PB_structure_melt$DNA_CONC <-  MGS_metadata$DNA_CONC[match(MGS_metadata$Short_sample_ID,PB_structure_melt$Short_sample_ID)]
PB_structure_melt$Clean_reads <-  MGS_metadata$Clean_reads[match(MGS_metadata$Short_sample_ID,PB_structure_melt$Short_sample_ID)]
PB_structure_melt$NEXT_ID <-  MGS_metadata$NEXT_ID[match(MGS_metadata$Short_sample_ID,PB_structure_melt$Short_sample_ID)]
PB_structure_melt$Age_months <-  MGS_metadata$Age_months[match(MGS_metadata$Short_sample_ID,PB_structure_melt$Short_sample_ID)]

result <- aov(value ~ variable, PB_structure_melt)
summary(result)

PB_mod0 <- lmer(value ~ DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = PB_structure_melt)
PB_mod1  = lmer(value ~ DNA_CONC + Clean_reads + variable + (1|NEXT_ID), REML = F, data = PB_structure_melt)

BIC(PB_mod0, PB_mod1)
summary(PB_mod1) # variableperc_TDB        2.710e-01  6.218e-02  5.310e+02   4.358 1.58e-05 ***
                 # variableperc_Singleton -2.864e-01  6.218e-02  5.310e+02  -4.606 5.15e-06 ***

# dependent on time?
PB_mod2 <- lmer(value ~ DNA_CONC + Clean_reads + variable + (1|NEXT_ID), REML = F, data = PB_structure_melt)
PB_mod3  = lmer(value ~ DNA_CONC + Clean_reads + Age_months + variable + (1|NEXT_ID), REML = F, data = PB_structure_melt)

BIC(PB_mod2, PB_mod3)
summary(PB_mod3) # variableperc_TDB        2.710e-01  6.218e-02  5.310e+02   4.358 1.58e-05 ***
                 # variableperc_Singleton -2.864e-01  6.218e-02  5.310e+02  -4.606 5.15e-06 ***


boxplot(PB_structure_melt$value~PB_structure_melt$variable)


PB_structure_time <- mixed_models_taxa(MGS_metadata, "Short_sample_ID", PB_structure_perc_CLR[,c("ab_PPB", "ab_TDB", "ab_Singletons")], 
                                        c("Age_months"), 
                                        'skip time as covariate')

PB_structure_pheno <- mixed_models_taxa(MGS_metadata, "Short_sample_ID", PB_structure_perc_CLR[,c("ab_PPB", "ab_TDB", "ab_Singletons")],
                                        c("infant_gestational_age",
                                          "infant_place_delivery", 
                                          "infant_ffq_feeding_mode_complex"), 
                                        'time_as_covariate')


PB_structure_solid_food <- linear_model_taxa(MGS_metadata, "Short_sample_ID", PB_structure_perc_CLR[,c("ab_PPB", "ab_TDB", "ab_Singletons")],
                                             "infant_food_solid_intro_M6")



###### OUTPUT #####

