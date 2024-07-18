### LIFELINES NEXT GUT MICROBIOME ANALYSIS TAXA INFANTS ####
### AUTHOR:  TRISHLA SINHA
### ORIGINAL SCRIPT: 3rd AUGUST, 2023
### LAST UPDATE: 16th November, 2023
### Updated by Cyrus on 30 April 2024 : adapted for HMO associations

#load packages 
library(tidyverse)
library(mmrm)
library(wesanderson)
library(reshape2)
library(here)

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

.run_likelihoodRatioTest = function(model1, model2) {
  class_checkpoint = ("mmrm" %in% class(model1) & "mmrm" %in% class(model2))
  reml_checkpoint = (model1$reml == F & model2$reml == F)
  if (class_checkpoint == FALSE) stop ("One or both models do not belong to MMRM class. Execution stopped!")
  if (reml_checkpoint == FALSE) stop ("One or both models use REML for calculation. The results will be unreliable. Execution stopped!")
  result.matrix = matrix(rep(NA,10),ncol = 5)
  colnames(result.matrix) <- c("#Df", "LogLik", "Df", "Chisq", "Pr(>Chisq)")
  rownames(result.matrix) <- c("Model1","Model2")
  result.matrix[1,1] <- length(model1$beta_est)
  result.matrix[2,1] <- length(model2$beta_est)
  result.matrix[1,2] <- model1$neg_log_lik
  result.matrix[2,2] = model2$neg_log_lik
  result.matrix[1,3] = 0
  result.matrix[2,3] = length(model2$beta_est) - length(model1$beta_est)
  result.matrix[1,4] = 0
  result.matrix[2,4] = 2*abs(result.matrix[2,2] - result.matrix[1,2])
  result.matrix[1,5] = NA
  result.matrix[2,5] = pchisq(result.matrix[2,4], round(abs(result.matrix[2,3])),lower.tail = F)
  result.matrix 
}

.run_mmrm_single = function(bacteria,phenotype,time,time_cat,NEXT_ID,covariates) {
  if(class(bacteria) !="numeric") stop("bacteria should be numeric vector!")
  if(class(time)!="numeric") stop("time should be 'double'!")
  if(class(time_cat)!='factor') stop("time_cat should be 'factor'")
  if(class(NEXT_ID)!='factor') stop("NEXT_ID should be 'factor'")
  
  formula.null.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),"+ Time1 + Time2 + Time3 + us(Time_cat|NEXT_ID)")
  formula.run.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),"+ Time1 + Time2 + Time3 + trait + trait:Time1 + us(Time_cat|NEXT_ID)")
  formula.null.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),"+ Time1 + Time2 + Time3 + cs(Time_cat|NEXT_ID)")
  formula.run.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),"+ Time1 + Time2 + Time3 + trait + trait:Time1 + cs(Time_cat|NEXT_ID)")
  
  
  data.fit = try(data.frame(bac = bacteria,trait = phenotype,Time = time,Time_cat = time_cat,NEXT_ID = NEXT_ID,covariates))
  if(class(data.fit)[1]!="data.frame") stop("check that your input data is syncronized, i.e. all inputs have the same samples in the same order")
  data.fit = data.fit[complete.cases(data.fit),]
  data.fit2 = data.frame(data.fit, poly(data.fit$Time,3))
  colnames(data.fit2)[(ncol(data.fit2)-2): ncol(data.fit2)] = c("Time1","Time2","Time3")
  #run ML models with US matrix
  nullModel.ml= try(mmrm(as.formula(formula.null.us), data= data.fit2,reml = F))
  runModel.ml = try(mmrm(as.formula(formula.run.us), data= data.fit2,reml = F))
  
  #run REML models with US matrix
  runModel.reml = try(mmrm(as.formula(formula.run.us), data= data.fit2,reml = T))
  #if(class(runModel.reml)[1] != "try-error") 
  summary.runModel.reml <-try(summary(runModel.reml))
  
  if (class(nullModel.ml)[1] == "try-error"|
      class(runModel.ml)[1] == "try-error"|
      class(runModel.reml)[1] == "try-error") {
    runModel.reml = try(mmrm(as.formula(formula.run.cs), data= data.fit2,reml = T))
    if(class(runModel.reml)[1] != "try-error") summary.runModel.reml <-try(summary(runModel.reml))
    
    nullModel.ml = try(mmrm(as.formula(formula.null.cs), data= data.fit2,reml = F))
    runModel.ml = try(mmrm(as.formula(formula.run.cs), data= data.fit2,reml = F))
    model.type = "CompoundSymmetry"
  } else {model.type = "Unstructured"}
  
  anova.results = try(.run_likelihoodRatioTest(nullModel.ml, runModel.ml))
  if(class(anova.results)[1] != "try-error") anova.results <- data.frame(type = "Success",
                                                                         Covar.type = model.type,
                                                                         bac = NA,
                                                                         trait = NA,
                                                                         N=nrow(data.fit2),
                                                                         Chisq = anova.results[2,4],
                                                                         P=anova.results[2,5])
  
  reml.results = summary.runModel.reml$coef[grep("trait",rownames(summary.runModel.reml$coef)),]
  reml.results = data.frame(type = "Success",
                            Covar.type = model.type,
                            bac = NA,
                            trait = NA,
                            N = nrow(data.fit2),
                            levels = rownames(reml.results),
                            reml.results)  
  if(class(data.fit2$trait)=="factor") reml.results$levels <- sub("trait","",reml.results$levels)
  list(nested = anova.results,
       byLevel = reml.results)
}


run_mmrm_analysis = function( bacteria, phenotypes, time, time_cat, NEXT_ID,covariates){
  output_names = c("nested","byLevel")
  output = sapply(output_names, function(x) NULL)
  for (i in 1:ncol(bacteria)) {
    print (i)
    for(j in 1:ncol(phenotypes)){
      print(j)
      print(paste("Bacterium: ", colnames(bacteria)[i]))
      print(paste("Phenotype: ", colnames(phenotypes)[j]))
      singleAssoc = .run_mmrm_single(bacteria[,i],phenotypes[,j],time = time,time_cat = time_cat, NEXT_ID = NEXT_ID,covariates = covariates)
      print(singleAssoc)
      singleAssoc$nested[,"bac"] = colnames(bacteria)[i]
      singleAssoc$nested[,"trait"] = colnames(phenotypes)[j]
      singleAssoc$byLevel[,"bac"] = colnames(bacteria)[i]
      singleAssoc$byLevel[,"trait"] = colnames(phenotypes)[j]
      output$nested = rbind(output$nested,singleAssoc$nested)
      output$byLevel = rbind(output$byLevel,singleAssoc$byLevel)
    }
    
  }
  
  output
  
}


# Load metadata and phenotypes
metadata<-read.delim(file = here("data_meta/LLNEXT_metadata_03_08_2023.txt"))
cross_phenotypes<-read.delim(file = here("data_phenotypes/masterfile_cross_sectional_2023_11_15.txt"))
cross_selection <-read.delim(file = here("data_phenotypes/masterfile_cross_sectional_selection_AFTER_correlation_&_correction_09_04_2024.txt"))


# Merging files and selecting infant relevant data
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"

# Selecting data relevant for infant 
infant_cross_selection<-cross_selection[cross_selection$prioritized_phenotypes==1,] #prioritized_phenotypes: make sure to select this column!
column_names <- as.vector(na.omit(infant_cross_selection[[1]])) # removed NAs 30 april 2024
infant_cross_phenotypes <- cross_phenotypes %>% select(all_of(column_names))
infant_metadata_cross_phenotypes<-left_join(metadata_infants, infant_cross_phenotypes)
row.names(infant_metadata_cross_phenotypes)<-infant_metadata_cross_phenotypes$NG_ID

# For this specific data remove duplicates for categorical timepoints 
infant_metadata_cross_phenotypes_1 <- infant_metadata_cross_phenotypes %>%
  distinct(SAMPLE_ID, .keep_all = TRUE) 

# Remove phenotypes with too many NA's 
infant_metadata_cross_phenotypes_2 <-drop_phenotypes(infant_metadata_cross_phenotypes_1, 200, 5)
names (infant_metadata_cross_phenotypes_2)[2]<-"NEXT_ID"
infant_metadata_cross_phenotypes_2$NEXT_ID=as.factor(infant_metadata_cross_phenotypes_2$NEXT_ID)

#Taxa
#taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")

# Analysis

# Select outcome 
#taxa_1 <- taxa[match(rownames(infant_metadata_cross_phenotypes_2), rownames(taxa)),]


# Select and normalize phenotypes 
metadata_columns <- c("FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                      "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_2",
                      "reads_lost_QC", "isolation_method", "NG_ID_short",
                      "exact_age_days_at_collection", "exact_age_months_at_collection",
                      "exact_age_years_at_collection", "SAMPLE_ID",
                      "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                      "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                      "shannon", "next_id_mother", "next_id_partner", "sibling_number", "timepoint") # I took out "NG_ID", "clean_reads_FQ_1","dna_conc", "BATCH_NUMBER", "NEXT_ID", "Timepoint_categorical" so that it stays.


phenos <- infant_metadata_cross_phenotypes_2[, !(colnames(infant_metadata_cross_phenotypes_2) %in% metadata_columns)] 
phenos <- run_invrank_dataFrame(phenos)


# # Select and normalize covariates 
# infant_metadata_cross_phenotypes_2$BATCH_NUMBER=as.factor(infant_metadata_cross_phenotypes_2$BATCH_NUMBER)
# covariates<-infant_metadata_cross_phenotypes_2[,c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER")]
# covariates_normalized<-run_invrank_dataFrame(covariates) 
# 
# results_1<-run_mmrm_analysis(taxa_1, phenos,time = infant_metadata_cross_phenotypes_2$exact_age_months_at_collection, time_cat = infant_metadata_cross_phenotypes_2$Timepoint_categorical, NEXT_ID = infant_metadata_cross_phenotypes_2$NEXT_ID, covariates =  covariates_normalized )
# 
# results_nested<-as.data.frame(results_1$nested)
# results_nested$FDR<-p.adjust (results_nested$P, method = "fdr")
# results_nested$trait_group <- sapply(strsplit(as.character(results_nested$trait), "_"), function(x) paste(x[1:2], collapse = "_"))
# 
# results_level<-as.data.frame(results_1$byLevel)
# results_level$FDR<-p.adjust (results_level$P, method = "fdr")
# results_level$trait_group <- sapply(strsplit(as.character(results_level$trait), "_"), function(x) paste(x[1:2], collapse = "_"))
# 
# setwd("~/Desktop/LLNEXT/Analysis/results/taxa/")
# 
# write.table(results_nested, "taxa_mmrm_nested_cross_sectional_infants.txt", sep = "\t", row.names = F)
# write.table(results_level, "taxa_mmrm_by_levels_cross_sectional_infants.txt", sep = "\t", row.names = F)
# 
# 
# 
# # Plots 
# infant_metadata_cross_phenotypes_2$Timepoint_categorical<-factor(infant_metadata_cross_phenotypes_2$Timepoint_categorical, levels = c("P12", "P28", "B", "W2" ,"M1", "M2", "M3", "M6", "M9", "M12"))
# 
# all<-merge(infant_metadata_cross_phenotypes_2, taxa_1, by="row.names" )
# row.names(all)<-all$Row.names
# 
# blood_group_genotype<-all %>%
#   filter(!is.na(mother_genetics_blood_group_genotype))
# 
# 
# ggplot(blood_group_genotype, aes(Timepoint_categorical, y=Actinomyces_urogenitalis.t__SGB15868, fill=mother_genetics_blood_group_genotype, color=mother_genetics_blood_group_genotype)) +
#   scale_fill_manual(values=c("#e60049", "#0bb4ff", "#50e991", "#f46a9b", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#b30000", "#7c1158", "#4421af"))+
#   scale_color_manual(values=c("#e60049", "#0bb4ff", "#50e991", "#f46a9b", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#b30000", "#7c1158", "#4421af"))+
#   geom_boxplot(alpha=0.4, outlier.colour = NA)+
#   geom_point(alpha=0.6,
#              position = position_jitterdodge(jitter.width = 0.01, jitter.height = 0))+
#   # geom_jitter() +
#   ggtitle("")+
#   theme_bw()+labs(x="", y = "CLR Actinomyces urogenitalis", fill="Maternal blood group", color="Maternal blood group")+
#   theme(
#     plot.title = element_text(color="black", size=22, face="bold"),
#     axis.title.x = element_text(color="black", size=22, face="bold"),
#     axis.title.y = element_text(color="black", size=18, face="bold"),
#     axis.text.y = element_text(face="bold", size=10),
#     axis.text.x = element_text(size=22, angle = 60, hjust = 1))
# #strip.text.x = element_text(size = 10))
# 
# 
# filtered_cols <- results_level %>%
#   filter(FDR < 0.05, trait == "birth_deliverybirthcard_mode_binary", levels=="traitVG") %>%
#   select(bac) %>%
#   unlist() %>%
#   as.character()
# 
# # Now, grep these column names in 'all'
# selected_columns <- grep(paste(filtered_cols, collapse = "|"), names(all), value = TRUE)
# selected_columns <- c(selected_columns, "birth_deliverybirthcard_mode_binary")
# 
# # Subset 'all' dataframe with these columns
# subset_all <- all[, selected_columns]
# 
# all_long <- melt(subset_all, id.vars = "birth_deliverybirthcard_mode_binary", 
#                  variable.name = "Species", value.name = "CLR_Abundance")
# all_long <- na.omit(all_long)
# all_long$Species <- sub("\\..*", "", all_long$Species)
# 
# 
# # PDF 10x 8
# ggplot(all_long, aes(Species, y=CLR_Abundance, fill=birth_deliverybirthcard_mode_binary, color=birth_deliverybirthcard_mode_binary)) +
#   scale_fill_manual(values=c("#e60049", "blue"))+
#   scale_color_manual(values=c("#e60049", "blue"))+
#   geom_boxplot(alpha=0.4, outlier.colour = NA)+
#   geom_point(alpha=0.6,
#              position = position_jitterdodge(jitter.width = 0.01, jitter.height = 0))+
#   # geom_jitter() +
#   ggtitle("")+
#   theme_bw()+labs(x="", y = "CLR_Abundance", fill="", color="")+
#   theme(
#     plot.title = element_text(color="black", size=18, face="bold"),
#     axis.title.x = element_text(color="black", size=18, face="bold"),
#     axis.title.y = element_text(color="black", size=18, face="bold"),
#     axis.text.y = element_text(face="bold", size=18),
#     axis.text.x = element_text(size=18, angle = 90, hjust = 1))+
#   coord_flip()
# #strip.text.x = element_text(size = 10))
# 
# 
# 
# 
# filtered_cols <- results_level %>%
#   filter(FDR < 0.002, trait == "mother_birthcard_parity", levels=="trait") %>%
#   select(bac) %>%
#   unlist() %>%
#   as.character()
# 
# # Now, grep these column names in 'all'
# selected_columns <- grep(paste(filtered_cols, collapse = "|"), names(all), value = TRUE)
# selected_columns <- c(selected_columns, "mother_birthcard_parity")
# 
# # Subset 'all' dataframe with these columns
# subset_all <- all[, selected_columns]
# 
# all_long <- melt(subset_all, id.vars = "mother_birthcard_parity", 
#                  variable.name = "Species", value.name = "CLR_Abundance")
# all_long <- na.omit(all_long)
# all_long$Species <- sub("\\..*", "", all_long$Species)
# 
# ggplot(all_long, aes(x = as.factor(mother_birthcard_parity), y = CLR_Abundance, color = Species)) +
#   geom_boxplot(alpha = 0.4, outlier.colour = NA) +  # Add boxplots
#   geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.01, jitter.height = 0)) +  # Add jittered points
#   scale_color_manual(values = c("#1F77B4", "#2CA02C", "#FF7F0E")) +  # Manually set colors for each species
#   labs(x = "Parity", y = "CLR Abundance", color = "Species") +
#   theme_minimal() +
#   ggtitle("Parity vs CLR Abundance by Species")
# 
