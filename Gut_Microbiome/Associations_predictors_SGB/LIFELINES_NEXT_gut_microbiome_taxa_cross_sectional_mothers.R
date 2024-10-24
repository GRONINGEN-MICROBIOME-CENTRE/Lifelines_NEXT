### LIFELINES NEXT GUT MICROBIOME ANALYSIS TAXA MOTHERS ####
### AUTHOR:  rishla Sinha & Alexander Kurilshikov
### ORIGINAL SCRIPT: 3rd AUGUST, 2023
### LAST UPDATE: 25th June, 2024

rm (list = ls())

#load packages 
library(tidyverse)
library(mmrm)

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
  
  data.fit = try(data.frame(bac = bacteria,trait = phenotype,Time = time,Time_cat = time_cat,NEXT_ID = NEXT_ID,covariates))
  if(class(data.fit)[1]!="data.frame") stop("check that your input data is syncronized, i.e. all inputs have the same samples in the same order")
  data.fit = data.fit[complete.cases(data.fit),]
  
  if(length(table(as.character(data.fit$Time_cat))) > 3) {
    formula.null.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + Time3 + trait + us(Time_cat|NEXT_ID)")
    formula.run.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + Time3 + trait + trait:Time1 + us(Time_cat|NEXT_ID)")
    formula.null.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + Time3 + trait + cs(Time_cat|NEXT_ID)")
    formula.run.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + Time3 + trait + trait:Time1 + cs(Time_cat|NEXT_ID)")
    
    data.fit2 = data.frame(data.fit, poly(data.fit$Time,3))
    colnames(data.fit2)[(ncol(data.fit2)-2): ncol(data.fit2)] = c("Time1","Time2","Time3")
  } else if (length(table(as.character(data.fit$Time_cat))) > 2) {
    formula.null.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + trait + us(Time_cat|NEXT_ID)")
    formula.run.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + trait + trait:Time1 + us(Time_cat|NEXT_ID)")
    formula.null.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + trait + cs(Time_cat|NEXT_ID)")
    formula.run.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + trait + trait:Time1 + cs(Time_cat|NEXT_ID)")
    
    data.fit2 = data.frame(data.fit, poly(data.fit$Time,2))
    colnames(data.fit2)[(ncol(data.fit2)-1): ncol(data.fit2)] = c("Time1","Time2")
  } else if (length(table(as.character(data.fit$Time_cat))) > 1){
    formula.null.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + trait + us(Time_cat|NEXT_ID)")
    formula.run.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + trait + trait:Time1 + us(Time_cat|NEXT_ID)")
    formula.null.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + trait + cs(Time_cat|NEXT_ID)")
    formula.run.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + trait + trait:Time1 + cs(Time_cat|NEXT_ID)")
    
    data.fit2 = data.frame(data.fit, poly(data.fit$Time,1))
    colnames(data.fit2)[ncol(data.fit2)] = c("Time1")
  } else {
    stop ("Data has only one timepoint! Swtich to lm() to run the models")
  }
  #run ML models with US matrix
  nullModel.reml = try(mmrm(as.formula(formula.null.us), data= data.fit2,reml = T))
  runModel.reml = try(mmrm(as.formula(formula.run.us), data= data.fit2,reml = T))
  
  if(class(runModel.reml)[1] != "try-error") summary.runModel.reml <-try(summary(runModel.reml))
  if(class(nullModel.reml)[1] != "try-error") summary.nullModel.reml <-try(summary(nullModel.reml))
  if(!exists("summary.runModel.reml")) {
    summary.runModel.reml = NA
    class(summary.runModel.reml) = "try-error"
  }
  if(!exists("summary.nullModel.reml")) {
    summary.nullModel.reml = NA
    class(summary.nullModel.reml) = "try-error"
  }
  
  if (class(nullModel.reml)[1] == "try-error"|
      class(runModel.reml)[1] == "try-error"|
      class(summary.nullModel.reml)[1] == "try-error"|
      class(summary.runModel.reml)[1] == "try-error") {
    runModel.reml = try(mmrm(as.formula(formula.run.cs), data= data.fit2,reml = T))
    nullModel.reml = try(mmrm(as.formula(formula.null.cs), data= data.fit2,reml = T))
    if(class(runModel.reml)[1] != "try-error") summary.runModel.reml <-try(summary(runModel.reml))
    if(class(nullModel.reml)[1] != "try-error") summary.nullModel.reml <-try(summary(nullModel.reml))
    model.type = "CompoundSymmetry"
  } else {model.type = "Unstructured"}
  
  if(class(nullModel.reml)[1]=="try-error"|
     class(runModel.reml)[1]=="try-error") {
    trait.results = data.frame(row.names = NULL,
                               type = "Failure",
                               Covar.type = model.type,
                               bac = NA,
                               trait = NA,
                               N = nrow(data.fit2),
                               levels = NA,
                               "Estimate" = NA,
                               "Std. Error" = NA,
                               "df" = NA,
                               "t value"= NA,
                               "Pr(>|t|)" = NA 
    )
    timeInt.results = data.frame(row.names = NULL,
                                 type = "Failure",
                                 Covar.type = model.type,
                                 bac = NA,
                                 trait = NA,
                                 N = nrow(data.fit2),
                                 levels = NA,
                                 "Estimate" = NA,
                                 "Std. Error" = NA,
                                 "df" = NA,
                                 "t value"= NA,
                                 "Pr(>|t|)" = NA
    )
  } else {
    trait.results = data.frame(row.names = NULL,
                               type = "Success",
                               Covar.type = model.type,
                               bac = NA,
                               trait = NA,
                               N = nrow(data.fit2),
                               levels = sub("trait",
                                            "",
                                            rownames(summary.nullModel.reml$coef)[
                                              grep("trait",rownames(summary.nullModel.reml$coef))]),
                               summary.nullModel.reml$coef[
                                 grep("trait",rownames(summary.nullModel.reml$coef)),,drop =F])
    timeInt.results = data.frame(row.names = NULL,
                                 type = "Success",
                                 Covar.type = model.type,
                                 bac = NA,
                                 trait = NA,
                                 N = nrow(data.fit2),
                                 levels = sub("trait",
                                              "",
                                              rownames(summary.runModel.reml$coef)[
                                                grep(":",rownames(summary.runModel.reml$coef))]),
                                 summary.runModel.reml$coef[grep(":",rownames(summary.runModel.reml$coef)),,drop = F])
  }
  
  list(trait = trait.results,
       time.trait = timeInt.results)
}

run_mmrm_analysis = function( bacteria, phenotypes, time, time_cat, NEXT_ID,covariates){
  output_names = c("trait","time.trait")
  output = sapply(output_names, function(x) NULL)
  for (i in 1:ncol(bacteria)) {
    print (i)
    for(j in 1:ncol(phenotypes)){
      print(j)
      singleAssoc = .run_mmrm_single(bacteria[,i],phenotypes[,j],time = time,time_cat = time_cat, NEXT_ID = NEXT_ID,covariates = covariates)
      singleAssoc$trait[,"bac"] = colnames(bacteria)[i]
      singleAssoc$trait[,"trait"] = colnames(phenotypes)[j]
      singleAssoc$time.trait[,"bac"] = colnames(bacteria)[i]
      singleAssoc$time.trait[,"trait"] = colnames(phenotypes)[j]
      output$trait = rbind(output$trait,singleAssoc$trait)
      output$time.trait = rbind(output$time.trait,singleAssoc$time.trait)
    }
    
  }
  colnames(output$trait)[8] = "SE"
  colnames(output$time.trait)[8] = "SE"
  colnames(output$time.trait)[11] = "P"
  colnames(output$trait)[11] = "P"
  
  
  output$trait$FDR = p.adjust(output$trait$P)
  output
  
}


# Load metadata and phenotypes
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
cross_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_2023_11_15.txt")
cross_selection <-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_selection_AFTER_correlation_&_correction_09_05_2024.txt")


# Merging files and selecting mother relevant data
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_mothers<-metadata[metadata$Type=="mother", ]
metadata_mothers <- subset(metadata_mothers, Timepoint_categorical != "M1" & Timepoint_categorical != "M2")
names (metadata_mothers)[2]<-"next_id_mother"

# Selecting data relevant for mother 
mother_cross_selection<-cross_selection[cross_selection$mother_microbiome_selection==1,]
column_names <- mother_cross_selection[[1]]
mother_cross_phenotypes <- cross_phenotypes %>% select(all_of(column_names))
mother_cross_phenotypes <- subset(mother_cross_phenotypes, twin_pair != "pair_2")
mother_cross_phenotypes$FAMILY=NULL
mother_cross_phenotypes$infant_relations=NULL

mother_metadata_cross_phenotypes<-left_join(metadata_mothers, mother_cross_phenotypes, by="next_id_mother")
row.names(mother_metadata_cross_phenotypes)<-mother_metadata_cross_phenotypes$NG_ID
mother_metadata_cross_phenotypes_1<-mother_metadata_cross_phenotypes
# Remove phenotypes with too many NA's 
mother_metadata_cross_phenotypes_2 <-drop_phenotypes(mother_metadata_cross_phenotypes_1, 200, 5)
names (mother_metadata_cross_phenotypes_2)[2]<-"NEXT_ID"
mother_metadata_cross_phenotypes_2$NEXT_ID=as.factor(mother_metadata_cross_phenotypes_2$NEXT_ID)


#Taxa
taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/NEXT_metaphlan_4_CLR_transformed_fil_30_percent_SGB_mothers_03_08_2023.txt")

# Analysis

# Select outcome 
taxa <- taxa[match(rownames(mother_metadata_cross_phenotypes_2), rownames(taxa)),]

# Select and normalize phenotypes 
metadata_columns <- c("NG_ID", "NEXT_ID","Modified_NEXT_ID_without_preg_number", "days_from_first_collection","FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                      "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                      "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                      "exact_age_days_at_collection", "exact_age_months_at_collection",
                      "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                      "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                      "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                      "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint")


phenos <- mother_metadata_cross_phenotypes_2[, !(colnames(mother_metadata_cross_phenotypes_2) %in% metadata_columns)] 

phenos<-run_invrank_dataFrame(phenos)

# Select and normalize covariates 
mother_metadata_cross_phenotypes_2$BATCH_NUMBER=as.factor(mother_metadata_cross_phenotypes_2$BATCH_NUMBER)
covariates<-mother_metadata_cross_phenotypes_2[,c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER")]
covariates_normalized<-run_invrank_dataFrame(covariates) 

phenos$next_id_infant=NULL
results_1<-run_mmrm_analysis(taxa, phenos,time = mother_metadata_cross_phenotypes_2$exact_age_months_at_collection, time_cat = mother_metadata_cross_phenotypes_2$Timepoint_categorical, NEXT_ID = mother_metadata_cross_phenotypes_2$NEXT_ID, covariates =  covariates_normalized )

results_trait<-as.data.frame(results_1$trait)
results_trait$FDR<-p.adjust (results_trait$P, method = "fdr")
results_trait$trait_group <- sapply(strsplit(as.character(results_trait$trait), "_"), function(x) paste(x[1:2], collapse = "_"))


results_trait_time_interaction<-as.data.frame(results_1$time.trait)
results_trait_time_interaction$FDR<-p.adjust (results_trait_time_interaction$P, method = "fdr")
results_trait_time_interaction$trait_group <- sapply(strsplit(as.character(results_trait_time_interaction$trait), "_"), function(x) paste(x[1:2], collapse = "_"))

setwd("~/Desktop/LLNEXT/Analysis/results/taxa/mothers/")


write.table(results_trait, "taxa_mmrm_trait_cross_sectional_mothers.txt", sep = "\t", row.names = F)
write.table(results_trait_time_interaction, "taxa_mmrm_trait_int_time_cross_sectional_mothers.txt", sep = "\t", row.names = F)

all<-merge(mother_metadata_cross_phenotypes_2, taxa, by="row.names")
row.names(all)<-all$Row.names

ggplot(all, aes(x = mother_ffq_group_coffee_ml_per_day_m2, y = `GGB3570_SGB4777.t__SGB4777`)) +
  geom_point(color = "#1F77B4") +  # Add the scatter plot points with a specific color
  geom_smooth(method = "lm", se = T, color = "#FF7F0E") +  # Add linear regression line with a different color
  labs(x = "Coffee Consumption (ml/day)", y = "CLR Abundance of SGB4777") +
  theme_minimal() +
  ggtitle("Coffee Consumption vs CLR Abundance of SGB4777")

