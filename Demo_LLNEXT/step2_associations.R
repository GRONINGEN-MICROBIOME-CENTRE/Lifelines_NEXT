############################ Microbiome-phenotype associations ###########
# Authors: Sergio Andreu Sanchez, Trishla Sinha, Alex Kurilshikov 
# Last update: 25th of Nov, 2024 
# This script is dedicated phenotypes-microbiome features associations
# using mixed models and MMRM

# loading packages, data 

library(mmrm)
library(lmerTest)
library(vegan)
library (mgcv) 

metaphlan = read.table("simulated_metaphlan.txt",header=T)
metadata = read.table("simulated_metadata.txt",header=T)

metadata$Timepoint = factor(metadata$Timepoint,levels = c("P12","P28","B","W2","M1","M3","M6","M9","M12"))
metadata$Covariate1 = factor(metadata$Covariate1)
metadata$Binary_phenotype_dynamic = factor(metadata$Binary_phenotype_dynamic)
metadata$Factor_phenotype_dynamic = factor(metadata$Factor_phenotype_dynamic)
metadata$Binary_phenotype_static = factor(metadata$Binary_phenotype_static)
metadata$Factor_phenotype_static = factor(metadata$Factor_phenotype_static)
metadata$Family = factor(metadata$Family)
metaphlan_infant = metaphlan[metadata$Type=="Infant",]
metadata_infant = metadata[metadata$Type=="Infant",]

metaphlan = decostand(metaphlan, method = "clr",pseudocount = min(metaphlan[metaphlan>0])/2)


# 2.1 Static phenotypes -------------------------------------------------------
# Adopted from https://github.com/GRONINGEN-MICROBIOME-CENTRE/Lifelines_NEXT/Gut_Microbiome/Associations_predictors_SGB/LIFELINES_NEXT_gut_microbiome_alpha_diversity_cross_sectional_infants.R

## Creating functions
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

## Example call
run_mmrm_analysis(metaphlan_infant[,1:5], 
                  metadata_infant[,c("Quant_phenotype_static","Quant_phenotype_static","Factor_phenotype_static"),drop =F ],
                  time = as.double(factor(metadata_infant$Timepoint,levels = c("W2","M1","M3","M6","M12"))),
                  time_cat = factor(metadata_infant$Timepoint,levels = c("W2","M1","M3","M6","M12")),
                  NEXT_ID = factor(metadata_infant$Family),
                  covariates = metadata_infant[,c("Covariate1","Covariate2")]
                  )

# 2.2 Dynamic phenotypes -------------------------------------------------------
# Adopted from https://github.com/GRONINGEN-MICROBIOME-CENTRE/Lifelines_NEXT/Gut_Microbiome/Associations_predictors_SGB/LIFELINES_NEXT_gut_microbiome_alpha_diversity_taxa_pathways_dynamic_infants.R

## main function
gam_function <- function(taxa, phenotypes, phenotype_list) {
   
  covariates <- phenotypes[, c("Covariate1","Covariate2")] 
  
  result = foreach(i = 1:ncol(taxa), .combine = rbind) %:% 
    foreach(j = 1:length(phenotype_list), .combine = rbind) %do% {
      print(paste(i, j))
      data.fit <- data.frame(bac = taxa[,i],
                             trait = phenotypes[,phenotype_list[j]],
                             Time = as.integer(phenotypes$Timepoint),
                             NEXT_ID = phenotypes$Family,
                             covariates)
      data.fit <- data.fit[complete.cases(data.fit),]
      #Timepoint data in the real analysis is given in days; thus we add noise here
      data.fit$Time = data.fit$Time + rnorm(nrow(data.fit),sd = 0.1)
      
      mod_gam1 <- tryCatch({
        bam(bac ~ trait +  Covariate1 + Covariate2 +s(Time) + s(NEXT_ID, bs = 're'), data = data.fit, discrete = TRUE)
      }, warning = function(w) {
        return(NULL)  
      }, error = function(e) {
        return(NULL)  
      })
      
      if (is.null(mod_gam1)) {
        converged = FALSE
        beta = statistics[,1]
        SE = statistics[,2]
        P = statistics[,4]
        N = nrow(data.fit)
      } else {
        summary.results <- summary(mod_gam1)
        statistics = summary.results$p.table[2:(which(rownames(summary.results$p.table) == 'Covariate1Covar1.2') - 1),, drop = FALSE]
        converged = TRUE
        beta = statistics[,1]
        SE = statistics[,2]
        P = statistics[,4]
        N = nrow(data.fit)
      }
      
      output = data.frame(
        species = colnames(taxa)[i],
        variable = phenotype_list[j],
        effect.level = if(!is.null(statistics)) sub("^trait", "", rownames(statistics)) else NA,
        beta = beta,
        SE = SE,
        P = P,
        N = N,
        converged = converged
      )
      
      output
    }
  result
}

## Example call
gam_function(metaphlan[metadata$Type == "Infant",1:5],
             metadata[metadata$Type=="Infant",],
             phenotype_list = c("Quant_phenotype_dynamic","Binary_phenotype_dynamic","Factor_phenotype_dynamic"))

