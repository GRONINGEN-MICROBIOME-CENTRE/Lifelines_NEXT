#### Functions used for associations with phenotypes ###
# Description: All functions used for associations of microbiome with phenotypes in manuscript: Mom Matters: Maternal influences on infant gut microbiome and health in 714 mother-infant pairs in the Lifelines NEXT cohort
# Authors: Trishla Sinha, Alexander Kurilshikov 

####### Load libraries #####

library(tidyverse)
library(data.table)
library(ggplot2)
library(dplyr)
library(lmerTest)
library(pheatmap)
library(vegan)
library(foreach)
library(mmrm)
library (mgcv) 

#### Functions ####

# Run inverse rank transformations on numeric data in variables 

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

# Centered log-ratio transformation 
do_clr_externalWeighting = function(interest_matrix, core_matrix){
  interest_matrix = interest_matrix + min(core_matrix[core_matrix>0])/2
  core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  
  # estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 1, gm_mean)
  tibble(ID = rownames(core_matrix) , Denominator = Gmean_core) %>% write_tsv("~/Desktop/Core_matrix.tsv")
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}

# Nullyfying zeros 

nullify_zeros = function(table,original){
  for (i in 1:ncol(table)){
    transformed = table[,i]
    original_record = original[,i]
    nonzero_records = transformed[which(original_record>0)]
    min_nonzero = min(nonzero_records)
    second_nonzero = min(nonzero_records[nonzero_records>min_nonzero])
    dist = abs(second_nonzero - min_nonzero)
    transformed[which(original_record == 0)] = min_nonzero - dist
    table[,i] = transformed
  }
  table
}

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

# Running mixed effects models correcting for dna concentration, read depth, batch and age #

mixed_models_cor_time <- function(metadata, ID, CLR_transformed_data, pheno_list) {
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
      Model0 = as.formula(paste( c(Bug2,  " ~ dna_conc + clean_reads_FQ_1 + BATCH_NUMBER+ exact_age_months_at_collection+(1|NEXT_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ dna_conc  + clean_reads_FQ_1 + BATCH_NUMBER+ exact_age_months_at_collection+",pheno2, "+ (1|NEXT_ID)"), collapse="" ))
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

mixed_models_without_time_correction <- function(metadata, ID, CLR_transformed_data, pheno_list) {
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
      Model0 = as.formula(paste( c(Bug2,  " ~ dna_conc + clean_reads_FQ_1 + BATCH_NUMBER+(1|NEXT_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ dna_conc  + clean_reads_FQ_1 + BATCH_NUMBER+",pheno2, "+ (1|NEXT_ID)"), collapse="" ))
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

# Run associations using with cross sectional phenotypes (which remain the same over time like delivery mode) using mmrm # 
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

.run_mmrm_single = function(outcome,phenotype,time,time_cat,NEXT_ID,covariates) {
  if(class(outcome) !="numeric") stop("outcome (bacteria or diversity matrix) should be numeric vector!")
  if(class(time)!="numeric") stop("time should be 'double'!")
  if(class(time_cat)!='factor') stop("time_cat should be 'factor'")
  if(class(NEXT_ID)!='factor') stop("NEXT_ID should be 'factor'")
  
  data.fit = try(data.frame(outcome = outcome,trait = phenotype,Time = time,Time_cat = time_cat,NEXT_ID = NEXT_ID,covariates))
  if(class(data.fit)[1]!="data.frame") stop("check that your input data is syncronized, i.e. all inputs have the same samples in the same order")
  data.fit = data.fit[complete.cases(data.fit),]
  
  if(length(table(as.character(data.fit$Time_cat))) > 3) {
    formula.null.us = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + Time3 + trait + us(Time_cat|NEXT_ID)")
    formula.run.us = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + Time3 + trait + trait:Time1 + us(Time_cat|NEXT_ID)")
    formula.null.cs = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + Time3 + trait + cs(Time_cat|NEXT_ID)")
    formula.run.cs = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + Time3 + trait + trait:Time1 + cs(Time_cat|NEXT_ID)")
    
    data.fit2 = data.frame(data.fit, poly(data.fit$Time,3))
    colnames(data.fit2)[(ncol(data.fit2)-2): ncol(data.fit2)] = c("Time1","Time2","Time3")
  } else if (length(table(as.character(data.fit$Time_cat))) > 2) {
    formula.null.us = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + trait + us(Time_cat|NEXT_ID)")
    formula.run.us = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + trait + trait:Time1 + us(Time_cat|NEXT_ID)")
    formula.null.cs = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + trait + cs(Time_cat|NEXT_ID)")
    formula.run.cs = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + trait + trait:Time1 + cs(Time_cat|NEXT_ID)")
    
    data.fit2 = data.frame(data.fit, poly(data.fit$Time,2))
    colnames(data.fit2)[(ncol(data.fit2)-1): ncol(data.fit2)] = c("Time1","Time2")
  } else if (length(table(as.character(data.fit$Time_cat))) > 1){
    formula.null.us = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + trait + us(Time_cat|NEXT_ID)")
    formula.run.us = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + trait + trait:Time1 + us(Time_cat|NEXT_ID)")
    formula.null.cs = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + trait + cs(Time_cat|NEXT_ID)")
    formula.run.cs = paste("outcome ~ ", paste(colnames(covariates),collapse = "+"),
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
                               model_sucess = "Failure",
                               Covar.type = model.type,
                               outcome = NA,
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
                                 model_sucess = "Failure",
                                 Covar.type = model.type,
                                 outcome = NA,
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
                               model_sucess = "Success",
                               Covar.type = model.type,
                               outcome = NA,
                               trait = NA,
                               N = nrow(data.fit2),
                               levels = sub("trait",
                                            "",
                                            rownames(summary.nullModel.reml$coef)[
                                              grep("trait",rownames(summary.nullModel.reml$coef))]),
                               summary.nullModel.reml$coef[
                                 grep("trait",rownames(summary.nullModel.reml$coef)),,drop =F])
    timeInt.results = data.frame(row.names = NULL,
                                 model_sucess = "Success",
                                 Covar.type = model.type,
                                 outcome = NA,
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

run_mmrm_analysis = function(outcome, phenotypes, time, time_cat, NEXT_ID,covariates){
  output_names = c("trait","time.trait")
  output = sapply(output_names, function(x) NULL)
  for (i in 1:ncol(outcome)) {
    print (i)
    for(j in 1:ncol(phenotypes)){
      print(j)
      singleAssoc = .run_mmrm_single(outcome[,i],phenotypes[,j],time = time,time_cat = time_cat, NEXT_ID = NEXT_ID,covariates = covariates)
      singleAssoc$trait[,"outcome"] = colnames(outcome)[i]
      singleAssoc$trait[,"trait"] = colnames(phenotypes)[j]
      singleAssoc$time.trait[,"outcome"] = colnames(outcome)[i]
      singleAssoc$time.trait[,"trait"] = colnames(phenotypes)[j]
      output$trait = rbind(output$trait,singleAssoc$trait)
      output$time.trait = rbind(output$time.trait,singleAssoc$time.trait)
    }
    
  }
  colnames(output$trait)[8] = "SE"
  colnames(output$time.trait)[8] = "SE"
  colnames(output$time.trait)[11] = "p"
  colnames(output$trait)[11] = "p"
  
  
  
  output
  
}

# Run associations with dynamic phenotypes using GAMs #

gam_function <- function(outcome, phenotypes, metadata_columns, covariate_columns) {
  # Remove metadata columns from phenotypes to get the variables to correlate
  phenotypes2correlate <- phenotypes[, !(colnames(phenotypes) %in% metadata_columns)]
  phenotypes2correlate[sapply(phenotypes2correlate, is.character)] <- lapply(phenotypes2correlate[sapply(phenotypes2correlate, is.character)], as.factor)
  
  # Extract covariates dynamically
  covariates <- phenotypes[, covariate_columns, drop = FALSE]
  
  # Inverse rank normalization of phenotypes and covariates
  phenotypes2correlate <- run_invrank_dataFrame(phenotypes2correlate)
  covariates <- run_invrank_dataFrame(covariates)
  
  result = foreach(i = 1:ncol(outcome), .combine = rbind) %:% 
    foreach(j = 1:ncol(phenotypes2correlate), .combine = rbind) %do% {
      outcome_name = colnames(outcome)[i]
      phenotype_name = colnames(phenotypes2correlate)[j]
      print(paste("Processing outcome:", outcome_name, "- Phenotype:", phenotype_name))
      
      data.fit <- data.frame(
        outcome= outcome[, i],
        trait = phenotypes2correlate[, j],
        Time = phenotypes$exact_age_months_at_collection,
        NEXT_ID = phenotypes$NEXT_ID
      )
      
      # Add covariates
      for (covariate in covariate_columns) {
        data.fit[[covariate]] <- covariates[[covariate]]
      }
      
      data.fit$NEXT_ID <- factor(data.fit$NEXT_ID)
      data.fit <- data.fit[complete.cases(data.fit),]
      
      # Run model
      formula_str <- paste("outcome ~ trait +", paste(covariate_columns, collapse = " + "), "+ s(Time) + s(NEXT_ID, bs = 're')")
      mod_gam1 <- tryCatch({
        bam(as.formula(formula_str), data = data.fit, discrete = TRUE)
      }, warning = function(w) {
        return(NULL)  
      }, error = function(e) {
        return(NULL)  
      })
      
      if (is.null(mod_gam1)) {
        statistics <- NULL
        converged <- FALSE
        beta <- NA
        SE <- NA
        P <- NA
        N <- nrow(data.fit)
      } else {
        summary.results <- summary(mod_gam1)
        statistics <- summary.results$p.table[2:(which(rownames(summary.results$p.table) == covariate_columns[1]) - 1), , drop = FALSE]
        converged <- TRUE
        beta <- statistics[, 1]
        SE <- statistics[, 2]
        P <- statistics[, 4]
        N <- nrow(data.fit)
      }
      
      output <- data.frame(
        outcome = colnames(outcome)[i],
        trait = colnames(phenotypes2correlate)[j],
        effect.level = if (!is.null(statistics)) sub("^trait", "", rownames(statistics)) else NA,
        Estimate = beta,
        SE = SE,
        p = P,
        N = N,
        converged = converged
      )
      
      output
    }
  result
}

linear_model_taxa_simple_per_timepoint <- function(metadata, ID, CLR_transformed_data, pheno_list, tp = "W2" ) {
  metadata %>% filter(Timepoint_categorical == tp) -> metadata
  
  df<-merge(metadata, CLR_transformed_data, by='row.names')
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
      
      Model2 = as.formula(paste( c(Bug2,  " ~ clean_reads_FQ_1 + dna_conc + BATCH_NUMBER+",pheno2), collapse="" ))
      lm(Model2, df_pheno) -> resultmodel2
      
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = Summ_simple$`Pr(>|t|)`, Outcome =Bug, Timepoint=tp, Correction="techincal") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  
  return(p)
}

linear_model_taxa_cor_feeding_kcal_per_timepoint <- function(metadata, ID, CLR_transformed_data, pheno_list, tp = "W2" ) {
  metadata %>% filter(Timepoint_categorical == tp) -> metadata
  
  df<-merge(metadata, CLR_transformed_data, by='row.names')
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
      
      Model2 = as.formula(paste( c(Bug2,  " ~ clean_reads_FQ_1 + dna_conc + BATCH_NUMBER+ infant_ffq_energy_kcalday_item_sum_recal+ infant_ffq_feeding_mode_complex+",pheno2), collapse="" ))
      lm(Model2, df_pheno) -> resultmodel2
      
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = Summ_simple$`Pr(>|t|)`, Outcome =Bug,  Timepoint=tp, Correction="techincal_feeding_kcal") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  
  return(p)
}
