############################ PGLMM ANALYSIS  ###########
# Authors: Alex Kurilshikov
# Last update: 26st of Apr, 2026 
# PGLMM analysis of association of phenotypes to strain phylogeny
# last revision: switched to MDMR-based analysis instead of brms/MCMCGLMM


# Load libraries ----------------------------------------------------------
library(ape)
library(phangorn)
library(foreach)
library(stringr)
library(vegan)
library(parallel)
library(MDMR)


# Supportive functions ---------------------------------------------------------------
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


pheno_filtering = function(tree_name, phenotype){
  tree_dist = all_distances[[tree_name]]
  input_data = data.frame(covariates, pheno = all_phenos[,phenotype])
  input_data = input_data[complete.cases(input_data),]
  #syncing tree and data
  same_samples = intersect(input_data$NG_ID,rownames(tree_dist))
  tree_dist_Input = tree_dist[same_samples,same_samples]
  input_data_Input = input_data[match(same_samples,input_data$NG_ID),]
  if(is.numeric(input_data_Input$pheno)) input_data_Input$pheno = invrank(input_data_Input$pheno)
  
  sort(table(input_data_Input$pheno), decreasing = T) [1] / sum(!is.na(input_data_Input$pheno))
  
}


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
      #print(paste(colnames(cleaned_phenotypes)[i], ": accepted"))
    }
  } 
  cleaned_phenotypes[, discarded_phenotypes] <- NULL #drop non-selected phenotypes
  return(cleaned_phenotypes)
}

# main functions for analysis ---------------------------------------------
perform_mdmr_analysis = function(tree_name, phenotype) {
  tree_dist = cophenetic.phylo(trees[[tree_name]])
  input_data = data.frame(covars, pheno = all_phenos[,phenotype])
  input_data = input_data[complete.cases(input_data),]
  #syncing tree and data
  same_samples = intersect(input_data$NG_ID,rownames(tree_dist))
  tree_dist_Input = tree_dist[same_samples,same_samples]
  input_data_Input = input_data[match(same_samples,input_data$NG_ID),]
  #pheno_filtering
  pheno_filter = pheno_filtering(tree_name,phenotype)
  N = nrow(input_data_Input)
  
  if(N < 100 | pheno_filter >=0.75 ){
    s1 = data.frame(bacterium = tree_name,phenotype = phenotype , pheno_redundancy = pheno_filter,N = N,
                    Statistic_TimeFixed=NA,
                    Numer.DF_TimeFixed=NA,
                    p.value_TimeFixed=NA,
                    Statistic_TimeRandom=NA,
                    Numer.DF_TimeRandom=NA,
                    p.value_TimeRandom=NA,
                    perm.P = NA)
  } else {
    if(is.numeric(input_data_Input$pheno)) input_data_Input$pheno = invrank(input_data_Input$pheno)
    if(length(unique(input_data_Input$Timepoint))>1){
      model1_mixed_acc <- mixed.mdmr( ~ clean_reads_FQ_1 + dna_conc + Timepoint +  pheno + (1|NEXT_ID),use.ssd=0.99,D = tree_dist_Input,data = input_data_Input)
      summary1 = summary(model1_mixed_acc)
    } else {
      model1_mixed_acc <- mdmr(input_data_Input[,c("clean_reads_FQ_1","dna_conc", "pheno")],D = tree_dist_Input,perm.p = F)
      summary1 = summary(model1_mixed_acc)[,-3]
      colnames(summary1) = c("Statistic","Numer.DF","p.value")
      summary2 = summary(model1_mixed_acc)[,-3]
      colnames(summary2) = c("Statistic","Numer.DF","p.value")
      
    }
    #  model1_mixed_acc2 <- mixed.mdmr( ~ clean_reads_FQ_1 + dna_conc +   pheno + (1|Timepoint),use.ssd=0.99,D = tree_dist_Input,data = input_data_Input)
    
    
    #summary2 = summary(model1_mixed_acc2)
    s1 = data.frame(bacterium = tree_name,phenotype = phenotype , pheno_redundancy = pheno_filter,N = nrow(input_data_Input),
                    Statistic_TimeFixed=summary1["pheno",1],
                    Numer.DF_TimeFixed=summary1["pheno",2],
                    p.value_TimeFixed=summary1["pheno",3],
                    Statistic_TimeRandom=summary2["pheno",1],
                    Numer.DF_TimeRandom=summary2["pheno",2],
                    p.value_TimeRandom=summary2["pheno",3],
                    perm.P = NA)
  }
  #if(perm == T){
  #  perms = foreach(i = 1:nperm,.combine = c) %do% {
  #    print(i)
  #    input_perm = input_data_Input[sample(1:nrow(input_data_Input),replace = F),]
  #    rownames(input_perm) = rownames(input_data_Input)
  #    model1_mixed_perm <- mixed.mdmr( ~ clean_reads_FQ_1 + dna_conc + Timepoint_categorical +  pheno + (1|NEXT_ID),use.ssd=0.99,D = tree_dist_Input,data = input_perm)
  #    summary(model1_mixed_perm)["pheno",1]
  #  }
  #  s1$perm.P = sum(s1$Statistic <=perms) / nperm
  #}
  
  s1
}

perform_mdmr_analysis_permut = function(tree_name, phenotype,nperm=1000) {
  tree_dist = all_distances[[tree_name]]
  input_data = data.frame(covariates, pheno = all_phenos[,phenotype])
  input_data = input_data[complete.cases(input_data),]
  #syncing tree and data
  same_samples = intersect(input_data$NG_ID,rownames(tree_dist))
  tree_dist_Input = tree_dist[same_samples,same_samples]
  input_data_Input = input_data[match(same_samples,input_data$NG_ID),]
  #pheno_filtering
  pheno_filter = pheno_filtering(tree_name,phenotype)
  N = nrow(input_data_Input)
  
  perms = rep(NA,nperm)
  for(i in 1:nperm) {
    print(i)
    input_data_Input_perm = input_data_Input[sample(1:nrow(input_data_Input)),]
    rownames(input_data_Input_perm) = rownames(input_data_Input)
    
    if(length(unique(input_data_Input$Timepoint))>1){
      model1_mixed_acc <- mixed.mdmr( ~ clean_reads_FQ_1 + dna_conc + Timepoint +  pheno + (1|NEXT_ID),use.ssd=0.99,D = tree_dist_Input,data = input_data_Input_perm)
      summary1 = summary(model1_mixed_acc)
    } else {
      model1_mixed_acc <- mdmr(input_data_Input[,c("clean_reads_FQ_1","dna_conc", "pheno")],D = input_data_Input_perm,perm.p = F)
      summary1 = summary(model1_mixed_acc)[,-3]
      colnames(summary1) = c("Statistic","Numer.DF","p.value")
      
    }
    perms[i] = summary1["pheno",3]
    
  }
  perms
  
}


# Data loading and syncing ------------------------------------------------


metadata = read.table("simulated_metadata.txt",header=T)
align = read.dna("alignment_example.aln",format = "fasta")

#generating ultrametric GTR tree
trees = list(tree1 = pml_bb(align,model = "GTR",method = "ultrametric",
              rearrangement="NNI")$tree)


# preparing dataset for loading. 
metadata_infant = metadata[metadata$Type=="Infant",]
metadata_infant$NG_ID = rownames(metadata_infant)
covars = metadata_infant[,c("Covariate1","Covariate2","NG_ID","Timepoint","Family")]


# main call ---------------------------------------------------------------

results_mdmr = do.call(rbind,lapply(c("Binary_phenotype_dynamic","Quant_phenotype_dynamic","Factor_phenotype_dynamic"),
                      \(x) {perform_mdmr_analysis(
                        "tree1",
                        x
                      )
                        }))
results_mdmr$FDR.fixed = p.adjust(results_mdmr$p.value_TimeFixed,method = "BH")
results_mdmr$FDR.random = p.adjust(results_mdmr$p.value_TimeRandom,method = "BH")

perms = lapply(1:nrow(results_mdmr),\(x)
  suppressMessages(suppressWarnings(
    perform_mdmr_analysis_permut(results_mdmr[x,1],results_mdmr[x,2],nperm = 20000))))

for(i in 1:nrow(results_mdmr)) {
  results_mdmr$perm.P[i] = (perms[[i]] < results_mdmr$p.value_TimeFixed[i])/20000
}
