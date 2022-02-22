setwd("~/Desktop/Projects_2018/Baby_pilot/01.Data_analysis/03.SCRIPTS/Analysis_for_paper/16.Phylum_test/")

##### LIBRARIES #####
library(purrr)
library(dplyr)
library(lme4)
library(lmerTest)
library(RLRsim)
library(tidyverse)
##### FUNCTIONS #####
taxonomy_level_extraction <- function(metaphlan_output, level){
  taxas=rownames(as.data.frame(metaphlan_output))
  list_taxa=list()
  for (i in taxas){
    if (length (unlist(strsplit(i, "\\|"))) == level){
      list_taxa=c(list_taxa,i)
    }
  }
  taxa_table=as.data.frame(metaphlan_output)[unlist(list_taxa),]
  print(taxa_table)
}

# CLR Normalization function 
do_clr_externalWeighting = function(interest_matrix, core_matrix) {
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  #estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 1, gm_mean)
  
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}
##### INPUT DATA #####
samples_metadata <- read.table('../INPUT/Microbiome_metadata_phenos.txt', sep='\t', header=T)
samples_metadata$Type <- as.factor(samples_metadata$Type)
samples_metadata$PSEUDOID_number <- as.factor(samples_metadata$PSEUDOID_number)
samples_metadata$Timepoint_continuous <- NA
samples_metadata[samples_metadata$Timepoint=="P3",]$Timepoint_continuous <- 0
samples_metadata[samples_metadata$Timepoint=="P7",]$Timepoint_continuous <- 4
samples_metadata[samples_metadata$Timepoint=="B",]$Timepoint_continuous <- 6
samples_metadata[samples_metadata$Timepoint=="M1",]$Timepoint_continuous <- 7
samples_metadata[samples_metadata$Timepoint=="M2",]$Timepoint_continuous <- 8
samples_metadata[samples_metadata$Timepoint=="M3",]$Timepoint_continuous <- 9
for_model <- samples_metadata[,c("Type", "PSEUDOID_number", "Timepoint_continuous")]

metaphlan_3<-read.table("~/Resilio Sync/NEXT Pilot/01. Data/02.MICROBIOME/METAPHLAN3/merged_abundance_table_metaphlan_TS_FEB_2021.txt",sep='\t', header=T, row.names = 1)
metaphlan_3 <- metaphlan_3[-1,-1]
colnames(metaphlan_3)<- substr(colnames(metaphlan_3), 0, 14) 
metaphlan_3 <- metaphlan_3[,colnames(metaphlan_3) %in% row.names(samples_metadata)]

phylum <- taxonomy_level_extraction(metaphlan_3, 2)
phylum_filtered <- phylum[(apply(phylum, 1, function(x) sum(x > 5) > 0.01 * ncol(phylum))),]
phylum_clr <- do_clr_externalWeighting(phylum_filtered, phylum_filtered)
phylum_clr <- as.data.frame(t(phylum_clr))
names(phylum_clr) <- substr(colnames(phylum_clr), 16, 100)

phylum_clr <- merge(phylum_clr, for_model, by='row.names')
row.names(phylum_clr) <- phylum_clr$Row.names
phylum_clr$Row.names <- NULL

##### PROCESSING #####
Overall_result =tibble()
for (Bug2 in names(phylum_clr[,c(1:5)])) {
  ModelR = as.formula(paste( c(Bug2,  " ~ Type + Timepoint_continuous  + (1|PSEUDOID_number)"), collapse="" ))
  Model0 = as.formula(paste( c(Bug2,  " ~ Type + Timepoint_continuous"), collapse="" ))
  resultmodelR <-lmer(ModelR, phylum_clr, REML = F)
  resultmodel0<-lm(Model0, phylum_clr)
  if(data.frame(VarCorr(resultmodelR),comp="Variance")[1,4] == 0){ N = 2
  } else {
    exactLRT(m = resultmodelR, m0 = resultmodel0 ) -> P_comparison # p values <0.05 in favour of the complex model (modelR)
    if (P_comparison$p.value < 0.05){ N = 1
    } else{ N = 2}
  }
  if (N == 1){
    M = "Mixed"
    as.data.frame(summary(resultmodelR)$coefficients)[2,5]->p_simp
    as.data.frame(summary(resultmodelR)$coefficients)[2,1:3] -> Summ_simple
  } else{ 
    M ="Fixed"
    base_model = resultmodel0
    as.data.frame(anova(resultmodel0, base_model))[2,6]->p_simp
    as.data.frame(summary(resultmodel0)$coefficients)[4,1:3] -> Summ_simple
  }
  
  
  Summ_simple %>%  as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug = Bug2) -> temp_output
  #Summ_inter  %>% as_tibble() %>% mutate(P = p_int, Model_choice = M, Bug =Bug2, Model ="interaction") -> temp_output2
  #rbind(temp_output,temp_output2) -> temp_output
  rbind(Overall_result, temp_output) -> Overall_result
  
}

p =as.data.frame(Overall_result)
results=p
results$FDR=p.adjust(results$P, method = "BH")


boxplot(phylum_clr$Actinobacteria ~ phylum_clr$Type)
boxplot(phylum_clr$Proteobacteria ~ phylum_clr$Type)

boxplot(phylum_clr$Bacteroidetes ~ phylum_clr$Type)
boxplot(phylum_clr$Firmicutes ~ phylum_clr$Type)



