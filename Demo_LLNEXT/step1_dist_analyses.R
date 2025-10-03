############################ Distance matrix-based analyses ###########
# Authors: Sergio Andreu Sanchez, Trishla Sinha, Alex Kurilshikov 
# Last update: 21st of Nov, 2024 
# This script is dedicated to per-timepoint PERMANOVA, sample-sample distance 
# comparisons, and species-time associations

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(vegan)
library(foreach)

metaphlan = read.table("simulated_metaphlan.txt",header=T)
metadata = read.table("simulated_metadata.txt",header=T)


# 1.1 Sample-sample distances & tests -----------------------------------------
# Adopted from https://github.com/GRONINGEN-MICROBIOME-CENTRE/Lifelines_NEXT/Gut_Microbiome/Time_vs_microbiome_composition/Lifelines_NEXT_permut_analysis.R

Get_Pvalues_Permutations = function(dists_cons ,Which_type=  "mother:mother", Which_comparison = "P12-P28", n_permutations = 2000 ){
  dists_cons %>% as_tibble() -> dists_cons
  dists_cons %>% filter(Type1 == Type2) %>% filter(Type == Which_type) %>% 
    filter(TimeComp == Which_comparison) -> DataAnalysis
  
  Get_values = function(DataAnalysis){
    DataAnalysis %>% group_by(Comp) %>% summarise(Average_dist=mean(Distance) , Var = var(Distance), N = n()) ->
      DataAnalysis_sum
    Differences = tibble(
      Same_vs_related =  DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "SameIndv"] - DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "Related"],
      Same_vs_unrelated = DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "SameIndv"] - DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "Unrelated"],
      related_vs_unrelated = DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "Related"] - DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "Unrelated"]
    )
    Differences_SE = tibble(
      Same_vs_related = sqrt(DataAnalysis_sum$Var[DataAnalysis_sum$Comp == "SameIndv"] + DataAnalysis_sum$Var[DataAnalysis_sum$Comp == "Related"] ) / 
        sqrt (DataAnalysis_sum$N[DataAnalysis_sum$Comp == "SameIndv"] + DataAnalysis_sum$N[DataAnalysis_sum$Comp == "Related"] ),
      Same_vs_unrelated = sqrt(DataAnalysis_sum$Var[DataAnalysis_sum$Comp == "SameIndv"] + DataAnalysis_sum$Var[DataAnalysis_sum$Comp == "Unrelated"]) /
        sqrt(DataAnalysis_sum$N[DataAnalysis_sum$Comp == "SameIndv"] + DataAnalysis_sum$N[DataAnalysis_sum$Comp == "Unrelated"]),
      related_vs_unrelated = sqrt(    DataAnalysis_sum$Var[DataAnalysis_sum$Comp == "Related"] + DataAnalysis_sum$Var[DataAnalysis_sum$Comp == "Unrelated"]) /
        sqrt(    DataAnalysis_sum$N[DataAnalysis_sum$Comp == "Related"] + DataAnalysis_sum$N[DataAnalysis_sum$Comp == "Unrelated"])
        
    )
    Differences_Tvalue = Differences/Differences_SE
    
    return(Differences_Tvalue)
  }
  Get_Pvalues = function(DataAnalysis){
    Comps =  levels(DataAnalysis$Comp)
    N_comp = length(Comps)
    Res = tibble()
    for (Level_n in 1:N_comp ){
      Level1 = Comps[Level_n]
      if (Level1 ==  N_comp) { break}
      for (Level_n2 in (Level_n+1):N_comp ){
        Level2 = Comps[Level_n2]
        if (is.na(Level2) | Level_n == Level_n2 ){ next }
        #print(paste0(Level1, "_vs_",Level2 ))
        wilcox.test( filter(DataAnalysis, Comp == Level1)$Distance,
                     filter(DataAnalysis, Comp == Level2)$Distance,  alternative = "less") -> Test
        P_test = Test$p.value
        Res %>% rbind(tibble( Test = paste0(Level1,":", Level2), P = P_test  )) -> Res
      } 
    }
    return(Res)
  }
  Get_values(DataAnalysis) -> Original
  #Get_Pvalues(DataAnalysis) -> Original_p
  
  #Do permutations and repeat
  Permutations = tibble()
  for (i in 1:n_permutations){
    print (paste0(i, "/", n_permutations) )
    DataAnalysis_perm = DataAnalysis
    sample(DataAnalysis$Comp, dim(DataAnalysis)[1], replace = F ) -> DataAnalysis_perm$Comp
    Get_values(DataAnalysis_perm) -> Permuted
    #Get_Pvalues(DataAnalysis_perm) -> Permuted
    Permutations %>% rbind(. , Permuted) -> Permutations 
  }
  #Get P-values
  
  Original %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column("Comparison")  %>% rename(Difference = V1) %>%
    as_tibble() -> Original_2
  sapply(Original_2$Comparison, function(Comparison_n){ 
    Original_value = filter(Original_2, Comparison == Comparison_n)$Difference
    Permutations[[Comparison_n]] %>% as_vector() -> Permut
    sum(Permut <= Original_value)/n_permutations %>% return()
  }  ) -> Res
  
  
  #sapply(Original_p$Test, function(Comparison_n){ 
  #  Original_value = filter(Original_p, Test == Comparison_n)$P
  #  Permutations %>% filter(Test == Comparison_n) %>% select(P) %>% as_vector() -> Permut
  #  sum(Permut <= Original_value)/n_permutations %>% return()
  #}  )
  tibble( Comparison = names(Res), Pvalue = as.vector(Res), Type = Which_type, TimeComparison = Which_comparison, N_perm = n_permutations ) %>% return()  
}

# @Sergio: make example dists_cons

Get_Pvalues_Permutations(dists_cons)


# 1.2 PERMANOVA with phenotypes at each timepoint separately ------------------
# Adopted from https://github.com/GRONINGEN-MICROBIOME-CENTRE/Lifelines_NEXT/Gut_Microbiome/Time_vs_microbiome_composition/LIFELINES_NEXT_PERMANOVA.R

Call_adonis = function(phenotypes,Phenotype_list, Distance, perm=1000, cores=2){
  Distance = as.matrix(Distance)
  adonis_results<- data.frame(matrix(ncol =5, nrow= ncol(phenotypes)))  
  colnames(adonis_results) <- c("Phenotype","Df", "F", "R2", "p-value")
  Phenotype_list = Phenotype_list[!is.na(Phenotype_list)]
  adon<-foreach(i=1:length(Phenotype_list),.combine=rbind)%do%{
    #This is the phentype vector
    Pheno = Phenotype_list[i]
    print(Pheno)
    as_vector(phenotypes[Pheno]) %>% as.vector() -> Values
    #Filter NAs
    r<-row.names(phenotypes)[!is.na(Values)]
    phenotypes[r,] -> phenos2
    #Match with Distance
    distmat_cleaned <- as.dist(Distance[r,r])
    #Run adonis
    FR = formula( paste0("distmat_cleaned ~ Covariate1 + Covariate2 +", Pheno) )
    ad1<-adonis2(FR , phenos2, permutations=perm,parallel=cores, na.action = na.fail,by='terms')
    
    adonis_results[i,] <- c(Pheno,ad1$Df[which(rownames(ad1)==Pheno)], 
                            ad1$F[which(rownames(ad1)==Pheno)],
                            ad1$R2[which(rownames(ad1)==Pheno)],
                            ad1$'Pr(>F)'[which(rownames(ad1)==Pheno)])
  }
  
  adonis_results %>% drop_na() %>% return()
}
# Example: infant, W2

Call_adonis(metadata[metadata$Timepoint=="W2",],
            c("Factor_phenotype_dynamic","Quant_phenotype_dynamic","Binary_phenotype_dynamic"),
            vegdist(metaphlan[metadata$Timepoint=="W2",],method = "bray"))






