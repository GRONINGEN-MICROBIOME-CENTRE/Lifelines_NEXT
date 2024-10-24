############################ PERMUTATIONAL ANALYSIS TO LOOP OVER DISTANCES ###########
# Authors: Sergio Andreu Sanchez, Trishla Sinha 
# Last update: 1st of August, 2024 
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/figure_1")
load("Figure_1_D.RData")


# Example input 
#Which_type=  "mother:mother"
#Which_comparison = "P12-P28"
#n_permutations = 2000

Get_Pvalues_Permutations = function(dists_cons ,Which_type=  "mother:mother", Which_comparison = "P12-P28", n_permutations = 2000 ){
  dists_cons %>% as_tibble() -> dists_cons
  dists_cons %>% filter(Type1 == Type2) %>% filter(Type == Which_type) %>% 
    filter(TimeComp == Which_comparison) -> DataAnalysis
  
  Get_values = function(DataAnalysis){
    DataAnalysis %>% group_by(Comp) %>% summarise(Average_dist=mean(Distance) , SD = sd(Distance), N = n(), SE = SD/N ) ->
    DataAnalysis_sum
    Differences = tibble(
    Same_vs_related =  DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "SameIndv"] - DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "Related"],
    Same_vs_unrelated = DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "SameIndv"] - DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "Unrelated"],
    related_vs_unrelated = DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "Related"] - DataAnalysis_sum$Average_dist[DataAnalysis_sum$Comp == "Unrelated"]
    )
    Differences_SE = tibble(
      Same_vs_related = sqrt(    DataAnalysis_sum$SE[DataAnalysis_sum$Comp == "SameIndv"]^2 + DataAnalysis_sum$SE[DataAnalysis_sum$Comp == "Related"]^2 ),
      Same_vs_unrelated = sqrt(    DataAnalysis_sum$SE[DataAnalysis_sum$Comp == "SameIndv"]^2 + DataAnalysis_sum$SE[DataAnalysis_sum$Comp == "Unrelated"]^2),
      related_vs_unrelated = sqrt(    DataAnalysis_sum$SE[DataAnalysis_sum$Comp == "Related"]^2 + DataAnalysis_sum$SE[DataAnalysis_sum$Comp == "Unrelated"]^2) 
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


Get_Pvalues_Permutations(dists_cons, Which_type="infant:infant", Which_comparison = "M2-M3", n_permutations = 300 )


unique_types <- unique(dists_cons$Type)
unique_time_comparisons <- unique(dists_cons$TimeComp)

results <- list()
for (type in unique_types) {
  if (type == "infant:mother" ){ next }
  dists_cons %>% filter(Type==type) -> Comparisons_to_do
  Comparisons_to_do = unique(Comparisons_to_do$TimeComp)
  for (time_comp in Comparisons_to_do) {
    print(paste0(type, " ", time_comp))
    result <- Get_Pvalues_Permutations(dists_cons, Which_type = type, Which_comparison = time_comp, n_permutations = 2000)
    results[[paste(type, time_comp, sep = "_")]] <- result
  }
}

final_results <- bind_rows(results, .id = "Type_TimeComparison")
# calculated the t statistics and to calculate the p value we permutate the labels 2000 times 


final_results$FDR<-p.adjust (final_results$Pvalue)

write.table(final_results, "aitchison_distances_permutation_test.txt", sep = "\t", row.names = F)  
  
  