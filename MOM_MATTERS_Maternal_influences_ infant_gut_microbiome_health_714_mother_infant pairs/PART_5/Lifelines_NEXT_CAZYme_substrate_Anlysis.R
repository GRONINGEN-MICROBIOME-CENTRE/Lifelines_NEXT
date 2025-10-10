########################## CAZYMES SUBSTRATE ENRICHMENT ###########################################
### AUTHOR:  TRISHLA SINHA, SERGIO ANDREU SANCHEZ
### ORIGINAL SCRIPT: 9th June, 2025
### LAST UPDATE: 1st October , 2025 

library(clusterProfiler)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)


Annotation_caz = read_csv('~/Desktop/LLNEXT/Analysis/cazymes/AnnotationCazymes_ST1_Durcamon_preprint.csv', skip = 1)
 
Annotation_caz = Annotation_caz %>% 
  select(Family, Subfamily, ORIGIN, FUNCTION_IN_ORIGIN, FUNCTION_AT_DESTINATION_1, 
         FUNCTION_AT_DESTINATION_2, FUNCTION_AT_DESTINATION_3, Glycan_annotation)
Columns_association = colnames(Annotation_caz %>% select(-c(Family, Subfamily) ))
List_substrate_link = list()
for  (C in Columns_association){
  List_substrate_link[[C]] = Annotation_caz %>% select('Subfamily', C) %>%  separate_rows( !!sym(C) , sep = ",") %>%
    distinct()
  
}


# Create function that does enrichment for a given 'substrate level'
asso_del_mode<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/cazymes/cazyme_cross_phenotypes_results_all_infant_29_09_2025.txt")
asso_feed_mode<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/cazymes/cazyme_dynamic_phenotypes_results_all_30_09_2025.txt")
asso_early_late<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes/association_early_late_timepoint_infants.txt")



do_enrichment_phenotypes = function(Substrate_level, Association_table, Pheno = 'birth_deliverybirthcard_mode_binary', Positive = T){
  if (Positive == T){
    Set = filter(Association_table, FDR<0.05 & Estimate_cor_base > 0 & trait == Pheno )$outcome 
  }else {
  Set = filter(Association_table, FDR<0.05 & Estimate_cor_base < 0 & trait == Pheno )$outcome 
  }  
  Test = clusterProfiler::enricher( 
    gene=Set ,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 2, 
    universe = Association_table$outcome %>% unique(),
    TERM2GENE = List_substrate_link[[Substrate_level]] %>%  select(Substrate_level, 'Subfamily')
    )
  res_tbl_pheno <- as_tibble(as.data.frame(Test))
  print(res_tbl_pheno)
}


do_enrichment_timepoint = function(Substrate_level, Association_table, Pheno = Pheno, Positive = T){
  if (Positive == T){
    Set = filter(Association_table, FDR<0.05 & Estimate > 0, Pheno == Pheno )$Bug 
  }else {
    Set = filter(Association_table, FDR<0.05 & Estimate < 0, Pheno == Pheno )$Bug
  }  
  Test = clusterProfiler::enricher( 
    gene=Set ,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 2, 
    universe = Association_table$outcome %>% unique(),
    TERM2GENE = List_substrate_link[[Substrate_level]] %>%  select(Substrate_level, 'Subfamily')
    
  )
  res_tbl_timepoint <- as_tibble(as.data.frame(Test))
  print(res_tbl_timepoint)
}

do_enrichment_feed = function(Substrate_level, Association_table, Pheno = '', Positive = T){
  if (Positive == T){
    Set = filter(Association_table, FDR<0.05 & Estimate_cor_base > 0 & trait == Pheno )$outcome 
  }else {
    Set = filter(Association_table, FDR<0.05 & Estimate_cor_base < 0 & trait == Pheno )$outcome 
  }  
  Test = clusterProfiler::enricher( 
    gene=Set ,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 2, 
    universe = Association_table$outcome %>% unique(),
    TERM2GENE = List_substrate_link[[Substrate_level]] %>%  select(Substrate_level, 'Subfamily')
  )
  res_tbl_pheno <- as_tibble(as.data.frame(Test))
  print(res_tbl_pheno)
}

#  Testing enrichment associated with mode of delivery 

# Enriched in VG
do_enrichment_phenotypes('ORIGIN',asso_del_mode,  'birth_deliverybirthcard_mode_binary', Positive = T)
do_enrichment_phenotypes('FUNCTION_AT_DESTINATION_1',asso_del_mode,  'birth_deliverybirthcard_mode_binary', Positive=T)
do_enrichment_phenotypes('FUNCTION_AT_DESTINATION_2',asso_del_mode,  'birth_deliverybirthcard_mode_binary', Positive=T)
do_enrichment_phenotypes('FUNCTION_AT_DESTINATION_3',asso_del_mode,  'birth_deliverybirthcard_mode_binary', Positive=T)
do_enrichment_phenotypes('Glycan_annotation',asso_del_mode,  'birth_deliverybirthcard_mode_binary', Positive=T)
# Enriched in CS
do_enrichment_phenotypes('ORIGIN',asso_del_mode,  'birth_deliverybirthcard_mode_binary', Positive = F)
do_enrichment_phenotypes('FUNCTION_AT_DESTINATION_1',asso_del_mode,  'birth_deliverybirthcard_mode_binary', Positive=F)
do_enrichment_phenotypes('FUNCTION_AT_DESTINATION_2',asso_del_mode,  'birth_deliverybirthcard_mode_binary', Positive=F)
do_enrichment_phenotypes('FUNCTION_AT_DESTINATION_3',asso_del_mode,  'birth_deliverybirthcard_mode_binary', Positive=F)
do_enrichment_phenotypes('Glycan_annotation',asso_del_mode,  'birth_deliverybirthcard_mode_binary', Positive=F)


#  Testing enrichment associated with timepoint in infants 
#LATE
do_enrichment_timepoint('ORIGIN',asso_early_late,  'Timepoint_early_late', Positive = T)
do_enrichment_timepoint('FUNCTION_AT_DESTINATION_1',asso_early_late,  'Timepoint_early_late', Positive=T)
do_enrichment_timepoint('FUNCTION_AT_DESTINATION_2',asso_early_late,  'Timepoint_early_late', Positive=T)
do_enrichment_timepoint('FUNCTION_AT_DESTINATION_3',asso_early_late,  'Timepoint_early_late', Positive=T)
do_enrichment_timepoint('Glycan_annotation',asso_early_late,  'Timepoint_early_late', Positive=T)

##EARLY 
do_enrichment_timepoint('ORIGIN',asso_early_late,  'Timepoint_early_late', Positive = F)
do_enrichment_timepoint('FUNCTION_AT_DESTINATION_1',asso_early_late,  'Timepoint_early_late', Positive=F)
do_enrichment_timepoint('FUNCTION_AT_DESTINATION_2',asso_early_late,  'Timepoint_early_late', Positive=F)
do_enrichment_timepoint('FUNCTION_AT_DESTINATION_3',asso_early_late,  'Timepoint_early_late', Positive=F)
do_enrichment_timepoint('Glycan_annotation',asso_early_late,  'Timepoint_early_late', Positive=F)

#  Testing enrichment associated with ever never breastfed in infants  
#FF
do_enrichment_feed ('ORIGIN',asso_feed_mode,  'infant_ffq_feeding_mode_simple', Positive = T)
do_enrichment_feed ('FUNCTION_AT_DESTINATION_1',asso_feed_mode,  'infant_ffq_feeding_mode_simple', Positive=T)
do_enrichment_feed ('FUNCTION_AT_DESTINATION_2',asso_feed_mode,  'infant_ffq_feeding_mode_simple', Positive=T)
do_enrichment_feed ('FUNCTION_AT_DESTINATION_3',asso_feed_mode,  'infant_ffq_feeding_mode_simple', Positive=T)
do_enrichment_feed ('Glycan_annotation',asso_feed_mode,  'infant_ffq_feeding_mode_simple', Positive=T)

# BF
do_enrichment_feed ('ORIGIN',asso_feed_mode,  'infant_ffq_feeding_mode_simple', Positive = F)
do_enrichment_feed ('FUNCTION_AT_DESTINATION_1',asso_feed_mode,  'infant_ffq_feeding_mode_simple', Positive=F)
do_enrichment_feed ('FUNCTION_AT_DESTINATION_2',asso_feed_mode,  'infant_ffq_feeding_mode_simple', Positive=F)
do_enrichment_feed ('FUNCTION_AT_DESTINATION_3',asso_feed_mode,  'infant_ffq_feeding_mode_simple', Positive=F)
do_enrichment_feed ('Glycan_annotation',asso_feed_mode,  'infant_ffq_feeding_mode_simple', Positive=F)


################# Plotting ##########################

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/cazymes")
plot_enrichment <- function(enrichment_df, plot_title, file_name) {
  enrichment_df <- enrichment_df[enrichment_df$ID != "PG", ]
  enrichment_df <- enrichment_df[enrichment_df$ID != "Other", ]
  fdr_significant <- enrichment_df %>%
    filter(FDR < 0.05 & Count>5) %>%
    mutate(log_q = -log10(FDR))
  
  pdf(file_name, width = 8, height = 5)
  
  print(
    ggplot(fdr_significant, aes(x = Direction, y = reorder(Description, Count))) +
      geom_point(aes(size = Count, color = log_q)) +
      scale_color_gradient(low = "blue", high = "red", name = "-log10(pvalue)") +
      #geom_point(data = filter(fdr_significant, FDR < 0.05), 
                 #aes(x = Direction, y = reorder(Description, Count)), 
                 #size = 5.4, shape = 21, fill = NA, color = "black", stroke = 1)+
      theme_minimal() +
      labs(title = plot_title, x = "", y = "Function")
  )
  
  dev.off()
}

run_and_label_general <- function(sub_level, assoc_table, pheno, enrichment_fun, positive_labels = c("Group1", "Group2")) {
  group1 <- enrichment_fun(sub_level, assoc_table, pheno, Positive = TRUE) %>%
    mutate(Source = sub_level, Direction = positive_labels[1])
  
  group2 <- enrichment_fun(sub_level, assoc_table, pheno, Positive = FALSE) %>%
    mutate(Source = sub_level, Direction = positive_labels[2])
  
  bind_rows(group1, group2)
}

all_timepoint_enrichment <- bind_rows(
  run_and_label_general('ORIGIN', asso_early_late, 'Timepoint_early_late', do_enrichment_timepoint, c("Late", "Early")),
  run_and_label_general('FUNCTION_AT_DESTINATION_1', asso_early_late, 'Timepoint_early_late', do_enrichment_timepoint, c("Late", "Early")),
  run_and_label_general('FUNCTION_AT_DESTINATION_2', asso_early_late, 'Timepoint_early_late', do_enrichment_timepoint, c("Late", "Early")),
  run_and_label_general('FUNCTION_AT_DESTINATION_3', asso_early_late, 'Timepoint_early_late', do_enrichment_timepoint, c("Late", "Early")),
  run_and_label_general('Glycan_annotation', asso_early_late, 'Timepoint_early_late', do_enrichment_timepoint, c("Late", "Early"))
)
all_timepoint_enrichment$FDR<-p.adjust(all_timepoint_enrichment$pvalue, method = "fdr")

plot_enrichment(all_timepoint_enrichment, "CAZyme enrichment by timepoint", "cazyme_substrate_infant_timepoints_early_late.pdf")

write.table(all_timepoint_enrichment, "all_timepoint_enrichment_caz_substrates.txt", sep = "\t", row.names = F)




all_delivery_enrichment <- bind_rows(
  run_and_label_general('ORIGIN', asso_del_mode, 'birth_deliverybirthcard_mode_binary', do_enrichment_phenotypes, c("Vaginal","C-section" )),
  run_and_label_general('FUNCTION_AT_DESTINATION_1', asso_del_mode, 'birth_deliverybirthcard_mode_binary', do_enrichment_phenotypes, c("Vaginal","C-section")),
  run_and_label_general('FUNCTION_AT_DESTINATION_2', asso_del_mode, 'birth_deliverybirthcard_mode_binary', do_enrichment_phenotypes, c("Vaginal","C-section")),
  run_and_label_general('FUNCTION_AT_DESTINATION_3', asso_del_mode, 'birth_deliverybirthcard_mode_binary', do_enrichment_phenotypes, c("Vaginal","C-section")),
  run_and_label_general('Glycan_annotation', asso_del_mode, 'birth_deliverybirthcard_mode_binary', do_enrichment_phenotypes, c("Vaginal","C-section"))
)
all_delivery_enrichment$FDR<-p.adjust(all_delivery_enrichment$pvalue, method = "fdr")

plot_enrichment(all_delivery_enrichment, "CAZyme enrichment by delivery mode", "cazyme_delivery_mode.pdf")
write.table(all_delivery_enrichment, "all_delivery_enrichment_caz_substrates.txt", sep = "\t", row.names = F)


all_feed_enrichment <- bind_rows(
  run_and_label_general('ORIGIN', asso_feed_mode, 'infant_ffq_feeding_mode_simple', do_enrichment_phenotypes, c("FF", "BF" )),
  run_and_label_general('FUNCTION_AT_DESTINATION_1', asso_feed_mode, 'infant_ffq_feeding_mode_simple', do_enrichment_phenotypes, c("FF", "BF")),
  run_and_label_general('FUNCTION_AT_DESTINATION_2', asso_feed_mode, 'infant_ffq_feeding_mode_simple', do_enrichment_phenotypes, c("FF", "BF")),
  run_and_label_general('FUNCTION_AT_DESTINATION_3', asso_feed_mode, 'infant_ffq_feeding_mode_simple', do_enrichment_phenotypes, c("FF", "BF")),
  run_and_label_general('Glycan_annotation', asso_feed_mode, 'infant_ffq_feeding_mode_simple', do_enrichment_phenotypes, c("FF", "BF"))
)
all_feed_enrichment$FDR<-p.adjust(all_feed_enrichment$pvalue, method = "fdr")

plot_enrichment(all_feed_enrichment, "CAZyme enrichment by feeding mode", "cazyme_feeding_mode.pdf")
write.table(all_feed_enrichment, "all_feed_enrichment_caz_substrates.txt", sep = "\t", row.names = F)

