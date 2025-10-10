# MERGING AND CLEANING INFANT PERMANOVA AND TCAM RESULTS ###################
# Author: T.Sinha 

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint/")

# Results permanova 
timepoint_base<-read.delim("result_adonis_per_timepoint_infants_BY_MARGIN_10000_perm_cor_base_28_09_2025.txt")
timepoint_cor_del<-read.delim("result_adonis_per_timepoint_infants_BY_MARGIN_10000_perm_cor_delivery_28_09_2025.txt")
timepoint_cor_del_feed<-read.delim("result_adonis_per_timepoint_infants_BY_MARGIN_10000_perm_cor_delivery_feeding_28_09_2025.txt")

rename_except_cross <- function(df, suffix) {
  cols_excluded <- c("Phenotype", "Timepoint")
  colnames(df) <- ifelse(colnames(df) %in% cols_excluded, colnames(df), paste0(colnames(df), suffix))
  return(df)
}
timepoint_base<- rename_except_cross (timepoint_base, "_cor_base")
timepoint_cor_del  <- rename_except_cross (timepoint_cor_del, "_cor_delivery")
timepoint_cor_del_feed <- rename_except_cross (timepoint_cor_del_feed , "_cor_delivery_breastfed")

combined_results_timepoint_adonis <- reduce(
  list(timepoint_base, timepoint_cor_del, timepoint_cor_del_feed),
  full_join,
  by = c("Phenotype", "Timepoint")
)

combined_results_timepoint_adonis<-na.omit(combined_results_timepoint_adonis)
combined_results_timepoint_adonis<-combined_results_timepoint_adonis[combined_results_timepoint_adonis$Df_cor_base<3,]

# Supplementary Table S25
write.table(combined_results_timepoint_adonis, "combined_results_timepoint_adonis.txt", sep = "\t", row.names = F)

# Results permanova TCAM
TCAM_base<-read.delim("result_adonis_TCAM_INFANTS_BY_MARGIN_10000_perm_no_cor_28_09_2025.txt")
TCAM_cor_del<-read.delim("result_adonis_TCAM_INFANTS_BY_MARGIN_10000_perm_cor_delivery_28_09_2025.txt")
TCAM_cor_del_feed<-read.delim("result_adonis_TCAM_INFANTS_BY_MARGIN_10000_perm_cor_delivery_feeding_28_09_2025.txt")


TCAM_base<- rename_except_cross (TCAM_base, "_cor_base")
TCAM_cor_del  <- rename_except_cross (TCAM_cor_del , "_cor_delivery")
TCAM_cor_del_feed <- rename_except_cross (TCAM_cor_del_feed , "_cor_delivery_breastfed")

combined_results_TCAM_adonis <- reduce(
  list(TCAM_base, TCAM_cor_del, TCAM_cor_del_feed),
  full_join,
  by = c("Phenotype")
)


combined_results_TCAM_adonis<-combined_results_TCAM_adonis[combined_results_TCAM_adonis$Df_cor_base<3,]
# Supplementary Table S26
write.table(combined_results_TCAM_adonis, "combined_results_TCAM_adonis.txt", sep = "\t", row.names = F)

