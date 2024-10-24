#############################################################
# DYNAMICS OF MICROBIOME : Timepoint Microbiome Associations
#############################################################

# Author: Trishla Sinha 
# Date last update: 1st of August, 2024

###############################
#Functions
##############################
mixed_models_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list) {
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
##############################
# Loading libraries
##############################
library(gplots)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(MetBrewer)
library(RLRsim)
library(lmerTest)
library(ggforce)
library(pheatmap)
##############################
# MOTHERS 
##############################

metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
row.names(metadata)<-metadata$NG_ID
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_mothers<-metadata[metadata$Type=="mother", ]
metadata_mothers <- subset(metadata_mothers, Timepoint_categorical != "M1" & Timepoint_categorical != "M2")
metadata_mothers$Timepoint_categorical<-factor(metadata_mothers$Timepoint_categorical, levels = c("P12", "P28", "B", "M3"))

taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/NEXT_metaphlan_4_CLR_transformed_fil_30_percent_SGB_mothers_03_08_2023.txt")


time_on_microbiome <- mixed_models_taxa(metadata_mothers, 
                                                 "NG_ID", 
                                                taxa, 
                                                 c("Timepoint_categorical"))

write.table(time_on_microbiome, "timepoint_mother_SGB's.txt", sep = "\t", row.names = F)
time_on_microbiome_sig<- time_on_microbiome %>% filter(FDR < 0.05)

wide_data <- dcast(time_on_microbiome_sig, Bug ~ Feature, value.var = "t value")
colnames(wide_data) <- gsub("Timepoint_categorical", "", colnames(wide_data))

filtered_all_wide <- wide_data  %>%
  select(Bug, P28, B, M3)

row.names(filtered_all_wide)<-filtered_all_wide$Bug
filtered_all_wide$Bug=NULL

# Convert the data frame to a matrix

heatmap_data <- as.matrix(t(filtered_all_wide))

# Plot clustered heatmap
mother<-pheatmap(heatmap_data,
         cluster_rows = F, 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = F,
         main = "Heatmap of t values of significant associations of time with maternal microbiome", 
         angle_col = 90,
         fontsize_col = 8)
mother 
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/supplementary")
pdf("supplementary_figure_1_A.pdf", 
       width = 12, height = 5)
dev.off()

#################################
# INFANTS
####################################
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
row.names(metadata)<-metadata$NG_ID
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
metadata_infants$Timepoint_categorical=factor(metadata_infants$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))

taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")


time_on_microbiome <- mixed_models_taxa(metadata_infants, 
                                        "NG_ID", 
                                        taxa, 
                                        c("Timepoint_categorical"))

write.table(time_on_microbiome, "timepoint_infant_SGB's.txt", sep = "\t", row.names = F)

time_on_microbiome_sig<- time_on_microbiome %>% filter(FDR < 0.05)

wide_data <- dcast(time_on_microbiome_sig, Bug ~ Feature, value.var = "t value")
colnames(wide_data) <- gsub("Timepoint_categorical", "", colnames(wide_data))

filtered_all_wide <- wide_data  %>%
  select(Bug, M1, M2, M3, M6, M9, M12)

row.names(filtered_all_wide)<-filtered_all_wide$Bug
filtered_all_wide$Bug=NULL

# Convert the data frame to a matrix

heatmap_data <- as.matrix(t(filtered_all_wide))

# Plot clustered heatmap

pdf ("supplementary_figure_1_B.pdf", 
       width = 12, height = 5)
infant<-pheatmap(heatmap_data,
         cluster_rows = F, 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = F,
         main = "Heatmap of t values of significant associations of time with infant microbiome", 
         angle_col = 90,
         fontsize_col = 8)
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/supplementary")
infant
dev.off()

# Week 2 taxa 
metadata_infants_w2<-metadata_infants[metadata_infants$Timepoint_categorical=="W2",]
taxa_m2<- taxa[rownames(taxa) %in% rownames(metadata_infants_w2), ]
column_means <- colMeans(taxa_m2, na.rm = TRUE)

# Convert the means to a dataframe
means_df <- data.frame(Column = names(column_means), Mean = column_means)

# Arrange the dataframe in ascending order of the means
means_df <- means_df[order(means_df$Mean), ]

