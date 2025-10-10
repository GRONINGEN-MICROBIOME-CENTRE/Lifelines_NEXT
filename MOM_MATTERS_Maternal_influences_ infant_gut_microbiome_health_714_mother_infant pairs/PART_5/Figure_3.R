############################### Figure 3 #########################
# Author: Trishla Sinha
# Date created:5th August, 2024 
# Last updated: 4th October, 2025

#load packages 
library(tidyverse)
library(mmrm)
library(wesanderson)
library(reshape2)
library(pheatmap)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(readxl)
library(reshape2)


source("Functions_associations_phenotypes_LLNEXT.R")

# Figure 1 A 

# Difference in Bacteroides, Parabacteroides, and Phocaeicola in VG versus CS 

metadata<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
metadata<- metadata %>%
  filter(!(Type == "mother" & Timepoint_categorical %in% c("M1", "M2")))
metadata$Timepoint_categorical=factor(metadata$Timepoint_categorical, levels = c("P12","P28","B", "W2", "M1", "M2", "M3", "M6", 'M9', "M12"))
metadata$Type=factor(metadata$Type, levels = c("mother", "infant"))
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"
metadata_infants$Timepoint_categorical=factor(metadata_infants$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))


# Loading phenotypes 
cross_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_2023_11_15.txt")
cross_selection <-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_selection_AFTER_correlation_&_correction_10_07_2024.txt")
infant_cross_selection<-cross_selection[cross_selection$infant_microbiome_selection==1,]
column_names <- infant_cross_selection[[1]]
infant_cross_phenotypes <- cross_phenotypes %>% select(all_of(column_names))
infant_metadata_cross_phenotypes<-left_join(metadata_infants, infant_cross_phenotypes)
row.names(infant_metadata_cross_phenotypes)<-infant_metadata_cross_phenotypes$NG_ID
infant_metadata_cross_phenotypes_1 <- infant_metadata_cross_phenotypes %>%
  distinct(SAMPLE_ID, .keep_all = TRUE) 
infant_metadata_cross_phenotypes_2 <-drop_phenotypes(infant_metadata_cross_phenotypes_1, 200, 5)
names (infant_metadata_cross_phenotypes_2)[2]<-"NEXT_ID"
infant_metadata_cross_phenotypes_2$NEXT_ID=as.factor(infant_metadata_cross_phenotypes_2$NEXT_ID)
# Load microbiome data all levels 
metaphlan <-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_all_levels_infants_20_07_2023.txt")

filtered_metaphlan <- metaphlan[, grepl("g__", colnames(metaphlan)) & 
                                  !grepl("s__", colnames(metaphlan)) & 
                                  !grepl("t__", colnames(metaphlan))]
filtered_metaphlan_subset <- filtered_metaphlan[, grepl("g__Bacteroides|g__Parabacteroides|g__Phocaeicola", colnames(filtered_metaphlan))]
filtered_subset_bifido <- filtered_metaphlan[, grepl("g__Bifidobacterium|g__GGB6606", colnames(filtered_metaphlan))]

filtered_metaphlan_subset $all_bacteroides <- rowSums(filtered_metaphlan_subset )


all<-merge(infant_metadata_cross_phenotypes_2, filtered_metaphlan_subset, by="row.names")
row.names(all)<-all$Row.names
all$Row.names=NULL


ggplot(all[!is.na(all$birth_deliverybirthcard_mode_binary),], aes(x = Timepoint_categorical, y = all_bacteroides, fill = birth_deliverybirthcard_mode_binary, color=birth_deliverybirthcard_mode_binary)) +
  geom_boxplot(aes(color = birth_deliverybirthcard_mode_binary), alpha = 0.6, outlier.colour = NA, width = 0.5, position = position_dodge(width = 0.9)) +
  geom_point(alpha=0.3,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0), size = 1) +
  scale_fill_manual(values = c("#ac2120", "#65a03f")) +
  scale_color_manual(values = c("#ac2120", "#65a03f")) +
  ggtitle("") +
  theme_bw() +
  labs(x = "", y = "Infant Relative Abundance Bacteroides", fill = "Mode of Delivery", color = "Mode of Delivery") +
  theme(
    plot.title = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 8),
    axis.text.x = element_text(size = 16, angle = 60, hjust = 1),
    legend.position = "top" # Move legend to the top
  )
ggsave(filename = "~/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_3_new/Figure_3_A_Bacteroides_mode_delivery.pdf", 
       width = 4.7, height = 4)



# Figure 1B 
# Bacteroides with vaginal delivery 
setwd("~/Desktop/LLNEXT/Analysis/results/taxa/")
vaginal_only_results<-read.delim("SGB_associations_VAGINAL_ONLY_27_09_2025.txt")
vaginal_only_results$FDR <-p.adjust(vaginal_only_results$p, method = "fdr")
#write.table(vaginal_only_results, "vaginal_deliveries_birth_factors_all_01_10_2025.txt", sep = "\t", row.names = F)
str (vaginal_only_results)
bacteroides<- vaginal_only_results %>%
  filter(grepl("Bacteroides|Parabacteroides|Phocaeicola", outcome))
bacteroides <- bacteroides %>%
  filter(!(trait %in% c("birth_birthcard_place_delivery_hospital", "birth_birthcard_place_delivery_complex", "mother_birthcard_del_rupture_grade", "mother_deliverybirthcard_del_induced_IV_drip")))
#bacteroides$FDR <-p.adjust(bacteroides$p, method = "fdr")
#write.table(bacteroides, "bacteroides_vaginal_deliveries_birth_factors.txt", sep = "\t", row.names = F)
bacteroides$phenotype<-paste0(bacteroides$trait,bacteroides$levels)
setDT(bacteroides)
wide_data <- dcast(bacteroides, outcome ~ phenotype, value.var = "t.value")
wide_data<-as.data.frame.array(wide_data)
row.names(wide_data) <- wide_data$outcome
wide_data$outcome <- NULL
#wide_data[is.na(wide_data)] <- 0

# Create a matrix of significance markers
significance_matrix <- dcast(bacteroides, outcome ~ phenotype, value.var = "p")
significance_matrix<-as.data.frame(significance_matrix)
row.names(significance_matrix) <- significance_matrix$outcome
significance_matrix$outcome <- NULL


# Convert p-values to significance markers
significance_matrix <- ifelse(significance_matrix < 0.05, "*", "")

# Convert the data frame to a matrix
heatmap_data <- as.matrix(t(wide_data))
custom_colors <- colorRampPalette(c("blue", "white", "red"))

# Plot clustered heatmap with significance markers
setwd("~/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_3_new/")
pdf("Figure_3_B.pdf", width = 8, height = 6)
bacteroides_plot <- pheatmap(heatmap_data,
                        cluster_rows = T, 
                        cluster_cols = T, 
                        clustering_distance_cols = "euclidean", 
                        clustering_method = "complete",
                        color = custom_colors(20),
                        display_numbers = t(significance_matrix),
                        main = "", 
                        angle_col = 90,
                        fontsize_col = 8,
                        fontsize_row = 8, 
                        fontsize_number = 12)

dev.off()


################ Figure 3C ###################

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/taxa/")
#feeding mode
cs_results<-read.delim("~/Desktop/LLNEXT/Analysis/results/taxa/taxa_SGB_cross_phenotypes_results_23_09.txt")
long_results <- read.delim("~/Desktop/LLNEXT/Analysis/results/taxa/taxa_SGB_dynamic_phenotypes_results_infants_24_09_2025.txt")

feeding_mode_cs_results <- cs_results %>%
  filter(trait %in% c("infant_ffq_ever_never_breastfed", "infant_birthcard_feeding_mode_after_birth"))

feeding_mode_cs_results$phenotype<-paste0(feeding_mode_cs_results$trait,feeding_mode_cs_results$levels)
feeding_mode_cs_results_plot<-feeding_mode_cs_results[feeding_mode_cs_results$p_cor_delivery<0.05,]

feeding_mode_long_results <- long_results %>%
  filter(trait %in% "infant_ffq_feeding_mode_complex")

feeding_mode_long_results$phenotype<-paste0(feeding_mode_long_results$trait,feeding_mode_long_results$effect.level)
feeding_mode_long_results_plot<-feeding_mode_long_results[feeding_mode_long_results$p_cor_delivery<0.05,]


feeding_mode_cs_results_plot <- feeding_mode_cs_results_plot %>%
  select(c(outcome,
           trait,
           t.value_cor_delivery,
           p_cor_delivery,
           phenotype))

#calculate t.value for longitudinal results
feeding_mode_long_results_plot <- feeding_mode_long_results_plot %>%
  mutate(t.value_cor_delivery=Estimate_cor_delivery/SE_cor_delivery)

feeding_mode_long_results_plot <- feeding_mode_long_results_plot %>%
  select(c(outcome,
           trait,
           t.value_cor_delivery,
           p_cor_delivery,
           phenotype))

feeding_mode_results <- bind_rows(feeding_mode_cs_results_plot, feeding_mode_long_results_plot)
setDT(feeding_mode_results) 
wide_data <- dcast(feeding_mode_results, outcome ~ phenotype, value.var = 't.value_cor_delivery')
wide_data<-as.data.frame(wide_data)
row.names(wide_data) <- wide_data$outcome
wide_data$outcome <- NULL
wide_data[is.na(wide_data)] <- 0

feeding_mode_results$FDR<-p.adjust(feeding_mode_results$p_cor_delivery, method = "fdr")
significance_matrix <- dcast(feeding_mode_results, outcome ~ phenotype, value.var = 'FDR')
row.names(significance_matrix) <- significance_matrix$outcome
significance_matrix$outcome <- NULL
significance_matrix <- ifelse(!is.na(significance_matrix) & significance_matrix < 0.05, '*', '')


# Convert the data frame to a matrix
heatmap_data <- as.matrix(wide_data)
custom_colors <- colorRampPalette(c("blue", "white", "red"))
# Plot clustered heatmap with significance markers


pdf('Figure3_C_feeding_mode_heatmap.pdf', width = 20, height = 6)

feeding_mode_results_plot <- pheatmap(t(heatmap_data),
                                      cluster_rows = TRUE,
                                      clustering_distance_rows = 'euclidean',
                                      clustering_method = 'complete',
                                      color = custom_colors(20),
                                      display_numbers = t(significance_matrix),  
                                      main = '',
                                      angle_col = 90,
                                      fontsize_col = 8,
                                      fontsize_row = 7,
                                      fontsize_number = 10)

feeding_mode_results_plot
dev.off()

#################### Figure 3D #################################

source ("Functions_associations_phenotypes_LLNEXT.R") 

# Load metadata and dynamic phenotypes
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/metadata_basic_phenotypes.txt")
dynamic_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_longitudinal_2023_09_29.txt")
dynamic_selection<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_longitudinal_selection_AFTER_correlation_&_correction_10_07_2024.txt")

# Merging files and selecting infant relevant data
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"
metadata_infants$BATCH_NUMBER<-as.factor(metadata_infants$BATCH_NUMBER)
metadata_infants$Timepoint_categorical=factor(metadata_infants$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))


# Selecting data relevant for infant 
infant_dynamic_selection<-dynamic_selection[dynamic_selection$infant_microbiome_selection==1,]
column_names <- infant_dynamic_selection[[1]]
infant_dynamic_phenotypes <- dynamic_phenotypes %>% select(all_of(column_names))
infant_metadata_dynamic_phenotypes<-left_join(metadata_infants, infant_dynamic_phenotypes)
row.names(infant_metadata_dynamic_phenotypes)<-infant_metadata_dynamic_phenotypes$NG_ID

# For this specific data remove duplicates for categorical timepoints 
infant_metadata_dynamic_phenotypes_1 <- infant_metadata_dynamic_phenotypes %>%
  distinct(SAMPLE_ID, .keep_all = TRUE) 

# Remove phenotypes with too many NA's 
infant_metadata_dynamic_phenotypes_2 <-drop_phenotypes(infant_metadata_dynamic_phenotypes_1, 200, 5)
names (infant_metadata_dynamic_phenotypes_2)[2]<-"NEXT_ID"
infant_metadata_dynamic_phenotypes_2$NEXT_ID=as.factor(infant_metadata_dynamic_phenotypes_2$NEXT_ID)

GBM<-read.delim("~/Desktop/LLNEXT/Analysis/GBM/NEXT_GBM_merged_clean.tsv")
row.names(GBM)<-GBM$ID
GBM$ID=NULL
GBM<- na.omit(GBM)

# Setting a prevalence cut-off 
infant_NEXT_ID<-infant_metadata_cross_phenotypes_2 %>%
  select(NEXT_ID)
infant_GBM_all<-merge(infant_NEXT_ID,GBM, by="row.names" )
row.names(infant_GBM_all)<-infant_GBM_all$Row.names
infant_GBM_all$Row.names=NULL
unique_counts <- sapply(infant_GBM_all, function(x) length(unique(infant_GBM_all$NEXT_ID[x >0.001]))) 
infant_GBM_all_filt <- infant_GBM_all[, unique_counts >= 0.3*length(unique(infant_GBM_all$NEXT_ID)) ] # Setting a 30% cut-off on prevalence 
infant_GBM_all_filt$NEXT_ID=NULL


## Loading microbiome data 
metaphlan = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
metaphlan_clr = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
metaphlan_infants = metaphlan[match(as.character(metadata_infants$NG_ID),rownames(metaphlan)),grep("t__",colnames(metaphlan))]
colnames(metaphlan_infants) = sub(".*s__","",colnames(metaphlan_infants))
metaphlan_infants = metaphlan_infants[,colnames(metaphlan_clr)] # Filtering to only the columns in the filtered file 


metaphlan_infants<- metaphlan_infants[rownames(metaphlan_infants) %in% rownames(infant_GBM_all_filt), ]
infant_GBM_all_filt<- infant_GBM_all_filt[rownames(infant_GBM_all_filt) %in% rownames(metaphlan_infants), ]
infant_metadata_cross_phenotypes_2<- infant_metadata_cross_phenotypes_2[rownames(infant_metadata_cross_phenotypes_2) %in% rownames(infant_GBM_all_filt), ]

# Data transformations 
GBM_alr <-do_clr_externalWeighting(infant_GBM_all_filt,metaphlan_infants)
GBM_alr<-as.data.frame(GBM_alr)
GBM_alr_null <-nullify_zeros(GBM_alr,infant_GBM_all_filt)

all<-merge(infant_metadata_dynamic_phenotypes_2,GBM_alr_null ,by="row.names" )
all$Row.names=NULL
row.names(all)<-all$NG_ID



clean_data <- all %>%
  filter(!is.na(infant_ffq_feeding_mode_simple)) %>%
  mutate(infant_ffq_feeding_mode_simple = factor(
    infant_ffq_feeding_mode_simple,
    levels = c("excl_BF", "excl_FF")  
  ))




acetate_synthesis<-ggplot(clean_data, aes(
  x = infant_ffq_feeding_mode_simple,
  y = MGB043,
  fill = infant_ffq_feeding_mode_simple
)) +
  
  geom_jitter(aes(color = infant_ffq_feeding_mode_simple),
              width = 0.2, alpha = 0.6, size = 2) +
  labs(
    title = "",
    x = "Feeding Mode",
    y = "Acetate synthesis I"
  ) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  facet_wrap(~ Timepoint_categorical, nrow = 1)
acetate_synthesis


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_3_new")
ggsave("acetate_synthesis_by_feeding_mode.pdf", plot = acetate_synthesis,
       width = 4, height = 4, units = "in")






################ Figure 3E ###################
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes")
infant_taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
results_adonis_infants <- read.delim("adonis_cazyme_taxa_by_timepoint_results.txt")
significant <- results_adonis_infants[results_adonis_infants$p_value < .0001, ]

# Clean taxon names
significant$Taxon <- as.character(significant$Taxon)
significant$Taxon <- gsub("\\.t__", " ", significant$Taxon)
significant$Taxon <- gsub("_", " ", significant$Taxon)

# Clean taxon order names
taxon_order <- names(infant_taxa)
clean_taxon_order <- gsub("\\.t__", " ", taxon_order)
clean_taxon_order <- gsub("_", " ", clean_taxon_order)

# Set factor levels
significant$Taxon <- factor(significant$Taxon, levels = clean_taxon_order, ordered = TRUE)

# Timepoint and R2 filtering
significant$Timepoint <- factor(significant$Timepoint,
                                levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12"),
                                ordered = TRUE)

significant <- significant %>%
  mutate(R2 = round(R2, 3)) %>%
  filter(R2 > 0.06)

# Segment data for lines
segments <- significant %>%
  mutate(Timepoint_num = as.numeric(Timepoint)) %>%
  group_by(Taxon) %>%
  summarise(Timepoint2 = min(Timepoint_num),
            xend = max(Timepoint_num)) %>%
  mutate(Timepoint = levels(significant$Timepoint)[Timepoint2],
         xend = levels(significant$Timepoint)[xend])
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_3_new/")
pdf("pcos_bacteria_cazyme_profile_infants.pdf", width = 8, height = 6)
ggplot() +
  geom_segment(data = segments, aes(x = Timepoint, y = Taxon,
                                    xend = xend, yend = Taxon),
               color = "black", linewidth = 0.2) +
  geom_point(data = significant,
             aes(x = Timepoint, y = Taxon,
                 size = R2, fill = Timepoint),
             shape = 21, color = "black") +
  xlab("Timepoint") +
  ylab("") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_size_continuous(range = c(1, 4)) +
  scale_fill_manual(values = c("#B5DD88", "#41C0B4", "#4397BB",
                               "#ECA4C9", "#cb4563", "#A42097", "#390962"))

dev.off()


################# Figure 3 F ############
asso_early_late<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes/association_early_late_timepoint_infants.txt")
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_3_new/")

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
      scale_color_gradient(low = "blue", high = "red", name = "-log10(FDR)") +
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


save.image(file = "Figure_3_for_Asier.RData")
