# Making plot heatmap cazymes mode of delivery 
### AUTHOR:  TRISHLA SINHA
### ORIGINAL SCRIPT: 9th June, 2025
### LAST UPDATE: 12th September, 2025 

# Supplementary figure S6

library(pheatmap)
library(RColorBrewer)
library(tidyverse)
source("/Users/trishlasinha/Desktop/LLNEXT/Analysis/scripts/new_scripts_2025/pheatmap_empty.r")

# Function to prepare cazyme data 
Prepare_cazymes = function(File = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes/NEXT_cayman_infant.txt", Prevalence_min = 0.4, Transform = "log" ){
  caz = read_tsv(File)
  caz %>% filter(! feature %in% c("total_reads", "filtered_reads")) -> caz
  #Get unique counts + ambigous counts normalized
  caz %>% select(ID, feature, combined_rpkm ) %>% spread(key = feature, value = combined_rpkm) -> caz
  #Change ID
  caz$ID = caz$ID 
  #Get summary
  summary_as_df <- function(x) {
    result <- summary(x)
    data.frame(stat = names(result), value = as.vector(result))
  }
  summary_df <- caz %>%
    select(-ID) %>%
    purrr::map(summary_as_df) %>%
    bind_rows(.id = "Cazyme") %>% spread(stat, value) %>% as_tibble()
  summary_df = summary_df %>% mutate( `NA's` = ifelse( is.na(`NA's`), 0,  `NA's` ), N = dim(caz)[1], Present= N - `NA's`)
  summary_df = summary_df %>% mutate(Present_perc = Present/N ) %>% mutate(Keep = ifelse(Present_perc>=Prevalence_min, T, F ) )
  #Make NA into 0
  caz %>% mutate_all(~ replace(., is.na(.), 0)) -> caz
  #Remove Cazymes with high NA prop
  caz %>% select(c("ID", summary_df$Cazyme[summary_df$Keep==T] )) -> caz_anal
  if (Transform == "log"){
    Log_transform(caz_anal)  ->caz_anal_tranf
  } else if (Transform == "clr"){
    CLR_transformation(caz_anal)  -> caz_anal_tranf
  } else { caz_anal_tranf = NULL}
  
  return(list("caz" = caz, "caz_filtered" = caz_anal, "caz_tranf"=caz_anal_tranf  ,"Summary_caz" = summary_df ) )
}


Log_transform = function(df){
  colnames(df)[1] = "ID"
  PS = Pseudocount(df)
  df %>% select(-ID) -> df2
  log10(df2 + PS) %>% as_tibble() %>% mutate(ID = df$ID, .before=1) %>% return()
}
Pseudocount = function(df_tara){
  df_tara %>% select(-ID) %>% as.matrix() -> PS
  min(PS[PS!=0])/2 %>% return()
}
CLR_transformation = function( df_tara ){
  df_2 = select(df_tara, -ID)
  df_2 = df_2 + Pseudocount(df_tara)
  matrix_data = as.matrix(df_2)
  
  clr_matrix <- apply(matrix_data, 1, function(row){
    exp(mean(log(row)) ) -> geometric_mean
    log(row / geometric_mean)
  } ) %>% t()
  clr_matrix %>% as_tibble() %>% mutate(ID =df_tara$ID  , .before=1) %>% return()
  
}
Filter_prevalence = function(df_tara, Prevalence = 0.4 ){
  df_tara %>% select(-ID) %>% apply(2, function(x){ mean(x!=0) }  ) -> Prevalences
  Prevalences[Prevalences>=Prevalence] %>% names() -> KEEP
  df_tara %>% select(c("ID",KEEP) ) %>% return()
}


# Load data 
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes")
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

# Loading cazymes 
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes")
caz_tables = Prepare_cazymes(File = "NEXT_cayman_infant.txt", Prevalence_min = 0.3, Transform = "NULL" )
caz_filtered = as.data.frame (caz_tables[["caz_filtered"]])
row.names(caz_filtered)<-caz_filtered$ID
caz_filtered$ID<-NULL
metaphlan = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
metaphlan_clr = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
metaphlan_infants = metaphlan[match(as.character(metadata_infants$NG_ID),rownames(metaphlan)),grep("t__",colnames(metaphlan))]
colnames(metaphlan_infants) = sub(".*s__","",colnames(metaphlan_infants))
metaphlan_infants = metaphlan_infants[,colnames(metaphlan_clr)] # Filtering to only the columns in the filtered file 
metaphlan_infants <- metaphlan_infants[rownames(metaphlan_infants) %in% rownames(caz_filtered), ]
cazyme_alr = do_clr_externalWeighting(caz_filtered,metaphlan_infants)
cazyme_alr<-as.data.frame(cazyme_alr)
cazyme_alr_null <-nullify_zeros(cazyme_alr,caz_filtered)



# Load microbiome data all levels 
metaphlan <-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_all_levels_infants_20_07_2023.txt")

filtered_metaphlan <- metaphlan[, grepl("g__", colnames(metaphlan)) & 
                                !grepl("s__", colnames(metaphlan)) & 
                                !grepl("t__", colnames(metaphlan))]
filtered_metaphlan_subset <- filtered_metaphlan[, grepl("g__Bacteroides|g__Parabacteroides|g__Phocaeicola", colnames(filtered_metaphlan))]
filtered_subset_bifido <- filtered_metaphlan[, grepl("g__Bifidobacterium|g__GGB6606", colnames(filtered_metaphlan))]

filtered_metaphlan_subset $all_bacteroides <- rowSums(filtered_metaphlan_subset )

filtered_metaphlan_subset$all_bacteroides_binary <- ifelse(filtered_metaphlan_subset$all_bacteroides > 0.01, "Yes", "No")

infant_metadata_cross_phenotypes_M3<-infant_metadata_cross_phenotypes_2[infant_metadata_cross_phenotypes_2$Timepoint_categorical=="M3",]

MOD<-infant_metadata_cross_phenotypes_M3[, c("infant_ffq_ever_never_breastfed", "birth_deliverybirthcard_mode_binary", "Timepoint_categorical")]

all<-merge(cazyme_alr_null, MOD, by="row.names")
row.names(all)<-all$Row.names
all$Row.names=NULL
all_new<-merge(all, filtered_metaphlan_subset, by="row.names")
row.names(all_new)<-all_new$Row.names
all_new$Row.names=NULL
names (all_new)

cazyme_columns <- names(cazyme_alr_null)
cazyme_matrix <- as.matrix(all_new[, cazyme_columns])


annotation_row <- data.frame(
  birth_mode = all_new$birth_deliverybirthcard_mode_binary,
  Bacteroides = all_new$all_bacteroides_binary,
  row.names = rownames(all_new)
)

# Supplementary Figure S6
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/supp_figures/Parts_of_figures/cazymes")
png("pheatmap_cazyme_profile_M3.png", width = 1000, height = 800, res = 150)  # open PNG device, adjust size and resolution
pdf("pheatmap_cazyme_profile_M3.pdf", width = 10, height = 8)


pheatmap(
  mat = cazyme_matrix,
  annotation_row = annotation_row,
  scale = "row",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize = 2,
  border_color = NA
)
dev.off()


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/supp_figures/Parts_of_figures/cazymes")
pdf("pheatmap_cazyme_profile_M3.pdf", width = 10, height = 8)

pheatmap2(
  mat = cazyme_matrix,
  annotation_row = annotation_row,
  scale = "row",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize = 2,
  border_color = NA
)
dev.off()





kmeans(cazyme_matrix, 2) -> Clusters
tibble( ID = names(Clusters$cluster), Kmean_cluster = Clusters$cluster)
clusters<-tibble( ID = names(Clusters$cluster), Kmean_cluster = Clusters$cluster)
clusters<-as.data.frame(clusters)
row.names(clusters)<-clusters$ID
clusters$ID=NULL
all_new_clusters<-merge(all_new, clusters, by="row.names")

all_new_clusters$Kmean_cluster <- ifelse(all_new_clusters$Kmean_cluster == 1, 0, 1)

prop.table(table(all_new_clusters$Kmean_cluster, all_new_clusters$all_bacteroides_binary), 1)
chisq.test(table(all_new_clusters$Kmean_cluster, all_new_clusters$all_bacteroides_binary))

fisher.test(table(all_new_clusters$Kmean_cluster, all_new_clusters$all_bacteroides_binary))
fisher.test(table(all_new_clusters$Kmean_cluster, all_new_clusters$birth_deliverybirthcard_mode_binary))
fisher.test(table(all_new_clusters$Kmean_cluster, all_new_clusters$infant_ffq_ever_never_breastfed))




plot_data <- all_new[!is.na(all_new$birth_deliverybirthcard_mode_binary), ]

ggplot(plot_data, aes(
  x = Timepoint_categorical,
  y = all_bacteroides,
  fill = birth_deliverybirthcard_mode_binary,
  color = birth_deliverybirthcard_mode_binary
)) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA, width = 0.5, position = position_dodge(0.9)) +
  geom_point(alpha = 0.3, position = position_jitterdodge(jitter.width = 0.3), size = 1) +
  scale_fill_manual(values = c("#ac2120", "#65a03f")) +
  scale_color_manual(values = c("#ac2120", "#65a03f")) +
  labs(x = "", y = "Bacteroides/Parabacteroides/Phocaeicola", fill = "Mode of Delivery", color = "Mode of Delivery") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top"
  )
