### Supplementary figure 1 
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/taxa/dynamics_pregnancy")
time_on_microbiome<-read.delim("mothers_dynamic_time_pregnancy.txt")
time_on_microbiome_sig<- time_on_microbiome %>% filter(FDR < 0.05)
wide_data <- dcast(time_on_microbiome_sig, Bug ~ Feature, value.var = "t.value")
colnames(wide_data) <- gsub("Timepoint_categorical", "", colnames(wide_data))

filtered_all_wide <- wide_data  %>%
  select(Bug, P28, B, M3)

row.names(filtered_all_wide)<-filtered_all_wide$Bug
filtered_all_wide$Bug=NULL

# Convert the data frame to a matrix

heatmap_data <- as.matrix(t(filtered_all_wide))

# Plot clustered heatmap
time_mother_SGB<-pheatmap(heatmap_data,
         cluster_rows = F, 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = F,
         main = "Heatmap of t values of significant associations of time with maternal microbiome", 
         angle_col = 90,
         fontsize_col = 8)
time_mother_SGB
ggsave("~/Desktop/LLNEXT/Analysis/results/figures/supplementary/supplementary_figure_1.pdf",time_mother_SGB,
      width = 12, height = 5 )





