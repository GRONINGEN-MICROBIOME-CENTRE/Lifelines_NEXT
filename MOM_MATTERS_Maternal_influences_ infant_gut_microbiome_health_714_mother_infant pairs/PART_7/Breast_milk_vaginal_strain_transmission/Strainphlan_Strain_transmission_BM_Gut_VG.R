############# Strain transmission analysis in BM, VG and gut samples #############
#Author: Trishla Sinha 
# Date updated: 2nd July, 2025

library(ape)
library(phangorn)
library(pbapply)
library(tidyverse)
library(readr)


# Part 1: Making distances from the tress for further strainsharing analysis 

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees")
files <- system("ls RAxML_bestTree.t__*.StrainPhlAn4.tre", intern = TRUE)

trees = pblapply(files, read.tree)
trees <- pblapply(trees, function(x) {
  x$tip.label <- sub("_.*", "", x$tip.label)
  x
})


distmats = pblapply(trees, cophenetic.phylo)

distmats = pblapply(distmats,as.matrix)

outputs = paste0("NEXT_RAxML_distmats/",sub("^NEXT_RAxML_trees/","",sub("/RAxML_bestTree.*","",files)), "_DistMat.txt")


for(i in 1:length(outputs)){
  write.table(distmats[[i]],file = outputs[i],sep="\t",quote = F)
}

# Part 2: Now checking if there are incidences of strain sharing between mother-infant pairs using earlier defined thresholds   

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees/NEXT_RAxML_distmats")

# Load metadata
metadata <- read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Metadata_BM_VG_Gut_20_05_2025.txt")
metadata$Modified_NEXT_ID_without_preg_number <- sub("_.*", "", metadata$NEXT_ID)

# Read thresholds table
thresholds <- read_tsv("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees/Thresholds/Distance_tables/thresholds_Sinha_et_al.tsv")

# List distance matrix files
file_list <- list.files(pattern = "DistMat.txt$")

# Directory to save outputs
output_dir <- "/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees/Processed_RDS"
dir.create(output_dir, showWarnings = FALSE)


# Function to process each species file
process_species <- function(file_name) {
  SGB_name <- sub("_DistMat\\.txt$", "", file_name)
  
  # Read and process matrix
  strain <- read.delim(file_name, row.names = 1)
  strain[upper.tri(strain)] <- NA
  strain$Sample_ID_1 <- rownames(strain)
  
  # Long format
  sp1 <- strain %>%
    pivot_longer(!Sample_ID_1, names_to = "Sample_ID_2", values_to = "dist") %>%
    filter(Sample_ID_1 != Sample_ID_2) %>%
    drop_na()
  
  # Trim IDs
  sp1$Sample_ID_1 <- substr(sp1$Sample_ID_1, 1, 12)
  sp1$Sample_ID_2 <- substr(sp1$Sample_ID_2, 1, 12)
  
  # Merge metadata
  sp3 <- left_join(sp1, metadata %>%
                     select(Sample_ID_1 = NG_ID, NEXT_ID_long_1 = NEXT_ID,
                            NEXT_ID_short_1 = Modified_NEXT_ID_without_preg_number,
                            FAMILY_1 = FAMILY, Timepoint_categorical_1 = Timepoint_categorical, Type_1 = Type, mother_infant_1=mother_infant))
  
  sp3 <- left_join(sp3, metadata %>%
                     select(Sample_ID_2 = NG_ID, NEXT_ID_2 = NEXT_ID,
                            NEXT_ID_short_2 = Modified_NEXT_ID_without_preg_number,
                            FAMILY_2 = FAMILY, Timepoint_categorical_2 = Timepoint_categorical, Type_2 = Type, mother_infant_2=mother_infant))
  
  sp3 <- sp3 %>% filter(Sample_ID_1 != Sample_ID_2)
  
  # Annotate relationships
  sp3 <- sp3 %>%
    mutate(
      same_individual = ifelse(NEXT_ID_short_1 == NEXT_ID_short_2, "same_individual", "different_individual"),
      related = ifelse(FAMILY_1 == FAMILY_2, "related", "unrelated"),
      mother_infant=ifelse (mother_infant_1==mother_infant_2, "no", "yes"),
      nGD = dist / max(dist, na.rm = TRUE)
    )
  
  # Match threshold
  threshold_row <- thresholds %>% filter(file_name == SGB_name)
  if (nrow(threshold_row) == 0) {
    warning(paste("Threshold not found for", SGB_name))
    return(NULL)
  }
  
  threshold_value <- threshold_row$final_threshold
  
  # Strain sharing classification
  sp3 <- sp3 %>%
    mutate(Strain_sharing = ifelse(nGD <= threshold_value, "yes", "no"))
  
  
  
  # Save output
  output_path <- file.path(output_dir, paste0(SGB_name, ".rds"))
  write_rds(sp3, output_path)
  message("Saved: ", output_path)
  
  return(sp3)
}

# Run the function on all files
results <- lapply(file_list, process_species)


# Part 3: Plotting and saving strain sharing results 
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees/Processed_RDS")

# List all .rds files
rds_files <- list.files(pattern = "^t__.*\\.rds$")

# Output directories
plot_dir_related <- "plots_BM_strain_sharing_related"
plot_dir_unrelated <- "plots_BM_strain_sharing_unrelated"
distances<-"distances"

results_list <- list()

for (file in rds_files) {
  
  species_name <- tools::file_path_sans_ext(file)
  check <- read_rds(file)
  
  distance_test_bm <- check %>%
    filter((Type_1 == "Milk" | Type_2 == "Milk") & mother_infant == "yes")
  
  distance_re_unrela_bm <- ggplot(distance_test_bm, aes(x = related, y = nGD, fill = related)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    geom_jitter(aes(color = related), width = 0.2, alpha = 0.7, size = 2) +
    scale_fill_manual(values = c("related" = "skyblue", "unrelated" = "salmon")) +
    scale_color_manual(values = c("related" = "skyblue", "unrelated" = "salmon")) +
    labs(
      title = paste0("Median Distances Between Related and Unrelated Mother-Infant Pairs\n", species_name),
      x = "Related Status",
      y = "Distance"
    ) +
    theme_minimal()
  
  ggsave(
    filename = file.path(distances, paste0(species_name, "_BM_distances_unrelated_related_plot.pdf")),
    plot = distance_re_unrela_bm,
    device = "pdf", width = 7, height = 5
  )
  
  ### Related mother-infant pairs ###
  breast_related <- check %>%
    filter((Type_1 == "Milk" | Type_2 == "Milk") &
             mother_infant == "yes" &
             related == "related")
  
  ### Unrelated mother-infant pairs ###
  breast_unrelated <- check %>%
    filter((Type_1 == "Milk" | Type_2 == "Milk") &
             mother_infant == "yes" &
             related == "unrelated")
  
  if (nrow(breast_related) > 0) {
    # Bar plot: related pairs
    strain_summary_related <- breast_related %>%
      count(Strain_sharing)
    
    p_related <- ggplot(strain_summary_related, aes(x = Strain_sharing, y = n, fill = Strain_sharing)) +
      geom_bar(stat = "identity", color = "black", width = 0.6) +
      labs(
        title = paste0("Strain Sharing in Related Milk Mother-Infant Pairs\n", species_name),
        x = "Strain Sharing", y = "Number of Pairs"
      ) +
      theme_minimal(base_size = 14) +
      scale_fill_manual(values = c("yes" = "blue", "no" = "grey")) +
      theme(legend.position = "none")
    
    ggsave(
      filename = file.path(plot_dir_related, paste0(species_name, "_BM_strain_sharing_related_plot.pdf")),
      plot = p_related,
      device = "pdf", width = 7, height = 5
    )
    
    ## Timepoint-based heatmap for related
    long_df <- breast_related %>%
      mutate(Pair_ID = row_number()) %>%
      select(
        Pair_ID,
        Sample_ID_1, Sample_ID_2,
        Timepoint_1 = Timepoint_categorical_1,
        Timepoint_2 = Timepoint_categorical_2,
        Type_1, Type_2,
        mother_infant_1, mother_infant_2,
        FAMILY = FAMILY_1,
        dist, Strain_sharing
      ) %>%
      pivot_longer(
        cols = starts_with("Sample_ID"),
        names_to = "Side",
        values_to = "Sample_ID"
      ) %>%
      mutate(
        Timepoint = ifelse(Side == "Sample_ID_1", Timepoint_1, Timepoint_2),
        Type = ifelse(Side == "Sample_ID_1", Type_1, Type_2),
        mother_infant = ifelse(Side == "Sample_ID_1", mother_infant_1, mother_infant_2)
      ) %>%
      select(Pair_ID, FAMILY, Sample_ID, Timepoint, Type, mother_infant, dist, Strain_sharing) %>%
      distinct() %>%
      mutate(Timepoint = factor(Timepoint, levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12")))
    
    long_df <- long_df %>%
      filter(!(Timepoint == "M1" & Type == "Milk") & !is.na(Timepoint))
    
    heat_data_binary <- long_df %>%
      group_by(FAMILY, Timepoint) %>%
      summarise(Strain_shared = ifelse(any(Strain_sharing == "yes"), "yes", "no"), .groups = "drop")
    
    p_heat <- ggplot(heat_data_binary, aes(x = Timepoint, y = FAMILY, fill = Strain_shared)) +
      geom_tile(color = "white") +
      scale_fill_manual(values = c("yes" = "blue", "no" = "grey90")) +
      labs(
        title = paste0("Strain Sharing BM (Yes/No) Across Families and Timepoints\n", species_name),
        x = "Timepoint", y = "Family"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )
    
    ggsave(
      filename = file.path(plot_dir_related, paste0(species_name, "_strain_sharing_BM_timepoint_heatmap.pdf")),
      plot = p_heat,
      device = "pdf", width = 7, height = 5
    )
    
  } else {
    message("Skipping ", species_name, " (no related Milk mother-infant pairs)")
  }
  
  if (nrow(breast_unrelated) > 0) {
    strain_summary_unrelated <- breast_unrelated %>%
      count(Strain_sharing)
    
    p_unrelated <- ggplot(strain_summary_unrelated, aes(x = Strain_sharing, y = n, fill = Strain_sharing)) +
      geom_bar(stat = "identity", color = "black", width = 0.6) +
      labs(
        title = paste0("Strain Sharing in Unrelated Milk Mother-Infant Pairs\n", species_name),
        x = "Strain Sharing", y = "Number of Pairs"
      ) +
      theme_minimal(base_size = 14) +
      scale_fill_manual(values = c("yes" = "blue", "no" = "grey")) +
      theme(legend.position = "none")
    
    ggsave(
      filename = file.path(plot_dir_unrelated, paste0(species_name, "_BM_strain_sharing_unrelated_plot.pdf")),
      plot = p_unrelated,
      device = "pdf", width = 7, height = 5
    )
  } else {
    message("Skipping ", species_name, " (no unrelated Milk mother-infant pairs)")
  }
  
  shared_related <- breast_related %>% filter(Strain_sharing == "yes")
  shared_unrelated <- breast_unrelated %>% filter(Strain_sharing == "yes")
  
  timepoints_shared_related <- if (nrow(shared_related) > 0) {
    paste(
      unique(c(shared_related$Timepoint_categorical_1,
               shared_related$Timepoint_categorical_2)) %>% na.omit(),
      collapse = ";"
    )
  } else {
    NA_character_
  }
  
  result_entry <- tibble(
    Species = species_name,
    N_Related = nrow(breast_related),
    N_Unrelated = nrow(breast_unrelated),
    N_Related_Shared = nrow(shared_related),
    N_Unrelated_Shared = nrow(shared_unrelated),
    Timepoints_Shared_Related = timepoints_shared_related
  )
  
  results_list[[length(results_list) + 1]] <- result_entry
  
}  

# After the loop: bind all results into one dataframe
results_all <- bind_rows(results_list)

write.table(results_all, "strain_sharing_BM_mother_infant.txt", row.names = F, sep = "\t")


###### Vaginal strain sharing 

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees/Processed_RDS")

rds_files <- list.files(pattern = "^t__.*\\.rds$")

file="t__SGB9269.rds"

# Changed output directories to VG
plot_dir_related <- "plots_VG_strain_sharing_related"     
plot_dir_unrelated <- "plots_VG_strain_sharing_unrelated" 
distances <- "distances"  

results_list <- list()

for (file in rds_files) {
  
  species_name <- tools::file_path_sans_ext(file)
  check <- read_rds(file)
  
  # Filter for VG instead of Milk
  distance_test_vg <- check %>% 
    filter((Type_1 == "VG" | Type_2 == "VG") & mother_infant == "yes")  
  
  distance_re_unrela_vg <- ggplot(distance_test_vg, aes(x = related, y = nGD, fill = related)) +  
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    geom_jitter(aes(color = related), width = 0.2, alpha = 0.7, size = 2) +
    scale_fill_manual(values = c("related" = "skyblue", "unrelated" = "salmon")) +
    scale_color_manual(values = c("related" = "skyblue", "unrelated" = "salmon")) +
    labs(
      title = paste0("Median Distances Between Related and Unrelated Mother-Infant Pairs\n", species_name),
      x = "Related Status",
      y = "Distance"
    ) +
    theme_minimal()
  
  # Save plot with VG filename
  ggsave(
    filename = file.path(distances, paste0(species_name, "_VG_distances_unrelated_related_plot.pdf")),  
    plot = distance_re_unrela_vg,
    device = "pdf", width = 7, height = 5
  )
  
  ### Related mother-infant pairs ###
  vaginal_related <- check %>%
    filter((Type_1 == "VG" | Type_2 == "VG") &   
             mother_infant == "yes" &
             related == "related")
  
  ### Unrelated mother-infant pairs ###
  vaginal_unrelated <- check %>%
    filter((Type_1 == "VG" | Type_2 == "VG") &  
             mother_infant == "yes" &
             related == "unrelated")
  
  if (nrow(vaginal_related) > 0) {
    strain_summary_related <- vaginal_related %>%
      count(Strain_sharing)
    
    p_related <- ggplot(strain_summary_related, aes(x = Strain_sharing, y = n, fill = Strain_sharing)) +
      geom_bar(stat = "identity", color = "black", width = 0.6) +
      labs(
        title = paste0("Strain Sharing in Related VG Mother-Infant Pairs\n", species_name),  
        x = "Strain Sharing", y = "Number of Pairs"
      ) +
      theme_minimal(base_size = 14) +
      scale_fill_manual(values = c("yes" = "blue", "no" = "grey")) +
      theme(legend.position = "none")
    
    ggsave(
      filename = file.path(plot_dir_related, paste0(species_name, "_VG_strain_sharing_related_plot.pdf")),  
      plot = p_related,
      device = "pdf", width = 7, height = 5
    )
    
    ## Timepoint-based heatmap for related
    long_df <- vaginal_related %>%
      mutate(Pair_ID = row_number()) %>%
      select(
        Pair_ID,
        Sample_ID_1, Sample_ID_2,
        Timepoint_1 = Timepoint_categorical_1,
        Timepoint_2 = Timepoint_categorical_2,
        Type_1, Type_2,
        mother_infant_1, mother_infant_2,
        FAMILY = FAMILY_1,
        dist, Strain_sharing
      ) %>%
      pivot_longer(
        cols = starts_with("Sample_ID"),
        names_to = "Side",
        values_to = "Sample_ID"
      ) %>%
      mutate(
        Timepoint = ifelse(Side == "Sample_ID_1", Timepoint_1, Timepoint_2),
        Type = ifelse(Side == "Sample_ID_1", Type_1, Type_2),
        mother_infant = ifelse(Side == "Sample_ID_1", mother_infant_1, mother_infant_2)
      ) %>%
      select(Pair_ID, FAMILY, Sample_ID, Timepoint, Type, mother_infant, dist, Strain_sharing) %>%
      distinct() %>%
      mutate(Timepoint = factor(Timepoint, levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12")))
    
    
    long_df <- long_df %>%
      filter (!is.na(Timepoint))  # CHANGED TO VG
    
    heat_data_binary <- long_df %>%
      group_by(FAMILY, Timepoint) %>%
      summarise(Strain_shared = ifelse(any(Strain_sharing == "yes"), "yes", "no"), .groups = "drop")
    
    p_heat <- ggplot(heat_data_binary, aes(x = Timepoint, y = FAMILY, fill = Strain_shared)) +
      geom_tile(color = "white") +
      scale_fill_manual(values = c("yes" = "blue", "no" = "grey90")) +
      labs(
        title = paste0("Strain Sharing VG (Yes/No) Across Families and Timepoints\n", species_name),  
        x = "Timepoint", y = "Family"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )
    
    ggsave(
      filename = file.path(plot_dir_related, paste0(species_name, "_strain_sharing_VG_timepoint_heatmap.pdf")),  
      plot = p_heat,
      device = "pdf", width = 7, height = 5
    )
    
  } else {
    message("Skipping ", species_name, " (no related VG mother-infant pairs)")  
  }
  
  if (nrow(vaginal_unrelated) > 0) {
    strain_summary_unrelated <- vaginal_unrelated %>%
      count(Strain_sharing)
    
    p_unrelated <- ggplot(strain_summary_unrelated, aes(x = Strain_sharing, y = n, fill = Strain_sharing)) +
      geom_bar(stat = "identity", color = "black", width = 0.6) +
      labs(
        title = paste0("Strain Sharing in Unrelated VG Mother-Infant Pairs\n", species_name),  
        x = "Strain Sharing", y = "Number of Pairs"
      ) +
      theme_minimal(base_size = 14) +
      scale_fill_manual(values = c("yes" = "blue", "no" = "grey")) +
      theme(legend.position = "none")
    
    ggsave(
      filename = file.path(plot_dir_unrelated, paste0(species_name, "_VG_strain_sharing_unrelated_plot.pdf")),  
      plot = p_unrelated,
      device = "pdf", width = 7, height = 5
    )
  } else {
    message("Skipping ", species_name, " (no unrelated VG mother-infant pairs)")  
  }
  
  shared_related <- vaginal_related %>% filter(Strain_sharing == "yes")
  shared_unrelated <- vaginal_unrelated %>% filter(Strain_sharing == "yes")
  
  timepoints_shared_related <- if (nrow(shared_related) > 0) {
    paste(
      unique(c(shared_related$Timepoint_categorical_1,
               shared_related$Timepoint_categorical_2)) %>% na.omit(),
      collapse = ";"
    )
  } else {
    NA_character_
  }
  
  result_entry <- tibble(
    Species = species_name,
    N_Related = nrow(vaginal_related),
    N_Unrelated = nrow(vaginal_unrelated),
    N_Related_Shared = nrow(shared_related),
    N_Unrelated_Shared = nrow(shared_unrelated),
    Timepoints_Shared_Related = timepoints_shared_related
  )
  
  results_list[[length(results_list) + 1]] <- result_entry
  
} # end of for loop

results_all <- bind_rows(results_list)

write.table(results_all, "strain_sharing_VG_mother_infant.txt", row.names = FALSE, sep = "\t")  



