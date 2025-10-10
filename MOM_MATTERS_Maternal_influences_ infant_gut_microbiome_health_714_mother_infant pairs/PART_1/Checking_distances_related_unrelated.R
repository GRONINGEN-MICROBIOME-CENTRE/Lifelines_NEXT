# Checking distances of related versus unrelated 
# Author: T. Sinha 
# Date: 14th of August, 2025


library(readr)
library(tools)
library(dplyr)
library(cowplot)


# Initialize results list
results_list <- list()

# Initialize results lists for Milk and VG
results_list_milk <- list()
results_list_vg <- list()

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees/Processed_RDS")

# List all .rds files
rds_files <- list.files(pattern = "^t__.*\\.rds$")

# Define function to perform permutation test and calculate stats
perform_test <- function(df, species, group_label) {
  if (nrow(df) == 0) {
    # If no data, return NA results
    return(data.frame(
      species = species,
      group = group_label,
      n_related = NA,
      n_unrelated = NA,
      mean_related = NA,
      mean_unrelated = NA,
      obs_diff_mean = NA,
      sd_related = NA,
      sd_unrelated = NA,
      median_related = NA,
      median_unrelated = NA,
      obs_diff_median = NA,
      p_value = NA
    ))
  }
  
  # Counts
  n_related <- sum(df$related == "related")
  n_unrelated <- sum(df$related == "unrelated")
  
  # Means
  mean_related   <- mean(df$dist[df$related == "related"])
  mean_unrelated <- mean(df$dist[df$related == "unrelated"])
  obs_diff_mean <- mean_related - mean_unrelated
  
  # SDs
  sd_related   <- sd(df$dist[df$related == "related"])
  sd_unrelated <- sd(df$dist[df$related == "unrelated"])
  
  # Medians
  median_related   <- median(df$dist[df$related == "related"])
  median_unrelated <- median(df$dist[df$related == "unrelated"])
  obs_diff_median <- median_related - median_unrelated
  
  # Permutation test
  n_perm <- 10000
  set.seed(42)
  perm_diffs <- numeric(n_perm)
  print (n_perm)
  
  for (i in 1:n_perm) {
    perm_labels <- sample(df$related)
    perm_mean_related   <- mean(df$dist[perm_labels == "related"])
    perm_mean_unrelated <- mean(df$dist[perm_labels == "unrelated"])
    perm_diffs[i] <- perm_mean_related - perm_mean_unrelated
  }
  
  p_value <- mean(abs(perm_diffs) >= abs(obs_diff_mean))
  
  data.frame(
    species = species,
    group = group_label,
    n_related = n_related,
    n_unrelated = n_unrelated,
    mean_related = mean_related,
    mean_unrelated = mean_unrelated,
    obs_diff_mean = obs_diff_mean,
    sd_related = sd_related,
    sd_unrelated = sd_unrelated,
    median_related = median_related,
    median_unrelated = median_unrelated,
    obs_diff_median = obs_diff_median,
    p_value = p_value
  )
}

# Loop through each .rds file
for (file in rds_files) {
  
  species_name <- file_path_sans_ext(basename(file))
  
  df_check <- read_rds(file)
  df_check <- na.omit(df_check)
  
  # Filter for Milk-related pairs
  distance_test_bm <- df_check %>%
    filter((Type_1 == "Milk" | Type_2 == "Milk") & mother_infant == "yes")
  
  # Filter for VG-related pairs
  distance_test_vg <- df_check %>%
    filter((Type_1 == "VG" | Type_2 == "VG") & mother_infant == "yes")
  
  result_milk <- perform_test(distance_test_bm, species_name, "Milk")
  result_vg <- perform_test(distance_test_vg, species_name, "VG")
  
  results_list_milk[[species_name]] <- result_milk
  results_list_vg[[species_name]] <- result_vg
}

# Combine results for Milk and VG separately
final_results_milk <- do.call(rbind, results_list_milk)
final_results_vg <- do.call(rbind, results_list_vg)

# Write to separate files
write.table(final_results_milk, "all_species_related_versus_unrelated_results_milk.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(final_results_vg, "all_species_related_versus_unrelated_results_vg.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# Get species names 

species<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/link_SGB_species_names.txt")
row.names(species)<-species$SGB


### First making plot for breastmilk 

final_results_milk_fil <- final_results_milk[!is.na(final_results_milk$n_related) &
                                               final_results_milk$n_related > 2, ]
final_results_milk_fil<-merge(species, final_results_milk_fil, by="row.names")
final_results_milk_fil$Row.names=NULL
final_results_milk_fil$species_SGB<-paste0(final_results_milk_fil$Species, "_", final_results_milk_fil$SGB)
final_results_milk_fil$species_SGB <- gsub("_t__", " ", final_results_milk_fil$species_SGB)
final_results_milk_fil$species_SGB <- gsub("_", " ", final_results_milk_fil$species_SGB)


# Order species by obs_diff
distances_BM <- final_results_milk_fil  %>%
  mutate(species = factor(species_SGB, levels = species_SGB[order(-p_value)]))

# Reshape to long format for plotting
df_plot <- distances_BM %>%
  pivot_longer(
    cols = c(mean_related, mean_unrelated),
    names_to = "Group",
    values_to = "Mean"
  )

# Prepare significance star positions
sig_labels <- distances_BM %>%
  mutate(
    label = ifelse(p_value < 0.05, "***", ""),
    x_pos = pmax(mean_related, mean_unrelated) * 1.05
  ) %>%
  select(species, label, x_pos)


p_milk<-ggplot(df_plot, aes(y = species_SGB, x = Mean, fill = Group)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c(
    "mean_related" = "#ff7f0e",   
    "mean_unrelated" = "#efb386"  
  )) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Mean Distance",
    y = "SGB"
  ) +
  geom_text(
    data = sig_labels,
    aes(x = x_pos, y = species, label = label),
    inherit.aes = FALSE,
    size = 6
  )


# Vaginal 


### First making plot for VG ###

# Filter VG results
final_results_vg_fil <- final_results_vg[!is.na(final_results_vg$n_related) &
                                           final_results_vg$n_related > 2, ]

# Merge with species info
final_results_vg_fil <- merge(species, final_results_vg_fil, by = "row.names")
final_results_vg_fil$Row.names <- NULL

# Create combined species_SGB column
final_results_vg_fil$species_SGB <- paste0(final_results_vg_fil$Species, "_", final_results_vg_fil$SGB)
final_results_vg_fil$species_SGB <- gsub("_t__", " ", final_results_vg_fil$species_SGB)
final_results_vg_fil$species_SGB <- gsub("_", " ", final_results_vg_fil$species_SGB)

# Order species by mean_related
distances_VG <- final_results_vg_fil %>%
  mutate(species = factor(species_SGB, levels = species_SGB[order(-p_value)]))

# Reshape to long format for plotting
df_plot_vg <- distances_VG %>%
  pivot_longer(
    cols = c(mean_related, mean_unrelated),
    names_to = "Group",
    values_to = "Mean"
  )

# Prepare significance star positions
sig_labels_vg <- distances_VG %>%
  mutate(
    label = ifelse(p_value < 0.05, "***", ""),
    x_pos = pmax(mean_related, mean_unrelated) * 1.05
  ) %>%
  select(species, label, x_pos)


p_vg <-ggplot(df_plot_vg, aes(y = species_SGB, x = Mean, fill = Group)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c(
    "mean_related" = "#1f77b4",   
    "mean_unrelated" = "#88c7ea"  
  )) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Mean Distance",
    y = "SGB"
  ) +
  geom_text(
    data = sig_labels_vg,
    aes(x = x_pos, y = species, label = label),
    inherit.aes = FALSE,
    size = 6
  )


plot_grid(p_milk, p_vg)
