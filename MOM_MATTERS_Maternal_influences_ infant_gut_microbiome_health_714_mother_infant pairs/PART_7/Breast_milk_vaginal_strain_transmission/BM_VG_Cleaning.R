# Creating of metadata file 
library(tidyverse)
library(stringr)
library(purrr)
library(ggplot)
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal")
metadata<-read.delim("metadata_VG_BM.txt")

# KN
kneaddata <-read.delim("LLNEXT_KN_STATS_173_09_05_2025.txt")
kneaddata$NG_ID <- str_sub(kneaddata$KN_ID, 1, 12)
kneaddata$NG_ID <- sub("_.*", "", kneaddata$NG_ID)


# Merge with metadata 
KN_meta<-left_join(metadata, kneaddata)


df_long <- KN_meta %>%
  pivot_longer(
    cols = c("Raw_reads_P1", "Raw_reads_P2",  
             "Contaminants_P1", "Contaminants_P2", 
             "Clean_reads_P1", "Clean_reads_P2"),
    names_to = "Step",
    values_to = "Value"
  ) %>%
  group_by(Step) %>%
  mutate(Step = reorder(Step, Value, FUN = median)) %>%
  ungroup()

# Define consistent color palette
custom_colors <- c("#e60049", "#0bb4ff", "#50e991", "#f46a9b", "#9b19f5", "#ffa300",
                   "#dc0ab4", "#b3d4ff", "#00bfa0", "#b30000", "#7c1158", "#4421af")


ggplot(df_long, aes(x = Step, y = Value, fill = Type, color = Type)) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  geom_boxplot(alpha = 0.4) +
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "Read Counts",
    fill = "Type",
    color = "Type"
  ) +
  theme(
    plot.title = element_text(color = "black", size = 22, face = "bold"),
    axis.title.x = element_text(color = "black", size = 22, face = "bold"),
    axis.title.y = element_text(color = "black", size = 18, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 22, angle = 60, hjust = 1)
 
  )

summary_stats_by_type <- KN_meta %>%
  select(Raw_reads_P1, Raw_reads_P2,  
         Contaminants_P1, Contaminants_P2, Contaminants_O1, Contaminants_O2,
         Clean_reads_P1, Clean_reads_P2, Type) %>%
  group_by(Type) %>%
  summarise(across(c(Raw_reads_P1, Raw_reads_P2,  
                     Contaminants_P1, Contaminants_P2, Contaminants_O1, Contaminants_O2,
                     Clean_reads_P1, Clean_reads_P2), 
                   list(
                     non_missing = ~sum(!is.na(.)),
                     missing = ~sum(is.na(.)),
                     min = ~min(., na.rm = TRUE),
                     median = ~median(., na.rm = TRUE),
                     mean = ~mean(., na.rm = TRUE),
                     max = ~max(., na.rm = TRUE)
                   ), .names = "{.col}_{.fn}"))

summary_stats_by_type %>%
  select(Type, Raw_reads_P1_mean) %>%
  print()

summary_stats_by_type %>%
  select(Type, Clean_reads_P1_mean) %>%
  print()

summary_stats_by_type %>%
  select(Type, Clean_reads_P1_median) %>%
  print()

ggplot(KN_meta, aes(x =Type , y = Clean_reads_P1, fill = Type, color = Type)) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  geom_boxplot(alpha = 0.4) +
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "Read Counts",
    fill = "Type",
    color = "Type"
  ) +
  theme(
    plot.title = element_text(color = "black", size = 22, face = "bold"),
    axis.title.x = element_text(color = "black", size = 22, face = "bold"),
    axis.title.y = element_text(color = "black", size = 18, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 22, angle = 60, hjust = 1)
    
  )

# Metaphlan data 
BM_KN<-read.delim("BM_VG_final_173_09_05_2025.txt")
row.names(BM_KN)<-BM_KN$clade_name
BM_KN$clade_name=NULL
tax<-as.data.frame(t(BM_KN)) 
row.names(tax)<- str_sub(row.names(tax), 1, 12)

summary (tax$UNCLASSIFIED)
tax$NG_ID<-row.names(tax)
tax_meta<-left_join(KN_meta, tax)


library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)


tax_meta_long <- tax_meta %>%
  pivot_longer(
    cols = 28:ncol(tax_meta),
    names_to = "Taxon",
    values_to = "Abundance"
  ) %>%
  mutate(
    Taxonomic_Level = case_when(
      str_detect(Taxon, "t__") ~ "SGB",
      str_detect(Taxon, "s__") ~ "Species",
      str_detect(Taxon, "g__") ~ "Genus",
      str_detect(Taxon, "f__") ~ "Family",
      str_detect(Taxon, "o__") ~ "Order",
      str_detect(Taxon, "c__") ~ "Class",
      str_detect(Taxon, "p__") ~ "Phylum",
      str_detect(Taxon, "k__") ~ "Kingdom",
      TRUE ~ "Unclassified"
    ),
    Clean_Taxon = case_when(
      Taxonomic_Level == "SGB" ~ str_extract(Taxon, "t__[^|]+"),
      Taxonomic_Level == "Species" ~ str_extract(Taxon, "s__[^|]+"),
      Taxonomic_Level == "Genus" ~ str_extract(Taxon, "g__[^|]+"),
      Taxonomic_Level == "Family" ~ str_extract(Taxon, "f__[^|]+"),
      Taxonomic_Level == "Order" ~ str_extract(Taxon, "o__[^|]+"),
      Taxonomic_Level == "Class" ~ str_extract(Taxon, "c__[^|]+"),
      Taxonomic_Level == "Phylum" ~ str_extract(Taxon, "p__[^|]+"),
      Taxonomic_Level == "Kingdom" ~ str_extract(Taxon, "k__[^|]+"),
      TRUE ~ "Unclassified"
    )
  )


tax_summary <- tax_meta_long %>%
  group_by(Type, Timepoint_categorical, Taxonomic_Level, Clean_Taxon) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")
tax_summary_filtered <- tax_summary %>%
  filter(Mean_Abundance >= 1)

tax_summary_split <- tax_summary_filtered %>%
  split(.$Taxonomic_Level)

plot_tax_level <- function(df, level) {
  ggplot(df, aes(x = Clean_Taxon, y = Mean_Abundance, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    # facet_wrap(~ Timepoint_categorical) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 10),  # changed from x to y since flipped
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    labs(
      title = paste("Mean Abundance at", level, "Level"),
      x = "Taxon",
      y = "Mean Abundance"
    ) +
    coord_flip()
}
tax_plots <- imap(tax_summary_split, plot_tax_level)
tax_plots$Genus
tax_plots$Species

tax_summary_split[[8]] %>%
  arrange(desc(Mean_Abundance)) %>%
  select(Clean_Taxon  , Mean_Abundance) %>%
  print(n = 31)                

##### Negative control check ######
abundances <- tax %>%
  filter(NG_ID == "DSMKVGco2_MK") %>%    # Replace 'NG_ID' with your sample ID column
  select(-NG_ID)                         # Remove the ID column to keep only abundance values
abundances_long <- abundances %>%
  pivot_longer(
    cols = everything(),
    names_to = "Taxon",
    values_to = "Abundance"
  )
high_abundance <- abundances_long %>%
  filter(Abundance > 0.1) %>%
  arrange(desc(Abundance))
ggplot(high_abundance, aes(x = reorder(Taxon, Abundance), y = Abundance)) +
  geom_bar(stat = "identity", fill = "#69b3a2") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "DSMKVGco2_MK",
    x = "Taxon",
    y = "Relative Abundance"
  )





vg_data <- tax_meta %>%
  filter(Type == "VG")

# Find the B. breve column name
breve_col <- grep("breve", colnames(vg_data), value = TRUE)[1]

only_breve <- vg_data[, breve_col, drop = FALSE]

ggplot(vg_data, aes(x = .data[[breve_col]])) +
  geom_histogram(aes(y = ..count..), fill = "skyblue", color = "black", bins = 30) +
  stat_bin(aes(y = ..count.., label = ..count..), 
           bins = 30, geom = "text", vjust = -0.5, size = 3) +
  labs(
    title = paste("", breve_col, "in Type == VG"),
    x = "Bifidobacterium breve Abundance",
    y = "Frequency"
  ) +
  theme_minimal()


tax_meta_long_clean <- tax_meta_long %>%
  filter(!grepl("Sphingomonas|Microbacterium", Taxon))

plot_taxonomic_level <- function(data, level = "Genus", abundance_threshold = 5, sample_type = "Milk") {
  # Filter the data based on input arguments
  filtered_data <- data %>%
    filter(
      Taxonomic_Level == level,
      Abundance >= abundance_threshold,
      Type == sample_type
    ) %>%
    mutate(Sample_ID = paste(NEXT_ID, Timepoint_categorical, sep = "-"))
  
  # Create the stacked bar plot
  p <- ggplot(filtered_data, aes(x = Sample_ID, y = Abundance, fill = Clean_Taxon)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    labs(
      x = "Sample ID (NEXT_ID-Timepoint)",
      y = "Relative Abundance",
      fill = "Clean Taxon",
      title = paste("Relative Abundance of", level, "in", sample_type, "Samples")
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(p)
}
# Genus Milk 
plot_taxonomic_level(tax_meta_long_clean, level = "Genus", sample_type = "Milk", abundance_threshold = 5)

# Species Milk 
plot_taxonomic_level(tax_meta_long_clean, level = "Species", sample_type = "Milk", abundance_threshold = 5)

# Species vaginal 
plot_taxonomic_level(tax_meta_long_clean, level = "Species", sample_type = "VG", abundance_threshold = 5)

# Making final metata 
write.table(KN_meta, "KN_meta_BM_VG_DS.txt", row.names = F, sep = "\t")

