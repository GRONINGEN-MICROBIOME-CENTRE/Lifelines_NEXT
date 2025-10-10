### Sharing mother-infant all sources #####################
library(stringr)
library(tidyverse)

taxa_all<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/taxa/ALL_METAPHLAN_BM_VG_GUT_metaphlan_12_05_2025.txt")
row.names(taxa_all)<-taxa_all$clade_name
taxa_all$clade_name=NULL
taxa_all<-as.data.frame(t(taxa_all))
taxa_all$NG_ID <-row.names(taxa_all)
taxa_all$NG_ID <- str_sub(taxa_all$NG_ID, 1, 12)
row.names(taxa_all)<-taxa_all$NG_ID
taxa_all$NG_ID=NULL

taxa_SGB<-taxa_all[,grep("t__",colnames(taxa_all))]
colnames(taxa_SGB)=gsub(".*s__","",colnames(taxa_SGB))
names (taxa_SGB)
taxa_SGB$`Sphingomonas_sp_FARSPH|t__SGB24644`= NULL
taxa_SGB$`Microbacterium_ginsengisoli|t__SGB16909`=NULL
taxa_SGB$`Microbacterium_sp_T32|t__SGB16891`=NULL

# Now add metadata to this

meta<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata/Metadata_BM_VG_Gut_20_05_2025.txt")
row.names(meta)<-meta$NG_ID
summary (meta$Type)

all<-merge(meta, taxa_SGB, by="row.names")
row.names(all)<-all$Row.names
all$Row.names=NULL

# Make a plot of the summary stats 
all_RD <- all %>%
  pivot_longer(cols = c(raw_reads, human_reads, clean_reads),
               names_to = "Read_Type",
               values_to = "Count")
ggplot(all_RD, aes(x = Type, y = Count, fill = Read_Type)) +
  geom_boxplot() +
  labs(title = "Boxplot of Read Counts by Type",
       x = "Type",
       y = "Read Count") +
  theme_minimal()

# Now make plot with absence presence of taxa 
taxa_cols <- grep("t__", colnames(all), value = TRUE)

# Dichotamizing taxa
all <- all %>%
  mutate(across(all_of(taxa_cols), ~ ifelse(. > 0.01, 1, 0)))


# Step 1: Prepare data with taxa_cols
data_for_sharing <- all %>%
  select(FAMILY, Type, mother_infant, all_of(taxa_cols))



# Step 2: Split into mother and infant, using mother's Type
mothers <- data_for_sharing %>%
  filter(mother_infant == "mother") %>%
  rename_with(~ paste0(., "_mother"), taxa_cols) %>%
  rename(Type_mother = Type)

infants <- data_for_sharing %>%
  filter(mother_infant == "infant") %>%
  rename_with(~ paste0(., "_infant"), taxa_cols)

# Step 3: Join on FAMILY
paired <- inner_join(mothers, infants, by = "FAMILY")

# Step 4: Calculate sharing for each taxon
for (taxon in taxa_cols) {
  m_col <- paste0(taxon, "_mother")
  i_col <- paste0(taxon, "_infant")
  share_col <- paste0("shared_", taxon)
  
  paired[[share_col]] <- ifelse(paired[[m_col]] == 1 & paired[[i_col]] == 1, 1, 0)
}

# Step 5: Reshape to long format
sharing_long <- paired %>%
  select(FAMILY, Type_mother, starts_with("shared_")) %>%
  pivot_longer(cols = starts_with("shared_"),
               names_to = "taxa",
               names_prefix = "shared_",
               values_to = "shared")

percent_sharing <- sharing_long %>%
  group_by(Type_mother, taxa) %>%
  summarise(
    total_families = n_distinct(FAMILY),
    shared_families = n_distinct(FAMILY[shared == 1]),
    percent_shared = (shared_families / total_families) * 100,
    .groups = "drop"
  )

distinct_families_by_type <- paired %>%
  group_by(Type) %>%
  summarise(n_distinct_families = n_distinct(FAMILY))

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/New_supplementary_files")
#write.table(percent_sharing, "Sharing_SGBs_mother_infant_all_timepoints.txt", sep = "\t", row.names = F)




top_shared <- percent_sharing %>%
  group_by(Type_mother) %>%
  slice_max(order_by = percent_shared, n = 10) %>%
  ungroup()


top_shared <- top_shared %>%
  mutate(taxa_grouped = paste(Type_mother, taxa, sep = "_")) %>%
  arrange(Type_mother, percent_shared) %>%
  mutate(taxa_grouped = fct_inorder(taxa_grouped))


type_colors <- c("Gut" = "#a8f0a5", "Milk" = "#ff7f0e", "VG" = "#004b9a")

ggplot(top_shared, aes(x = taxa_grouped, y = percent_shared, fill = Type_mother)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ Type_mother, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = type_colors) +
  labs(
    title = "",
    x = "SGB",
    y = "Percent Shared"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

# Selecting strains for analysis 

Strainphlan_selection <- percent_sharing %>%
  filter(Type_mother %in% c("VG", "Milk"), shared_families> 2)

write.table(Strainphlan_selection, "BM_GUT_VG_SGB_for_Strainphlan.txt", sep = "\t", row.names = F)

Strainphlan_selection$taxa <- str_extract(Strainphlan_selection$taxa, "t__.*")
writeLines(Strainphlan_selection$taxa, "extracted_taxa.txt")

meta <- meta %>% rename(label = NG_ID)


# Building trees 
tree_dir <- "/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees"
tree_files <- list.files(tree_dir, pattern = "bestTree.*\\.tre$", full.names = TRUE)


type_colors <- c("Gut" = "#a8f0a5", "Milk" = "#ff7f0e", "VG" = "#004b9a")
mother_infant_shapes <- c("mother" = 21, "infant" = 24)
mother_infant_fills <- c("mother" = "white", "infant" = "black")

out_dir <- file.path(tree_dir, "PDF_PLOTS")
dir.create(out_dir, showWarnings = FALSE)

# Loop through each tree file and plot
for (file in tree_files) {
  tree <- read.tree(file)
  tree$tip.label <- sub("_.*", "", tree$tip.label)
  
  p <- ggtree(tree, layout = "circular") %<+% meta +
    geom_tiplab(aes(color = Type), size = 0.5, fontface = "bold") +
    geom_tippoint(aes(shape = mother_infant, fill = mother_infant, color = Type), size = 3, stroke = 1) +
    scale_color_manual(values = type_colors) +
    scale_fill_manual(values = mother_infant_fills) +
    scale_shape_manual(values = mother_infant_shapes) +
    ggtitle(basename(file)) +
    theme(legend.position = "right")
  
  pdf_name <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(file)), ".pdf"))
  ggsave(pdf_name, plot = p, width = 8, height = 8)
}
