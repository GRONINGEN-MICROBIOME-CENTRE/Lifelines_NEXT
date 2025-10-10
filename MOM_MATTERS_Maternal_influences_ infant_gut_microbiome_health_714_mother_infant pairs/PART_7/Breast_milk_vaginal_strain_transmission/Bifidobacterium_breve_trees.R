######  Making B. breve zoomed in trees ###########################
library(ggtree)
library(tidyverse)
library(treeio)
library(phytools)
library(ggnewscale)

# Load metadata and set factor
meta <- read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata/Metadata_BM_VG_Gut_20_05_2025.txt")
meta$Timepoint_categorical[meta$mother_infant == "mother"] <- "mother"       
meta$Type <- as.factor(meta$Type)
meta <- na.omit(meta)
meta$Timepoint_categorical=factor(meta$Timepoint_categorical, levels = c("mother", "W2", "M1", "M2", "M3", "M6", 'M9', "M12"))

meta$Timepoint_early_late <- dplyr::case_when(
  meta$Timepoint_categorical %in% c("W2", "M1", "M2", "M3") ~ "early",
  meta$Timepoint_categorical %in% c("M6", "M9", "M12") ~ "late",
  meta$Timepoint_categorical == "mother" ~ "mother",
  TRUE ~ NA_character_
)


early_late_colors <- c(
  "mother"   = "white", # soft green
  "early"  = "#4397bb", # light-medium green
  "late" = "#390962"  # very deep green
)


# Identify families with at least one VG sample
vg_BM_families <- meta %>% 
  filter(Type == "VG" | Type == "Milk") %>% 
  pull(FAMILY) %>% 
  unique()


# Filter metadata to only those families
filtered_metadata <- meta %>% 
  filter(FAMILY %in% vg_BM_families)
length(unique(meta$FAMILY))
length(unique(filtered_metadata $FAMILY))

# Color and shape mappings
type_colors <- c(
  "Gut"  = "#b5dd88",  # very light mint green # Colors same as the colors for the supplementary table 
  "Milk" = "#e6550d",  # strong orange
  "VG"   = "#004b9a"   # deep blue
)



mother_infant_shapes <- c("mother" = 21, "infant" = 24)
mother_infant_fills <- c("mother" = "white", "infant" = "black")

# Generate distinct colors for families
family_list <- unique(filtered_metadata$FAMILY)
family_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(family_list)),
  family_list
)

# Load tree and prepare
tree_path <- "/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees/RAxML_bestTree.t__SGB17247.StrainPhlAn4.tre"
tree <- read.tree(tree_path)
tree <- midpoint.root(tree)

# Simplify tip labels to IDs only (remove suffix after underscore)
tree$tip.label <- sub("_.*", "", tree$tip.label)

# Filter tree to only tips present in filtered metadata
tree_filtered <- drop.tip(tree, setdiff(tree$tip.label, filtered_metadata$NG_ID))

# Plot 1: Full tree with all tips and points
all <- ggtree(tree, layout = "dendrogram") %<+% meta +
  geom_tiplab(aes(label = FAMILY, color = FAMILY), size = 1, show.legend = FALSE)+
  geom_tippoint(aes(shape = mother_infant, fill = Type), size = 5, stroke = 0.5) +
  scale_color_manual(values = family_colors, guide = "none") +
  scale_fill_manual(values = type_colors) +
  scale_shape_manual(values = mother_infant_shapes) +
  theme(legend.position = "none")

print(all)



# Plot 2: Tree filtered to only families with VG samples
only_vg_BM <- ggtree(tree_filtered, layout = "dendrogram") %<+% filtered_metadata +
  geom_tippoint(
    aes(
      shape = mother_infant,
      fill  = Timepoint_early_late,
      color = Type , 
    ),
    size = 4,
    stroke = 2
  ) +
  scale_color_manual(values = type_colors) +   
  scale_fill_manual(values = early_late_colors) +         
  scale_shape_manual(values = mother_infant_shapes) +
  new_scale_color() +
  #geom_tiplab(aes(label = FAMILY, color = FAMILY), size = 2, show.legend = FALSE)+
  theme(legend.position = "right")  

print(only_vg_BM)




setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/supp_figures/Parts_of_figures/Breve_tree_sharing")

ggsave(
  filename = "only_vg_BM.pdf",
  plot = only_vg_BM,
  device = "pdf",
  width = 16,
  height = 10,
  units = "in"
)



# Making zoomed in subtrees of desired families (Consisting of a BM and VG sample)

Do_subtree = function(Family, meta, tree){
  meta %>% filter(FAMILY == Family) -> submeta
  submeta %>% filter(NG_ID %in%  tree$tip.label) -> submeta
  keep.tip(tree, submeta$NG_ID) -> subtree
  return(subtree)  
  
}

target_families <- c("FAM0373", "FAM0362", "FAM0358", 
                     "FAM0121", "FAM0112", "FAM0554", "FAM0468")


for(fam in target_families){
  
  # make the zoomed subtree
  p <- Do_subtree(fam, meta, tree) %>%
    ggtree(layout = "dendrogram") %<+% meta +
    geom_tippoint(
      aes(shape = mother_infant, fill = Timepoint_early_late, color = Type ),
      size = 10,
      stroke = 4
    ) +
    geom_tiplab(
      aes(label = Timepoint_categorical),
      size = 4,
      fontface = "bold",
      show.legend = FALSE,
      vjust = 3,
      angle = 0
    ) +
    scale_color_manual(values = type_colors) + 
    scale_fill_manual(values = early_late_colors) +
    scale_shape_manual(values = mother_infant_shapes) +
    theme(legend.position = "none")
  
  # save as PDF with the family name
  ggsave(
    filename = paste0(fam, "_raw.pdf"),
    plot = p,
    device = "pdf",
    width = 5,
    height = 4,
    units = "in"
  )
}


