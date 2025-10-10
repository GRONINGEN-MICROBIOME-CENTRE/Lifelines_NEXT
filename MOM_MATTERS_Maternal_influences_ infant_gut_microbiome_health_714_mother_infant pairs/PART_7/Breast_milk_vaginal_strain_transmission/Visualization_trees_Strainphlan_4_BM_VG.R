library(ggtree)
library(tidyverse)
library(treeio)
library(phytools)

library(ggtree)
library(tidyverse)
library(ape)
library(RColorBrewer)
library(ggsci)

# Load metadata and set factor
meta <- read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata/Metadata_BM_VG_Gut_20_05_2025.txt")
meta$Type <- as.factor(meta$Type)

# Identify families with at least one VG sample
vg_families <- meta %>% 
  filter(Type == "VG") %>% 
  pull(FAMILY) %>% 
  unique()

# Filter metadata to only those families
filtered_metadata <- meta %>% 
  filter(FAMILY %in% vg_families)

# Color and shape mappings
type_colors <- c("Gut" = "#2ca02c", "Milk" = "#ff7f0e", "VG" = "#1f77b4")
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
  geom_tiplab(aes(color = FAMILY), size = 2, fontface = "bold", show.legend = FALSE) +
  geom_tippoint(aes(shape = mother_infant, fill = Type), size = 2, stroke = 1) +
  scale_color_manual(values = family_colors, guide = "none") +
  scale_fill_manual(values = type_colors) +
  scale_shape_manual(values = mother_infant_shapes) +
  theme(legend.position = "none")

print(all)

# Plot 2: Tree filtered to only families with VG samples
only_vg <- ggtree(tree_filtered, layout = "dendrogram") %<+% filtered_metadata +
  #geom_tiplab(aes(color = FAMILY), size = 2, fontface = "bold", show.legend = FALSE) +
  geom_tippoint(aes(shape = mother_infant, fill = Type), size = 4, stroke = 0.5) +
  scale_color_manual(values = family_colors, guide = "none") +
  scale_fill_manual(values = type_colors) +
  scale_shape_manual(values = mother_infant_shapes) +
  theme(legend.position = "none")

print(only_vg)



##### All PDF's #####################

library(ggtree)
library(tidyverse)
library(treeio)
library(phytools)
library(RColorBrewer)

# Read metadata
meta <- read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata/Metadata_BM_VG_Gut_20_05_2025.txt")
meta$Type <- as.factor(meta$Type)

# Get families with Milk samples
bm_families <- meta %>%
  filter(Type == "Milk") %>%
  pull(FAMILY) %>%
  unique()

# Filter metadata for those families
filtered_metadata <- meta %>%
  filter(FAMILY %in% bm_families)

# Tree directory and files
tree_dir <- "/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees"
tree_files <- list.files(tree_dir, pattern = "bestTree.*\\.tre$", full.names = TRUE)

# Colors
type_colors <- c("Gut" = "#2ca02c", "Milk" = "#ff7f0e", "VG" = "#1f77b4")
mother_infant_shapes <- c("mother" = 21, "infant" = 24)
mother_infant_fills <- c("mother" = "white", "infant" = "black")
family_list <- unique(filtered_metadata$FAMILY)
family_colors <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(family_list)), family_list)

# Output PDF
pdf("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees/all_trees_breastmilk_only.pdf", width = 8, height = 12)

for (tree_file in tree_files) {
  # Read and process tree
  tree <- read.tree(tree_file)
  tree <- midpoint.root(tree)
  tree$tip.label <- sub("_.*", "", tree$tip.label)
  
  tree_filtered <- drop.tip(tree, setdiff(tree$tip.label, filtered_metadata$NG_ID))
  
  # Plot
  p <- ggtree(tree_filtered, layout = "rectangular") %<+% filtered_metadata +
    geom_tippoint(aes(shape = mother_infant, fill = Type), size = 2, stroke = 1) +
    geom_tiplab(aes(color = FAMILY), size = 2, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = family_colors, guide = "none") +
    scale_fill_manual(values = type_colors) +
    scale_shape_manual(values = mother_infant_shapes) +
    ggtitle(basename(tree_file)) +
    theme(legend.position = "right")
  
  print(p)
}

dev.off()




library(ggtree)
library(tidyverse)
library(treeio)
library(phytools)

meta<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata/Metadata_BM_VG_Gut_20_05_2025.txt")
meta$Type<-as.factor(meta$Type)
summary (meta$Type)
vg_families <- meta %>%
  filter(Type == "VG") %>%
  pull(FAMILY) %>%
  unique()
filtered_metadata <- meta %>%
  filter(FAMILY %in% vg_families & mother_infant=="mother")

# Directory with trees
tree_dir <- "/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees"
tree_files <- list.files(tree_dir, pattern = "bestTree.*\\.tre$", full.names = TRUE)


type_colors <- c("Gut" = "#2ca02c", "Milk" = "#ff7f0e", "VG" = "#1f77b4")
mother_infant_shapes <- c("mother" = 21, "infant" = 24)
mother_infant_fills <- c("mother" = "white", "infant" = "black")

library(ggsci)
family_list <- unique(filtered_metadata$FAMILY)
family_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(family_list)),
                          family_list)

tree <- read.tree("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees/RAxML_bestTree.t__SGB17247.StrainPhlAn4.tre")
tree <- midpoint.root(tree)
tree$tip.label <- sub("_.*", "", tree$tip.label)
tree_filtered <- drop.tip(tree, setdiff(tree$tip.label, filtered_metadata$NG_ID))



p <- ggtree(tree_filtered, layout = "dendrogram") %<+% filtered_metadata +
  geom_tiplab(aes(color = FAMILY), size = 2, fontface = "bold", show.legend = FALSE) +  
  geom_tippoint(aes(shape = mother_infant, fill = Type), size = 2, stroke = 1) +
  scale_color_manual(values = family_colors, guide = "none") +  
  scale_fill_manual(values = type_colors) +
  scale_shape_manual(values = mother_infant_shapes) +
  theme(legend.position = "none")


p




##### All PDF's #####################

library(ggtree)
library(tidyverse)
library(treeio)
library(phytools)
library(RColorBrewer)

# Read metadata
meta <- read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata/Metadata_BM_VG_Gut_20_05_2025.txt")
meta$Type <- as.factor(meta$Type)

# Get families with Milk samples
bm_families <- meta %>%
  filter(Type == "Milk") %>%
  pull(FAMILY) %>%
  unique()

# Filter metadata for those families
filtered_metadata <- meta %>%
  filter(FAMILY %in% bm_families)

# Tree directory and files
tree_dir <- "/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees"
tree_files <- list.files(tree_dir, pattern = "bestTree.*\\.tre$", full.names = TRUE)

# Colors
type_colors <- c("Gut" = "#2ca02c", "Milk" = "#ff7f0e", "VG" = "#1f77b4")
mother_infant_shapes <- c("mother" = 21, "infant" = 24)
mother_infant_fills <- c("mother" = "white", "infant" = "black")
family_list <- unique(filtered_metadata$FAMILY)
family_colors <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(family_list)), family_list)

# Output PDF
pdf("/Users/trishlasinha/Desktop/LLNEXT/Analysis/Breastmilk_vaginal/Strainphlan_trees/all_trees_breastmilk_only.pdf", width = 8, height = 12)

for (tree_file in tree_files) {
  # Read and process tree
  tree <- read.tree(tree_file)
  tree <- midpoint.root(tree)
  tree$tip.label <- sub("_.*", "", tree$tip.label)
  
  tree_filtered <- drop.tip(tree, setdiff(tree$tip.label, filtered_metadata$NG_ID))
  
  # Plot
  p <- ggtree(tree_filtered, layout = "rectangular") %<+% filtered_metadata +
    geom_tippoint(aes(shape = mother_infant, fill = Type), size = 2, stroke = 1) +
    geom_tiplab(aes(color = FAMILY), size = 2, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = family_colors, guide = "none") +
    scale_fill_manual(values = type_colors) +
    scale_shape_manual(values = mother_infant_shapes) +
    ggtitle(basename(tree_file)) +
    theme(legend.position = "right")
  
  print(p)
}

dev.off()


