setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Phylogenetic tree
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################


library(ape)
library(ggplot2)
library(ggtree)
library(dplyr)
library(tidytree)

library(MetBrewer)

##############################
# Input data
##############################
metadata <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/metadata_combined_for_exp.txt', sep='\t', header=T)
metadata$Alex_ID <- paste0(metadata$FAM_ID, '_', metadata$Type, '_', substr(metadata$Short_sample_ID, 1,1), '_', metadata$source, '_', metadata$Timepoint)
row.names(metadata) <- metadata$Alex_ID
metadata$label <- metadata$Alex_ID

RPKM_combined_0.95 <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/RPKM_counts_combined_0.95_UPD_final_for_exp.txt', sep='\t', header=T)

folder = "02.CLEAN_DATA/dist.matrices/"  
file_list = list.files(path=folder, pattern="*.txt")  
virus=lapply(paste0(folder, file_list), function(x) read.table(x, sep='\t', header=T))
names(virus) <- gsub(".dist.txt", '', file_list)

MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata$source <- "MGS"
MGS_metadata$Alex_ID <- paste0(MGS_metadata$FAM_ID, '_', MGS_metadata$Type, '_', substr(MGS_metadata$Short_sample_ID_bact, 1,1), '_', MGS_metadata$source, '_', MGS_metadata$Timepoint)

folder = "01.RAW_DATA/strainphlan_4_distance_matrix/"  
file_list = list.files(path=folder, pattern="*.csv")  
bacterium=lapply(paste0(folder, file_list), function(x) read.csv(x, header=T))
names(bacterium) <- gsub("_dmat_.csv", '', file_list)

## since my parsing scripts have certain ID format: 

bacterium <- lapply(bacterium, function(x) {
  row.names(x) <- substr(row.names(x), 0, 12) # first remove data processing artifacts
  colnames(x) <- substr(colnames(x), 0, 12) 
  
  row.names(x) <- MGS_metadata$Alex_ID[match(row.names(x), MGS_metadata$NG_ID)]
  colnames(x) <- MGS_metadata$Alex_ID[match(colnames(x), MGS_metadata$NG_ID)] # then change the sequencing IDs to IDs for parsing
  x
})

bacterium_to_test <- list()
for (n in 1:NROW(bacterium)) {
  
  bacteriumN <- bacterium[[n]]
  bacteriumName <- names(bacterium[n])
  
  if ( length(grep('Infant', colnames(bacteriumN)))!=0 & length(grep('Mother', colnames(bacteriumN)))!=0) {
    bacterium_to_test[[bacteriumName]] <- bacterium[[bacteriumName]]
  }
  
}

Cross <- met.brewer("Cross")

metaphlan <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)
metaphlan_species <- metaphlan[grep('s__', row.names(metaphlan)),]

virus_host_transmission <- read.table("03a.RESULTS/Bacterial_hosts_transmission_for_transmitted_viruses.txt", sep='\t', header=T)
virus_host_transmission_reconstructed <- virus_host_transmission[!is.na(virus_host_transmission$easy_name),]


virus_metadata <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/Viruses_shared_min_5_families_UPD_final.txt', sep='\t', header=T)

##############################
# ANALYSIS
##############################

#### VIRUS
dist_85k <- virus[["LN_6A08_VL_306_NODE_3_length_85266_cov_2453.209609"]]
tree <- nj(as.matrix(dist_85k))
plot(tree, cex=0.5)


p <- ggtree(tree, options(ignore.negative.edge = T))

p[["data"]]$Type <- metadata$Type[match(p[["data"]]$label, metadata$Alex_ID)]
p[["data"]]$FAMILY <- metadata$FAM_ID[match(p[["data"]]$label, metadata$Alex_ID)]
p[["data"]]$source <- metadata$source[match(p[["data"]]$label, metadata$Alex_ID)]
p[["data"]]$Timepoint <- metadata$Timepoint[match(p[["data"]]$label, metadata$Alex_ID)]


a <- unique(p[["data"]]$FAMILY)
a <- a[!is.na(a)]

pdf('./04.PLOTS/L_85266_LS0_tree.pdf', width=15/2.54, height=12/2.54)
p+ geom_tippoint(aes(shape = Type, fill = FAMILY, color = source), size = 2, stroke=1) +
  scale_shape_manual(values=c("Mother" = 21, "Infant"=24)) + 
  scale_fill_manual(values=setNames(Cross[8:2], a)) + 
  scale_color_manual(values = c("VLP"="black", "MGS"="brown")) + 
  guides(shape = guide_legend(override.aes = list(color = NULL)),
         fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  geom_tiplab(aes(label=Timepoint), linetype = NULL, offset=.001) + 
  ggtitle("L_85266_LS0")
dev.off()


#### BACTERIAL HOST
dist_85k_host <- bacterium_to_test[["SGB1836_group"]]
# only samples from the same families:
dist_85k_host <- dist_85k_host[grep(paste(a,collapse="|"), 
                                    colnames(dist_85k_host)),
                               grep(paste(a,collapse="|"), 
                                    colnames(dist_85k_host))]



bactree <- nj(as.matrix(dist_85k_host))

bp <- ggtree(bactree, options(ignore.negative.edge = T))

bp[["data"]]$Type <- metadata$Type[match(bp[["data"]]$label, metadata$Alex_ID)]
bp[["data"]]$FAMILY <- metadata$FAM_ID[match(bp[["data"]]$label, metadata$Alex_ID)]
bp[["data"]]$source <- metadata$source[match(bp[["data"]]$label, metadata$Alex_ID)]
bp[["data"]]$Timepoint <- metadata$Timepoint[match(bp[["data"]]$label, metadata$Alex_ID)]


pdf('./04.PLOTS/Host_L_85266_LS0_tree.pdf', width=15/2.54, height=18/2.54)
bp+ geom_tippoint(aes(shape = Type, fill = FAMILY, color = source), size = 2, stroke=1) +
  scale_shape_manual(values=c("Mother" = 21, "Infant"=24)) + 
  scale_fill_manual(values=setNames(Cross[8:2], a)) + 
  scale_color_manual(values = c("VLP"="black", "MGS"="brown")) + 
  guides(shape = guide_legend(override.aes = list(color = NULL)),
         fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  geom_tiplab(aes(label=Timepoint), linetype = NULL, offset=.4) + 
  ggtitle("Bacteroides uniformis")
dev.off()










##############################
# OUTPUT
##############################


