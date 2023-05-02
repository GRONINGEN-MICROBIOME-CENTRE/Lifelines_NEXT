setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the differences in bacteriome and virome 
# communitues between infant and maternal samples at the 
# species and vOTUs levels USING BRAY-CURTIS DISTANCE
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################
library(vegan)

library(ggplot2)
library(ggExtra)
library(ggrepel)
##############################
# Input data
##############################
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_with_phenos.txt', sep='\t', header=T)

MGS_metadata$Type <- factor(MGS_metadata$Type, levels=c('Infant', 'Mother'), ordered=T)

microbiome <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_with_phenos.txt', sep='\t', header=T)

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header = T)
VLP_metadata$Type <- factor(VLP_metadata$Type, levels=c('Infant', 'Mother'), ordered=T)
##############################
# ANALYSIS
##############################

#### NMDS for bacteriome

# preapring phenotypes
for_bacplot <- MGS_metadata
row.names(for_bacplot) <- MGS_metadata$Short_sample_ID_bact
for_bacplot <- for_bacplot[colnames(microbiome),]
#CHOOSING TECHNICAL PHENOTYPES, SAMPLE TYPE & BACTERIAL ALPHA DIVERSITY
for_bacplot <- for_bacplot[,c("Type", "Clean_reads", "DNA_CONC", "bacterial_alpha_diversity", "metaphlan_unknown_perc")]

ord <- metaMDS(t(microbiome), distance = "bray", k=2)
en = envfit(ord, for_bacplot, permutations = 999, na.rm = TRUE)
en$factors #r2 0.7416; p-value < 0.001
en$vectors

centroids <- as.data.frame(scores(en, "factors"))
centroids$Type <- c(gsub('Type', '', row.names(centroids)))

spp.scrs <- as.data.frame(scores(en, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Type <- 'Infant'
spp.scrs$Species <- c('N clean reads', 'DNA concentration', 'Alpha diversity (bacteria)', '% unknown in Metaphlan4')

### for plotting:
data.scores = as.data.frame(scores(ord, "sites"))
data.scores <- merge(data.scores, for_bacplot, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL

pdf('./04.PLOTS/Bacterial_spp_Bray_NMDS.pdf', width=10/2.54, height=10/2.54)
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Type)) + 
  geom_point(size = 2, alpha=0.8) + 
  geom_point(data=centroids, aes(fill=Type),shape=23, size=4, color='black', ) + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Type, color=Type), linetype = 2)+
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_label_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', size = 2) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        legend.position = "bottom") + 
  scale_color_manual(values=c("#F8766D", "#00BFC4"))

ggMarginal(gg, type="densigram", groupFill=T)
dev.off()

### calculating the significance of separation:
if (identical(colnames(microbiome), MGS_metadata$Short_sample_ID_bact)) {
  bacteria_permanova <- adonis2(t(microbiome) ~ Type +  bacterial_alpha_diversity + Clean_reads + DNA_CONC + metaphlan_unknown_perc, 
                                data=MGS_metadata,  
                                permutations = 999,
                                method="bray",
                                by = "margin")
}

#### NMDS for bacteriome

# preapring phenotypes
for_virplot <- VLP_metadata
row.names(for_virplot) <- VLP_metadata$Short_sample_ID
for_virplot <- for_virplot[colnames(RPKM_counts_VLP),]
#CHOOSING TECHNICAL PHENOTYPES, SAMPLE TYPE & VIRAL ALPHA DIVERSITY
for_virplot <- for_virplot[,c("Type", "Clean_reads", "DNA_CONC", "viral_alpha_diversity", "bacterial_contamination_perc_reads")]

#### yes, I use the same object names and don't make a function out of it
ord <- metaMDS(t(RPKM_counts_VLP), distance = "bray", k=2)
en = envfit(ord, for_virplot, permutations = 999, na.rm = TRUE)
en$factors #r2 0.33; p-value < 0.001
en$vectors

centroids <- as.data.frame(scores(en, "factors"))
centroids$Type <- c(gsub('Type', '', row.names(centroids)))

spp.scrs <- as.data.frame(scores(en, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Type <- 'Infant'
spp.scrs <- spp.scrs[spp.scrs$Species!='Clean_reads',] # because it was not significant in envfit
spp.scrs$Species <- c('DNA concentration', 'Alpha diversity (viruses)', '% bacterial contamination')

### for plotting:
data.scores = as.data.frame(scores(ord, "sites"))
data.scores <- merge(data.scores, for_virplot, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL

pdf('./04.PLOTS/Viral_vOTUs_Bray_NMDS.pdf', width=10/2.54, height=10/2.54)
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Type)) + 
  geom_point(size = 2, alpha=0.8) + 
  geom_point(data=centroids, aes(fill=Type),shape=23, size=4, color='black', ) + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Type, color=Type), linetype = 2)+
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_label_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', alpha = 0.8, size = 2) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        legend.position = "bottom") + 
  scale_color_manual(values=c("#F8766D", "#00BFC4"))

ggMarginal(gg, type="densigram", groupFill=T)
dev.off()


### calculating the significance of separation:
RPKM_counts_VLP <- RPKM_counts_VLP[,VLP_metadata$Short_sample_ID]
if (identical(colnames(RPKM_counts_VLP), VLP_metadata$Short_sample_ID)) {
  viruses_permanova <- adonis2(t(RPKM_counts_VLP) ~ Type +viral_alpha_diversity + Clean_reads + DNA_CONC + bacterial_contamination_perc_reads, 
                                data=VLP_metadata,  
                                permutations = 999,
                                method="bray",
                                by = "margin")
}



###### OUTPUT #####
write.table(bacteria_permanova, '03a.RESULTS/Bacterial_spp_bray_PERMANOVA.txt', sep='\t', quote=F)
write.table(bacteria_permanova, '03a.RESULTS/Viral_vOTUs_bray_PERMANOVA.txt', sep='\t', quote=F)
