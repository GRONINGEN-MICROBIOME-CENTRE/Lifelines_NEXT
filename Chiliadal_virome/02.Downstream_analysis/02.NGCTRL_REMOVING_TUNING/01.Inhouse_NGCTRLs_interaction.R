setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

##########################################
# Exploring sharedness of vOTUs from 
# negative controls and own samples
##########################################

##############################
# Loading libraries
##############################
library(readr)
library(vegan)

library(ggplot2)
library(ggExtra)
library(ggrepel)

library(patchwork)

library(parallel)

library(ggforce)
library(ggsignif)
##############################
# Functions
##############################
my_theme <-   theme_bw() + 
  theme(axis.title = element_text(size=12,face="bold"),
        axis.text = element_text(size=10,face="bold"),
        legend.title = element_text(size=12,face="bold"),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-3,-10))

# RPKM_table function requires three dfs of the certain format: 
# counts: N of reads aligned to the sequence (contig)
# coverage: N of non-zero bases of the sequence length
# Rows are contigs and columns are samples.
# contigs_metadta: to extract the length of contigs
# RPKM_table returns a table of RPKM counts. 
# Be aware that obtained RPKM counts will be representing counts 
# normalized on the number of ALIGNED reads (to the database), 
# not the initial size of quality-trimmed library. 

RPKM_table <- function(counts, coverage, contigs_metadata){
  counts <- counts[order( row.names(counts) ),]
  
  coverage <- coverage[order( row.names(coverage) ),]
  
  if (identical(colnames(counts), colnames(coverage))) {
    
    RPM <- as.data.frame( t(t(counts) / (colSums(counts)/1000000)) )
    
    RPM$length <- contigs_metadata$POST_CHV_length[match(row.names(RPM), contigs_metadata$New_CID)]
    
    RPKM <- RPM / (RPM$length / 1000) 
    
    RPKM$length <- NULL
    
    RPKM[coverage<=0.75] <-0
  } else {
    print("Counts and coverage tables have different sets of samples")
  }
  
    return(RPKM)
}
##############################
# Input data
##############################
# for some reason, there was a warning when attempting to read VLP_metadata saved without quotes
VLP_metadata <- read.table('06.CLEAN_DATA/Intermediate/VLP_metadata_temporary.txt', sep='\t', header=T)
VLP_metadata <- VLP_metadata[order(VLP_metadata$Sequencing_ID_VLP),]
VLP_metadata$Timepoint_VLP <- factor(VLP_metadata$Timepoint_VLP, levels=c("P12", "P28", "B", "M1", "M3", "M6", "M12", "BLANK"), ordered = T)

VLP_counts <- read_delim('04.RAW_DATA/05.Abundance_tables/w_neg_der95/count_table.txt')
names(VLP_counts)[1] <- "vOTUs"
VLP_counts <- as.data.frame(VLP_counts)
row.names(VLP_counts) <- VLP_counts$vOTUs
VLP_counts$vOTUs <- NULL

VLP_contig_coverage <- read_delim('04.RAW_DATA/05.Abundance_tables/w_neg_der95/coverage_table.txt')
names(VLP_contig_coverage)[1] <- "vOTUs"
VLP_contig_coverage <- as.data.frame(VLP_contig_coverage)
row.names(VLP_contig_coverage) <- VLP_contig_coverage$vOTUs
VLP_contig_coverage$vOTUs <- NULL

VLP_contigs_metadata <- read_delim('04.RAW_DATA/04.Virome_discovery/Extended_table_of_origin')
VLP_contigs_metadata <- VLP_contigs_metadata[VLP_contigs_metadata$New_CID %in% row.names(VLP_counts), ]

RPKM_counts_VLP <- RPKM_table(VLP_counts, VLP_contig_coverage, VLP_contigs_metadata)

# Adding new columns to metadata (for the time being):
identical(VLP_metadata$Sequencing_ID_VLP, colnames(RPKM_counts_VLP))
VLP_metadata$Virus_richness <- colSums(RPKM_counts_VLP>0)

VLP_metadata$shared_with_NGCTRL1 <- sapply(colnames(RPKM_counts_VLP), function(x) sum( rowSums( RPKM_counts_VLP[, c('CHV199906F12', x)] > 0) == 2 ))
VLP_metadata$shared_with_NGCTRL2 <- sapply(colnames(RPKM_counts_VLP), function(x) sum( rowSums( RPKM_counts_VLP[, c('CHV200013F12', x)] > 0) == 2 ))

VLP_metadata$perc_shared_NGCTRL1 <- VLP_metadata$shared_with_NGCTRL1/VLP_metadata$Virus_richness*100
VLP_metadata$perc_shared_NGCTRL2 <- VLP_metadata$shared_with_NGCTRL2/VLP_metadata$Virus_richness*100
summary(VLP_metadata$perc_shared_NGCTRL1)
summary(VLP_metadata$perc_shared_NGCTRL2)

keep_RPKM_safe <- RPKM_counts_VLP

RPKM_counts_VLP["NEXT_V1036_N310_L4097_cov_36.631234_E0_P0_F0",] <- 0

VLP_metadata$abundance_NGCTRL1 <- colSums(RPKM_counts_VLP[RPKM_counts_VLP$CHV199906F12 >0, ]) / colSums(RPKM_counts_VLP) * 100
VLP_metadata$abundance_NGCTRL2 <- colSums(RPKM_counts_VLP[RPKM_counts_VLP$CHV200013F12 >0, ]) / colSums(RPKM_counts_VLP) * 100
summary(VLP_metadata$abundance_NGCTRL1)
sd(VLP_metadata$abundance_NGCTRL1)
summary(VLP_metadata$abundance_NGCTRL2)
sd(VLP_metadata$abundance_NGCTRL2)

View(VLP_metadata[VLP_metadata$BATCH_ISO=="IB47",])

VLP_metadata[VLP_metadata$BATCH_ISO=="IB47","Sequencing_ID_VLP"]

for (i in VLP_metadata[VLP_metadata$BATCH_ISO=="IB47","Sequencing_ID_VLP"]) {
  
  if (i != "CHV200013F12") {
    
    VLP_metadata[,paste0("shared_with_", i)] <- sapply(colnames(RPKM_counts_VLP), function(x) sum( rowSums( RPKM_counts_VLP[, c(i, x)] > 0) == 2 ))
    
  }
  
}

for (i in VLP_metadata[VLP_metadata$BATCH_ISO=="IB47","Sequencing_ID_VLP"]) {
  
  if (i != "CHV200013F12") {
    
    VLP_metadata[,paste0("abundance_", i)] <- colSums(RPKM_counts_VLP[RPKM_counts_VLP[,i] >0, ]) / colSums(RPKM_counts_VLP) * 100
    
  }
  
}

##############################
# ANALYSIS
##############################

#### NMDS

# preapring phenotypes
for_virplot <- VLP_metadata
row.names(for_virplot) <- VLP_metadata$Sequencing_ID_VLP
for_virplot <- for_virplot[ colnames(RPKM_counts_VLP), ]
for_virplot[!is.na(for_virplot$Type_VLP) & for_virplot$Type_VLP=="M",]$Timepoint_VLP <- "Mother"
for_virplot$Timepoint_VLP <- factor(for_virplot$Timepoint_VLP, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)
#CHOOSING TECHNICAL PHENOTYPES, SAMPLE TYPE & VIRAL ALPHA DIVERSITY
for_virplot <- for_virplot[,c("Type_VLP", "Reads_HQ", "DNA_CONC_VLP", "total.enrichmnet.score", "N_disc_viruses", "Timepoint_VLP")]

ord <- metaMDS(t(RPKM_counts_VLP), distance = "bray", k=2, parallel=8)
en = envfit(ord, for_virplot, permutations = 999, na.rm = TRUE)
en$factors #r2 0.36; p-value < 0.001
en$vectors

centroids <- as.data.frame(scores(en, "factors"))
centroids$Type <- c(gsub('Type_VLP|Timepoint_VLP', '', row.names(centroids)))

spp.scrs <- as.data.frame(scores(en, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Type <- 'Infant'
spp.scrs$Species <- c('Clean_reads','DNA concentration', 'virus enrichment', 'N discovered viruses')
spp.scrs <- spp.scrs[spp.scrs$Species!="virus enrichment",] # because it was not significant in envfit

### for plotting:
data.scores = as.data.frame(scores(ord, "sites"))
data.scores <- merge(data.scores, for_virplot, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL
data.scores$no_ngctrl_Type <- data.scores$Type_VLP
data.scores[!is.na(data.scores$no_ngctrl_Type) & data.scores$no_ngctrl_Type=="NC","no_ngctrl_Type"] <- NA
data.scores$Timepoint <- VLP_metadata$Timepoint_VLP[match(VLP_metadata$Sequencing_ID_VLP, row.names(data.scores))]
data.scores[!is.na(data.scores$Timepoint) & data.scores$Type_VLP=="M",]$Timepoint <- "Mother"
data.scores$Timepoint <- factor(data.scores$Timepoint, levels = c("M1", "M3", "M6", "M12", "Mother"), ordered=T)

pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/Viral_vOTUs_Bray_NMDS.pdf', width=14/2.54, height=10/2.54)
gg = ggplot(data = data.scores[!is.na(data.scores$Type_VLP),], aes(x = NMDS1, y = NMDS2, color=Type_VLP)) + 
  geom_point(size = 1.5, alpha=0.8) + 
  geom_point(data=centroids, aes(fill=Type),shape=23, size=4, color='black') + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Type_VLP, color=Type_VLP), linetype = 2)+
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
  scale_color_manual(values=c("#F8766D", "#00BFC4", "black"))

ggMarginal(gg, type="boxplot", groupFill=T)
dev.off()

gg = ggplot(data = data.scores[!is.na(data.scores$Timepoint),], aes(x = NMDS1, y = NMDS2, color=Timepoint)) + 
  stat_ellipse(geom = "polygon", alpha = 0.3, aes(group = Timepoint, color=Timepoint, fill=Timepoint), linetype = 2)+
  geom_point(size = 1.2, alpha=0.6) + 
  geom_point(data=centroids[-c(1,2),], aes(fill=Type),shape=23, size=4, color='black') +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        legend.position = "bottom") 

xAxisBoxPlot_v <- ggplot(data = data.scores[!is.na(data.scores$Timepoint),], aes(x = Timepoint, y = NMDS1)) +
  geom_boxplot(aes(fill=Timepoint), alpha=0.7, show.legend = FALSE) + 
  coord_flip() +
  theme_void() +
  theme(plot.title = element_text(hjust=0.0, size=7), plot.tag = element_text(face="bold", vjust=-7, size=7))

NMDS_by_time <- xAxisBoxPlot_v / gg + 
  plot_layout(heights = c(2,8)) & guides(color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         alpha="none")

pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/Viral_vOTUs_Bray_NMDS_by_time.pdf', width=14/2.54, height=10/2.54)
NMDS_by_time
dev.off()

###### EXPLORING OUTLIERS IN NMDS ######

###### EXPLORING NGCTRLs vOTUs sharedness #####

summary(lm(shared_with_NGCTRL1 ~ as.factor(BATCH_ISO), VLP_metadata))
tmp <- TukeyHSD(aov(shared_with_NGCTRL1 ~ as.factor(BATCH_ISO), VLP_metadata))
View(tmp[["as.factor(BATCH_ISO)"]])

summary(lm(shared_with_NGCTRL2 ~ as.factor(BATCH_ISO), VLP_metadata))
tmp2 <- TukeyHSD(aov(shared_with_NGCTRL2 ~ as.factor(BATCH_ISO), VLP_metadata))
View(tmp2[["as.factor(BATCH_ISO)"]])


#### ADD PAIR-WISE SIGNIFICANCE!

#pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/N_vOTUs_shared_NGCTRL1.pdf', width=27/2.54, height=10/2.54)
pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/Log10N_vOTUs_shared_NGCTRL1.pdf', width=27/2.54, height=10/2.54)
ggplot(VLP_metadata, aes(as.factor(BATCH_ISO), (shared_with_NGCTRL1 + 1) ) ) +
  geom_boxplot(aes(fill=BATCH_SEQ, color=BATCH_SEQ), alpha=0.5, outlier.shape = NA) +
  geom_sina(aes(color=BATCH_SEQ), size=0.7)+
  theme_bw() +
  #labs (y="N of shared contigs", x="Isolation batch") +
  labs (y="Log10 N of shared contigs", x="Isolation batch") + 
  scale_x_discrete(labels=gsub('IB', "", sort(unique(VLP_metadata$BATCH_ISO)))) +
  scale_y_log10() +
  # geom_signif(stat="identity",
  #             data=data.frame(x=c(1, 25), xend=c(1.125, 2.125),
  #                             y=c(5.8, 8.5), annotation=c("**", "NS")),
  #             aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) + 
  my_theme + 
  ggtitle("Number of NGCTRL1 viruses shared with batches")
dev.off()

#pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/N_vOTUs_shared_NGCTRL2.pdf', width=27/2.54, height=10/2.54)
pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/Log10N_vOTUs_shared_NGCTRL2.pdf', width=27/2.54, height=10/2.54)
ggplot(VLP_metadata, aes(as.factor(BATCH_ISO), (shared_with_NGCTRL2 + 1))) +
  geom_boxplot(aes(fill=BATCH_SEQ, color=BATCH_SEQ), alpha=0.5, outlier.shape = NA) +
  geom_sina(aes(color=BATCH_SEQ), size=0.7)+
  theme_bw() +
  #labs (y="N of shared contigs", x="Isolation batch") +
  labs (y="Log10 N of shared contigs", x="Isolation batch") + 
  scale_x_discrete(labels=gsub('IB', "", sort(unique(VLP_metadata$BATCH_ISO)))) +
  scale_y_log10() +
  my_theme + 
  ggtitle("Number of NGCTRL2 viruses shared with batches")
dev.off()


wilcox.test(VLP_metadata$shared_with_NGCTRL1 ~ VLP_metadata$BATCH_SEQ)
pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/Per_seq_batch_vOTUs_shared_NGCTRL1.pdf', width=14/2.54, height=10/2.54)
ggplot(VLP_metadata, aes(as.factor(BATCH_SEQ), (shared_with_NGCTRL1 + 1))) +
  geom_boxplot(aes(fill=BATCH_SEQ, color=BATCH_SEQ), alpha=0.5, outlier.shape = NA) +
  geom_sina(aes(color=BATCH_SEQ), size=0.7)+
  theme_bw() +
  #labs (y="N of shared contigs", x="Isolation batch") +
  labs (y="Log10 N of shared contigs", x="Sequencing batches") + 
  scale_y_log10() +
  my_theme + 
  ggtitle("Number of NGCTRL1 viruses shared with seq batches") + 
  geom_signif(comparisons=list(c("SB01", "SB02")), annotations="p-value=1.2e-04",
              y_position = 4.5, tip_length = 0)
dev.off()

wilcox.test(VLP_metadata$shared_with_NGCTRL2 ~ VLP_metadata$BATCH_SEQ)
pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/Per_seq_batch_vOTUs_shared_NGCTRL2.pdf', width=14/2.54, height=10/2.54)
ggplot(VLP_metadata, aes(as.factor(BATCH_SEQ), (shared_with_NGCTRL2 + 1))) +
  geom_boxplot(aes(fill=BATCH_SEQ, color=BATCH_SEQ), alpha=0.5, outlier.shape = NA) +
  geom_sina(aes(color=BATCH_SEQ), size=0.7)+
  theme_bw() +
  #labs (y="N of shared contigs", x="Isolation batch") +
  labs (y="Log10 N of shared contigs", x="Sequencing batches") + 
  scale_y_log10() +
  my_theme + 
  ggtitle("Number of NGCTRL2 viruses shared with seq batches") + 
  geom_signif(comparisons=list(c("SB01", "SB02")), annotations="p-value=4.7e-04",
              y_position = 4.4, tip_length = 0)
dev.off()


pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/Per_Type_perc_shared_NGCTRL1.pdf', width=14/2.54, height=10/2.54)
#pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/Per_Type_perc_shared_NGCTRL2.pdf', width=14/2.54, height=10/2.54)
ggplot(VLP_metadata, aes(Timepoint_VLP, perc_shared_NGCTRL1)) +
#ggplot(VLP_metadata, aes(Timepoint_VLP, perc_shared_NGCTRL2)) +
  geom_boxplot(aes(fill=Type_VLP, color=Type_VLP), alpha=0.5, outlier.shape = NA, show.legend = F) +
  geom_sina(aes(color=Type_VLP), size=0.7, show.legend = F) +
  facet_grid(~Type_VLP, scales = "free") +
  theme_bw() +
  labs (y="vOTU overlap with NGCTRl1 ratio", x="Type") +
  #labs (y="vOTU overlap with NGCTRl2 ratio", x="Type") +
  my_theme
dev.off()

#pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/Per_Type_abundance_vOTUs_shared_NGCTRL1.pdf', width=14/2.54, height=10/2.54)
pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/No_NEXT_V1036_N310_L4097_Per_Type_abundance_vOTUs_shared_NGCTRL2.pdf', width=14/2.54, height=10/2.54)
#ggplot(VLP_metadata, aes(Timepoint_VLP, abundance_NGCTRL1)) +
ggplot(VLP_metadata, aes(Timepoint_VLP, abundance_NGCTRL2)) +
  geom_boxplot(aes(fill=Type_VLP, color=Type_VLP), alpha=0.5, outlier.shape = NA, show.legend = F) +
  geom_sina(aes(color=Type_VLP), size=0.7, show.legend = F) +
  facet_grid(~Type_VLP, scales = "free") +
  theme_bw() +
  #labs (y="Abundance of vOTUs shared with NGCTRL1", x="Type") +
  labs (y="Abundance of vOTUs shared with NGCTRL2", x="Type") +
  my_theme
dev.off()


new_tmp <- reshape2::melt(VLP_metadata[,c(14,15,85, 90:114)])

new_tmp <- reshape2::melt(VLP_metadata[,c(14,15,89, 115:139)])

pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/Log10N_vOTUs_shared_Batch47.pdf', width=21/2.54, height=90/2.54)
ggplot(new_tmp, aes(as.factor(BATCH_ISO), (value + 1))) +
  geom_boxplot(aes(fill=BATCH_SEQ, color=BATCH_SEQ), alpha=0.5, outlier.shape = NA) +
  geom_sina(aes(color=BATCH_SEQ), size=0.7)+
  theme_bw() +
  facet_wrap(~variable, ncol = 1) +
  #labs (y="N of shared contigs", x="Isolation batch") +
  labs (y="Log10 N of shared contigs", x="Isolation batch") + 
  scale_x_discrete(labels=gsub('IB', "", sort(unique(VLP_metadata$BATCH_ISO)))) +
  scale_y_log10() +
  my_theme + 
  ggtitle("Number of ISO47 viruses shared with batches")
dev.off()


for (i in unique(new_tmp$variable)) {
  print(paste0(i, ' ', median(VLP_metadata[,i])))
}

pdf('05.PLOTS/02.NGCTRL_REMOVING_TUNING/Abundance_vOTUs_shared_Batch47.pdf', width=21/2.54, height=90/2.54)
ggplot(new_tmp, aes(as.factor(BATCH_ISO), value)) +
  geom_boxplot(aes(fill=BATCH_SEQ, color=BATCH_SEQ), alpha=0.5, outlier.shape = NA) +
  geom_sina(aes(color=BATCH_SEQ), size=0.7)+
  theme_bw() +
  facet_wrap(~variable, ncol = 1) +
  labs (y="Shared vOTU abundance", x="Isolation batch") + 
  scale_x_discrete(labels=gsub('IB', "", sort(unique(VLP_metadata$BATCH_ISO)))) +
  my_theme + 
  ggtitle("Abundance of ISO47 viruses shared with batches")
dev.off()

##############################
# OUTPUT
##############################
write.table(RPKM_counts_VLP, "06.CLEAN_DATA/Intermediate/RPKM_counts_VLP_1102.txt", sep='\t')


