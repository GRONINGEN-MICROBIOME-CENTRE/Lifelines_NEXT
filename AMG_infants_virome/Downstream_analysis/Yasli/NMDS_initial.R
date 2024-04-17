setwd('~/Desktop/Projects_2024/AMG_paper/')
##########################################
# NMDS
##########################################

##############################
# Loading libraries
##############################
library(vegan)
library(ggplot2)
library(patchwork)
library(MetBrewer)
library(ggforce)
##############################
# Functions
##############################

##############################
# Input data
##############################

metadata <- read.table('metadata_with_qc_v1.tsv', sep='\t', header=T)

# CHECK YOU METADATA, YOU HAVE NAs in Sample_name & DOUBLE LINE CONTAINING LN_7C08_VL_405
metadata <- metadata[-grep("LN_7C08_VL_405", metadata$Sample_ID),]
# ASSIGNING TYPE TO NEGATIVE CONTROL FROM PILOT
metadata[metadata$Sample_name=="LN_7C08_VL_405",]$Type <- 'Neg_ctrl'

RPKM_counts_VLP <- read.table("RPKM_counts_VLP_v1.txt", sep='\t', header=T)

# HOW MANY vOTU & SAMPLES?
dim(RPKM_counts_VLP) 
# [1] 373,301   2583

# HOW MANY vOTU DETECTED PER SAMPLE?
N_detected_viruses <- as.data.frame(colSums(RPKM_counts_VLP>0))
colnames(N_detected_viruses) <- 'N_detected_viruses'

# ADDING N_DETECTED vOTU PER SAMPLE TO METADATA
metadata$N_detected_viruses <- N_detected_viruses$N_detected_viruses[match(metadata$Sample_name, row.names(N_detected_viruses))]

# REMOVING vOTU THAT DID NOT PASS PRESENCE CUT-OFF & SAMPLES THAT DID NOT HAVE VIRUSES DETECTED
#RPKM_counts_VLP <- RPKM_counts_VLP[rowSums(RPKM_counts_VLP)>0,colSums(RPKM_counts_VLP)>0]

# HOW MANY vOTU & SAMPLES SURVIVED?
#dim(RPKM_counts_VLP)
# [1] 371,696   2509

for_virplot <- metadata[metadata$Sample_name %in% colnames(RPKM_counts_VLP),] 
row.names(for_virplot) <- for_virplot$Sample_name

# SORTING FOR ENV FIT
for_virplot <- for_virplot[ colnames(RPKM_counts_VLP), ]

data.scores <- read.table('data.scores.txt', sep='\t', header=T) # IT WILL BE OVERWRITTEN BY PILOT AND THEN READ AGAIN LATER

# COLOR PALETTE:
Peru <- met.brewer('Peru1')
##############################
# ANALYSIS
##############################
# 1ST QUESTION: IS EVERYTHING OK WITH THE PIPELINE? -> CHECKING WITH PILOT
PILOT <- RPKM_counts_VLP[,colnames(RPKM_counts_VLP) %in% metadata[metadata$cohort=='garmaeva',]$Sample_name] 
# REMOVING UNDETECTED IN PILOT VIRUSES
PILOT <- PILOT[rowSums(PILOT)>0,]
# HOW MANY vOTUs DETECTED IN PILOT SAMPLES?
dim(PILOT)
# [1] 209348    206
##### ->> 56.32% OF ALL DETECTED vOTUs ARE PRESENT IN PILOT DATASET 

# CALCULATING NMDS WITH 2 DIMENSIONS USING BRAY-CURTIS DISSIMILARITY AND 8 CPUs
ord <- metaMDS(t(PILOT), distance = "bray", k=2, parallel=8)

# PULLING NMDS1 AND NMDS2 FROM ORDINATION OBJECT
data.scores = as.data.frame(scores(ord, "sites"))
# MERGING NMDS DATA & METADATA 
data.scores <- merge(data.scores, for_virplot, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL
# FACTORIZING & ORDERING TIMEPOINT FOR PILOT DATASET ONLY:
data.scores[data.scores$Type=='Mother',]$Timepoint <- 'Mother'
data.scores[data.scores$Type=='Neg_ctrl',]$Timepoint <- 'NGCTRL'
data.scores$Timepoint <- factor(data.scores$Timepoint, levels = c('M1', 'M2', 'M3', 'M6', 'M12', 'Mother', 'NGCTRL'), ordered = T)

# HOW DOES NMDS LOOK LIKE FOR PILOT?

pv <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Timepoint)) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = Timepoint, fill=Timepoint), linetype = 2) +
  geom_point(size = 1.5, alpha=0.8) + 
  #xlim(-0.8, 0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=5),
        axis.title = element_text(size = 7), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         alpha="none")

xAxisBoxPlot_v <- ggplot(data = data.scores, aes(x = Timepoint, y = NMDS1)) +
  geom_boxplot(aes(fill=Timepoint), alpha=0.7, show.legend = FALSE) + 
  #ylim(-0.8,0.8)+
  coord_flip() +
  ggtitle("Shift in infant virome composition over time") +
  #scale_fill_manual(values = c("#3E0751", "#423A7F", "#3F678B", "#468E8B","#9FD55C", "#F9E855") ) +
  theme_void() +
  theme(plot.title = element_text(hjust=0.0, size=7), plot.tag = element_text(face="bold", vjust=-7, size=7))

Fig1B <- xAxisBoxPlot_v / pv + 
  plot_layout(heights = c(2,8)) & guides(color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         alpha="none")

pdf('PILOT_NMDS.pdf', width=12/2.54, height=10/2.54)
Fig1B
dev.off()

##### ->> NMDS FOR PILOT ALONE LOOKS FINE AND VERY SIMILAR TO THE ONE IN THE ORIGINAL PAPER
##### ->> ALL PIPELINE ALTERATIONS (SOFTWARE UPDATE; DARK MATTER CRITERION EXCLUSION ETC) DID NOT INFLUENCE THE OVERALL COMPOSITIONAL DISTRIBUTION

# HOW DOES NMDS LOOK LIKE FOR THE WHOLE DATASET?
data.scores <- read.table('data.scores.txt', sep='\t', header=T)
# MERGING NMDS DATA & METADATA 
data.scores <- merge(data.scores, for_virplot, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL
# COLLAPSING TIMEPOINT
data.scores[!is.na(data.scores$Type) & data.scores$Type=="Mother",]$Timepoint <- "Mother"

data.scores$Timepoint <- factor(data.scores$Timepoint, levels=c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16",
                                                                "M17", "M18", "M19", "M20", "M21", "M22", "M23", "M24", "M25", "M26", "M27", "M28", "M29", "M30", "M31", 
                                                                "M32", "M33", "M34", "M35", "M36", "M37", "M39", "Y2-5", "Mother", NA), ordered = T)

# NMDS COLORING BY COHORT
pdf('OVERALL_TRIMMED_NMDS.pdf', width=12/2.54, height=12/2.54)
ggplot(data = data.scores[!is.na(data.scores$Type),], aes(x = NMDS1, y = NMDS2, color=cohort)) + 
  geom_point(size = 1.5, alpha=0.6) + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = cohort, color=cohort), linetype = 2)+
  xlim(0.0005,0.003)+
  ylim(0.0045,0.0062) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        legend.position = "bottom") + 
  scale_color_manual(values = Peru) +
  guides(color=guide_legend(nrow=2,byrow=TRUE, title.position = 'top', title.hjust = 0.5))
dev.off()

# NMDS COLORING BY TIMEPOINT
pdf('OVERALL_TRIMMED_NMDS_TIMEPOINT.pdf', width=16/2.54, height=16/2.54)
ggplot(data = data.scores[!is.na(data.scores$Type),], aes(x = NMDS1, y = NMDS2, color=Timepoint)) + 
  geom_point(size = 1.5, alpha=0.6) + 
  #stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Timepoint, color=Timepoint), linetype = 2)+
  xlim(0.0005,0.003)+
  ylim(0.0045,0.0062) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        legend.position = "bottom") + 
  #scale_color_manual(values = Peru) +
  guides(color=guide_legend(nrow=6,byrow=TRUE, title.position = 'top', title.hjust = 0.5))
dev.off()

# AVERAGE NUMBER OF DETECTED vOTUs PER COHORT PER SAMPLE TYPE
# FACTORIZING COHORT:
pdf('RICHNESS_PER_TYPE_AND_COHORT.pdf', width=16/2.54, height=10/2.54)
data.scores$cohort <- factor(data.scores$cohort, levels = c('garmaeva', 'shah', 'walters', 'maqsood', 'liang'), ordered = T)
ggplot(data = data.scores[!is.na(data.scores$Type),], aes(x = Type, y = log10(N_detected_viruses), color=cohort)) +
  geom_sina(alpha=0.6) +
  geom_boxplot(alpha=0.3) + 
  facet_wrap(~cohort, scales = 'free') +
  labs(x='Sample Type', y='vOTU Richness') +
  theme_bw() +
  theme(axis.text=element_text(size=5),
        axis.title = element_text(size = 7, face = "bold", colour = "grey30"), 
        strip.text = element_text(size=5),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 5, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 7, colour = "grey30"),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5))
dev.off()

