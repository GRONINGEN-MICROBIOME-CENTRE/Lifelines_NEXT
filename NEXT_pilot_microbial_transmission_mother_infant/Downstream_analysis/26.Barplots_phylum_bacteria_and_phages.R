setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the bacterial phylum-level composition
# of infant and maternal gut vs vOTUs aggregated at host
# phylum level
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################
library(dplyr)
library(reshape2)
library(ggplot2)
library(MetBrewer)
library(patchwork)
##############################
# Input data
##############################
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)

bac_phylum <- read.table('02.CLEAN_DATA/Microbiome_phyla_unfiltred.txt', sep='\t', header=T)
bac_phylum["gap",] <- 100-colSums(bac_phylum) # due to 'Unclassified' in metaphlan4, not all colSums will be equal to 100
bac_phylum_filt <- bac_phylum[rowSums(bac_phylum!=0)>0.10*ncol(bac_phylum),]
bac_phylum_filt["Other",] <- 100-colSums(bac_phylum_filt) # the sum abundance of rarw phyla are here
bac_phylum_filt <- bac_phylum_filt[row.names(bac_phylum_filt)!='gap',]
bac_phylum_filt$Phylum <- row.names(bac_phylum_filt)

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)

host_assignment <- read.table('02.CLEAN_DATA/Host_prediction_to_genus_m90_refined_taxonomy_no_generalists.txt', sep='\t', header=T)

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)
RPKM_counts_VLP$Virus <- row.names(RPKM_counts_VLP)

host_assignment_combined <- merge(host_assignment[,c("Virus", "Phylum")], RPKM_counts_VLP, by='Virus', all.y = T)
host_assignment_combined[is.na(host_assignment_combined$Phylum),"Phylum"] <- 'Unassigned'
host_assignment_combined$Virus <- NULL

RPKM_counts_BacPhyl <- aggregate(.~Phylum, host_assignment_combined, sum)
row.names(RPKM_counts_BacPhyl) <- RPKM_counts_BacPhyl$Phylum
RPKM_counts_BacPhyl$Phylum <- NULL
RPKM_counts_BacPhyl_RA <- as.data.frame( t(t( RPKM_counts_BacPhyl )/colSums(RPKM_counts_BacPhyl))*100  )

#RPKM_counts_BacPhyl_RA_filt <- RPKM_counts_BacPhyl_RA[rowSums(RPKM_counts_BacPhyl_RA!=0)>0.10*ncol(RPKM_counts_BacPhyl_RA),]
#RPKM_counts_BacPhyl_RA_filt["Other",] <- 100 - round(colSums(RPKM_counts_BacPhyl_RA_filt), 2)
#RPKM_counts_BacPhyl_RA_filt$Phylum <- row.names(RPKM_counts_BacPhyl_RA_filt)

RPKM_counts_MGS <- read.table('02.CLEAN_DATA/RPKM_counts_MGS.txt', sep='\t', header=T)
RPKM_counts_MGS$Virus <- row.names(RPKM_counts_MGS)

host_assignment_combined_MGS <- merge(host_assignment[,c("Virus", "Phylum")], RPKM_counts_MGS, by='Virus', all.y = T)
host_assignment_combined_MGS[is.na(host_assignment_combined_MGS$Phylum),"Phylum"] <- 'Unassigned'
host_assignment_combined_MGS$Virus <- NULL

RPKM_counts_MGS_BacPhyl <- aggregate(.~Phylum, host_assignment_combined_MGS, sum)
row.names(RPKM_counts_MGS_BacPhyl) <- RPKM_counts_MGS_BacPhyl$Phylum
RPKM_counts_MGS_BacPhyl$Phylum <- NULL
RPKM_counts_MGS_BacPhyl_RA <- as.data.frame( t(t( RPKM_counts_MGS_BacPhyl )/colSums(RPKM_counts_MGS_BacPhyl))*100  )
# getting only those that are absent in VLP dataset
RPKM_counts_MGS_BacPhyl_RA <- RPKM_counts_MGS_BacPhyl_RA[,MGS_metadata[ !(MGS_metadata$Universal_fecal_ID %in% VLP_metadata$Universal_fecal_ID),  ]$Short_sample_ID]


RPKM_counts_BacPhyl_all <- merge(RPKM_counts_BacPhyl_RA, RPKM_counts_MGS_BacPhyl_RA, by='row.names', all=T)
row.names(RPKM_counts_BacPhyl_all) <- RPKM_counts_BacPhyl_all$Row.names
RPKM_counts_BacPhyl_all$Row.names <- NULL
RPKM_counts_BacPhyl_all[is.na(RPKM_counts_BacPhyl_all)] <- 0

RPKM_counts_BacPhyl_all_filt <- RPKM_counts_BacPhyl_all[rowSums(RPKM_counts_BacPhyl_all!=0)>0.10*ncol(RPKM_counts_BacPhyl_all),]
RPKM_counts_BacPhyl_all_filt["Other",] <- 100 - round(colSums(RPKM_counts_BacPhyl_all_filt), 2)
RPKM_counts_BacPhyl_all_filt[c(phyla_palette[!(phyla_palette$Phylum %in% row.names(RPKM_counts_BacPhyl_all_filt)),]$Phylum),] <- 0
RPKM_counts_BacPhyl_all_filt$Phylum <- row.names(RPKM_counts_BacPhyl_all_filt)
##############################
# ANALYSIS
##############################
Renoir <- met.brewer("Renoir")

Kandinsky <- met.brewer("Kandinsky")


phyla_palette <- data.frame(unique(c(levels(bac_phylum_plot$Phylum), levels(vOTUs_bac_phylum_plot$Phylum))))
colnames(phyla_palette) <- 'Phylum'
phyla_palette$Phylum <- phyla_palette[order(phyla_palette$Phylum),]
phyla_palette$Phylum <- phyla_palette[c(1:8,10:12,14,9,13),]
phyla_palette$Color <- c(Kandinsky[1:2], Renoir[1:10], Kandinsky[3:4])


###### BACTERIAL PHYLA
bac_phylum_plot <- melt(bac_phylum_filt, id.vars="Phylum")
bac_phylum_plot$Kingdom <- gsub('.*k__','', sapply(strsplit(bac_phylum_plot$Phylum, '\\|'), "[", 1))
bac_phylum_plot$Phylum <- gsub('.*p__','',bac_phylum_plot$Phylum)
bac_phylum_plot$Phylum <- factor(bac_phylum_plot$Phylum, levels=unique(bac_phylum_plot$Phylum), ordered=T)
bac_phylum_plot$Type <- MGS_metadata$Type[ match( bac_phylum_plot$variable, MGS_metadata$Short_sample_ID_bact)  ]
bac_phylum_plot$Timepoint <- MGS_metadata$Timepoint[ match( bac_phylum_plot$variable, MGS_metadata$Short_sample_ID_bact)  ]
bac_phylum_plot$Timepoint <- factor(bac_phylum_plot$Timepoint, levels = c("P3", "P7", "B",
                                                                          "M1", "M2", "M3",
                                                                       "M6", "M9", "M12"), ordered = T)
bac_phylum_plot$Individual_ID <- MGS_metadata$Individual_ID[ match( bac_phylum_plot$variable, MGS_metadata$Short_sample_ID_bact)  ]

# ordering the future graphs in by the abundance of Bacteroidetes phylum at the first sampled timepoint
ordering_df <- bac_phylum_plot[bac_phylum_plot$Phylum=='Actinobacteria' & bac_phylum_plot$Type=='Infant' & bac_phylum_plot$Timepoint=='M1',c('value', 'Individual_ID')]
ordering_df <- ordering_df[order(ordering_df$value,decreasing = T),]
ordering_df_m <- bac_phylum_plot[bac_phylum_plot$Phylum=='Actinobacteria' & bac_phylum_plot$Type=='Mother' & bac_phylum_plot$Timepoint=='P3',c('value', 'Individual_ID')]
ordering_df_m <- ordering_df_m[order(ordering_df_m$value, decreasing = T),]

bac_phylum_plot$Individual_ID <- factor(bac_phylum_plot$Individual_ID, levels = c(ordering_df$Individual_ID, 
                                                                                setdiff(unique(bac_phylum_plot[bac_phylum_plot$Type=='Infant',]$Individual_ID), ordering_df$Individual_ID),
                                                                                                    ordering_df_m$Individual_ID,
                                                                                                    setdiff(unique(bac_phylum_plot[bac_phylum_plot$Type=='Mother',]$Individual_ID), ordering_df_m$Individual_ID)), ordered=T )

p1 <- ggplot(bac_phylum_plot[bac_phylum_plot$Type=='Mother',], aes(Individual_ID, value, fill=Phylum)) + 
  geom_bar(position = "fill", stat = "identity",width=1) + 
  facet_wrap(~Timepoint, scales = 'free_y', nrow = 6) + 
  xlab(label = "Mothers microbiome") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Phyla') +
  theme_bw() + 
  scale_fill_manual(values=phyla_palette[phyla_palette$Phylum %in% unique(bac_phylum_plot$Phylum),]$Color) +
  theme(axis.text.y=element_text(size=16), 
        axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12), 
        legend.position = "none")

p2 <- ggplot(bac_phylum_plot[bac_phylum_plot$Type=='Infant',], aes(Individual_ID, value, fill=Phylum)) + 
  geom_bar(position = "fill", stat = "identity",width=1) + 
  facet_wrap(~Timepoint, scales = 'free_y', nrow = 6) + 
  xlab(label = "Infants microbiome") + 
  ylab(label = "") + 
  labs(fill='Phyla') +
  theme_bw() + 
  scale_fill_manual(values=phyla_palette[phyla_palette$Phylum %in% unique(bac_phylum_plot$Phylum),]$Color) +
  theme(axis.text.y=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x =element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.position = "none")

###### vOTUs at BACTERIAL PHYLA level
vOTUs_bac_phylum_plot <- melt(RPKM_counts_BacPhyl_all_filt, id.vars="Phylum")

vOTUs_bac_phylum_plot$Phylum <- factor(vOTUs_bac_phylum_plot$Phylum, levels=c(phyla_palette$Phylum), ordered=T)
vOTUs_bac_phylum_plot$Source <- substr(vOTUs_bac_phylum_plot$variable, 8,8)
vOTUs_bac_phylum_plot$Type <- NA
vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='V',]$Type <- VLP_metadata$Type[ match( vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='V',]$variable, VLP_metadata$Short_sample_ID)  ]
vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='M',]$Type <- MGS_metadata$Type[ match( vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='M',]$variable, MGS_metadata$Short_sample_ID)  ]

vOTUs_bac_phylum_plot$Timepoint <- NA
vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='V',]$Timepoint <- VLP_metadata$Timepoint[ match( vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='V',]$variable, VLP_metadata$Short_sample_ID)  ]
vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='M',]$Timepoint <- MGS_metadata$Timepoint[ match( vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='M',]$variable, MGS_metadata$Short_sample_ID)  ]
vOTUs_bac_phylum_plot$Timepoint <- factor(vOTUs_bac_phylum_plot$Timepoint, levels = c("P3", "P7", "B",
                                                                                      "M1", "M2", "M3",
                                                                                      "M6", "M9", "M12"), ordered = T)
vOTUs_bac_phylum_plot$Individual_ID <- NA
vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='V',]$Individual_ID <- VLP_metadata$Individual_ID[ match( vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='V',]$variable, VLP_metadata$Short_sample_ID)  ]
vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='M',]$Individual_ID <- MGS_metadata$Individual_ID[ match( vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='M',]$variable, MGS_metadata$Short_sample_ID)  ]

# ordering the future graphs in by the abundance of Bacteroidetes phylum at the first sampled timepoint
vOTUs_ordering_df <- vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Phylum=='Actinobacteria' & vOTUs_bac_phylum_plot$Type=='Infant' & vOTUs_bac_phylum_plot$Timepoint=='M1',c('value', 'Individual_ID')]
vOTUs_ordering_df <- vOTUs_ordering_df[order(vOTUs_ordering_df$value,decreasing = T),]
vOTUs_ordering_df_m <- vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Phylum=='Actinobacteria' & vOTUs_bac_phylum_plot$Type=='Mother' & vOTUs_bac_phylum_plot$Timepoint=='P3',c('value', 'Individual_ID')]
vOTUs_ordering_df_m <- vOTUs_ordering_df_m[order(vOTUs_ordering_df_m$value, decreasing = T),]

vOTUs_bac_phylum_plot$Individual_ID <- factor(vOTUs_bac_phylum_plot$Individual_ID, levels = c(vOTUs_ordering_df$Individual_ID, 
                                                                                  setdiff(unique(vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Type=='Infant',]$Individual_ID), vOTUs_ordering_df$Individual_ID),
                                                                                  vOTUs_ordering_df_m$Individual_ID,
                                                                                  setdiff(unique(vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Type=='Mother',]$Individual_ID), vOTUs_ordering_df_m$Individual_ID)), ordered=T )
vOTUs_bac_phylum_plot$Sequencing <- vOTUs_bac_phylum_plot$Source
vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='M',]$Sequencing <- 'MGS'
vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Source=='V',]$Sequencing <- 'VLP'

p3 <- ggplot(vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Type=='Mother',], aes(Individual_ID, value, fill=Phylum, alpha=Sequencing)) + 
  geom_bar(position = "fill", stat = "identity",width=1) + 
  facet_wrap(~Timepoint, scales = 'free_y', nrow = 6) + 
  xlab(label = "Mothers virome") + 
  ylab(label = "") + 
  labs(fill='Phyla', alpha='Sequencing') +
  theme_bw() + 
  scale_fill_manual(values=phyla_palette[phyla_palette$Phylum %in% unique(RPKM_counts_BacPhyl_all_filt$Phylum),]$Color) +
  scale_alpha_manual(values = c(0.7, 1)) +
  theme(axis.text.y=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x =element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.position = "none")

p4 <- ggplot(vOTUs_bac_phylum_plot[vOTUs_bac_phylum_plot$Type=='Infant',], aes(Individual_ID, value, fill=Phylum, alpha=Sequencing)) + 
  geom_bar(position = "fill", stat = "identity",width=1) + 
  facet_wrap(~Timepoint, scales = 'free_y', nrow = 6) + 
  xlab(label = "Infants virome") + 
  ylab(label = "") + 
  labs(fill='Host phyla',alpha='Sequencing') +
  theme_bw() + 
  scale_fill_manual(values=phyla_palette[phyla_palette$Phylum %in% unique(RPKM_counts_BacPhyl_all_filt$Phylum),]$Color) +
  scale_alpha_manual(values = c(0.7, 1)) +
  theme(axis.text.y=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x =element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12))

combined_plot <- p1 + p2 + p3 + p4

combined_plot <- combined_plot +
  plot_layout(ncol = 4, guides = "collect") + 
  plot_annotation(title = "")

# Print the combined plot
pdf('./04.PLOTS/Phyla_barplots_bacteria_vOTUs.pdf', width=30/2.54, height=22/2.54)
combined_plot
dev.off()

#!!!! NOT ALL Phyla are reflected!


###### OUTPUT #####

