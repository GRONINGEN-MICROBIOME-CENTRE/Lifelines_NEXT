setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the bacterial species-level stability
# of infant and maternal gut 
#############################################################

##############################
# Functions
##############################
source("03.SCRIPTS/NEXT_pilot_FUP_downstream/stability_functions.R")

##############################
# Loading libraries
##############################
library(gplots)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(MetBrewer)
##############################
# Input data
##############################

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels=c("P7", "B",
                                                                  "M1", "M2", "M3",
                                                                  "M6", "M12"), ordered=T)

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)

vOTUs_infants <- RPKM_counts_VLP[,VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID]
vOTUs_infants <- vOTUs_infants[rowSums(vOTUs_infants)!=0,]

vOTUs_mothers <- RPKM_counts_VLP[,VLP_metadata[VLP_metadata$Type=='Mother',]$Short_sample_ID]
vOTUs_mothers <- vOTUs_mothers[rowSums(vOTUs_mothers)!=0,]

Demuth <- met.brewer("Demuth")
Monet <- met.brewer('Monet')

##############################
# ANALYSIS
##############################
# vOTUs:
Infant_timepoints <- c('M1', 'M2', 'M3', 'M6', 'M12')
Mother_timepoints <- c("P7", "B", "M1", "M2", "M3")

### Infant vOTUs ###
# N bacterial species retained from M1 over different timepoints
virstability_M1 <- stability_initial(VLP_metadata[VLP_metadata$Type=='Infant',], vOTUs_infants, Infant_timepoints, 'Short_sample_ID', 'Richness')
virstability_M1_stat <- summary_stat_bootstrap(virstability_M1, Infant_timepoints)
virstability_M1_stat$Condition <- factor(virstability_M1_stat$Condition, levels = c('Retained', 'Richness', 'Not_retained'), ordered=T)
virstability_M1_stat <- virstability_M1_stat[row.names(virstability_M1_stat)!="Richness_M1",]

virstability_M1_abundance <- stability_initial(VLP_metadata[VLP_metadata$Type=='Infant',], vOTUs_infants, Infant_timepoints, 'Short_sample_ID', 'abundance')
virstability_M1_abundance_stat <- summary_stat_bootstrap(virstability_M1_abundance, Infant_timepoints)
virstability_M1_abundance_stat <- virstability_M1_abundance_stat[row.names(virstability_M1_abundance_stat)!="Total_space_M1",]
virstability_M1_abundance_stat$Condition <- factor(virstability_M1_abundance_stat$Condition, levels=c('Retained', 'Total_space', 'Not_retained', ordered=T))

### this plot contains only 11 infants that had M1 VLP sequenced
pdf('./04.PLOTS/vOTUs_retained_from_M1_over_time.pdf', width=7/2.54, height=10/2.54)
ggplot(virstability_M1_stat[virstability_M1_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_errorbar(data=virstability_M1_stat[virstability_M1_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "N detected vOTUs") + 
  labs(fill='Present vOTUs', color='', alpha='Transparency') + 
  theme_bw() + 
  theme(axis.title =element_text(size=9,face="bold"),
        legend.position = "bottom",
        legend.title = element_text(size=8,face="bold"),
        legend.text = element_text(size=7)) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from M1', 'Not present at M1')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  scale_alpha_manual(values=c(0.9,0.3)) +
  guides(color="none", 
         fill=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'), 
         alpha=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'))
dev.off()

pdf('./04.PLOTS/vOTUs_retained_from_M1_over_time_abundance.pdf', width=7/2.54, height=10/2.54)
ggplot(virstability_M1_abundance_stat[virstability_M1_abundance_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position="identity") +
  geom_errorbar(data=virstability_M1_abundance_stat[virstability_M1_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Present vOTUs', color='',alpha='Transparency') + 
  theme_bw() + 
  theme(axis.title =element_text(size=9,face="bold"),
        legend.position = "bottom",
        legend.title = element_text(size=8,face="bold"),
        legend.text = element_text(size=7)) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from M1', 'Not present at M1')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained', 'Total space')) +
  guides(color="none",
         fill=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'), 
         alpha=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'))
dev.off()


virpresence_from_M1 <- persistent_taxa(VLP_metadata[VLP_metadata$Type=='Infant',], vOTUs_infants, Infant_timepoints, "Short_sample_ID")
virpresence_from_M1_prevalent <- virpresence_from_M1[virpresence_from_M1$M1 > 0.5*length((VLP_metadata[VLP_metadata$Type=='Infant' & VLP_metadata$Timepoint=='M1',]$Short_sample_ID)),]
row.names(virpresence_from_M1_prevalent) <- gsub('.*_length','Length', row.names(virpresence_from_M1_prevalent) )
row.names(virpresence_from_M1_prevalent) <- gsub("\\...\\K\\d+", "", row.names(virpresence_from_M1_prevalent), perl = TRUE)
Okeeffe1_continuous <- met.brewer("OKeeffe1", n=9)
Okeeffe1_continuous <- Okeeffe1_continuous[9:1]

pdf('./04.PLOTS/vOTUs_persistent_over_time.pdf', width=20/2.54, height=20/2.54)
heatmap.2(as.matrix(virpresence_from_M1_prevalent), col=as.vector(Okeeffe1_continuous),
          #dendrogram ="none", 
          trace="none",
          #Rowv=FALSE, 
          Colv=FALSE,  
          margins=c(18,18)
)
dev.off()

## Personal virome fractions
# had to consider only 14 infants to exclude those that have less than 3 timepoints
p_vir_frac_infants <- personal_biome_fraction(VLP_metadata[VLP_metadata$Type=='Infant' & VLP_metadata$N_timepoints>2,], vOTUs_infants, 'Short_sample_ID')

PPB_vir_freq_infants <- as.data.frame(table(unlist(unname( p_vir_frac_infants[["PPB_bacteria"]] ))))

p_vir_frac_infants_individual <- p_vir_frac_infants[["personal_biome_fractions"]]
p_vir_frac_infants_individual$Infant <- row.names(p_vir_frac_infants_individual)
p_vir_frac_infants_individual <- melt(p_vir_frac_infants_individual[,c('PPB', 'TDB', 'Singletons', 'Infant')], id.vars = 'Infant')
p_vir_frac_infants_individual_df <- p_vir_frac_infants_individual[p_vir_frac_infants_individual$variable=='PPB',]
p_vir_frac_infants_individual_df <- p_vir_frac_infants_individual_df[order(p_vir_frac_infants_individual_df$value,decreasing = T),]
p_vir_frac_infants_individual$Infant <- factor(p_vir_frac_infants_individual$Infant, levels=p_vir_frac_infants_individual_df$Infant, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction.pdf', width=10/2.54, height=10/2.54)
ggplot(p_vir_frac_infants_individual, aes(Infant, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity") + 
  xlab(label = "Infant") + 
  ylab(label = "% of infant personal virome") + 
  labs(fill='Fraction') +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c("PPV", "TDV", "Singletons")) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.position = "bottom")
dev.off()

p_vir_frac_infants_per_sample <- p_vir_frac_infants[["personal_biome_per_sample"]]
p_vir_frac_infants_per_sample$Sample <- row.names(p_vir_frac_infants_per_sample)

p_vir_frac_infants_per_sample_ab <- melt(p_vir_frac_infants_per_sample[,c(colnames(p_vir_frac_infants_per_sample)[grep('ab_', colnames(p_vir_frac_infants_per_sample))], "Sample")], id='Sample')
p_vir_frac_infants_per_sample_ab$Infant <- VLP_metadata$Individual_ID[match(p_vir_frac_infants_per_sample_ab$Sample, VLP_metadata$Short_sample_ID)] 
p_vir_frac_infants_per_sample_ab$Timepoint <- VLP_metadata$Timepoint[match(p_vir_frac_infants_per_sample_ab$Sample, VLP_metadata$Short_sample_ID)]
p_vir_frac_infants_per_sample_ab$Infant <- factor(p_vir_frac_infants_per_sample_ab$Infant, levels=p_vir_frac_infants_individual_df$Infant, ordered = T)
p_vir_frac_infants_per_sample_ab$Timepoint <- factor(p_vir_frac_infants_per_sample_ab$Timepoint, levels=Infant_timepoints, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction_abundance.pdf', width=10/2.54, height=7/2.54)
ggplot(p_vir_frac_infants_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Infant samples") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Fraction') +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPV', 'TDV', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()


plot_vir_frac_infants_per_sample_N <- melt(p_vir_frac_infants_per_sample[,c(colnames(p_vir_frac_infants_per_sample)[grep('N_', colnames(p_vir_frac_infants_per_sample))], "Sample")], id='Sample')
plot_vir_frac_infants_per_sample_N$Infant <- VLP_metadata$Individual_ID[match(plot_vir_frac_infants_per_sample_N$Sample, VLP_metadata$Short_sample_ID)] 
plot_vir_frac_infants_per_sample_N$Timepoint <- VLP_metadata$Timepoint[match(plot_vir_frac_infants_per_sample_N$Sample, VLP_metadata$Short_sample_ID)]
plot_vir_frac_infants_per_sample_N$Infant <- factor(plot_vir_frac_infants_per_sample_N$Infant, levels=p_vir_frac_infants_individual_df$Infant, ordered = T)
plot_vir_frac_infants_per_sample_N$Timepoint <- factor(plot_vir_frac_infants_per_sample_N$Timepoint, levels=Infant_timepoints, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction_N.pdf', width=10/2.54, height=7/2.54)
ggplot(plot_vir_frac_infants_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Infant samples") + 
  ylab(label = "% of sample richness") + 
  labs(fill='Fraction') +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPV', 'TDV', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()

### Mother vOTUs ###
virstability_P7 <- stability_initial(VLP_metadata[VLP_metadata$Type=='Mother',], vOTUs_mothers, Mother_timepoints, 'Short_sample_ID', 'Richness')
virstability_P7_stat <- summary_stat_bootstrap(virstability_P7, Mother_timepoints)
virstability_P7_stat$Condition <- factor(virstability_P7_stat$Condition, levels = c('Retained', 'Richness', 'Not_retained'), ordered=T)
virstability_P7_stat <- virstability_P7_stat[row.names(virstability_P7_stat)!="Richness_P7",]



virstability_P7_abundance <- stability_initial(VLP_metadata[VLP_metadata$Type=='Mother',], vOTUs_mothers, Mother_timepoints, 'Short_sample_ID', 'abundance')
virstability_P7_abundance_stat <- summary_stat_bootstrap(virstability_P7_abundance, Mother_timepoints)
virstability_P7_abundance_stat$Condition <- factor(virstability_P7_abundance_stat$Condition, levels = c('Retained', 'Total_space', 'Not_retained'), ordered=T)
virstability_P7_abundance_stat <- virstability_P7_abundance_stat[row.names(virstability_P7_abundance_stat)!="Total_space_P7",]

pdf('./04.PLOTS/vOTUs_retained_from_P7_over_time.pdf', width=7/2.54, height=10/2.54)
ggplot(virstability_P7_stat[virstability_P7_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_errorbar(data=virstability_P7_stat[virstability_P7_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap,
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "N detected vOTUs") + 
  labs(fill='Present vOTUs', color='') + 
  theme_bw() + 
  theme(axis.title =element_text(size=9,face="bold"),
        legend.position = "bottom",
        legend.title = element_text(size=8,face="bold"),
        legend.text = element_text(size=7)) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from P7', 'Not present at P7')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  scale_alpha_manual(values=c(0.9,0.3)) +
  guides(color="none", 
         fill=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'), 
         alpha=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'))
dev.off()


pdf('./04.PLOTS/vOTUs_retained_from_P7_over_time_abundance.pdf', width=7/2.54, height=10/2.54)
ggplot(virstability_P7_abundance_stat[virstability_P7_abundance_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_errorbar(data=virstability_P7_abundance_stat[virstability_P7_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap,
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Present vOTUs', color='') + 
  theme_bw() + 
  theme(axis.title =element_text(size=9,face="bold"),
        legend.position = "bottom",
        legend.title = element_text(size=8,face="bold"),
        legend.text = element_text(size=7)) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from P7', 'Not present at P7')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained', 'Total space')) +
  guides(color="none", 
         fill=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'), 
         alpha=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'))
dev.off()

virpresence_from_P7 <- persistent_taxa(VLP_metadata[VLP_metadata$Type=='Mother',], vOTUs_mothers, Mother_timepoints, "Short_sample_ID")
virpresence_from_P7_prevalent <- virpresence_from_P7[virpresence_from_P7$P7 > 0.5*length((VLP_metadata[VLP_metadata$Type=='Mother' & VLP_metadata$Timepoint=='P7',]$Short_sample_ID)),]
row.names(virpresence_from_P7_prevalent) <- gsub('.*_length','Length', row.names(virpresence_from_P7_prevalent) )
row.names(virpresence_from_P7_prevalent) <- gsub("\\...\\K\\d+", "", row.names(virpresence_from_P7_prevalent), perl = TRUE)

Okeeffe1_continuous <- met.brewer("OKeeffe1", n=19)
Okeeffe1_continuous <- Okeeffe1_continuous[19:1]

pdf('./04.PLOTS/vOTUs_persistent_over_time_mothers.pdf', width=20/2.54, height=20/2.54)
heatmap.2(as.matrix(virpresence_from_P7_prevalent), col=as.vector(Okeeffe1_continuous),
          #dendrogram ="none", 
          trace="none",
          #Rowv=FALSE, 
          Colv=FALSE,  
          margins=c(18,18)
)
dev.off()



p_vir_frac_mothers <- personal_biome_fraction(VLP_metadata[VLP_metadata$Type=='Mother' & VLP_metadata$N_timepoints>2,], vOTUs_mothers, 'Short_sample_ID')

PPB_vir_freq_mothers <- as.data.frame(table(unlist(unname( p_vir_frac_mothers[["PPB_bacteria"]] ))))

p_vir_frac_mothers_individual <- p_vir_frac_mothers[["personal_biome_fractions"]]
p_vir_frac_mothers_individual$Mother <- row.names(p_vir_frac_mothers_individual)
p_vir_frac_mothers_individual <- melt(p_vir_frac_mothers_individual[,c('PPB', 'TDB', 'Singletons', 'Mother')], id.vars = 'Mother')
p_vir_frac_mothers_individual_df <- p_vir_frac_mothers_individual[p_vir_frac_mothers_individual$variable=='PPB',]
p_vir_frac_mothers_individual_df <- p_vir_frac_mothers_individual_df[order(p_vir_frac_mothers_individual_df$value,decreasing = T),]
p_vir_frac_mothers_individual$Mother <- factor(p_vir_frac_mothers_individual$Mother, levels=p_vir_frac_mothers_individual_df$Mother, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction_mothers.pdf', width=10/2.54, height=10/2.54)
ggplot(p_vir_frac_mothers_individual, aes(Mother, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity") + 
  xlab(label = "Mother") + 
  ylab(label = "% of maternal personal virome") + 
  labs(fill='Fraction') +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c("PPV", "TDV", "Singletons")) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.position = "bottom")
dev.off()

p_vir_frac_mothers_per_sample <- p_vir_frac_mothers[["personal_biome_per_sample"]]
p_vir_frac_mothers_per_sample$Sample <- row.names(p_vir_frac_mothers_per_sample)

p_vir_frac_mothers_per_sample_ab <- melt(p_vir_frac_mothers_per_sample[,c(colnames(p_vir_frac_mothers_per_sample)[grep('ab_', colnames(p_vir_frac_mothers_per_sample))], "Sample")], id='Sample')
p_vir_frac_mothers_per_sample_ab$Mother <- VLP_metadata$Individual_ID[match(p_vir_frac_mothers_per_sample_ab$Sample, VLP_metadata$Short_sample_ID)] 
p_vir_frac_mothers_per_sample_ab$Timepoint <- VLP_metadata$Timepoint[match(p_vir_frac_mothers_per_sample_ab$Sample, VLP_metadata$Short_sample_ID)]
p_vir_frac_mothers_per_sample_ab$Mother <- factor(p_vir_frac_mothers_per_sample_ab$Mother, levels=p_vir_frac_mothers_individual_df$Mother, ordered = T)
p_vir_frac_mothers_per_sample_ab$Timepoint <- factor(p_vir_frac_mothers_per_sample_ab$Timepoint, levels=Mother_timepoints, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction_abundance_mothers.pdf', width=25/2.54, height=7/2.54)
ggplot(p_vir_frac_mothers_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Maternal samples") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Fraction') +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPV', 'TDV', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()


plot_vir_frac_mothers_per_sample_N <- melt(p_vir_frac_mothers_per_sample[,c(colnames(p_vir_frac_mothers_per_sample)[grep('N_', colnames(p_vir_frac_mothers_per_sample))], "Sample")], id='Sample')
plot_vir_frac_mothers_per_sample_N$Mother <- VLP_metadata$Individual_ID[match(plot_vir_frac_mothers_per_sample_N$Sample, VLP_metadata$Short_sample_ID)] 
plot_vir_frac_mothers_per_sample_N$Timepoint <- VLP_metadata$Timepoint[match(plot_vir_frac_mothers_per_sample_N$Sample, VLP_metadata$Short_sample_ID)]
plot_vir_frac_mothers_per_sample_N$Mother <- factor(plot_vir_frac_mothers_per_sample_N$Mother, levels=p_vir_frac_mothers_individual_df$Mother, ordered = T)
plot_vir_frac_mothers_per_sample_N$Timepoint <- factor(plot_vir_frac_mothers_per_sample_N$Timepoint, levels=Mother_timepoints, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction_N_mothers.pdf', width=25/2.54, height=7/2.54)
ggplot(plot_vir_frac_mothers_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Maternal samples") + 
  ylab(label = "% of sample richness") + 
  labs(fill='Fraction') +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPV', 'TDV', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()

##### OUTPUT #####
write.table(virstability_M1, '03a.RESULTS/N_retained_vOTUs_from_M1_infants.txt', sep='\t', quote=F)
write.table(virstability_M1_abundance, '03a.RESULTS/Abundance_retained_vOTUs_from_M1_infants.txt', sep='\t', quote=F)
write.table(virpresence_from_M1, '03a.RESULTS/vOTUs_retained_from_M1_over_time_infants.txt', sep='\t', quote=F)
write.table(PPB_vir_freq_infants, '03a.RESULTS/vOTUs_members_of_PPV_infants_freq.txt', sep='\t', quote=F)
write.table(p_vir_frac_infants_per_sample, '03a.RESULTS/Fractions_of_personal_infant_virome_N_and_abundance.txt', sep='\t', quote=F)

write.table(virstability_P7, '03a.RESULTS/N_retained_vOTUs_from_P7_mothers.txt', sep='\t', quote=F)
write.table(virstability_P7_abundance, '03a.RESULTS/Abundance_retained_vOTUs_from_P7_mothers.txt', sep='\t', quote=F)
write.table(virpresence_from_P7, '03a.RESULTS/vOTUs_retained_from_P7_over_time_mothers.txt', sep='\t', quote=F)
write.table(PPB_vir_freq_mothers, '03a.RESULTS/vOTUs_members_of_PPV_mothers_freq.txt', sep='\t', quote=F)
write.table(p_vir_frac_mothers_per_sample, '03a.RESULTS/Fractions_of_personal_maternal_virome_N_and_abundance.txt', sep='\t', quote=F)


