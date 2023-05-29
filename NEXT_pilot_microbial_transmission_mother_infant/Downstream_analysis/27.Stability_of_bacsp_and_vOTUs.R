setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the bacterial species-level stability
# of infant and maternal gut vs vOTUs aggregated at host
# species level stability
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
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels=c("P3", "P7", "B",
                                                                  "M1", "M2", "M3",
                                                                  "M6", "M9", "M12"), ordered=T)
#MGS_metadata <- MGS_metadata[MGS_metadata$Type=='Infant',]
species <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)

species_infants <- species[,MGS_metadata[MGS_metadata$Type=='Infant',]$Short_sample_ID_bact]
species_infants <- species_infants[rowSums(species_infants)!=0,]

species_mothers <- species[,MGS_metadata[MGS_metadata$Type=='Mother',]$Short_sample_ID_bact]
species_mothers <- species_mothers[rowSums(species_mothers)!=0,]

Demuth <- met.brewer("Demuth")
Monet <- met.brewer('Monet')
##############################
# ANALYSIS
##############################
# bacteria:
Infant_timepoints <- c('M1', 'M2', 'M3', 'M6', 'M9', 'M12')
Mother_timepoints <- c("P3", "P7", "B", "M1", "M2", "M3")
### Infant Bacteria ###



# N bacterial species retained from M1 over different timepoints
bacstability_M1 <- stability_initial(MGS_metadata[MGS_metadata$Type=='Infant',], species_infants, Infant_timepoints, 'Short_sample_ID_bact', 'Richness')
bacstability_M1_stat <- summary_stat(bacstability_M1, Infant_timepoints)

bacstability_M1_abundance <- stability_initial(MGS_metadata[MGS_metadata$Type=='Infant',], species_infants, Infant_timepoints, 'Short_sample_ID_bact', 'abundance')
bacstability_M1_abundance_stat <- summary_stat(bacstability_M1_abundance, Infant_timepoints)


pdf('./04.PLOTS/BacSp_retained_from_M1_over_time.pdf', width=10/2.54, height=7/2.54)
ggplot(bacstability_M1_stat[bacstability_M1_stat$Condition!='Richness',], aes(Timepoint, mean_value, fill=Condition)) + 
  geom_bar(stat = "identity", position="stack") +
  geom_errorbar(data=bacstability_M1_stat[bacstability_M1_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=mean_value - sd_value , 
                    ymax=mean_value + sd_value , 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "N detected species") + 
  labs(fill='Present species', color='') + 
  theme_bw() + 
  theme(axis.title =element_text(size=10,face="bold")) + 
  scale_fill_manual(values = c(Demuth[6:5], "#FFFFFF"), labels=c('Not present at M1', 'Retained from M1', '')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  guides(color="none")
dev.off()

pdf('./04.PLOTS/BacSp_retained_from_M1_over_time_abundance.pdf', width=10/2.54, height=7/2.54)
ggplot(bacstability_M1_abundance_stat[bacstability_M1_abundance_stat$Condition!='Total_space',], aes(Timepoint, mean_value, fill=Condition)) + 
  geom_bar(stat = "identity", position="stack") +
  geom_errorbar(data=bacstability_M1_abundance_stat[bacstability_M1_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=mean_value - sd_value, 
                    ymax=mean_value + sd_value, 
                colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Present species', color='') + 
  theme_bw() + 
  theme(axis.title =element_text(size=10,face="bold")) + 
  scale_fill_manual(values = c(Demuth[6:5], "#FFFFFF"), labels=c('Not present at M1', 'Retained from M1', '')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  guides(color="none")
dev.off()


######## will not be used probabl y########
stability <- as.data.frame(matrix(NA, nrow=length( unique(MGS_metadata$Individual_ID) ), ncol=6))
row.names(stability) <- unique(MGS_metadata$Individual_ID)
colnames(stability) <- c(paste0('Richness_', Infant_timepoints))

timepoint_pairs <- Map(function(x, y) c(x, y), Infant_timepoints[-length(Infant_timepoints)], Infant_timepoints[-1])

for (i in row.names(stability)) {
 
  for (timepoint in Infant_timepoints) {
    if ( length(MGS_metadata[MGS_metadata$Individual_ID == i & MGS_metadata$Timepoint==timepoint,]$Short_sample_ID_bact) != 0 ) {
      stability[i,paste0('Richness_', timepoint)] <- sum(species_infants[, MGS_metadata[MGS_metadata$Individual_ID == i & MGS_metadata$Timepoint==timepoint,]$Short_sample_ID_bact ] != 0)
    }
  }
  
  
  for (tp_pair in timepoint_pairs) {
    A <- MGS_metadata[MGS_metadata$Individual_ID == i & MGS_metadata$Timepoint == tp_pair[1], ]$Short_sample_ID_bact
    B <- MGS_metadata[MGS_metadata$Individual_ID == i & MGS_metadata$Timepoint == tp_pair[2], ]$Short_sample_ID_bact
    
    if (length(A) != 0 && length(B) != 0) {
      columns <- c(A, B)
      presence <- rowSums(species_infants[, columns] != 0) == 2
      stability[i, paste0(tp_pair[1], "_to_", tp_pair[2])] <- sum(presence)
      
    }
  }
}

# relocating M1_to_M2 for better representation:
stability <- stability[,c(  1:grep('M12', colnames(stability))[1], 
                            grep('M1_to_M2', colnames(stability)),  
                            (grep('M12', colnames(stability))[1] + 1):(grep('M1_to_M2', colnames(stability))-1)) ]

# for visualization:
for (i in Infant_timepoints) {
  
  
  if (i!='M1') {
    stability[,paste0('Not_retained_', i)] <- stability_M1[,paste0('Richness_', i)] - stability[,colnames(stability)[grep(i, colnames(stability))[2] ] ]
  } else {
    stability[,paste0('Not_retained_', i)] <- 0 
  }
  
  
}

stability_stat <- data.frame( colMeans(stability, na.rm = T), apply(stability, 2, sd, na.rm=T) )

colnames(stability_stat) <- c("mean_value", "sd_value")
stability_stat$Timepoint <- gsub(".*_", "", row.names(stability_stat))
stability_stat$Timepoint <- factor(stability_stat$Timepoint, levels=Infant_timepoints, ordered = T)
stability_stat["M1_to_M1",] <- stability_stat[grep('Richness_M1', row.names(stability_stat))[1],]
stability_stat$Condition <- gsub("_M.*", "", row.names(stability_stat))
stability_stat[grep('to', stability_stat$Condition),]$Condition <- 'Retained'

pdf('./04.PLOTS/BacSp_retained_from_previous_over_time_N.pdf', width=7/2.54, height=9/2.54)
ggplot(stability_stat[stability_stat$Condition!='Richness',], aes(Timepoint, mean_value, fill=Condition)) + 
  geom_bar(stat = "identity", position="stack") +
  geom_errorbar(data=stability_stat[stability_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=mean_value - sd_value, 
                    ymax=mean_value + sd_value, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "N detected species") + 
  labs(fill='Present species', color='') + 
  theme_bw() + 
  theme(axis.title =element_text(size=10,face="bold"),
        legend.position = 'bottom', 
        legend.text = element_text(size=5),
        legend.title = element_text(size=6)) + 
  scale_fill_manual(values = c(Demuth[6:5], "#FFFFFF"), labels=c('Not retained', 'Retained from previous', '')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  guides(color="none")
dev.off()

######## will not be used probably ########

bacpresence_from_M1 <- persistent_taxa(MGS_metadata[MGS_metadata$Type=='Infant',], species_infants, Infant_timepoints, "Short_sample_ID_bact")
bacpresence_from_M1_prevalent <- bacpresence_from_M1[bacpresence_from_M1$M1 > 0.5*length((MGS_metadata[MGS_metadata$Type=='Infant' & MGS_metadata$Timepoint=='M1',]$Short_sample_ID_bact)),]
row.names(bacpresence_from_M1_prevalent) <- gsub('.*s__','', row.names(bacpresence_from_M1_prevalent) )
row.names(bacpresence_from_M1_prevalent) <- gsub('_', ' ', row.names(bacpresence_from_M1_prevalent))
Okeeffe1_continuous <- met.brewer("OKeeffe1", n=26)
Okeeffe1_continuous <- Okeeffe1_continuous[26:1]

pdf('./04.PLOTS/BacSp_persistent_over_time.pdf', width=20/2.54, height=20/2.54)
heatmap.2(as.matrix(bacpresence_from_M1_prevalent), col=as.vector(Okeeffe1_continuous),
          #dendrogram ="none", 
          trace="none",
          #Rowv=FALSE, 
          Colv=FALSE,  
          margins=c(18,18)
          )
dev.off()


## Personal microbiome fractions

p_bac_frac_infants <- personal_biome_fraction(MGS_metadata[MGS_metadata$Type=='Infant',], species_infants, 'Short_sample_ID_bact')

PPB_bac_freq_infants <- as.data.frame(table(unlist(unname( p_bac_frac_infants[["PPB_bacteria"]] ))))

p_bac_frac_infants_individual <- p_bac_frac_infants[["personal_biome_fractions"]]
p_bac_frac_infants_individual$Infant <- row.names(p_bac_frac_infants_individual)
p_bac_frac_infants_individual <- melt(p_bac_frac_infants_individual[,c('PPB', 'TDB', 'Singletons', 'Infant')], id.vars = 'Infant')
p_bac_frac_infants_individual_df <- p_bac_frac_infants_individual[p_bac_frac_infants_individual$variable=='PPB',]
p_bac_frac_infants_individual_df <- p_bac_frac_infants_individual_df[order(p_bac_frac_infants_individual_df$value,decreasing = T),]
p_bac_frac_infants_individual$Infant <- factor(p_bac_frac_infants_individual$Infant, levels=p_bac_frac_infants_individual_df$Infant, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction.pdf', width=10/2.54, height=10/2.54)
ggplot(p_bac_frac_infants_individual, aes(Infant, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity") + 
  xlab(label = "Infant") + 
  ylab(label = "% of infant personal bacteriome") + 
  labs(fill='Fraction') +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5]) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.position = "bottom")
dev.off()

p_bac_frac_infants_per_sample <- p_bac_frac_infants[["personal_biome_per_sample"]]
p_bac_frac_infants_per_sample$Sample <- row.names(p_bac_frac_infants_per_sample)

p_bac_frac_infants_per_sample_ab <- melt(p_bac_frac_infants_per_sample[,c(colnames(p_bac_frac_infants_per_sample)[grep('ab_', colnames(p_bac_frac_infants_per_sample))], "Sample")], id='Sample')
p_bac_frac_infants_per_sample_ab$Infant <- MGS_metadata$Individual_ID[match(p_bac_frac_infants_per_sample_ab$Sample, MGS_metadata$Short_sample_ID_bact)] 
p_bac_frac_infants_per_sample_ab$Timepoint <- MGS_metadata$Timepoint[match(p_bac_frac_infants_per_sample_ab$Sample, MGS_metadata$Short_sample_ID_bact)]
p_bac_frac_infants_per_sample_ab$Infant <- factor(p_bac_frac_infants_per_sample_ab$Infant, levels=p_bac_frac_infants_individual_df$Infant, ordered = T)
p_bac_frac_infants_per_sample_ab$Timepoint <- factor(p_bac_frac_infants_per_sample_ab$Timepoint, levels=Infant_timepoints, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction_abundance.pdf', width=25/2.54, height=7/2.54)
ggplot(p_bac_frac_infants_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Infant samples") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Fraction') +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPB', 'TDB', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()


plot_bac_frac_infants_per_sample_N <- melt(p_bac_frac_infants_per_sample[,c(colnames(p_bac_frac_infants_per_sample)[grep('N_', colnames(p_bac_frac_infants_per_sample))], "Sample")], id='Sample')
plot_bac_frac_infants_per_sample_N$Infant <- MGS_metadata$Individual_ID[match(plot_bac_frac_infants_per_sample_N$Sample, MGS_metadata$Short_sample_ID_bact)] 
plot_bac_frac_infants_per_sample_N$Timepoint <- MGS_metadata$Timepoint[match(plot_bac_frac_infants_per_sample_N$Sample, MGS_metadata$Short_sample_ID_bact)]
plot_bac_frac_infants_per_sample_N$Infant <- factor(plot_bac_frac_infants_per_sample_N$Infant, levels=p_bac_frac_infants_individual_df$Infant, ordered = T)
plot_bac_frac_infants_per_sample_N$Timepoint <- factor(plot_bac_frac_infants_per_sample_N$Timepoint, levels=Infant_timepoints, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction_N.pdf', width=25/2.54, height=7/2.54)
ggplot(plot_bac_frac_infants_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Infant samples") + 
  ylab(label = "% of sample richness") + 
  labs(fill='Fraction') +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPB', 'TDB', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()

### Mother Bacteria ###
bacstability_P3 <- stability_initial(MGS_metadata[MGS_metadata$Type=='Mother',], species_mothers, Mother_timepoints, 'Short_sample_ID_bact', 'Richness')
bacstability_P3_stat <- summary_stat(bacstability_P3, Mother_timepoints)


bacstability_P3_abundance <- stability_initial(MGS_metadata[MGS_metadata$Type=='Mother',], species_mothers, Mother_timepoints, 'Short_sample_ID_bact', 'abundance')
bacstability_P3_abundance_stat <- summary_stat(bacstability_P3_abundance, Mother_timepoints)


pdf('./04.PLOTS/BacSp_retained_from_P3_over_time.pdf', width=10/2.54, height=7/2.54)
ggplot(bacstability_P3_stat[bacstability_P3_stat$Condition!='Richness',], aes(Timepoint, mean_value, fill=Condition)) + 
  geom_bar(stat = "identity", position="stack") +
  geom_errorbar(data=bacstability_P3_stat[bacstability_P3_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=mean_value - sd_value , 
                    ymax=mean_value + sd_value , 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "N detected species") + 
  labs(fill='Present species', color='') + 
  theme_bw() + 
  theme(axis.title =element_text(size=10,face="bold")) + 
  scale_fill_manual(values = c(Demuth[6:5], "#FFFFFF"), labels=c('Not present at P3', 'Retained from P3', '')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  guides(color="none")
dev.off()

pdf('./04.PLOTS/BacSp_retained_from_P3_over_time_abundance.pdf', width=10/2.54, height=7/2.54)
ggplot(bacstability_P3_abundance_stat[bacstability_P3_abundance_stat$Condition!='Total_space',], aes(Timepoint, mean_value, fill=Condition)) + 
  geom_bar(stat = "identity", position="stack") +
  geom_errorbar(data=bacstability_P3_abundance_stat[bacstability_P3_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=mean_value - sd_value, 
                    ymax=mean_value + sd_value, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Present species', color='') + 
  theme_bw() + 
  theme(axis.title =element_text(size=10,face="bold")) + 
  scale_fill_manual(values = c(Demuth[6:5], "#FFFFFF"), labels=c('Not present at P3', 'Retained from P3', '')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  guides(color="none")
dev.off()

bacpresence_from_P3 <- persistent_taxa(MGS_metadata[MGS_metadata$Type=='Mother',], species_mothers, Mother_timepoints, "Short_sample_ID_bact")
bacpresence_from_P3_prevalent <- bacpresence_from_P3[bacpresence_from_P3$P3 > 0.5*length((MGS_metadata[MGS_metadata$Type=='Mother' & MGS_metadata$Timepoint=='P3',]$Short_sample_ID_bact)),]
row.names(bacpresence_from_P3_prevalent) <- gsub('.*s__','', row.names(bacpresence_from_P3_prevalent) )
row.names(bacpresence_from_P3_prevalent) <- gsub('_', ' ', row.names(bacpresence_from_P3_prevalent))
Okeeffe1_continuous <- met.brewer("OKeeffe1", n=21)
Okeeffe1_continuous <- Okeeffe1_continuous[21:1]

pdf('./04.PLOTS/BacSp_persistent_over_time_mothers.pdf', width=20/2.54, height=20/2.54)
heatmap.2(as.matrix(bacpresence_from_P3_prevalent), col=as.vector(Okeeffe1_continuous),
          #dendrogram ="none", 
          trace="none",
          #Rowv=FALSE, 
          Colv=FALSE,  
          margins=c(18,18)
)
dev.off()

p_bac_frac_mothers <- personal_biome_fraction(MGS_metadata[MGS_metadata$Type=='Mother',], species_mothers, 'Short_sample_ID_bact')

PPB_bac_freq_mothers <- as.data.frame(table(unlist(unname( p_bac_frac_mothers[["PPB_bacteria"]] ))))

p_bac_frac_mothers_individual <- p_bac_frac_mothers[["personal_biome_fractions"]]
# since it is impossible to judge about the stability in 2 timepoints
p_bac_frac_mothers_individual <- p_bac_frac_mothers_individual[p_bac_frac_mothers_individual$N_timepoints > 2,]
p_bac_frac_mothers_individual$Mother <- row.names(p_bac_frac_mothers_individual)
p_bac_frac_mothers_individual <- melt(p_bac_frac_mothers_individual[,c('PPB', 'TDB', 'Singletons', 'Mother')], id.vars = 'Mother')
p_bac_frac_mothers_individual_df <- p_bac_frac_mothers_individual[p_bac_frac_mothers_individual$variable=='PPB',]
p_bac_frac_mothers_individual_df <- p_bac_frac_mothers_individual_df[order(p_bac_frac_mothers_individual_df$value,decreasing = T),]
p_bac_frac_mothers_individual$Mother <- factor(p_bac_frac_mothers_individual$Mother, levels=p_bac_frac_mothers_individual_df$Mother, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction_mothers.pdf', width=10/2.54, height=10/2.54)
ggplot(p_bac_frac_mothers_individual, aes(Mother, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity") + 
  xlab(label = "Mother") + 
  ylab(label = "% of maternal personal bacteriome") + 
  labs(fill='Fraction') +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5]) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.position = "bottom")
dev.off()

p_bac_frac_mothers_per_sample <- p_bac_frac_mothers[["personal_biome_per_sample"]]
p_bac_frac_mothers_per_sample$Sample <- row.names(p_bac_frac_mothers_per_sample)

p_bac_frac_mothers_per_sample_ab <- melt(p_bac_frac_mothers_per_sample[,c(colnames(p_bac_frac_mothers_per_sample)[grep('ab_', colnames(p_bac_frac_mothers_per_sample))], "Sample")], id='Sample')
p_bac_frac_mothers_per_sample_ab$Mother <- MGS_metadata$Individual_ID[match(p_bac_frac_mothers_per_sample_ab$Sample, MGS_metadata$Short_sample_ID_bact)] 
# excluding those that have only 2 timepoints: 
p_bac_frac_mothers_per_sample_ab <- p_bac_frac_mothers_per_sample_ab[p_bac_frac_mothers_per_sample_ab$Mother %in% p_bac_frac_mothers_individual$Mother,]
p_bac_frac_mothers_per_sample_ab$Timepoint <- MGS_metadata$Timepoint[match(p_bac_frac_mothers_per_sample_ab$Sample, MGS_metadata$Short_sample_ID_bact)]
p_bac_frac_mothers_per_sample_ab$Mother <- factor(p_bac_frac_mothers_per_sample_ab$Mother, levels=p_bac_frac_mothers_individual_df$Mother, ordered = T)
p_bac_frac_mothers_per_sample_ab$Timepoint <- factor(p_bac_frac_mothers_per_sample_ab$Timepoint, levels=Mother_timepoints, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction_abundance_mothers.pdf', width=25/2.54, height=7/2.54)
ggplot(p_bac_frac_mothers_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Maternal samples") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Fraction') +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPB', 'TDB', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()


plot_bac_frac_mothers_per_sample_N <- melt(p_bac_frac_mothers_per_sample[,c(colnames(p_bac_frac_mothers_per_sample)[grep('N_', colnames(p_bac_frac_mothers_per_sample))], "Sample")], id='Sample')
plot_bac_frac_mothers_per_sample_N$Mother <- MGS_metadata$Individual_ID[match(plot_bac_frac_mothers_per_sample_N$Sample, MGS_metadata$Short_sample_ID_bact)] 
# excluding those that have only 2 timepoints: 
plot_bac_frac_mothers_per_sample_N <- plot_bac_frac_mothers_per_sample_N[plot_bac_frac_mothers_per_sample_N$Mother %in% p_bac_frac_mothers_individual$Mother,]
plot_bac_frac_mothers_per_sample_N$Timepoint <- MGS_metadata$Timepoint[match(plot_bac_frac_mothers_per_sample_N$Sample, MGS_metadata$Short_sample_ID_bact)]
plot_bac_frac_mothers_per_sample_N$Mother <- factor(plot_bac_frac_mothers_per_sample_N$Mother, levels=p_bac_frac_mothers_individual_df$Mother, ordered = T)
plot_bac_frac_mothers_per_sample_N$Timepoint <- factor(plot_bac_frac_mothers_per_sample_N$Timepoint, levels=Mother_timepoints, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction_N_mothers.pdf', width=25/2.54, height=7/2.54)
ggplot(plot_bac_frac_mothers_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Maternal samples") + 
  ylab(label = "% of sample richness") + 
  labs(fill='Fraction') +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPB', 'TDB', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()

###### OUTPUT #####