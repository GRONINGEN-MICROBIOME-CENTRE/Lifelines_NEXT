setwd("~/Desktop/Projects_2018/Baby_pilot/01.Data_analysis/03.SCRIPTS/Clean_code/04.2.Virome_imputed_stat/")
load('.RData')

##### LIBRARIES #####
library(tidyverse)
library(ggplot2)
library(vegan)
library(ggExtra)
library(MASS)
library(extrafont)
font_import()
loadfonts(device="pdf")
library(ggpubr)
##### FUNCTIONS #####

#BC_dist function takes counts_table in the format rows==bacteria, columns==samples and
#samples metadata. BC_dist returns a list of two objects: distance matrix that can be accessed
#my_list[[1]] and pcoa data merged with samples metadata, my_list[[2]]
BC_dist <- function(counts_table, metadata){
  dist.bray = vegdist(as.data.frame(t(counts_table)), method = "bray")
  dist_matrix <- as.matrix(dist.bray)
  my_pcoa <- as.data.frame(cmdscale(dist.bray, k=5))
  my_pcoa <- merge(my_pcoa, metadata, by='row.names')
  row.names(my_pcoa) <- my_pcoa$Row.names
  my_pcoa$Row.names <- NULL
  output <- list(dist_matrix, my_pcoa)
  return(output)
}

##### INPUT DATA #####
viruses_tm_shared <- read.table('../INTERMEDIATE_OUTPUT/baby_shared_imputed_data.txt', sep='\t', header=T)
viruses_tm_shared$n_not_sh_own_m <- viruses_tm_shared$n_total - viruses_tm_shared$n_sh_own_m
viruses_tm_shared$perc_shared <- viruses_tm_shared$n_sh_own_m/viruses_tm_shared$n_total
viruses_tm_shared$perc_not_shared <- viruses_tm_shared$n_not_sh_own_m/viruses_tm_shared$n_total

viruses_vm_shared <- read.table('../../Analysis_for_paper/INTERMEDIATE_OUTPUT/baby_shared_viruses_data.txt', sep='\t',header = T)
row.names(viruses_vm_shared) <- viruses_vm_shared$new_sample_id
viruses_vm_shared$n_not_sh_own_m <- viruses_vm_shared$n_total - viruses_vm_shared$n_sh_own_m
viruses_vm_shared$perc_shared <- viruses_vm_shared$n_sh_own_m/viruses_vm_shared$n_total
viruses_vm_shared$perc_not_shared <- viruses_vm_shared$n_not_sh_own_m/viruses_vm_shared$n_total

bacteria_shared <- read.table('../../Analysis_for_paper/02.Microbiome_shared/baby_shared_bactSP_data_with_filters.txt', sep='\t', header=T)
row.names(bacteria_shared) <- bacteria_shared$new_sample_id
bacteria_shared$n_not_sh_own_m <- bacteria_shared$n_total - bacteria_shared$n_sh_own_m
bacteria_shared$perc_shared <- bacteria_shared$n_sh_own_m/bacteria_shared$n_total
bacteria_shared$perc_not_shared <- bacteria_shared$n_not_sh_own_m/bacteria_shared$n_total

virome_imputed <- read.table('../INPUT/MGS_RPKM_counts_LLNEXT_pilot_updated.txt', sep='\t', header=T)
virome_imputed <- virome_imputed[,colnames(virome_imputed) %in% row.names(samples_metadata)]
virome_imputed <- virome_imputed[,order(colnames(virome_imputed))]
samples_metadata <- read.table('../INPUT/Microbiome_metadata_phenos.txt', sep='\t', header=T)
virome_imputed <- virome_imputed[rowSums(virome_imputed)>0,]

##### PROCESSING #####

wilcox.test(viruses_tm_shared[row.names(viruses_vm_shared),]$perc_shared, viruses_vm_shared$perc_shared, paired=T)
#V = 189, p-value = 7.629e-06

cor.test(bacteria_shared[row.names(viruses_vm_shared),]$perc_shared,viruses_vm_shared$perc_shared)
# Pearson's product-moment correlation: cor=0.2879098, p-value = 0.232

cor.test(viruses_tm_shared$perc_shared,bacteria_shared$perc_shared)
#Pearson's product-moment correlation, cor=0.4718302, p-value = 0.006405

sharedness <- cbind(bacteria_shared[,"perc_shared", drop=F], viruses_tm_shared[,"perc_shared", drop=F])
colnames(sharedness) <- c('perc_shared_bacteria','perc_shared_viruses_in_tm')
sharedness <- merge(sharedness, viruses_vm_shared[,'perc_shared', drop=F], all.x = T, by='row.names')
row.names(sharedness) <- sharedness$Row.names
sharedness$Row.names <- NULL
colnames(sharedness)[3] <- 'perc_shared_viruses_in_vm'

ggplot(sharedness, aes(perc_shared_bacteria, perc_shared_viruses_in_vm)) + 
  geom_point(shape=1, color="#00BFC4", size=3, stroke=0.8) + 
  geom_smooth(aes(perc_shared_bacteria, perc_shared_viruses_in_vm), method = 'lm', color="#09009B") + 
  labs(x="Shared bacterial species in mother-infant pairs (%)", y="Shared viruses in mother-infant pairs (%), VM") +
  theme_classic() + theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + 
  annotate("text", x=0.45, y=0.6, label="R=0.29\np-value=0.23", size=5)

ggplot(sharedness, aes(perc_shared_bacteria, perc_shared_viruses_in_tm)) + 
  geom_point(shape=1, color="#00BFC4", size=3, stroke=0.8) + 
  geom_smooth(aes(perc_shared_bacteria, perc_shared_viruses_in_tm), method = 'lm', color="#09009B") + 
  labs(x="Shared bacterial species in mother-infant pairs (%)", y="Shared viruses from TM in mother-infant pairs (%), TM") +
  theme_classic() + theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + 
  annotate("text", x=0.45, y=0.6, label="R=0.47\np-value=0.006", size=5)

# PCoA for imputed virome

dist.virome.bray <- BC_dist(virome_imputed, samples_metadata)

beta.cmd=cmdscale(as.matrix(dist.virome.bray[[1]]),k=10,eig = T) 
PCoA1=beta.cmd$eig[1]/sum(beta.cmd$eig)*100 
PCoA2=beta.cmd$eig[2]/sum(beta.cmd$eig)*100 
PCoA3=beta.cmd$eig[3]/sum(beta.cmd$eig)*100 
PCoA4=beta.cmd$eig[4]/sum(beta.cmd$eig)*100 
PCoA5=beta.cmd$eig[5]/sum(beta.cmd$eig)*100 
beta.cmd=data.frame(beta.cmd$points)
bc_virome_plot <- dist.virome.bray[[2]]

p <- ggplot(bc_virome_plot, aes(V1,V2, color=as.factor(Type))) + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Type, color=Type), linetype = 2)+
  geom_point(size=3) + 
  theme_bw() + 
  xlab(label = "PC1 (16.2% of total variation)") + 
  ylab(label = "PC2 (5.8% of total variation)") + 
  labs(color='Type') +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=16,face="bold"), 
        legend.title = element_text(size=14,face="bold"), 
        legend.text=element_text(size=14),
        legend.position = "bottom")
S4A <- ggMarginal(p, type="densigram", groupFill=T)

samples_metadata_plot <- samples_metadata[,"Type", drop=F]
samples_metadata_plot <- samples_metadata_plot[order(row.names(samples_metadata_plot)), , drop=F]
identical(row.names(samples_metadata_plot), row.names(dist.virome.bray[[1]]))
ord <- metaMDS(dist.virome.bray[[1]], distance = "bray", k=5)
en = envfit(ord, samples_metadata_plot, permutations = 999, na.rm = TRUE)
en$factors #Type 0.69  0.001 ***; 

# Heatmap for imputed virome based on BC


pdf("FigureS4A.pdf", width=21/2.54, height=10/2.54, family = "Arial")

dev.off()

##### OUTPUT #####
write.table(sharedness, '../INTERMEDIATE_OUTPUT/Figure_5DE_data.txt', sep='\t', quote=F)

