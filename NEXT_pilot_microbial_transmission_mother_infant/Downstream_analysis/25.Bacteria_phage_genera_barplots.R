setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the top prevalent genera and their phages
# aggregated at the level of genera of the host
#############################################################
# Authors: Nataliia Kuzub & Sanzhima Garmaeva

##############################
# Functions
##############################
taxonomy_abundance <- function(taxonomy_table) {
  my_results=matrix(ncol = 5, nrow=ncol(taxonomy_table))
  rownames(my_results) = colnames(taxonomy_table)
  colnames(my_results) = c("Mean","Prevalence", "nz_mean", "N_of_0", "perc_missing") 
  ##Function to calculate mean excluding 0 values
  nzmean <- function(a){
    mean(a[a!=0])
  }
  ##Function to calculate nº of 0
  zsum <- function(a){
    sum (a==0)
  }
  ##Function to calculate nº of non-0
  nsum <- function(a){
    sum (a!=0)
  }
  ## Loop for each column (taxonomy) in the taxonomy table
  for (i in 1:ncol(taxonomy_table)) {
    #Calculate mean for each column
    aa = mean(taxonomy_table[,i])
    #Calculate number of non-zeros (individuals)
    bb = nsum(taxonomy_table[,i])
    #Calculate mean without taking into account the 0
    cc = nzmean(taxonomy_table[,i])
    #Calculate number of zeros 
    dd = zsum(taxonomy_table[,i])
    ee= (dd/(dd+bb))*100
    my_results[i,1] = aa
    my_results[i,2] = bb
    my_results[i,3] = cc
    my_results[i,4] = dd
    my_results[i,5] = ee
  }
  return(my_results)
}
##############################
# Loading libraries
##############################
library(dplyr)

library(ggplot2)
library(forcats)

library(patchwork)
library(data.table)

##############################
# Input data
##############################
top20_bacteria <- read.table('03a.RESULTS/Top20_prevalent_bacgenera_by_Type.txt', sep='\t', header=T)
top20_bacteria$Type <- factor(top20_bacteria$Type, levels=c('Mother', 'Infant'), ordered = T)

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)
VLP_metadata$Type <- factor(VLP_metadata$Type, levels=c('Infant', 'Mother'), ordered=T)

contigs_metadata <- read.table('02.CLEAN_DATA/VLP_viral_contigs_metadata.txt', sep='\t', header=T)
colnames(contigs_metadata)[1] <- 'Virus'

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)
RPKM_binary_VLP <- RPKM_counts_VLP
RPKM_binary_VLP[RPKM_binary_VLP!=0] <- 1
RPKM_binary_VLP$Virus = row.names(RPKM_binary_VLP)

RPKM_counts_VLP$Virus = row.names(RPKM_counts_VLP)

microbiome_all <- read.table('01.RAW_DATA/Metaphlan4_all_samples/LLNEXT_metaphlan_4_complete_10_02_2023.txt', sep='\t', header=T)
taxa_all <- data.frame(microbiome_all$clade_name)
taxa_all <- data.frame(taxa_all[grep('g__', taxa_all$microbiome_all.clade_name),])
taxa_all <- data.frame(taxa_all[grep('s__', taxa_all$taxa_all.grep..g__...taxa_all.microbiome_all.clade_name...., invert=T),])
colnames(taxa_all)[1] <- 'Host.genus'
taxa_all$Genus <- gsub('.*g__','', sapply(strsplit(taxa_all$Host.genus, '\\|'), "[", 6))

host_assignment <- read.csv('01.RAW_DATA/Host_assignment/MERGED_Host_prediction_to_genus_m90.csv')
# taxonomy refinement to match bacteriome taxonomy:
host_assignment$Host.genus <- gsub('_([[:upper:]]);', ';',host_assignment$Host.genus, perl=T)
host_assignment$Host.genus <- gsub('_[[:upper:]]$', '', host_assignment$Host.genus, perl = TRUE)
# manual work: if bacterium is not classified, adding '_unclassified' to the lowest known taxon
# order:
host_assignment[grep('o__;', host_assignment$Host.genus),'Host.genus'] <- 'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacilli_unclassified;f__Bacilli_unclassified;g__Bacilli_unclassified'
# family:
host_assignment[grep('f__;g__$', host_assignment$Host.genus),'Host.genus'] <- gsub('o__([^;]+);f__;g__', 'o__\\1;f__\\1_unclassified;g__\\1_unclassified', host_assignment[grep('f__;g__$', host_assignment$Host.genus),]$Host.genus)
# (exception with known genus but not family)
host_assignment[grep('o__UMGS1883;f__;', host_assignment$Host.genus),'Host.genus'] <- 'd__Bacteria;p__Firmicutes;c__Clostridia;o__UMGS1883;f__UMGS1883_unclassified;g__UMGS1540'
# genus:
host_assignment[grep('g__$', host_assignment$Host.genus),'Host.genus'] <- gsub('f__([^;]+);g__', 'f__\\1;g__\\1_unclassified', host_assignment[grep('g__$', host_assignment$Host.genus),]$Host.genus)
# manual work: renaming taxons that have different naming in MetaPhlan4 taxonomy
# genus:
host_assignment$Host.genus <- gsub('g__Paraclostridium', 'g__Paeniclostridium', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('g__ER4', 'g__Oscillospiraceae_unclassified', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('g__CAG-83', 'g__Oscillospiraceae_unclassified', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('g__Agathobacter', 'g__Lachnospiraceae_unclassified', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('g__F0422', 'g__Veillonella', host_assignment$Host.genus)
# phylum:
host_assignment$Host.genus <- gsub('Bacteroidota', 'Bacteroidetes', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('Verrucomicrobiota', 'Verrucomicrobia', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('Actinobacteriota', 'Actinobacteria', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('Synergistota', 'Synergistetes', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('Fusobacteriota', 'Fusobacteria', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('Spirochaetota', 'Spirochaetes', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('Planctomycetota', 'Planctomycetes', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('Elusimicrobiota', 'Elusimicrobia', host_assignment$Host.genus)
host_assignment$Host.genus <- gsub('Thermoplasmatota', 'Candidatus_Thermoplasmatota', host_assignment$Host.genus)

# in metaphlan, Euryarchaeota is a phylum:
host_assignment$Host.genus <- gsub('Methanobacteriota', 'Euryarchaeota', host_assignment$Host.genus)
# in metaphlan, Campylobacterota is Proteobacteria:
host_assignment$Host.genus <- gsub('Campylobacterota', 'Proteobacteria', host_assignment$Host.genus)
# since Cyanobacteria and Deinococcus_Thermus are present in the metaphlan database, but not in our samples:
host_assignment <- host_assignment[-grep('Cyanobacteria', host_assignment$Host.genus),]
host_assignment <- host_assignment[-grep('Deinococcota', host_assignment$Host.genus),]
# single predictions that are not in metaphlan:
host_assignment <- host_assignment[-grep('Thermotogota', host_assignment$Host.genus),]
host_assignment <- host_assignment[-grep('Dependentiae', host_assignment$Host.genus),]
host_assignment <- host_assignment[-grep('Aquificota', host_assignment$Host.genus),]
host_assignment <- host_assignment[-grep('Myxococcota', host_assignment$Host.genus),]

host_assignment$Phylum <- gsub('.*p__','', sapply(strsplit(host_assignment$Host.genus, '\\;'), "[", 2))
host_assignment$Genus <- gsub('.*g__','', sapply(strsplit(host_assignment$Host.genus, '\\;'), "[", 6))
host_assignment <- host_assignment[!duplicated(host_assignment[, c("Virus", "Genus")]), ]

# no generalists: 
host_assignment_unique <- host_assignment[order(host_assignment$Virus, host_assignment$Confidence.score),]
host_assignment_unique  <- host_assignment_unique [ !duplicated(host_assignment_unique$Virus), ]

##############################
# ANALYSIS
##############################

host_assignment_combined <- merge(host_assignment[, c("Virus", "Host.genus")], RPKM_binary_VLP, by="Virus")

host_infection_freq <- as.data.frame(table(host_assignment$Virus))
colnames(host_infection_freq) <- c("Virus", "N_host_genera")
host_infection_freq <- host_infection_freq %>% mutate(generalist = ifelse(N_host_genera > 1, 1, 0))
host_infection_freq <- merge(host_infection_freq, contigs_metadata, all.y = T, by='Virus')
# "generalist" if infects at least 2 genera
host_infection_freq[is.na(host_infection_freq[,c("N_host_genera")]), c("N_host_genera", "generalist")] <- 0

# keeping the length column to filter for length in case it's necessary
host_assignment_combined <- merge(host_infection_freq[,c("Virus", "N_host_genera", "generalist", "temperate", "length")], 
                                  host_assignment_combined, by='Virus')

# Mothers: 
mothers <- host_assignment_combined[,c("generalist", "temperate", "Host.genus", VLP_metadata[VLP_metadata$Type=='Mother',]$Short_sample_ID)]
mothers <- mothers[rowSums(mothers[,VLP_metadata[VLP_metadata$Type=='Mother',]$Short_sample_ID])>0,]
mothers$N_virus <- 1
mothers_genera <- aggregate(.~Host.genus, mothers, sum)

mothers_genera$Prevalence <- rowSums(mothers_genera[,VLP_metadata[VLP_metadata$Type=='Mother',]$Short_sample_ID]!=0)
mothers_genera$generalist_perc <- mothers_genera$generalist/mothers_genera$N_virus*100
mothers_genera$temperate_perc <- mothers_genera$temperate/mothers_genera$N_virus*100


top20_mothers <- mothers_genera[,c("Host.genus", "Prevalence", "N_virus", "generalist_perc", "temperate_perc")]
top20_mothers <- top20_mothers[order(top20_mothers$Prevalence, top20_mothers$N_virus, decreasing = F),]
top20_mothers <- top20_mothers[(nrow(top20_mothers) - 19):nrow(top20_mothers),]
top20_mothers$Type <- 'Mother'

# Infants: 
infants <- host_assignment_combined[,c("generalist", "temperate", "Host.genus", VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID)]
infants <- infants[rowSums(infants[,VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID])>0,]
infants$N_virus <- 1
infants_genera <- aggregate(.~Host.genus, infants, sum)

infants_genera$Prevalence <- rowSums(infants_genera[,VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID]!=0)
infants_genera$generalist_perc <- infants_genera$generalist/infants_genera$N_virus*100
infants_genera$temperate_perc <- infants_genera$temperate/infants_genera$N_virus*100


top20_infants <- infants_genera[,c("Host.genus", "Prevalence", "N_virus", "generalist_perc", "temperate_perc")]
top20_infants <- top20_infants[order(top20_infants$Prevalence, top20_infants$N_virus, decreasing = F),]
top20_infants <- top20_infants[(nrow(top20_infants) - 19):nrow(top20_infants),]
top20_infants$Type <- 'Infant'

top20 <- rbind(top20_mothers, top20_infants)
top20$Type <- factor(top20$Type, levels=c('Mother','Infant'), ordered=T)
top20$genus <- gsub('.*g__', '', top20$Host.genus)
top20$phylum <- gsub('.*p__','', sapply(strsplit(top20$Host.genus, '\\;'), "[", 2))
top20$N_virus <- log10(top20$N_virus)




### Bacteria
top20_bacteria$Prevalence_virus <- NA 
top20_bacteria[top20_bacteria$Type=='Mother',]$Prevalence_virus <- top20[top20$Type=='Mother',]$ord[match( top20_bacteria[top20_bacteria$Type=='Mother',]$genus, top20[top20$Type=='Mother',]$genus)]
top20_bacteria[top20_bacteria$Type=='Infant',]$Prevalence_virus <- top20[top20$Type=='Infant',]$ord[match(top20_bacteria[top20_bacteria$Type=='Infant',]$genus, top20[top20$Type=='Infant',]$genus)] 
top20_bacteria[is.na(top20_bacteria$Prevalence_virus),]$Prevalence_virus <- 1 #random


top20_bacteria <- as.data.table(top20_bacteria)
top20_bacteria[, ord := sprintf("%02i", frank(top20_bacteria, Type, Prevalence, Prevalence_virus, ties.method = "first"))]

p0 <- ggplot(top20_bacteria, aes(x = Prevalence, y = ord, fill = phylum)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y', drop = TRUE) +
  labs(y = "", x = "Detected in\nN samples", title = "Bacterial genera") +
  scale_y_discrete(labels = setNames(as.character(top20_bacteria$genus), top20_bacteria$ord)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = 'none',
        title =element_text(size=10)) +
  scale_fill_manual(values = c('#FF9C0E', '#1F77B4', '#ED3419', '#B446B3'))

### Viruses
top20 <- as.data.table(top20)
top20[, ord := sprintf("%02i", frank(top20, Type, Prevalence, ties.method = "first"))]

p1 <- ggplot(top20, aes(x = Prevalence, y = ord, fill = phylum)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y', drop = TRUE) +
  labs(y = "", x = "Detected in\nN samples", title = "vOTUs aggregated\nat bacterial genus\nlevel") +
  scale_y_discrete(labels = setNames(as.character(top20$genus), top20$ord)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = 'none',
        title =element_text(size=10)) +
  scale_fill_manual(values = c('#FF9C0E', '#1F77B4', '#ED3419', '#B446B3'))

p2 <- ggplot(top20, aes(x = N_virus, y = ord, fill = phylum)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y', drop = TRUE) +
  labs(y = "", x = "Number vOTUs\n assigned (log10)") +
  scale_y_discrete(labels = setNames(as.character(top20$genus), top20$ord)) +
  theme_bw() +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = 'none') +
  scale_fill_manual(values = c('#FF9C0E', '#1F77B4', '#ED3419', '#B446B3'))

p3 <- ggplot(top20, aes(x = generalist_perc, y = ord, fill = phylum)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y', drop = TRUE) +
  labs(y = "", x = "% generalists") +
  scale_y_discrete(labels = setNames(as.character(top20$genus), top20$ord)) +
  theme_bw() +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = 'none') +
  scale_fill_manual(values = c('#FF9C0E', '#1F77B4', '#ED3419', '#B446B3'))

p4 <- ggplot(top20, aes(x = temperate_perc, y = ord, fill = phylum)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Type, ncol = 1, scales = 'free_y', drop = TRUE, strip.position = 'right') +
  labs(y = "", x = "% temperate") +
  scale_y_discrete(labels = setNames(as.character(top20$genus), top20$ord)) +
  theme_bw() +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank()) +
  scale_fill_manual(values = c('#FF9C0E', '#1F77B4', '#ED3419', '#B446B3'))

# Combine the plots using patchwork
combined_plot <- p0 + p1 + p2 + p3 + p4

combined_plot <- combined_plot +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "")

# Print the combined plot
pdf('./04.PLOTS/Bacterial_genera_and_vOTUs_aggregated_at_bacgen_level_prevalence.pdf', width=40/2.54, height=22/2.54)
combined_plot
dev.off()


#### RPKM table no generalists:
host_assignment_unique_combined <- merge(host_assignment_unique[, c("Virus", "Host.genus")], RPKM_counts_VLP, by="Virus", all.y=T)
host_assignment_unique_combined[is.na(host_assignment_unique_combined$Host.genus),]$Host.genus <- 'Unassigned'
host_assignment_unique_combined <- host_assignment_unique_combined[,colnames(host_assignment_unique_combined)!='Virus']
RPKM_counts_VLP_by_BacGen <- aggregate(.~Host.genus, host_assignment_unique_combined, sum)

##### think about the fact that this prevalence table does not solve the problem with counting generalists and temperates per genus in infants and mothers separately

###### OUTPUT #####
write.table(host_assignment, '02.CLEAN_DATA/Host_prediction_to_genus_m90_refined_taxonomy.txt', sep='\t', quote=F, row.names=F)
write.table(host_assignment_unique, '02.CLEAN_DATA/Host_prediction_to_genus_m90_refined_taxonomy_no_generalists.txt', sep='\t', quote=F, row.names=F)
write.table(RPKM_counts_VLP_by_BacGen, '02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_genus.txt', sep='\t', quote=F, row.names=F)
