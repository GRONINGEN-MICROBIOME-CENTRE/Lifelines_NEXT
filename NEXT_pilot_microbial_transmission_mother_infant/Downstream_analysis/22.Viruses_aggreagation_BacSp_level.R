setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore differential abundance of gut bacterial
# species and viral vOTUs between maternal and infant gut 
#############################################################


##############################
# Functions
##############################

##############################
# Loading libraries
##############################
library(vegan)
library(dplyr)
##############################
# Input data
##############################

host_assignment_species <- read.csv('01.RAW_DATA/Host_assignment/MERGED_Host_prediction_to_genome_m90.csv', header=T)

# pre-refining:

# order:
host_assignment_species[grep('o__;f__;g__;s__$', host_assignment_species$Host.taxonomy),'Host.taxonomy'] <- gsub('c__([^;]+);o__;f__;g__;s__$', 
                                                                                                                 'c__\\1;o__\\1_unclassified;f__\\1_unclassified;g__\\1_unclassified;s__\\1 unclassified', 
                                                                                                                 host_assignment_species[grep('o__;f__;g__;s__$', host_assignment_species$Host.taxonomy),]$Host.taxonomy)
host_assignment_species[grep('o__;', host_assignment_species$Host.taxonomy),'Host.taxonomy'] <- gsub('c__([^;]+);o__;', 
                                                                                                     'c__\\1;o__\\1_unclassified;', 
                                                                                                     host_assignment_species[grep('o__;', host_assignment_species$Host.taxonomy),]$Host.taxonomy)

# family:
host_assignment_species[grep('f__;g__;', host_assignment_species$Host.taxonomy),'Host.taxonomy'] <- gsub('o__([^;]+);f__;g__;s__', 
                                                                                                         'o__\\1;f__\\1_unclassified;g__\\1_unclassified;s__\\1 unclassified', 
                                                                                                         host_assignment_species[grep('f__;g__;', host_assignment_species$Host.taxonomy),]$Host.taxonomy)
host_assignment_species[grep('f__;', host_assignment_species$Host.taxonomy),'Host.taxonomy'] <- gsub('o__([^;]+);f__;', 
                                                                                                     'o__\\1;f__\\1_unclassified;', 
                                                                                                     host_assignment_species[grep('f__;g__', host_assignment_species$Host.taxonomy),]$Host.taxonomy)
# genus:
host_assignment_species[grep('g__;s__$', host_assignment_species$Host.taxonomy),'Host.taxonomy'] <- gsub('f__([^;]+);g__;s__$', 
                                                                                                         'f__\\1;g__\\1_unclassified;s__\\1 unclassified', 
                                                                                                         host_assignment_species[grep('g__;s__$', host_assignment_species$Host.taxonomy),]$Host.taxonomy)
host_assignment_species[grep('g__;s__', host_assignment_species$Host.taxonomy),'Host.taxonomy'] <- gsub('f__([^;]+);g__;', 
                                                                                                        'f__\\1;g__\\1_unclassified;', 
                                                                                                        host_assignment_species[grep('g__;s__', host_assignment_species$Host.taxonomy),]$Host.taxonomy)
# species: 
host_assignment_species[grep('s__$', host_assignment_species$Host.taxonomy),'Host.taxonomy'] <- gsub('g__([^;]+);s__$', 
                                                                                                     'g__\\1;s__\\1 unclassified', 
                                                                                                     host_assignment_species[grep('s__$', host_assignment_species$Host.taxonomy),]$Host.taxonomy)


host_assignment_species$Kingdom <- gsub('.*d__','', sapply(strsplit(host_assignment_species$Host.taxonomy, '\\;'), "[", 1))
host_assignment_species$Phylum <- gsub('.*p__','', sapply(strsplit(host_assignment_species$Host.taxonomy, '\\;'), "[", 2))
host_assignment_species$Class <- gsub('.*c__','', sapply(strsplit(host_assignment_species$Host.taxonomy, '\\;'), "[", 3))
host_assignment_species$Order <- gsub('.*o__','', sapply(strsplit(host_assignment_species$Host.taxonomy, '\\;'), "[", 4))
host_assignment_species$Family <- gsub('.*f__','', sapply(strsplit(host_assignment_species$Host.taxonomy, '\\;'), "[", 5))
host_assignment_species$Genus <- gsub('.*g__','', sapply(strsplit(host_assignment_species$Host.taxonomy, '\\;'), "[", 6))
host_assignment_species$Species <- gsub('.*s__','', sapply(strsplit(host_assignment_species$Host.taxonomy, '\\;'), "[", 7))

full_taxonomy <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

for (i in full_taxonomy) {
  host_assignment_species[,i] <- gsub('_([[:upper:]])', '',host_assignment_species[,i], perl=T)
}

host_assignment_species$Host.taxonomy.refined <- paste0('d__', host_assignment_species$Kingdom, ';',
                                                        'p__', host_assignment_species$Phylum, ';',
                                                        'c__', host_assignment_species$Class, ';',
                                                        'o__', host_assignment_species$Order, ';',
                                                        'f__', host_assignment_species$Family, ';',
                                                        'g__', host_assignment_species$Genus, ';',
                                                        's__', host_assignment_species$Species)

host_assignment_species <- host_assignment_species[order(host_assignment_species$Virus, host_assignment_species$Confidence.score ), ]
host_assignment_species_unique <- host_assignment_species[ !duplicated(host_assignment_species$Virus), ]

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep = '\t', header=T)
RPKM_counts_VLP$Virus <- row.names(RPKM_counts_VLP)

host_assignment_species_unique_combined <- merge(host_assignment_species_unique, RPKM_counts_VLP, by='Virus', all.y = T)
host_assignment_species_unique_combined[is.na(host_assignment_species_unique_combined$Host.taxonomy.refined),]$Host.taxonomy.refined <- 'Unassigned'

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)

host_assignment_species_unique_combined <- host_assignment_species_unique_combined[,c('Host.taxonomy.refined',VLP_metadata$Short_sample_ID)]

##############################
# ANALYSIS
##############################
## Creating species_guild based table:
RPKM_VLP_by_BacSp <- aggregate(.~Host.taxonomy.refined, host_assignment_species_unique_combined, sum)
row.names(RPKM_VLP_by_BacSp) <- RPKM_VLP_by_BacSp$Host.taxonomy.refined
RPKM_VLP_by_BacSp$Host.taxonomy.refined <- NULL

# filtering
RPKM_VLP_by_BacSp_filt <- RPKM_VLP_by_BacSp[(rowSums(RPKM_VLP_by_BacSp!=0) > 0.05*ncol(RPKM_VLP_by_BacSp)),  ]

# CLR transformation:
RPKM_VLP_by_BacSp_filt_t <- as.data.frame(t(RPKM_VLP_by_BacSp_filt))
my_pseudocount_normal=min(RPKM_VLP_by_BacSp_filt_t[RPKM_VLP_by_BacSp_filt_t!=0])/2
RPKM_VLP_by_BacSp_filt_CLR<-decostand(RPKM_VLP_by_BacSp_filt_t, "clr", pseudocount=my_pseudocount_normal)


##############################
# OUTPUT
##############################
write.table(host_assignment_species_unique, '02.CLEAN_DATA/Host_prediction_to_genome_m90_refined_taxonomy_no_generalists.txt', sep='\t', quote=F, row.names=F)
write.table(RPKM_VLP_by_BacSp, '02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_species.txt', sep='\t', quote = F)
write.table(RPKM_VLP_by_BacSp_filt, '02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_species_filt.txt', sep='\t', quote=F)
write.table(RPKM_VLP_by_BacSp_filt_CLR, '02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_species_filtered_CLR_transformed.txt', sep='\t', quote=F)

