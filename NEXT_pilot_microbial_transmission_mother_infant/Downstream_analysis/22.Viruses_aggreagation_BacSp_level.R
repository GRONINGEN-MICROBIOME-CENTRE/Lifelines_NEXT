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
##############################
# Input data
##############################

host_assignment_species <- read.csv('01.RAW_DATA/Host_assignment/MERGED_Host_prediction_to_genome_m90.csv', header=T)
# pre-cleaning:
host_assignment_species$Genus <- gsub('.*g__','', sapply(strsplit(host_assignment_species$Host.taxonomy, '\\;'), "[", 6))
host_assignment_species$Species <- gsub('.*s__','', sapply(strsplit(host_assignment_species$Host.taxonomy, '\\;'), "[", 7))
host_assignment_species <- host_assignment_species[host_assignment_species$Species!='',]
host_assignment_species <- host_assignment_species[order(host_assignment_species$Virus, host_assignment_species$Confidence.score ), ]
host_assignment_species_unique <- host_assignment_species[ !duplicated(host_assignment_species$Virus), ]
host_assignment_species_unique$Species_new <- gsub(' ', '_', host_assignment_species_unique$Species)

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep = '\t', header=T)


microbiome <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header = T)
microbiome_species <- data.frame(row.names(microbiome))
colnames(microbiome_species)[1] <- 'Taxonomy'
microbiome_species$Genus <- gsub('.*g__','', sapply(strsplit(microbiome_species$Taxonomy, '\\|'), "[", 6))
microbiome_species$Species <- gsub('.*s__','', sapply(strsplit(microbiome_species$Taxonomy, '\\|'), "[", 7))
##############################
# ANALYSIS
##############################
## Creating species_guild based table:
RPKM_VLP_by_BacSp <- data.frame(matrix(ncol = ncol(RPKM_counts_VLP), nrow=length(unique(host_assignment_species_unique$Species_new)) ) )
colnames(RPKM_VLP_by_BacSp) <- colnames(RPKM_counts_VLP)
row.names(RPKM_VLP_by_BacSp) <- unique(host_assignment_species_unique$Species_new)

for (i in colnames(RPKM_VLP_by_BacSp) ) {
  
  for (j in row.names(RPKM_VLP_by_BacSp) ) {
    
    RPKM_VLP_by_BacSp[row.names(RPKM_VLP_by_BacSp)==j,i] <- sum( RPKM_counts_VLP[ c( host_assignment_species_unique[ host_assignment_species_unique$Species_new==j, ]$Virus ) , i] )   
    
  }
  
}
RPKM_VLP_by_BacSp_raw  <- RPKM_VLP_by_BacSp


# filtering
RPKM_VLP_by_BacSp_filt <- RPKM_VLP_by_BacSp[(rowSums(RPKM_VLP_by_BacSp!=0) > 0.05*ncol(RPKM_VLP_by_BacSp)),  ]

#### refining host assignment taxa based on most prevalent species and microbiome notation
host_assignment_species_unique <- host_assignment_species_unique[host_assignment_species_unique$Species_new %in% row.names(RPKM_VLP_by_BacSp_filt),]
host_assignment_species_unique$Species_new <- gsub('_[[:upper:]]', '',host_assignment_species_unique$Species_new)
# there is a typo in Metaphlan:
host_assignment_species_unique$Species_new <- sub('Hydrogeniiclostridium', 'Hydrogeniiclostidium', host_assignment_species_unique$Species_new)

### remaking the species_guild for phages:
RPKM_VLP_by_BacSp <- data.frame(matrix(ncol = ncol(RPKM_counts_VLP), nrow=length(unique(host_assignment_species_unique$Species_new)) ) )
colnames(RPKM_VLP_by_BacSp) <- colnames(RPKM_counts_VLP)
row.names(RPKM_VLP_by_BacSp) <- unique(host_assignment_species_unique$Species_new)

for (i in colnames(RPKM_VLP_by_BacSp) ) {
  
  for (j in row.names(RPKM_VLP_by_BacSp) ) {
    
    RPKM_VLP_by_BacSp[row.names(RPKM_VLP_by_BacSp)==j,i] <- sum( RPKM_counts_VLP[ c( host_assignment_species_unique[ host_assignment_species_unique$Species_new==j, ]$Virus ) , i] )   
    
  }
  
}

# no need to filter at 0.05 prevalence, as they have been already filtered

# CLR transformation:
RPKM_VLP_by_BacSp_filt <- as.data.frame(RPKM_VLP_by_BacSp)
my_pseudocount_normal=min(RPKM_VLP_by_BacSp_filt[RPKM_VLP_by_BacSp_filt!=0])/2
RPKM_VLP_by_BacSp_filt_CLR<-decostand(RPKM_VLP_by_BacSp_filt, "clr", pseudocount=my_pseudocount_normal)


##############################
# OUTPUT
##############################
write.table(RPKM_VLP_by_BacSp_raw, '02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_raw.txt', sep='\t', quote = F)
write.table(RPKM_VLP_by_BacSp, '02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_filt.txt', sep='\t', quote=F)
write.table(RPKM_VLP_by_BacSp_filt_CLR, '02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_filtered_CLR_transformed.txt', sep='\t', quote=F)

