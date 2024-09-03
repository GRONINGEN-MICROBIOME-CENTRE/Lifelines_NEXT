##########################################
# Getting only virus contigs with at least 
# Medium-quality virus contigs to increase
# the speed of dereplication 
##########################################
args = commandArgs(trailingOnly=TRUE)
##############################
# Loading libraries
##############################

##############################
# Functions
##############################

##############################
# Input data
##############################

##############################
# ANALYSIS
##############################
Extended_TOF <- read.table(paste0(args[1]), sep='\t', header=T)

Extended_TOF <- Extended_TOF[Extended_TOF$POST_CHV_length >= 1000 &
                            !(grepl('longer than expected genome length', Extended_TOF$warnings)) & 
                            (Extended_TOF$viral_genes >= Extended_TOF$host_genes) & 
			    Extended_TOF$plasmid=="No"
                                       ,]

# Simulating ETOF for Viral RefSeq 
VREF_ETOF <- read.table('/scratch/p282752/databases/viral_refseq_apr_24/Extended_TOF', sep='\t', header=T)
VREF_ETOF <- VREF_ETOF[!duplicated(VREF_ETOF$Original_CID),]
VREF_ETOF$New_CID <- VREF_ETOF$Original_CID
VREF_ETOF$provirus <- "No"
VREF_ETOF$POST_CHV_length <- VREF_ETOF$Original_length
VREF_ETOF$GND_pruned <- "No"
VREF_ETOF$CHV_pruned <- "No"
##############################
# OUTPUT
##############################
write.table(Extended_TOF[,"New_CID"], paste0(dirname(args[1]), '/QualFilt_virus_contigs_IDs'), sep='\t', row.names=F, col.names=F, quote=F)
write.table(Extended_TOF, paste0(dirname(args[1]), '/Extended_TOF_filtered'), sep='\t', row.names=F, col.names=T, quote=F)
write.table(VREF_ETOF, "/scratch/p282752/databases/viral_refseq_apr_24/Extended_TOF_simulated", sep='\t', row.names=F, col.names=T, quote=F)
