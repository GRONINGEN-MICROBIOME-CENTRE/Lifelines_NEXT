##########################################
# Designing new names and metadata 
# for extended and pruned control contigs
##########################################
args = commandArgs(trailingOnly=TRUE)
##############################
# Loading libraries
##############################
library(stringr)
##############################
# Functions
##############################

##############################
# Input data
##############################
contigs <- read.table(paste0(args[1], 'contigs_for_dereplication'), col.names = "contig_id")

contigs_length <- read.table(paste0(args[1], 'contigs_for_dereplication_length'), col.names = c("contig_id", "POST_CBR_length"))

VLP_contigs_PD_metadata <- merge(contigs, contigs_length, by="contig_id")

##############################
# ANALYSIS
##############################

##### new (provisional) contigs names: 

# Old name: CHV200013F12_NODE_28562_length_1013_cov_0.285553_extended_circular
# Changes:
# 0) add source (NEXT)
# 1) keep only VXXX from the SAMPLEID and change V to N (N2000) (V stands for virome, N stands for negative control)
# 3) merge NODE & its number (N28562)
# 4) merge (final after extension & pruning length?)	(L748)
# 5) COBRA status: extended or untouched? (E0)
# 6) CheckV: pruned or not?	(P0), always P0 here as no pruning has been performed on these contigs
# 7) fragment_number: all F0 since unpruned;

# New name: NEXT_N2000_N28562_L1013_cov_0.285553_E1_P0_F0


# this line prunes the post COBRA (so all potential _extended and prophage coordinates characters are trimmed away)
# it is very dependant on the number of underscores in the contig id, so if the sample name contains underscores, change the number from 6 to smth else
VLP_contigs_PD_metadata$Original_CID <- sub("(([^_]*_){6}[^_]*).*", "\\1", VLP_contigs_PD_metadata$contig_id)

# here, the original contig length is extracted from its id, so it is also dependant on "_" in contig & sample id
VLP_contigs_PD_metadata$Original_length <- as.numeric(sapply(str_split(VLP_contigs_PD_metadata$Original_CID, "_"), `[[`, 5))

#  post-COBRA contig ids
VLP_contigs_PD_metadata$POST_CBR_CID <- VLP_contigs_PD_metadata$contig_id

colnames(VLP_contigs_PD_metadata)[grep('contig_id', colnames(VLP_contigs_PD_metadata))] <- "POST_CHV_name"

# simulating POST_CHV_length from POST_CBR_length
VLP_contigs_PD_metadata$POST_CHV_length <- VLP_contigs_PD_metadata$POST_CBR_length

VLP_contigs_PD_metadata$CHV_pruned <- "No"

VLP_contigs_PD_metadata[,c("region_types", "region_lengths", "region_coords_bp")] <- NA

# adding cohort name & removing box location (make sure to include it for samples themselves in the sample metadata)
# this line is very much dependant on the ID structure; make sure to adjust it

if ( length(grep('CHV', VLP_contigs_PD_metadata$Original_CID)) != 0 ) {
  
  # in-house negative control:
  VLP_contigs_PD_metadata$New_CID <- gsub("CHV(\\d{4})\\d*[_A-Z]*\\d*", "NEXT_N\\1", VLP_contigs_PD_metadata$Original_CID)
  
} else {
  
  # BaseClear sequencing negative control
  if ( length(grep('BCB', VLP_contigs_PD_metadata$Original_CID)) != 0 ) {
    
    VLP_contigs_PD_metadata$New_CID <- gsub("BCBPB(\\d{4})", "BSCL_N\\1", VLP_contigs_PD_metadata$Original_CID)
    
  } else {
    
    # BaseClear sequenucing positive control
    VLP_contigs_PD_metadata$New_CID <- gsub("BCPPB(\\d{4})", "BSCL_P\\1", VLP_contigs_PD_metadata$Original_CID)
    
    # manual magic:
    VLP_contigs_PD_metadata[grep("BCPPB36444", VLP_contigs_PD_metadata$Original_CID),]$New_CID <- gsub("BCPPB36444", "BSCL_P3644", VLP_contigs_PD_metadata[grep("BCPPB36444", VLP_contigs_PD_metadata$Original_CID),]$Original_CID)
    VLP_contigs_PD_metadata[grep("BCPPlate14", VLP_contigs_PD_metadata$Original_CID),]$New_CID <- gsub("BCPPlate14", "BSCL_P3614", VLP_contigs_PD_metadata[grep("BCPPlate14", VLP_contigs_PD_metadata$Original_CID),]$Original_CID) #randomly assigned by me
    
  }
  
}

# shorten NODE info
VLP_contigs_PD_metadata$New_CID <- gsub("ODE_", "", VLP_contigs_PD_metadata$New_CID)

# change the length in the name to the new length (postCheckV length)
VLP_contigs_PD_metadata$New_CID <- mapply(function(CID, postCHV_length) gsub("length_[0-9]+", paste0("L", postCHV_length), CID), 
                                          VLP_contigs_PD_metadata$New_CID, 
                                          VLP_contigs_PD_metadata$POST_CHV_length)

# add COBRA status: extended or untouched? (E0 or E1)
VLP_contigs_PD_metadata$COB_status <- "E0"
if (length(VLP_contigs_PD_metadata[grep('extended',VLP_contigs_PD_metadata$POST_CBR_CID),]$COB_status)!=0) {
  VLP_contigs_PD_metadata[grep('extended',VLP_contigs_PD_metadata$POST_CBR_CID),]$COB_status <- "E1"
} else {
  print("No contigs were extended")
}

# add CheckV status: pruned or untouched? (P0 or P1)
VLP_contigs_PD_metadata$CHV_status <- "P0"

# add Fragment number after CHV:
VLP_contigs_PD_metadata$fragment_N <- 0

VLP_contigs_PD_metadata$New_CID <- paste0(VLP_contigs_PD_metadata$New_CID, '_', 
                                          VLP_contigs_PD_metadata$COB_status, '_', 
                                          VLP_contigs_PD_metadata$CHV_status, '_',
                                          paste0('F', str_pad(VLP_contigs_PD_metadata$fragment_N, 1, pad = "0")))

VLP_contigs_PD_metadata[,c("CenoteTaker3", "DeepVirFinder", "geNomad", "VIBRANT", "VirSorter2")] <- NA

VLP_contigs_PD_metadata[,c("provirus", "proviral_length", 
                           "gene_count", "viral_genes", "host_genes",
                           "checkv_quality", "miuvig_quality", "completeness",
                           "completeness_method", "contamination", "kmer_freq",
                           "warnings")] <- NA

VLP_contigs_PD_metadata <- VLP_contigs_PD_metadata[,c("Original_CID", "Original_length", "CenoteTaker3",
                                                      "DeepVirFinder", "geNomad", "VIBRANT",
                                                      "VirSorter2", "POST_CBR_CID", "COB_status", "POST_CBR_length",
                                                      "CHV_pruned", "region_types", "region_lengths", 
                                                      "region_coords_bp", "POST_CHV_name", "POST_CHV_length", 
                                                      "New_CID", "provirus", "proviral_length", 
                                                      "gene_count", "viral_genes", "host_genes",
                                                      "checkv_quality", "miuvig_quality", "completeness",
                                                      "completeness_method", "contamination", "kmer_freq",
                                                      "warnings")]

write.table(VLP_contigs_PD_metadata, paste0(args[1], 'Extended_TOF'), sep='\t', row.names=F, col.names=T, quote=F)

