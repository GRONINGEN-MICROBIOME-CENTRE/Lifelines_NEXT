##########################################
# Designing new names and metadata 
# for extended and pruned NEXT virus
# contigs
##########################################
args = commandArgs(trailingOnly=TRUE)
##############################
# Loading libraries
##############################
library(stringr)
library(reshape2)
##############################
# Functions
##############################

##############################
# Input data
##############################

##############################
# ANALYSIS
##############################
prophage_pruning <- read.table(paste0(args[1], 'CheckV_pruning/','contamination.tsv'), sep='\t', header=T)

quality_summary <- read.table(paste0(args[1], 'CheckV_pruning/','quality_summary.tsv'), sep='\t', header=T)

virus_discovery <- read.table(paste0(args[1], args[2], '_table_of_origin'), sep='\t', header=F)

table_of_origin <- dcast(virus_discovery, V1~V2)
table_of_origin[-1] <- as.integer(table_of_origin[-1] != 0)
table_of_origin[is.na(table_of_origin)] <- 0
colnames(table_of_origin)[grep('V1', colnames(table_of_origin))] <- "Original_CID"

# what is the max number of fragments in pruned contigs?
if (sum(prophage_pruning$provirus=="Yes")!=0) {
  # extract the maximum fragment number from the pruned contigs (they have a specific matching pattern "_Nfrag_startcoord-endcoord_length", where coord are given for the prophage region)
  max_frag <- max(sapply(str_split(sub("_([0-9]+)-([0-9]+)_([0-9]+)$", "", quality_summary[grep('-', quality_summary$contig_id),]$contig_id), "_"), function(x) x[[length(x)]]))
  print(paste0("Max number of fragments: ", max_frag))
} else {
  print("No prophages identified") 
} 

##### new (provisional) contigs names: 

# Old name: CHV200013F12_NODE_28562_length_1013_cov_0.285553_1_266-1013_1013
# Changes:
# 0) add source (NEXT)
# 1) keep only VXXX from the SAMPLEID	(V2000)
# 3) merge NODE & its number (N28562)
# 4) merge (final after extension & pruning length?)	(L748)
# 5) COBRA status: extended or untouched? (E0)
# 6) CheckV: pruned or not?	(P1)
# 7) fragment_number: F0 if unpruned, FX if pruned, where X - number of the fragment (F1);
# the number of leading zeroes to be deicded when all samples will be run through CheckV, but
# for now I baldly decide that it cannot exceed 10 fragments
# New name: NEXT_V2000_N28562_L1013_cov_0.285553_E0_P1_F01	

##### Contigs (postdiscovery) metadata: 

VLP_contigs_PD_metadata <- quality_summary

# this line prunes the post COBRA & CheckV IDs (so all potential _extended and prophage coordinates characters are trimmed away)
# it is very dependant on the number of underscores in the contig id, so if the sample name contains underscores, change the number from 6 to smth else
VLP_contigs_PD_metadata$Original_CID <- sub("(([^_]*_){6}[^_]*).*", "\\1", VLP_contigs_PD_metadata$contig_id)

# here, the original contig length is extracted from its id, so it is also dependant on "_" in contig & sample id
VLP_contigs_PD_metadata$Original_length <- as.numeric(sapply(str_split(VLP_contigs_PD_metadata$Original_CID, "_"), `[[`, 5))

# trim away prophage pruning coordinates to extract post-COBRA contig ids
VLP_contigs_PD_metadata$POST_CBR_CID <- gsub("_[0-9]+$","", sub("_([0-9]+)-([0-9]+)_([0-9]+)$", "", VLP_contigs_PD_metadata$contig_id))

colnames(VLP_contigs_PD_metadata)[grep('contig_id', colnames(VLP_contigs_PD_metadata))] <- "POST_CHV_name"

colnames(VLP_contigs_PD_metadata)[grep('contig_length', colnames(VLP_contigs_PD_metadata))] <- "POST_CHV_length"

colnames(prophage_pruning)[grep('contig_id', colnames(prophage_pruning))] <- "POST_CBR_CID"

colnames(prophage_pruning)[grep('contig_length', colnames(prophage_pruning))] <- "POST_CBR_length"

colnames(prophage_pruning)[grep('provirus', colnames(prophage_pruning))] <- "CHV_pruned"

VLP_contigs_PD_metadata <- merge(VLP_contigs_PD_metadata, 
             prophage_pruning[,c("POST_CBR_CID", "POST_CBR_length", "CHV_pruned", "region_types", "region_lengths", "region_coords_bp")], 
             by='POST_CBR_CID')

# adding cohort name & removing box location (make sure to include it for samples themselves in the sample metadata)
# this line is very much dependant on the ID structure; make sure to adjust it
VLP_contigs_PD_metadata$New_CID <- gsub("CHV(\\d{4})\\d*[_A-Z]*\\d*", "NEXT_V\\1", VLP_contigs_PD_metadata$Original_CID)

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
if (sum(VLP_contigs_PD_metadata$CHV_pruned=="Yes")!=0) {
  VLP_contigs_PD_metadata[VLP_contigs_PD_metadata$CHV_pruned=="Yes",]$CHV_status <- "P1"
}

# add Fragment number after CHV:
VLP_contigs_PD_metadata$fragment_N <- 0
VLP_contigs_PD_metadata[VLP_contigs_PD_metadata$CHV_pruned=="Yes",]$fragment_N <- as.numeric( sapply (str_split(sub("_([0-9]+)-([0-9]+)_([0-9]+)$", "", VLP_contigs_PD_metadata[VLP_contigs_PD_metadata$CHV_pruned=="Yes",]$POST_CHV_name), "_"), function(x) x[[length(x)]]))

VLP_contigs_PD_metadata$New_CID <- paste0(VLP_contigs_PD_metadata$New_CID, '_', 
                                          VLP_contigs_PD_metadata$COB_status, '_', 
                                          VLP_contigs_PD_metadata$CHV_status, '_',
                                          paste0('F', str_pad(VLP_contigs_PD_metadata$fragment_N, 1, pad = "0")))

VLP_contigs_PD_metadata <- merge(VLP_contigs_PD_metadata, table_of_origin, by="Original_CID")

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

##############################
# OUTPUT
##############################
write.table(VLP_contigs_PD_metadata, paste0(args[1], 'Extended_TOF'), sep='\t', row.names=F, col.names=T, quote=F)
