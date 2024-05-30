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
POST_CBR_LENGTH <- read.table(paste0(args[1], 'Prophage_pruning/','POST_CBR_LENGTH'), sep='\t', col.names = c('POST_CBR_CID', "POST_CBR_length"))

plasmids <- read.table(paste0(args[1], 'Prophage_pruning/', args[2], '_extended_viral_plasmid_summary.tsv'), sep='\t', header=T)

prophage_prepruning <- read.table(paste0(args[1], 'Prophage_pruning/', args[2], '_extended_viral_virus_summary.tsv'), sep='\t', header=T)

prophage_pruning <- read.table(paste0(args[1], 'Prophage_pruning/','contamination.tsv'), sep='\t', header=T)

quality_summary <- read.table(paste0(args[1], 'Prophage_pruning/', 'quality_summary.tsv'), sep='\t', header=T)

virus_discovery <- read.table(paste0(args[1], args[2], '_table_of_origin'), sep='\t', header=F)

# transforming virus_discovery stat into wide format
table_of_origin <- dcast(virus_discovery, V1~V2)
table_of_origin[-1] <- as.integer(table_of_origin[-1] != 0)
table_of_origin[is.na(table_of_origin)] <- 0
colnames(table_of_origin)[grep('V1', colnames(table_of_origin))] <- "Original_CID"
##############################
# ANALYSIS
##############################

##### new contigs IDs: 

# Old name: CHV002101E03_NODE_183_length_31652_cov_6.159651_extended_partial|provirus_122336_143812

# Changes:
# 1) add source cohort (NEXT)
# 2) keep only VXXX from the SAMPLE_ID	(V0021)
# 3) shorten NODE_ to N and merge the node & its number (N183)
# 4) shorten length_ to L and merge it with the final length of the fragment after extension & pruning (L30294)
# 5) shorten cov_x.xxxx to K and merge with its number rounded to 1 digit after .
# 6) COBRA status: extended or untouched? (E1)
# 7) pruning status: pruned or not pruned? (P1); two sources: geNomad and CheckV
# 8) fragment_number: F0 if unpruned, FX if pruned, where X - number of the fragment (F1);
# the number of leading zeroes is empirically deduced from previous runs
# New name: NEXT_V0021_N183_L21477_K6.2_E1_P1_F1

##### Contigs (postdiscovery) metadata: 

# Quality of final virus contigs:
VLP_contigs_PD_metadata <- quality_summary

# this line prunes the post COBRA & geNomad & CheckV IDs (so all potential _extended and prophage coordinates characters are trimmed away)
# it is very dependant on the number of underscores in the contig id, so if the sample name contains underscores, change the number from 6 to smth else
VLP_contigs_PD_metadata$Original_CID <- sub("(([^_]*_){6}[^_]*).*", "\\1", sub("\\|.*", "", VLP_contigs_PD_metadata$contig_id))

# here, the original contig length is extracted from its id, so it is also dependant on "_" in contig & sample id
VLP_contigs_PD_metadata$Original_length <- as.numeric(sapply(str_split(VLP_contigs_PD_metadata$Original_CID, "_"), `[[`, 5))

# trim away CheckV prophage pruning coordinates to extract post-geNomad contig ids
VLP_contigs_PD_metadata$POST_GND_CID <- sub("_[1-9]_[0-9]+-[0-9]+_[0-9]+$", "", VLP_contigs_PD_metadata$contig_id)

colnames(VLP_contigs_PD_metadata)[grep('contig_id', colnames(VLP_contigs_PD_metadata))] <- "POST_CHV_CID"

colnames(VLP_contigs_PD_metadata)[grep('contig_length', colnames(VLP_contigs_PD_metadata))] <- "POST_CHV_length"

# CheckV pruning information of post-geNomad contigs:
colnames(prophage_pruning)[grep('contig_id', colnames(prophage_pruning))] <- "POST_GND_CID"

colnames(prophage_pruning)[grep('contig_length', colnames(prophage_pruning))] <- "POST_GND_length"

colnames(prophage_pruning)[grep('provirus', colnames(prophage_pruning))] <- "CHV_pruned"

# Combining CheckV quality assessment and CheckV prophage pruning:
VLP_contigs_PD_metadata <- merge(VLP_contigs_PD_metadata, 
                                 prophage_pruning[,c("POST_GND_CID", "POST_GND_length", "CHV_pruned", 
                                                     "region_types", "region_lengths", "region_coords_bp")], 
                                                  by='POST_GND_CID', all=T)

# Since I decided that sometimes CheckV overtrims some virus contigs and did not run 
# contigs containing DTR and ITRs identified by geNomad through CheckV pruning:
VLP_contigs_PD_metadata[is.na(VLP_contigs_PD_metadata$POST_GND_length),"CHV_pruned"] <- "Untouched"
VLP_contigs_PD_metadata[is.na(VLP_contigs_PD_metadata$POST_GND_length),"POST_GND_length"] <- VLP_contigs_PD_metadata[is.na(VLP_contigs_PD_metadata$POST_GND_length),"POST_CHV_length"]

# geNomad pruning information:
prophage_prepruning$GND_pruned <- "No"

prophage_prepruning[grep('\\|provirus', prophage_prepruning$seq_name),"GND_pruned"] <- "Yes"

colnames(prophage_prepruning)[grep('seq_name', colnames(prophage_prepruning))] <- "POST_GND_CID"

colnames(prophage_prepruning)[grep('topology', colnames(prophage_prepruning))] <- "GND_topology"

colnames(prophage_prepruning)[grep('coordinates', colnames(prophage_prepruning))] <- "GND_coordinates"

# Combining CheckV quality assessment and geNomad prophage prepruning:
VLP_contigs_PD_metadata <- merge(VLP_contigs_PD_metadata, 
                                 prophage_prepruning[,c("POST_GND_CID", "GND_topology", "GND_coordinates", "taxonomy", "GND_pruned")], 
                                 by="POST_GND_CID", all=T)

# since not all contigs were recognized by geNomad as viral:
VLP_contigs_PD_metadata[is.na(VLP_contigs_PD_metadata$GND_pruned), ]$GND_pruned <- "No"

VLP_contigs_PD_metadata$POST_CBR_CID <- sub("\\|.*", "", VLP_contigs_PD_metadata$POST_GND_CID)   

VLP_contigs_PD_metadata <- merge(VLP_contigs_PD_metadata, POST_CBR_LENGTH, by='POST_CBR_CID', all=T)

# adding cohort name & removing box location (make sure to include it for samples themselves in the sample metadata)
# this line is very much dependant on the ID structure; make sure to adjust it
VLP_contigs_PD_metadata$New_CID <- gsub("CHV(\\d{4})\\d*[_A-Z]*\\d*", "NEXT_V\\1", VLP_contigs_PD_metadata$Original_CID)

# shorten NODE info
VLP_contigs_PD_metadata$New_CID <- gsub("ODE_", "", VLP_contigs_PD_metadata$New_CID)

# change the length in the name to the new length (postCheckV length)
VLP_contigs_PD_metadata$New_CID <- mapply(function(CID, postCHV_length) gsub("length_[0-9]+", paste0("L", postCHV_length), CID), 
                                          VLP_contigs_PD_metadata$New_CID, 
                                          VLP_contigs_PD_metadata$POST_CHV_length)

#### rounding the k-mer coverage:
VLP_contigs_PD_metadata$kmer_cov <- round(as.numeric(gsub(".*_", "", VLP_contigs_PD_metadata$Original_CID)), 1)

VLP_contigs_PD_metadata$New_CID <- mapply(function(CID, kmer_cov) gsub("cov_[0-9]+.*", paste0("K", kmer_cov), CID), 
                                          VLP_contigs_PD_metadata$New_CID, 
                                          VLP_contigs_PD_metadata$kmer_cov)

print(paste0("N unique contigs: ", length(VLP_contigs_PD_metadata$POST_CHV_CID)))
print(paste0("N unique New_CIDs: ", length(VLP_contigs_PD_metadata$New_CID)))

# add COBRA (extension) status: extended or untouched? (E0 or E1)
VLP_contigs_PD_metadata$COB_status <- "E0"

if (length(VLP_contigs_PD_metadata[grep('extended',VLP_contigs_PD_metadata$POST_CBR_CID),]$COB_status)!=0) {
  
  VLP_contigs_PD_metadata[grep('extended',VLP_contigs_PD_metadata$POST_CBR_CID),]$COB_status <- "E1"
  
} else {
  print("No contigs were extended")
}

# add pruning status: pruned or not pruned? (P0 or P1)
# combines both pre-pruning and pruning
VLP_contigs_PD_metadata$PRU_status <- "P0"

if (sum(VLP_contigs_PD_metadata$CHV_pruned=="Yes")!=0 | sum(VLP_contigs_PD_metadata$GND_pruned=="Yes")!=0) {
  
  VLP_contigs_PD_metadata[(VLP_contigs_PD_metadata$GND_pruned=="Yes" | VLP_contigs_PD_metadata$CHV_pruned=="Yes"),]$PRU_status <- "P1"

  } else {
  print("No contigs were pruned")
}

# Calculate the fragment number

# 1. Get the start coordinates of geNomad and CheckV pruning:
VLP_contigs_PD_metadata$GND_start <- as.numeric(gsub('-.*', '', VLP_contigs_PD_metadata$GND_coordinates))
VLP_contigs_PD_metadata$CHV_start <- NA
VLP_contigs_PD_metadata[VLP_contigs_PD_metadata$CHV_pruned=='Yes',]$CHV_start <- as.numeric(sub(".*_([0-9]+)-[0-9]+.*$", "\\1", VLP_contigs_PD_metadata[VLP_contigs_PD_metadata$CHV_pruned=='Yes',]$POST_CHV_CID))

# 2. Get the consecutive number of the trimmed sequence:

# sorting by Original contig id, then by the start of geNomad, and then by CheckV start coordinate
VLP_contigs_PD_metadata <- VLP_contigs_PD_metadata[order(VLP_contigs_PD_metadata$Original_CID, 
                                                         VLP_contigs_PD_metadata$GND_start,
                                                         VLP_contigs_PD_metadata$CHV_start),]

VLP_contigs_PD_metadata$fragment_N <- as.numeric(ave(VLP_contigs_PD_metadata$Original_CID, 
                                                     VLP_contigs_PD_metadata$Original_CID, 
                                                     FUN = seq_along))

VLP_contigs_PD_metadata[VLP_contigs_PD_metadata$PRU_status=="P0","fragment_N"] <- 0


VLP_contigs_PD_metadata$New_CID <- paste0(VLP_contigs_PD_metadata$New_CID, '_', 
                                          VLP_contigs_PD_metadata$COB_status, '_', 
                                          VLP_contigs_PD_metadata$PRU_status, '_',
                                          paste0('F', str_pad(VLP_contigs_PD_metadata$fragment_N, 1, pad = "0")))

VLP_contigs_PD_metadata <- merge(VLP_contigs_PD_metadata, table_of_origin, by="Original_CID")

VLP_contigs_PD_metadata <- VLP_contigs_PD_metadata[,c("New_CID", "provirus", "POST_CHV_length", 
                                                      "proviral_length", "gene_count", "viral_genes", 
                                                      "host_genes", "checkv_quality", "miuvig_quality", 
                                                      "completeness", "completeness_method", "contamination", 
                                                      "kmer_freq", "warnings", "taxonomy", "Original_CID", 
                                                      "Original_length", "CenoteTaker3", "DeepVirFinder", 
                                                      "geNomad", "VIBRANT", "VirSorter2", "POST_CBR_CID", 
                                                      "COB_status", "POST_CBR_length", "POST_GND_CID", 
                                                      "GND_pruned", "POST_GND_length", "GND_topology", "GND_coordinates",
                                                      "CHV_pruned", "region_types", "region_lengths", 
                                                      "region_coords_bp", "POST_CHV_CID", 
                                                      "PRU_status", "fragment_N")]

# getting plasmid info:
VLP_contigs_PD_metadata$plasmid <- "No"

if (dim(plasmids)[1] != 0) {
  
  VLP_contigs_PD_metadata[VLP_contigs_PD_metadata$POST_GND_CID %in% plasmids$seq_name,]$plasmid <- "Yes"
  VLP_contigs_PD_metadata[VLP_contigs_PD_metadata$POST_GND_CID %in% plasmids$seq_name,]$GND_topology <- plasmids$topology[match(VLP_contigs_PD_metadata[VLP_contigs_PD_metadata$POST_GND_CID %in% plasmids$seq_name,]$POST_GND_CID, plasmids$seq_name)]
  VLP_contigs_PD_metadata$conjugation_genes <- plasmids$conjugation_genes[match(VLP_contigs_PD_metadata$POST_GND_CID, plasmids$seq_name)]
  VLP_contigs_PD_metadata$amr_genes <- plasmids$amr_genes[match(VLP_contigs_PD_metadata$POST_GND_CID, plasmids$seq_name)]
  VLP_contigs_PD_metadata$plasmid_score <- plasmids$plasmid_score[match(VLP_contigs_PD_metadata$POST_GND_CID, plasmids$seq_name)]
  VLP_contigs_PD_metadata$plasmid_fdr <- plasmids$fdr[match(VLP_contigs_PD_metadata$POST_GND_CID, plasmids$seq_name)]
} else {
  
  VLP_contigs_PD_metadata$conjugation_genes <- NA
  VLP_contigs_PD_metadata$amr_genes <- NA
  VLP_contigs_PD_metadata$plasmid_score <- NA
  VLP_contigs_PD_metadata$plasmid_fdr <- NA
  
}
##############################
# OUTPUT
##############################
write.table(VLP_contigs_PD_metadata, paste0(args[1], 'Prophage_pruning/', 'Extended_TOF'), sep='\t', row.names=F, col.names=T, quote=F)
