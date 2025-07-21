##########################################
# Designing new names and metadata 
# for pruned DB virus contigs
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
POST_CBR_LENGTH <- read.table(paste0(args[1], '/Prophage_pruning/','POST_CBR_LENGTH'), sep='\t', col.names = c('POST_CBR_CID', "POST_CBR_length"))

plasmids <- read.table(paste0(args[1], '/Prophage_pruning/', args[2], '_plasmid_summary.tsv'), sep='\t', header=T)

prophage_prepruning <- read.table(paste0(args[1], '/Prophage_pruning/', args[2], '_virus_summary.tsv'), sep='\t', header=T)

prophage_pruning <- read.table(paste0(args[1], '/Prophage_pruning/','contamination.tsv'), sep='\t', header=T)

quality_summary <- read.table(paste0(args[1], '/Prophage_pruning/', 'quality_summary.tsv'), sep='\t', header=T)

##############################
# ANALYSIS
##############################

##### new contigs IDs: 

# Old names:
# Guerin_CrAss_cdzh01002743_1
# Yutin_CrAss_phage_2
# NCBI_CrAss_OQ383781.1
# NL_crAss000001
# Benler_CYGL01000011.1
# COPSAC_OTU_1062
# VREF_NC_029549.1
# GVD_Broecker_Sample_ID29-1191-0_NODE_1524_length_8276_cov_8.926043
# GPD_ivig_1
# MGV-GENOME-0364295
# IMGVR_UViG_2504643025_000001|2504643025|2504645545|381-5974

# Changes:
# 1) add DB source
# 2) keep the SAMPLE_ID if available
# 3) Simulate node: N0
# 4) Add new length (L7564)
# 5) Simulate cov: K0.0
# 6) COBRA status: E0
# 7) pruning status: pruned or not pruned? (P1); two sources: geNomad and CheckV
# 8) fragment_number: F0 if unpruned, FX if pruned, where X - number of the fragment (F1);
# the number of leading zeroes is empirically deduced from previous runs
# New name: Guerin_CrAss_cs_ms_24_1_L91105_K0_K0.0_E0_P1_F1

##### Contigs (postdiscovery) metadata: 

# Quality of final virus contigs:
VLP_contigs_PD_metadata <- quality_summary

# this line prunes the post geNomad & CheckV IDs (so all potential _extended and prophage coordinates characters are trimmed away)
# it is very dependant on the number of underscores in the contig id, so if the sample name contains underscores, change the number from 6 to smth else
VLP_contigs_PD_metadata$Original_CID <- sub("_[1-9]_[0-9]+-[0-9]+_[0-9]+$", "\\1", sub("\\|.*", "", VLP_contigs_PD_metadata$contig_id))

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
if ( sum(is.na(VLP_contigs_PD_metadata$GND_pruned))!=0 ) {
  VLP_contigs_PD_metadata[is.na(VLP_contigs_PD_metadata$GND_pruned), ]$GND_pruned <- "No"
}

VLP_contigs_PD_metadata$POST_CBR_CID <- sub("\\|.*", "", VLP_contigs_PD_metadata$POST_GND_CID)

VLP_contigs_PD_metadata <- merge(VLP_contigs_PD_metadata, POST_CBR_LENGTH, by='POST_CBR_CID', all=T)

# adding cohort name & removing box location (make sure to include it for samples themselves in the sample metadata)
# this line is very much dependant on the ID structure; make sure to adjust it
VLP_contigs_PD_metadata$New_CID <- VLP_contigs_PD_metadata$Original_CID

# shorten NODE info
if (length(grep('_NODE_', VLP_contigs_PD_metadata$New_CID)) == nrow(VLP_contigs_PD_metadata)) {
  VLP_contigs_PD_metadata$New_CID <- gsub('NODE_', 'N', VLP_contigs_PD_metadata$New_CID)
} else {
  VLP_contigs_PD_metadata$New_CID <- paste0(VLP_contigs_PD_metadata$New_CID, '_N0')
}

# change the length in the name to the new length (postCheckV length)
if (length(grep('_length_', VLP_contigs_PD_metadata$New_CID)) == nrow(VLP_contigs_PD_metadata)) {
  VLP_contigs_PD_metadata$New_CID <- mapply(function(CID, postCHV_length) gsub("length_[0-9]+", paste0("L", postCHV_length), CID), 
                                            VLP_contigs_PD_metadata$New_CID, 
                                            VLP_contigs_PD_metadata$POST_CHV_length)
} else {
  VLP_contigs_PD_metadata$New_CID <- paste0(VLP_contigs_PD_metadata$New_CID,"_L", VLP_contigs_PD_metadata$POST_CHV_length)
}

#### rounding the k-mer coverage:


if (length(grep('_cov_', VLP_contigs_PD_metadata$New_CID)) == nrow(VLP_contigs_PD_metadata)) {
  VLP_contigs_PD_metadata$kmer_cov <- round(as.numeric(gsub(".*_", "", VLP_contigs_PD_metadata$Original_CID)), 1)
  VLP_contigs_PD_metadata$New_CID <- mapply(function(CID, kmer_cov) gsub("cov_[0-9]+.*", paste0("K", kmer_cov), CID), 
                                            VLP_contigs_PD_metadata$New_CID, 
                                            VLP_contigs_PD_metadata$kmer_cov)
} else {
  VLP_contigs_PD_metadata$New_CID <- paste0(VLP_contigs_PD_metadata$New_CID, '_K0.0')
}

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
VLP_contigs_PD_metadata[VLP_contigs_PD_metadata$CHV_pruned=='Yes',]$CHV_start <- as.numeric(sub(".*_[0-9]+_([0-9]+)-[0-9]+.*$", "\\1", VLP_contigs_PD_metadata[VLP_contigs_PD_metadata$CHV_pruned=='Yes',]$POST_CHV_CID))

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

# add table of origin
VLP_contigs_PD_metadata[,c("CenoteTaker3", "DeepVirFinder", 
                           "geNomad", "VIBRANT", "VirSorter2")] <- NA

VLP_contigs_PD_metadata$Original_length <- VLP_contigs_PD_metadata$POST_CBR_length

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
write.table(VLP_contigs_PD_metadata, paste0(args[1], '/Prophage_pruning/', 'Extended_TOF'), sep='\t', row.names=F, col.names=T, quote=F)