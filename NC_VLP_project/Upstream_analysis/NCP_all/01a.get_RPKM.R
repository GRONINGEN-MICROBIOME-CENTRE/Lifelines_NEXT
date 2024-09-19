##############################
# Loading libraries
##############################
.libPaths("/home1/p309176/R_libraries")
library(readr)
##############################
# Functions
##############################
RPKM_table <- function(counts, coverage, contigs_metadata){
  counts <- counts[order( row.names(counts) ),]
  
  coverage <- coverage[order( row.names(coverage) ),]
  
  if (identical(colnames(counts), colnames(coverage))) {
    
    RPM <- as.data.frame( t(t(counts) / (colSums(counts)/1000000)) )
    
    RPM$length <- contigs_metadata$POST_CHV_length[match(row.names(RPM), contigs_metadata$New_CID)]
    
    RPKM <- RPM / (RPM$length / 1000) 
    
    RPKM$length <- NULL
    
    RPKM[coverage<=0.75] <-0
  } else {
    print("Counts and coverage tables have different sets of samples")
  }
  
  return(RPKM)
}

##############################
# Input data
##############################

setwd('/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R')

VLP_counts <- as.data.frame(read_tsv('../VIR_DB/mapping/VLP_to_w_neg_der95_NCP/coverage_NCP/count_table.txt'))
names(VLP_counts)[1] <- "vOTUs"
row.names(VLP_counts) <- VLP_counts$vOTUs
VLP_counts$vOTUs <- NULL

VLP_contig_coverage <- as.data.frame(read_tsv('../VIR_DB/mapping/VLP_to_w_neg_der95_NCP/coverage_NCP/coverage_table.txt'))
names(VLP_contig_coverage)[1] <- "vOTUs"
row.names(VLP_contig_coverage) <- VLP_contig_coverage$vOTUs
VLP_contig_coverage$vOTUs <- NULL

VLP_contigs_metadata <- read.delim('../VIR_DB/virus_contigs/MERGED_Extended_TOF_NCP')

VLP_contigs_metadata <- VLP_contigs_metadata[VLP_contigs_metadata$New_CID %in% row.names(VLP_counts), ]

RPKM_counts_VLP <- RPKM_table(VLP_counts, VLP_contig_coverage, VLP_contigs_metadata)

RPKM_counts_VLP_filt <- RPKM_counts_VLP[rowSums(RPKM_counts_VLP) > 0, colSums(RPKM_counts_VLP) > 0]

##############################
# OUTPUT
##############################
write.table(RPKM_counts_VLP_filt, "RPKM_counts_VLP_NCP.txt", sep='\t')
