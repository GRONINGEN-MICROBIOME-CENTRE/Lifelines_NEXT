##########################################
# creating RPKM table
##########################################
args = commandArgs(trailingOnly=TRUE)
##############################
# Loading libraries
##############################
##############################
# Functions
##############################
# RPKM_table function requires three dfs of the certain format: 
# counts: N of reads aligned to the sequence (contig)
# coverage: N of non-zero bases of the sequence length
# Rows are contigs and columns are samples.
# contigs_metadta: to extract the length of contigs
# RPKM_table returns a table of RPKM counts. 
# Be aware that obtained RPKM counts will be representing counts 
# normalized on the number of ALIGNED reads (to the database), 
# not the initial size of quality-trimmed library. 

RPKM_table <- function(counts, coverage, contigs_metadata){
  counts <- counts[order( row.names(counts) ),]
  
  coverage <- coverage[order( row.names(coverage) ),]
  
  if (identical(colnames(counts), colnames(coverage))) {
    
    RPM <- as.data.frame( t(t(counts) / (colSums(counts)/1000000)) )
    
    RPM$length <- contigs_metadata$V3[match(row.names(RPM), contigs_metadata$V1)]
    
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
overlapping <- read.table("/scratch/p282752/ANALYSIS_CHILIADAL/MGS_VLP_samples_full_overlap.txt", sep='\t', header=T)
dim(overlapping)
controls <- read.table("/scratch/p282752/ANALYSIS_CHILIADAL/controls.list", sep='\t', header=F)

print(colnames(controls))

samples_to_select <- c(overlapping$VLP, overlapping$MGS, controls$V1)
length(samples_to_select)

VLP_counts <- read.table(paste0(args[1], '/count_table.txt'), sep='\t', header=T, row.names=1)
dim(VLP_counts)
VLP_counts <- VLP_counts[,colnames(VLP_counts) %in% samples_to_select]
dim(VLP_counts)

VLP_contig_coverage <- read.table(paste0(args[1], '/coverage_table.txt'), sep='\t', header=T, row.names=1)
dim(VLP_contig_coverage)
VLP_contig_coverage <- VLP_contig_coverage[, colnames(VLP_contig_coverage) %in% samples_to_select]
dim(VLP_contig_coverage)

VLP_contigs_metadata <- read.table(paste0(args[2]), sep='\t', header=F)
dim(VLP_contigs_metadata)

VLP_contigs_metadata <- VLP_contigs_metadata[VLP_contigs_metadata$V1 %in% row.names(VLP_counts),]
VLP_contigs_metadata$V3 <- as.numeric(VLP_contigs_metadata$V3)
##############################
# ANALYSIS
##############################
RPKM_counts_VLP <- RPKM_table(VLP_counts, VLP_contig_coverage, VLP_contigs_metadata)
RPKM_counts_VLP <- RPKM_counts_VLP[rowSums(RPKM_counts_VLP)>0,]
##############################
# OUTPUT
##############################
write.table(RPKM_counts_VLP, paste0(args[1], '/RPKM_counts_VLP.txt'), sep='\t', quote=F)
