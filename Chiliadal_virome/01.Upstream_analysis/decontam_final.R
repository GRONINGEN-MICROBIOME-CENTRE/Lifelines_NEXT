##########################################
# cleaning RPKM table
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
RPKM <- read.table("/scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/mapping/VLP_MGS_overlap/RPKM_counts_VLP.txt", sep='\t', header=T)
dim(RPKM)

DF <- read.table("/scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/VLP_MGS_decontamination/instrain_compare_ALL/output/instrain_compare_ALL_comparisonsTable.tsv", sep='\t', header=T)
dim(DF)

controls <- read.table("/scratch/p282752/ANALYSIS_CHILIADAL/controls.list", sep='\t', header=F)
##############################
# ANALYSIS
##############################
DF$name1 <- gsub("_filtered.sorted.bam", "", DF$name1)
DF$name2 <- gsub("_filtered.sorted.bam", "", DF$name2)

swap_indices <- DF$name1 %in% controls$V1 & !DF$name2 %in% controls$V1

DF[swap_indices, c("name1", "name2")] <- DF[swap_indices, c("name2", "name1")]

DF <- DF[DF$name2 %in% controls$V1 & !DF$name1 %in% controls$V1, ]

DF <- DF[DF$popANI>=0.99 ,]

clean_RPKM <- as.matrix(RPKM)

row_idx <- match(DF$scaffold, rownames(clean_RPKM))
col_idx <- match(DF$name1,   colnames(clean_RPKM))

keep <- !is.na(row_idx) & !is.na(col_idx)

clean_RPKM[cbind(row_idx[keep], col_idx[keep])] <- 0

clean_RPKM <- as.data.frame(clean_RPKM, check.names = FALSE)

##############################
# OUTPUT
##############################
write.table(clean_RPKM, "/scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/VLP_MGS_decontamination/clean_RPKM_UPD.txt", sep='\t', quote=F)
