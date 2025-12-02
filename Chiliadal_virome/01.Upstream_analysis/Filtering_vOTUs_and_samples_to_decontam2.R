##########################################
# Preparing files for strain-aware decon-
# tamination
##########################################

##############################
# Loading libraries
##############################
library(readr)
##############################
# Functions
##############################

##############################
# INPUT
##############################
# LIST OF CONTROL SAMPLES:
controls <- read.table("/scratch/p282752/ANALYSIS_CHILIADAL/controls.list", sep='\t', header=F)
dim(controls)
rpkm <- read.delim("../VIR_DB/mapping/VLP_MGS_overlap/RPKM_counts_VLP.txt")
print(dim(rpkm))

genomes_stb <- data.frame(row.names(rpkm), row.names(rpkm)) 
##############################
# ANALYSIS
##############################
# KEEPING ONLY vOTUs SHARED TO CONTROLS:
rpkm_CTRL <- rpkm[, colnames(rpkm) %in% controls$V1]
rpkm_CTRL <- rpkm_CTRL[rowSums(rpkm_CTRL) >0,]
print(dim(rpkm_CTRL))

# KEEPING ONLY vOTUs PRESENT IN CONTROLS
control_vOTUs <- row.names(rpkm_CTRL[rowSums(rpkm_CTRL) > 0,])
rpkm <- rpkm[row.names(rpkm) %in% control_vOTUs,]

# KEEPING ONLY vOTUs PRESENT IN > 1 CONTROLS OR SAMPLES
rpkm <- rpkm[rowSums(rpkm>0) > 1,]
rpkm <- rpkm[,colSums(rpkm) >0]

print(dim(rpkm))

##############################
# OUTPUT
##############################
for (i in colnames(rpkm)) {
  votus_keep <- row.names(rpkm[rpkm[,i] > 0,])
  
  write.table(votus_keep, paste0("../VIR_DB/VLP_MGS_decontamination/keep_vOTUs_per_sample/", i, "_vOTUs_to_keep"), sep='\t', row.names=F, col.names=F, quote=F)
}

write.table(colnames(rpkm), "../VIR_DB/VLP_MGS_decontamination/SAMPLES_and_NCs_run_decon", sep = "\t", row.names = F, col.names = F, quote=F)

write.table(row.names(rpkm), "../VIR_DB/VLP_MGS_decontamination/vOTUs_shared_by_samples_and_ncs", sep = "\t", row.names = F, col.names = F, quote=F)

write.table(genomes_stb, "../VIR_DB/VLP_MGS_decontamination/genomes.stb", sep = "\t", row.names = F, col.names = F, quote=F)
