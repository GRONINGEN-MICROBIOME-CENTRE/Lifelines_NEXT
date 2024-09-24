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
sample_metadata = read.table("../METADATA/Chiliadal_metadata_ver_04_03042024.txt",header=T,sep="\t", quote = "")
print(dim(sample_metadata))

# CHANGING FROM NEW IDs TO VNP:
prechili <- read.table('../METADATA/BC_PILOT_new_ids_chiliadal_format.txt', sep='\t')
prechili <- prechili[,c('V1', 'V2')]

sample_metadata[grep('CHV3', sample_metadata$Sequencing_ID_VLP),]$Sequencing_ID_VLP <- prechili$V1[match(sample_metadata[grep('CHV3', sample_metadata$Sequencing_ID_VLP),]$Sequencing_ID_VLP, prechili$V2)]

# LIST OF CONTROL SAMPLES:
controls = sample_metadata[sample_metadata$Type_VLP %in% c("BLANK", "POSITIVE"), "Sequencing_ID_VLP"]

rpkm <- read.delim("../VIR_DB/mapping/VLP_N_vir_genes/RPKM_counts_VLP.txt")
print(dim(rpkm))

##############################
# ANALYSIS
##############################
# KEEPING ONLY vOTUs SHARED TO CONTROLS:
rpkm_CTRL <- rpkm[, colnames(rpkm) %in% controls]
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
  
  write.table(votus_keep, paste0("../VIR_DB/decontamination/keep_vOTUs_per_sample/", i, "_vOTUs_to_keep"), sep='\t', row.names=F, col.names=F, quote=F)
}

write.table(colnames(rpkm), "../VIR_DB/decontamination/SAMPLES_and_NCs_run_decon", sep = "\t", row.names = F, col.names = F, quote=F)

write.table(row.names(rpkm), "../VIR_DB/decontamination/vOTUs_shared_by_samples_and_ncs", sep = "\t", row.names = F, col.names = F, quote=F)

