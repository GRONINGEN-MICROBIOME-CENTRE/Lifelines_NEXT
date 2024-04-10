##########################################
# Splitting Extended_TOF to increase
# speed of renaming of contigs from 
# in-house negative controls
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

##############################
# ANALYSIS
##############################
Extended_TOF <- read.table('../Extended_TOF', sep='\t', header=T)

frag_list <- read.table(paste0(args[1], '_contigs'), sep='\t', header=F)

Extended_TOF <- Extended_TOF[Extended_TOF$POST_CBR_CID %in% frag_list$V1,]


##############################
# OUTPUT
##############################
write.table(Extended_TOF, paste0(args[1], '_Extended_TOF'), sep='\t', row.names=F, col.names=T, quote=F)
