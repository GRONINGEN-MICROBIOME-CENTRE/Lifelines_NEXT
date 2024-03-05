##########################################
# Combining per sample table of origin  
# and converting them from long to wide 
# format
##########################################
args = commandArgs(trailingOnly=TRUE)
##############################
# Loading libraries
##############################
library(reshape2)
library(dplyr)

##############################
# Functions
##############################

##############################
# Input data
##############################

##############################
# ANALYSIS
##############################
temp = list.files(path=args[1], pattern=args[2], full.names=T)

myfiles = lapply(temp, read.table)

# coverting table of origin from long to wide format
import_data=lapply(myfiles, function(x) {
                   out <- list(dcast(x, V1~V2))
                   return(out)
                   })

# concatenating per sample tables of origin
table_of_origin_redundant<-bind_rows(import_data)

# converting tool names and NAs to boolean 
table_of_origin_redundant[-1] <- as.integer(table_of_origin_redundant[-1] != 0) 
table_of_origin_redundant[is.na(table_of_origin_redundant)] <- 0

##############################
# OUTPUT
##############################
write.table(table_of_origin_redundant, paste0(args[1], 'table_of_origin'), sep='\t', quote=F, row.names=F)
