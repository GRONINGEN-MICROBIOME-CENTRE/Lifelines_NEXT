library(reshape2)
library(dplyr)

temp = list.files(path="/scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/clean_viral/table_of_origin",pattern="*_table_of_origin",full.names=T)

myfiles = lapply(temp, read.table)

import_data=lapply(myfiles, function(x) {
                   out <- list(dcast(x, V1~V2))
                   return(out)
                   })
table_of_origin_redundant<-bind_rows(import_data)
table_of_origin_redundant[-1] <- as.integer(table_of_origin_redundant[-1] != 0) 

table_of_origin_redundant[is.na(table_of_origin_redundant)] <- 0

all_viral_noneg405_99_der95  <- read.table('../clean_viral/dereplication/all_viral_noneg405_99_der95_IDs', sep='\t')

table_of_origin <- merge(table_of_origin_redundant, all_viral_noneg405_99_der95, by="V1")

pVOGs_stat <- read.table('../clean_viral/table_of_origin/pVOG_stat_all', sep='\t', header=F)
colnames(pVOGs_stat)[4] <- "N_pVOGs_per_10kb"

table_of_origin_complete <-merge(table_of_origin, pVOGs_stat[,c(1,4)], by="V1", all.x=T)

all_pvogs_detected <- read.table('../clean_viral/table_of_origin/all_pVOGs_detected', sep='', header=F)
all_pvogs_detected$V1 <- sub("_[^_]+$", "", all_pvogs_detected$V1)
list_temperate_vogs <- read.table('/data/umcg-sgarmaeva/pvogs/temperate/list_VOGs')
list_temperate_vogs_function <- read.table('/data/umcg-sgarmaeva/pvogs/temperate/temperate_pVOGs_uniq.txt', sep='', header=T)
list_repressor_vogs <- read.table('/data/umcg-sgarmaeva/pvogs/temperate/list_VOGs_repressor_halmarks.txt')

temperate_vogs_in_data <- all_pvogs_detected[all_pvogs_detected$V3 %in% list_temperate_vogs$V1,]

temperate_vogs_in_data$integrase <- NA
temperate_vogs_in_data[temperate_vogs_in_data$V3 %in% list_temperate_vogs_function[grep('integrase', list_temperate_vogs_function$type),]$pVOG,]$integrase <- 1
temperate_vogs_in_data[is.na(temperate_vogs_in_data$integrase),]$integrase <- 0

temperate_vogs_in_data$recombinase <- NA
temperate_vogs_in_data[temperate_vogs_in_data$V3 %in% list_temperate_vogs_function[list_temperate_vogs_function$type=='SSR',]$pVOG,]$recombinase <- 1
temperate_vogs_in_data[is.na(temperate_vogs_in_data$recombinase),]$recombinase <- 0

temperate_bacteriophages_integrase <- aggregate(integrase ~ V1, temperate_vogs_in_data, sum)
temperate_bacteriophages_recombinase <- aggregate(recombinase ~ V1, temperate_vogs_in_data, sum)
temperate_bacteriophages <- merge(temperate_bacteriophages_integrase, temperate_bacteriophages_recombinase, by='V1', all=T)
temperate_bacteriophages$keep <- F
temperate_bacteriophages[temperate_bacteriophages$integrase!=0,]$keep <- T
temperate_bacteriophages[temperate_bacteriophages$integrase==0 & 
                           temperate_bacteriophages$recombinase!=0 &
                           (temperate_bacteriophages$V1 %in% all_pvogs_detected[all_pvogs_detected$V3 %in% list_repressor_vogs$V1,]$V1),]$keep <- T

table_of_origin_complete$temperate <- NA
table_of_origin_complete[table_of_origin_complete$V1 %in% unique( temperate_bacteriophages[temperate_bacteriophages$keep==T,]$V1 ),]$temperate <- 1

table_of_origin_complete[is.na(table_of_origin_complete)] <- 0

write.table(table_of_origin_complete, '../clean_viral/table_of_origin/table_of_origin_noneg405_99_der95', sep='\t', quote=F, row.names=F)
