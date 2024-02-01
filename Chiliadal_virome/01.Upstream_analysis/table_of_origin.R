library(reshape2)
library(dplyr)

temp = list.files(path="../VIR_DB/table_of_origin",pattern="*_table_of_origin",full.names=T)

myfiles = lapply(temp, read.table)

import_data=lapply(myfiles, function(x) {
                   out <- list(dcast(x, V1~V2))
                   return(out)
                   })
table_of_origin_redundant<-bind_rows(import_data)
table_of_origin_redundant[-1] <- as.integer(table_of_origin_redundant[-1] != 0) 

table_of_origin_redundant[is.na(table_of_origin_redundant)] <- 0

write.table(table_of_origin_redundant, '../VIR_DB/table_of_origin/table_of_origin', sep='\t', quote=F, row.names=F)
