## This chunk of code is now a part of the qc_assembly_summary_stats.R

library(readr)
library(dplyr)
library(openxlsx)
library(stringr)
library(UpSetR)

setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\data")

to_garmaeva <- as.data.frame(read_tsv("table_of_origin_garmaeva"))
to_shah <- as.data.frame(read_tsv("table_of_origin_shah"))
to_walters <- as.data.frame(read_tsv("table_of_origin_walters"))
to_liang <- as.data.frame(read_tsv("table_of_origin_liang"))
to_maqsood <- as.data.frame(read_tsv("table_of_origin_maqsood"))
to_full <- rbind(to_garmaeva, to_shah, to_walters, to_liang, to_maqsood)


to_full_per_sample <- to_full
to_full_per_sample$V1 <- sub("_NODE.*", "", to_full_per_sample$V1)
to_full_per_sample_test <- summarise_all(group_by(to_full_per_sample, V1), sum)


summary_per_tool <- as.data.frame(matrix(NA, nrow=4, ncol=2))
summary_per_tool$V1 <- colnames(to_full[,c(2:5)])
for (i in summary_per_tool$V1) {
  summary_per_tool[summary_per_tool$V1==i,]$V2 <- sum(to_full[,i])
}



listInput <- list(geNomad=to_full[to_full$geNomad==1,]$V1,
                  VIBRANT=to_full[to_full$VIBRANT==1,]$V1,
                  DeepVirFinder=to_full[to_full$DeepVirFinder==1,]$V1,
                  VirSorter2=to_full[to_full$VirSorter2==1,]$V1)

pdf('C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\plots\\VD_redundant_tools_overlap.pdf', width=18/2.54, height=10/2.54)
upset(fromList(listInput), order.by = "freq", sets.bar.color = "#C00000", 
      number.angles = 20,
      sets.x.label = "N detected virus contigs", scale.sets = "identity",
      text.scale = c(1, 1, 1, 1, 1, 0.75))
dev.off()
