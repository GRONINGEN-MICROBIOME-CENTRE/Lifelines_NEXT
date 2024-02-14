library(readr)

setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\data")

files_rename <- as.data.frame(read_tsv("filereport_read_run_PRJEB46943.txt"))
df_split <- data.frame(do.call('rbind', strsplit(as.character(files_rename$fastq_ftp), ';')))

files_rename <- cbind(df_split, files_rename[c(2,3)])
files_rename$X1 <- sub(".*/", "", files_rename$X1)
files_rename$X2 <- sub(".*/", "", files_rename$X2)

pairedfwd_files_rename <- files_rename[grepl("_1", files_rename$X1), ]
pairedrev_files_rename <- files_rename[grepl("_1", files_rename$X1), ]
unmatched_files_rename <- files_rename[!grepl("_", files_rename$X1), ]

pairedfwd_files_rename$sample_alias <- paste(pairedfwd_files_rename$sample_alias, "_1.fasq.gz", sep = "")
pairedfwd_files_rename <- pairedfwd_files_rename[c(1,3)]
colnames(pairedfwd_files_rename) <- c("before", "after")

pairedrev_files_rename$sample_alias <- paste(pairedrev_files_rename$sample_alias, "_2.fasq.gz", sep = "")
pairedrev_files_rename <- pairedrev_files_rename[c(2,3)]
colnames(pairedrev_files_rename) <- c("before", "after")

unmatched_files_rename$sample_alias <- paste(unmatched_files_rename$sample_alias, ".fasq.gz", sep = "")
unmatched_files_rename <- unmatched_files_rename[c(2,3)]
colnames(unmatched_files_rename) <- c("before", "after")

final_rename <- rbind(pairedfwd_files_rename, pairedrev_files_rename, unmatched_files_rename)

write.table(final_rename, file = "rename.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
