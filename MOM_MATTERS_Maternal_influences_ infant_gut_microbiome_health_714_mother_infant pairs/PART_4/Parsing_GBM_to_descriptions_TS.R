#Parsing GMB to descriptions # 

raw_lines <- readLines("GBM.inputfile.txt")  # Replace with actual file name

# Keep only lines that start with MGB and have a description
mgb_lines <- raw_lines[grepl("^MGB\\d+\\t", raw_lines)]

# Extract MGB ID and description
mgb_df <- as.data.frame(do.call(rbind, strsplit(mgb_lines, "\t")), stringsAsFactors = FALSE)
colnames(mgb_df) <- c("MGB_ID", "Pathway")

write.table(mgb_df, "MGBM_pathway_descriptions.txt", sep = "\t", row.names = F)

check<-read.delim("NEXT_GBM_merged_clean.tsv")
