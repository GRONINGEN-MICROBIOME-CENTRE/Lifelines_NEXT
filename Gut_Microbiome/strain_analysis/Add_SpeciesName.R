library(tidyverse)

Info = read_tsv("metadata/link_SGB_species_names.txt")
Files_to_add = c("LMERTEST_merged_output.tsv", "SGB_perc_sharing_mothertimepoint.tsv","SGB_perc_sharing.tsv")
for (File in Files_to_add){
	read_tsv(paste0("Results/", File) ) -> File_read
	if ("Species" %in% colnames(File_read) ) { next }
	File_read %>% left_join(Info, by="SGB") %>% write_tsv(. , paste0("Results/", File) )


}

