library(tidyverse)


#Find twins
metadata<-read.delim("metadata/LLNEXT_metadata_03_01_2024.txt")
metadata %>% as_tibble() %>% filter(Type == "infant") %>% group_by(FAMILY) %>% summarise(distinct_next_ids = n_distinct(NEXT_ID)) %>% filter(distinct_next_ids >= 2 ) -> families
metadata %>% filter(Type == "infant" & FAMILY %in% families$FAMILY) %>% as_tibble() -> metadata


Read_and_proces = function(File = "Distance_tables/t__SGB1815.rds", SGB_name ="t__SGB1815"){

	read_rds(File) -> sp4
	Size = dim(sp4)[1]

	sp5<-sp4 %>% filter(Type_1 == "infant" & Type_2 == "infant") %>% filter(NEXT_ID_short_1 != NEXT_ID_short_2) %>% filter(FAMILY_1 == FAMILY_2) %>% filter(FAMILY_1 %in% metadata$FAMILY)
	Size_twins = dim(sp5)[1]
	if ( Size_twins == 0) { return(NA) }

	sp5 %>% mutate(SGB =  SGB_name) %>% select(SGB,  NEXT_ID_short_1, NEXT_ID_short_2, FAMILY_1, Timepoint_categorical_1, Timepoint_categorical_2, dist, Strain_sharing ) %>% mutate(Strain_sharing = ifelse(Strain_sharing=="yes", 1, ifelse(Strain_sharing == "no", 0, NA))) %>% return()

}

Read_all = function(){

        file_list <- list.files("Distance_tables", pattern = "*.rds", full.names = TRUE)
        SGB_list = c()
        for (Entry in file_list){
                SGB_list = c(SGB_list, str_remove(basename(Entry), ".rds") )
        }
        return(list(file_list, SGB_list))
}

DO_ALL = T
#Read_and_proces("Distance_tables/t__SGB14947.rds", SGB_name = "t__SGB14947") -> Res
if (DO_ALL == T){
	Read_all() -> Info
	Merged_table = tibble()
	file_list = Info[[1]]
	SGB_list = Info[[2]]
	for(i in 1:length(file_list)){
		SGB = SGB_list[[i]]
		File_name = file_list[[i]]
		Read_and_proces(File_name,SGB) -> Res
		C = class(Res) ; if(length(C) > 1){ C = C[1] }
		if  (C == "logical") {  next }
		rbind(Merged_table, Res) -> Merged_table
	}
}

write_rds(Merged_table, "Results/SharingTwins.rds")

