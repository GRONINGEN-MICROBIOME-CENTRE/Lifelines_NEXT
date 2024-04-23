library(tidyverse)
source("Common_function.R")


Read_all = function(){

	file_list <- list.files("Distance_tables", pattern = "*.rds", full.names = TRUE)
	SGB_list = c()
	for (Entry in file_list){
		SGB_list = c(SGB_list, str_remove(basename(Entry), ".rds") )
	}
	return(list(file_list, SGB_list))
}


Attempt_all = function(){
	Input = Read_all()
	file_list = Input[[1]]
	SGB_list = Input[[2]]

	All = tibble()
	for (i in 1:length(file_list)){
		SGB = SGB_list[i]
		sp4 = read_rds(file_list[i])
		Mother_infant_formatting(sp4) ->  subset_mother_infant_pairs
		print(SGB)
		if (dim(subset_mother_infant_pairs)[1] == 0){ next }
		#If there is only one timepoint per mother, then if there are >= 3 entries per baby, it means it has at least three available timepoints
		#subset_mother_infant_pairs %>% group_by(next_id_infant) %>% summarise(N = n()) %>% group_by(N) %>% summarise(N_babies_with_X_timepoints = n()) -> Babies_info
		#write_tsv(Babies_info, paste0("NBabies_with_XTimepoints/", SGB, ".tsv") )
		
		#Babies_info %>% filter(N>=3) -> Babies_with_3
		#if (dim(Babies_with_3)[1] == 0){ Babies_with_3 = 0 
		#} else { Babies_with_3 = sum(Babies_with_3$N_babies_with_X_timepoints) }
		dim(subset_mother_infant_pairs)[1] -> N_samples
		#tibble(SGB = SGB, NBabies_with_3_timpeoints = Babies_with_3) %>% rbind(All, . ) -> All
		tibble(SGB = SGB, N = N_samples) %>% rbind(All, . ) -> All
	}
	write_tsv(All, "SGB_babies.tsv")
}
Attempt_one = function(SGB="t__SGB1862"){
	RDS = paste0("Distance_tables/",SGB, ".rds")
	sp4 = read_rds(RDS)
	print(sp4)
	Mother_infant_formatting(sp4) ->  subset_mother_infant_pairs
	print(subset_mother_infant_pairs)
	if (dim(subset_mother_infant_pairs)[1] == 0){ print("Empty!") ; q() }
	subset_mother_infant_pairs %>% group_by(next_id_infant) %>% summarise(N = n()) %>% group_by(N) %>% summarise(N_babies_with_X_timepoints = n()) -> Babies_info
	print(Babies_info)
	write_tsv(Babies_info, paste0("NBabies_with_XTimepoints/", SGB, ".tsv") )
}

Attempt_missing = function(){
	Input = Read_all()
	file_list = Input[[1]]
	SGB_list = Input[[2]]
	for (i in 1:length(file_list)){
		SGB = SGB_list[i]
		Output = paste0("NBabies_with_XTimepoints/", SGB, ".tsv")
		if (file.exists(Output)){ next }
		print(SGB)
		sp4 = read_rds(file_list[i]);Mother_infant_formatting(sp4) ->  subset_mother_infant_pairs
		if (dim(subset_mother_infant_pairs)[1] == 0){ next }
		subset_mother_infant_pairs %>% group_by(next_id_infant) %>% summarise(N = n()) %>% group_by(N) %>% summarise(N_babies_with_X_timepoints = n()) -> Babies_info
		write_tsv(Babies_info,Output)
	}
}

Merge_tables =  function(){
	file_list <- list.files("NBabies_with_XTimepoints", pattern = "*.tsv", full.names = TRUE)
	SGB_list = c()
	for (Entry in file_list){ SGB_list = c(SGB_list, str_remove(basename(Entry), ".tsv") ) }
	All = tibble()
	for (i in 1:length(file_list)){
		SGB = SGB_list[i]
		Babies_info = read_tsv(file_list[i])
		Babies_info %>% filter(N>=3) -> Babies_with_3
		if (dim(Babies_with_3)[1] == 0){ Babies_with_3 = 0
		} else {Babies_with_3 = sum(Babies_with_3$N_babies_with_X_timepoints) }
		tibble(SGB = SGB, NBabies_with_3_timpeoints = Babies_with_3) %>% rbind(All, . ) -> All		
	}
	write_tsv(All, "SGB_babieswith3.tsv")

}

Attempt_all()

#Attempt_one("t__SGB1863")
#Attempt_missing()
#Merge_tables()
