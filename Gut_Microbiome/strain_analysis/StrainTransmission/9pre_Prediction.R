library(tidyverse)
source("Common_function.R")


#Prepare a DF wtih all layers of information: mother-infant distances, Mother in birth or previous timepoint, infant in first two timpoints. Abundance mother, abundance infant, cross-sectional phenotypes.
#Next script: prediction of transfer not/transfer


Read_all = function(){

        file_list <- list.files("Distance_tables", pattern = "*.rds", full.names = TRUE)
        SGB_list = c()
        for (Entry in file_list){
                SGB_list = c(SGB_list, str_remove(basename(Entry), ".rds") )
        }
        return(list(file_list, SGB_list))
}



Individual_SGB = function(SGB, File){
	sp4 = read_rds(File)
	Mother_infant_formatting(sp4) ->  subset_mother_infant_pairs
	if (dim(subset_mother_infant_pairs)[1] == 0){ return() }

	subset_mother_infant_pairs %>% filter(infant_timepoint_categorical %in% c("W2", "M1") ) %>% arrange(infant_timepoint_categorical) %>% distinct(next_id_infant, .keep_all=T) -> subset_mother_infant_pairs
	dim(subset_mother_infant_pairs)[1] -> N_samples
	if ( N_samples < 100) { return() }
	#print(colnames(subset_mother_infant_pairs))
	
	subset_mother_infant_pairs %>% filter(Type_2 == "infant" )  -> Reversed
	Remove_col =  colnames(Reversed)[grepl("_1", colnames(Reversed))]
	Reversed %>% select(-Remove_col) -> Reversed
	colnames(Reversed) = str_replace(colnames(Reversed), "_2", "")

	subset_mother_infant_pairs %>% filter(Type_1 == "infant" )  -> Forward
	Remove_col =  colnames(Forward)[grepl("_2", colnames(Forward))]
	Forward %>% select(-Remove_col) -> Forward
	colnames(Forward) = str_replace(colnames(Forward), "_1", "")
	
	rbind(Forward %>% rename(NEXT_ID=NEXT_ID_long ) , Reversed) -> motherinfant

	#Add abundance infant
	Abundance = read.table("metadata/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt") %>% rownames_to_column("NG_ID") %>% as_tibble()
	Abundance %>% select(-NG_ID) %>% colnames() %>% sapply( function(x){ paste0("InfantAbundance_", x )   } ) -> NewCol
	colnames(Abundance) = c("NG_ID", NewCol)
	motherinfant = left_join(motherinfant, Abundance, by=c("Sample_ID" = "NG_ID") )
	print(motherinfant)
	#Add abundance mother
	Abundance_mother = read.table("metadata/NEXT_metaphlan_4_CLR_transformed_fil_SGB_mothers_03_08_2023.txt") %>% rownames_to_column("NG_ID") %>% as_tibble()
	Abundance_mother %>% select(-NG_ID) %>% colnames() %>% sapply( function(x){ paste0("MotherAbundance_", x )   } ) -> NewCol
	colnames(Abundance_mother) = c("NG_ID", NewCol)
	left_join(metadata %>% filter(NG_ID %in% c(subset_mother_infant_pairs$Sample_ID_1, subset_mother_infant_pairs$Sample_ID_2)  ) %>%  select(NG_ID, FAMILY), Abundance_mother) %>% drop_na() %>% select(-NG_ID) -> Abundance_mother
	
	motherinfant = left_join(motherinfant, Abundance_mother, by="FAMILY" ) -> motherinfant
	
	write_tsv(motherinfant, paste0("DataFrame_allInfo/", SGB, ".tsv") )
	
}


Attempt_all = function(){
        Input = Read_all()
        file_list = Input[[1]]
        SGB_list = Input[[2]]

        All = tibble()
        for (i in 1:length(file_list)){
                nSGB = SGB_list[i]
		nFile = file_list[i]
		print(paste(nSGB, nFile, collapse="\t"))
		Individual_SGB(nSGB, nFile)
        }
}


Attempt_all()
#Individual_SGB("t__SGB14546_group", "Distance_tables/t__SGB14546_group.rds")
