library(tidyverse)
library(vegan)
library(geosphere)



Latitude_longitude_info = read_csv("/groups/umcg-llnext/tmp01/NEXT_metadata/geo_location_NEXT/LLNEXT_geo_location_municipality.csv")

RunAnalysis = function(sp, SGB_name){

	rownames(sp) = substr(rownames(sp), 1, 12)
	colnames(sp) = rownames(sp)
	
	mothers <-  as_tibble(metadata) %>% filter(Type == "mother") %>% select(NEXT_ID, FAMILY) %>% distinct(FAMILY, .keep_all=T) %>% arrange(FAMILY)

	metadata %>% as_tibble() %>% left_join(mothers, by = "FAMILY", suffix = c("", "_Mother")) %>% rename(next_id_mother= NEXT_ID_Mother) %>%
  select(NG_ID,NEXT_ID, Type, FAMILY,next_id_mother  ) %>% metadata_motherIDs
	
	metadata_motherIDs %>% drop_na() -> metadata_motherIDs

	metadata_motherIDs %>% left_join(Latitude_longitude_info) -> metadata_motherIDs
	metadata_motherIDs %>% filter(NG_ID %in% rownames(sp) ) -> metadata_motherIDs
	rownames(sp)[! is.na((match(rownames(sp), metadata_motherIDs$NG_ID ) )) ] -> Keep
	sp[Keep, Keep] -> sp
	match(metadata_motherIDs$NG_ID, rownames(sp) ) -> Order  
	sp[Order, Order] -> sp2

	#Create covariate matrix
	# Get unique NG_IDs and families
	unique_ng_ids <- unique(metadata_motherIDs$NG_ID)
	unique_families <- unique(metadata_motherIDs$FAMILY)
	# Initialize the matrix with zeros
	n <- length(unique_ng_ids)
	matrix_covariate <- matrix(0, nrow = n, ncol = n, dimnames = list(unique_ng_ids, unique_ng_ids))
	# Populate the matrix based on family similarity
	for (i in 1:n) {
		for (j in 1:n) {
			if (metadata_motherIDs$FAMILY[metadata_motherIDs$NG_ID == unique_ng_ids[i]] ==  metadata_motherIDs$FAMILY[metadata_motherIDs$NG_ID == unique_ng_ids[j]]) {
      				matrix_covariate[i, j] <- 1
    			}
  		}
	}


	metadata_motherIDs %>% select(lon, lat) %>% as.matrix()  %>% distm(. ,  fun = distVincentySphere) -> geographical_distance

	Result = mantel.partial(sp2, geographical_distance, matrix_covariate, method = "pearson", permutations = 30000)
	tibble(SGB = SGB_name, Rho = Result$statistic, P = Result$signif, Permutations = Result$permutations, N = dim(sp2)[1], N_uniq = length(unique_ng_ids), N_families = length(unique_families)    ) -> Res

	return(Res)
}


args <- commandArgs(trailingOnly = TRUE)
file = args[1]
print(paste0("Initializing: ", file))


file_new = str_remove(basename(file), "\\.txt") 
SGB_name = str_remove(file_new, "\\_DistMat") 
print(SGB_name)

print("Reading file and metadata")
strain = read.delim(file)
metadata<-read.delim("metadata/LLNEXT_metadata_03_01_2024.txt")

RunAnalysis(strain, SGB_name) -> Out
print(Out)
paste0("Mantel_distance/", SGB_name, ".tsv") -> O
write_tsv(Out, O)


