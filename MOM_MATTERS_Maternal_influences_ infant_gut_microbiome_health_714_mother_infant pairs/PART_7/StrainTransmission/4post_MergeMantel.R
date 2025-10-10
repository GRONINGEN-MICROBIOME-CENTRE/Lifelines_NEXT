.libPaths(c("/groups/umcg-llnext/tmp01/umcg-sandreusanchez/Strain_transmission/Packages", .libPaths()))

library(tidyverse)
library(cluster)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ape)
library(phytools)

if (!file.exists("Mantel_distance/SummaryStatsMantel.csv") ){
	file_list <- list.files("Mantel_distance", pattern = "\\.tsv$", full.names = TRUE)
	# Read and bind rows
	merged_data <- map_dfr(file_list, read_tsv)
	# View merged data
	merged_data %>% arrange(P) %>% mutate(FDR = p.adjust(P)) -> merged_data
	print(merged_data)
	write_csv(merged_data, "Mantel_distance/SummaryStatsMantel.csv")
}

#Try to visualize, tree annotated with provicnes /  families?
Latitude_longitude_info = read_csv("/groups/umcg-llnext/tmp01/NEXT_metadata/geo_location_NEXT/LLNEXT_geo_location_municipality.csv")
Sign_SGB = c("t__SGB14242", "t__SGB1587",  "t__SGB1822", "t__SGB1897", "t__SGB2299", "t__SGB6584")

metadata<-read.delim("metadata/LLNEXT_metadata_03_01_2024.txt")

#assess the ideal number of clusters
#select(drop_na(Latitude_longitude_info_t), c(lon, lat) ) %>% clusGap(., FUNcluster = kmeans, K.max = 15) -> gap.stat
#factoextra::fviz_gap_stat(gap.stat)
select(drop_na(Latitude_longitude_info), c(lon, lat) ) %>% kmeans(4) -> Clusters
Latitude_longitude_info %>% mutate(Cluster=as.factor(Clusters$cluster) ) -> Latitude_longitude_info
#Latitude_longitude_info_t %>% mutate(Cluster=as.factor(Clusters$cluster)) %>% ggplot(aes(x=lon,y=lat)) + geom_label(aes(label=FAMILY, col=municipality ), size=2.3) + scale_color_manual(values=c25) + new_scale_color() + theme_bw() + stat_ellipse(aes(col=Cluster), type = "t") +  scale_color_manual(values=c25)

c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown", "purple", "red" )

Do_plot = function( SGB ){
	#1. Get tree
	print("Get and clean tree")
	read.tree( paste0("Trees/RAxML_bestTree.", SGB,".StrainPhlAn4.tre")) -> tree
	#2. Clean
	str_replace(tree$tip.label, "_kneaddata_cleaned_pair_metaphlan4", "") -> tree$tip.label
	substr(tree$tip.label, 1, 12) -> tree$tip.label

	#Prepare metadata
	print("Prepare metadata")
	mothers <-  as_tibble(metadata) %>% filter(Type == "mother") %>% select(NEXT_ID, FAMILY) %>% distinct(FAMILY, .keep_all=T) %>% arrange(FAMILY)
metadata %>% as_tibble() %>% left_join(mothers, by = "FAMILY", suffix = c("", "_Mother")) %>% rename(next_id_mother= NEXT_ID_Mother) %>%
  select(NG_ID,NEXT_ID, Type, FAMILY,next_id_mother  ) %>% left_join(Latitude_longitude_info) %>% filter(!is.na(lon)) -> Latitude_longitude_info_t
	#match
	Latitude_longitude_info_t %>% filter(NG_ID %in% tree$tip.label) -> Latitude_longitude_info_t
	keep.tip(tree, Latitude_longitude_info_t$NG_ID) -> tree
	Latitude_longitude_info_t[ match(tree$tip.label, Latitude_longitude_info_t$NG_ID),  ] -> Latitude_longitude_info_t
	print("Plotting tree")
	tree = midpoint.root(tree)

	ggtree(tree, layout="fan", open.angle=15, size=0.1) %<+% Latitude_longitude_info_t -> p
	p + geom_tiplab(aes(label=FAMILY, col=Type)) + new_scale_color() +  geom_tippoint(aes(col=Cluster)) +  scale_color_manual(values=c25) -> p
	
	Out = paste0("plots_associations/", SGB, "_TreeANDgeography.pdf")
	print(Out)
	ggsave(Out, p)
}

for (SGB in Sign_SGB){
	Do_plot(SGB)
}










