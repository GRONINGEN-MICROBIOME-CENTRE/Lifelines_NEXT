library(tidyverse)
source("Common_function.R")

Persistent_formatting = function(sp4){
  #Keep only comparisons among the same sample
  sp5<-sp4 %>% filter(related == "related")
  infant_pairs <- sp5 %>% filter( Type_1 == "infant" & Type_2 == "infant") 
  infant_pairs %>% filter(NEXT_ID_short_1 == NEXT_ID_short_2) -> infant_pairs
  #Keep only early and late timepoints
  Early = c("W2", "M1"); Late = c("M9", "M12")
  infant_pairs %>% filter(  ( (Timepoint_categorical_1 %in% Early) & (Timepoint_categorical_2 %in% Late) ) | ( (Timepoint_categorical_2 %in% Early) & (Timepoint_categorical_1 %in% Late) )  ) -> infant_pairs
  #If multiple comparisons within the same infant, annotate which one is the maximum time difference
infant_pairs %>% group_by(NEXT_ID_short_1) %>% mutate(Maximum_distance = ifelse(time_day_diff == max(time_day_diff), T, F ))  %>% ungroup() -> infant_pairs
  #Prepare output: Table per NEXT_ID. If there is strain sharing by the maximum time difference, then we consider it Persistant
	infant_pairs %>% select(NEXT_ID_short_1, time_day_diff, Maximum_distance, Strain_sharing,Timepoint_categorical_1, Timepoint_categorical_2) %>% rename(NEXT_ID =NEXT_ID_short_1 ) %>%
  group_by(NEXT_ID) %>% mutate(Persistence = any(Maximum_distance == TRUE & Strain_sharing == "yes")) %>% ungroup() -> Persistence_info
	return(Persistence_info)
} 


Do_one = function(file){
	print(file)
	Run_analysis(file)
}
Do_multiple = function(){
	rds_files <- list.files(path = "Distance_tables", pattern = "\\.rds$", full.names = TRUE)
	results = tibble()
	for (file in rds_files){
		print(file)
		Res = Run_analysis(file); Res = Res[[1]]	
		if (dim(Res)[1] == 0) { next }
		results %>% rbind(Res, . ) -> results
	}
	write_tsv(results, "Results/Persistance/AllPersistanceAss.tsv")
}
Run_analysis = function(file){
	SGB_name = str_remove(basename(file), "\\.rds") 
	sp4 = read_rds(file)
	Persistent_formatting(sp4) -> Persistance_table
	print(Persistance_table)
	#Continue only if we output something
	if (dim(Persistance_table)[1] < 2 ){ return(tibble()) }
	write_tsv(Persistance_table, paste0("StrainPersistence/", SGB_name, ".tsv"))
	#Next check if the strain present in the infant is transmitted from mother      
	Mother_infant_formatting(sp4) ->  subset_mother_infant_pairs
	if (dim(subset_mother_infant_pairs)[1] < 10){ return(list()) }
	subset_mother_infant_pairs %>% select(sample_id_infant,next_id_infant, Strain_sharing, infant_timepoint_categorical) %>% rename(NEXT_ID = next_id_infant) -> subset_mother_infant_pairs
	subset_mother_infant_pairs %>% left_join(Persistance_table %>% select(-Strain_sharing) %>% filter(Maximum_distance==T)  ) %>% drop_na() -> InfoAnalysis
	InfoAnalysis %>% filter(infant_timepoint_categorical %in% c(Timepoint_categorical_1, Timepoint_categorical_2) ) %>% filter(infant_timepoint_categorical %in% c("M1", "W2") )-> InfoAnalysis
	if (dim(InfoAnalysis)[1] < 20 ){ return(list()) }
	InfoAnalysis %>% group_by(Strain_sharing, Persistence) %>% summarise(count = n()) -> contingency_table
	#contingency_table_wide <- pivot_wider(contingency_table, names_from = Strain_sharing, values_from = count, values_fill = 0)
	#colnames(contingency_table_wide) <- c("Persistence", "No Strain Sharing", "Strain Sharing")
	#contingency_table_wide %>% as.data.frame() %>% column_to_rownames("Persistence") -> I
	contingency_table_wide2 <- pivot_wider(contingency_table, names_from = Persistence, values_from = count, values_fill = 0)
	contingency_table_wide2 %>% select(Strain_sharing,`TRUE`, `FALSE`) -> contingency_table_wide2
	contingency_table_wide2 %>% arrange(desc(Strain_sharing)) -> contingency_table_wide2
	print(contingency_table_wide2)
	contingency_table_wide2 %>% as.data.frame() %>% column_to_rownames("Strain_sharing") -> I
	I %>% fisher.test(.,alternative = "greater") -> Test
	print(Test)
	res = tibble(SGB = SGB_name, P=Test$p.value, Odds_ratio=Test$estimate, N_PersistentShared=I[1,1], N_PersisntentnotShared=I[2,1], N_SharedNoPersistent=I[1,2], N_NoSharedNoPersistent=I[2,2]  )
	return(list(res,InfoAnalysis ) )
}
#Do_multiple()
Do_one("Distance_tables/t__SGB17244.rds") -> Adolescentis
print(Adolescentis[[2]])
Adolescentis[[2]] %>% mutate(Strain_sharing = ifelse(Strain_sharing=="no","No early mother-infant\nstrain sharing", "Early mother-infant\nstrain sharing" ) )  %>% ggplot(aes(x=Strain_sharing, fill=Persistence )) + geom_bar(position = "fill") + theme_bw() + xlab("Early timepoint strain sharing") + ylab("Infant number") + geom_text(aes(label = ..count.., group = Persistence), stat = "count",  position = position_fill(vjust = 0.5),vjust = -0.5,color = "white") + labs(fill = "Late timepoint\nstrain consistency") + scale_fill_manual(values = c("FALSE" = "#dce0d9", "TRUE" = "#31081f")) + ylab("SGB17244/B.adolescentis\nproportion strain persistance") + xlab(NULL) -> Plot
ggsave("StrainPersistence/Plots/StrainConsistency_and_Sharing_t__SGB17244.pdf", Plot, width=4, height=2 )

Do_one("Distance_tables/t__SGB1861.rds") -> Theta
print(Theta[[2]])
Theta[[2]] %>% mutate(Strain_sharing = ifelse(Strain_sharing=="no","No early mother-infant\nstrain sharing", "Early mother-infant\nstrain sharing"     ) ) %>% ggplot(aes(x=Strain_sharing, fill=Persistence )) + geom_bar(position = "fill") + theme_bw() + xlab("Early timepoint strain sharing") + ylab("Infant number") + geom_text(aes(label = ..count.., group = Persistence), stat = "count",  position = position_fill(vjust = 0.5),vjust = -0.5,color = "white")  + labs(fill = "Late timepoint\nstrain consistency") + scale_fill_manual(values = c("FALSE" = "#dce0d9", "TRUE" = "#31081f"))  + ylab("SGB17244/B. thetaiotamicron\nproportion strain persistance") + xlab(NULL)  -> Plot2
ggsave("StrainPersistence/Plots/StrainConsistency_and_Sharing_t__SGB1861.pdf", Plot2, width=4, height=2 ) 


if (file.exists("Results/Persistance/AllPersistanceAss.tsv")){
	Results = read_tsv("Results/Persistance/AllPersistanceAss.tsv")
	if (! "Species" %in% colnames(Results) ) { 
		read_tsv("metadata/link_SGB_species_names.txt") -> link	
		left_join(Results, link) -> Results
		Results %>% mutate(FDR = p.adjust(P, "fdr")) %>% arrange(P) -> Results
		print(Results)
		write_tsv(Results,"Results/Persistance/AllPersistanceAss.tsv")
	}

}


