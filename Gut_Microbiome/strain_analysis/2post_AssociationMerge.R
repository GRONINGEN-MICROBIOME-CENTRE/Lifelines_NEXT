library(tidyverse)

folder ="Results"
Output = paste0(folder, "/LMERTEST_merged_output.tsv")

if (! file.exists(Output)){
	file_list = list.files(path=folder,pattern=".*LMERTEST.*\\.tsv$") # Names of all the species you have distance matrices for 
	Summary_statistics = tibble()
	for (i in file_list){
		if ( grepl("OverallSharing", i) ){ next }
		read_tsv(paste0(folder, "/", i)) -> d
		if (dim(d)[1] == 1 ){ 
			if (is.na(d$Estimate)) { next }
		 }  
		d %>% rbind(Summary_statistics, .) -> Summary_statistics	
	}
write_tsv(Summary_statistics, Output)
} else {
 Summary_statistics = read_tsv(Output)
}

Summary_statistics %>% mutate(FDR_all = p.adjust(`Pr(>|z|)`, "fdr")) -> Summary_statistics
Summary_statistics %>% arrange(FDR_all) -> Summary_statistics

Info = read_tsv("metadata/link_SGB_species_names.txt")
Summary_statistics %>% left_join(Info, by="SGB") -> Summary_statistics

write_tsv(Summary_statistics, Output) 

q()
################
#Make some plots#
################
source("Common_function.R")
library(ggforce)

pPlot_stuff = function(SGB, Phenotype, Plot_type="Bar" ){
	print(paste0("Preparing for plot: ", SGB, " with ", Phenotype))
	Input = paste0("Distance_tables/",SGB,".rds")
	out_name = paste0("plots_associations/sharing_",SGB,"_vs_",Phenotype,"_",Plot_type,".pdf")
	sp4 = read_rds(Input) ; Mother_infant_formatting(sp4) ->  subset_mother_infant_pairs
	
        Plot_stuff(subset_mother_infant_pairs, out_name)
        #if (Plot_type =="Box"){
	#	if ( class( as_vector(subset_mother_infant_pairs[,Phenotype])) == "integer") { subset_mother_infant_pairs[,Phenotype] = as.factor(as_vector(subset_mother_infant_pairs[,Phenotype])) } 
	#	subset_mother_infant_pairs %>% ggplot(aes_string(y="nGD", x=Phenotype)) +  geom_boxplot(outlier.shape = NA) + geom_sina(alpha=0.5) + theme_bw() -> Plot
	#}else if( Plot_type == "Bar"){
	#	if ( class( as_vector(subset_mother_infant_pairs[,Phenotype])) == "integer") { subset_mother_infant_pairs[,Phenotype] = as.factor(as_vector(subset_mother_infant_pairs[,Phenotype])) }

	#	counts <- subset_mother_infant_pairs %>%
	#		count(!!sym(Phenotype), Strain_sharing) %>%
  	#		group_by(!!sym(Phenotype)) %>%
  	#		mutate(proportion = n / sum(n)) %>% drop_na()

	#	ggplot(counts, aes_string(x = Phenotype, y = "proportion", fill = "Strain_sharing" )) +
  	#	geom_bar(stat = "identity", position = "stack") +
  	#	geom_text(aes(label = n, y = proportion), col="white", vjust = -0.5, size = 5) +
  	#	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values= c("#243b74", "#cb343f")) + theme_bw() -> Plot
	#}
	#ggsave(out_name, Plot)
}

Summary_statistics %>% select(SGB, Overall_name) %>% distinct(.keep_all=T) -> ForPlot
for (Entry in  seq(1, dim(ForPlot)[1])){
	E = ForPlot[Entry,]
	SGB = E$SGB
	Pheno = E$Overall_name
	pPlot_stuff(SGB, Pheno)
}

