library(tidyverse)
library(ggforce)
source("Common_function.R")
library(ggrepel)
#Make plot with timepoint vs SGB shating rate




#per RDS: Get mother-baby pairs, bet the % of transmission
#Get_average = 

Read_all = function(){

        file_list <- list.files("Distance_tables", pattern = "*.rds", full.names = TRUE)
        SGB_list = c()
        for (Entry in file_list){
                SGB_list = c(SGB_list, str_remove(basename(Entry), ".rds") )
        }
        return(list(file_list, SGB_list))
}

D = Read_all()
SGBs = D[[2]]
Files = D[[1]]

if (! file.exists("Results/SGB_perc_sharing.tsv") ) { 

Avr_sharing_table = tibble()
for (i in seq(1, length(SGBs))){
	SGB_name = SGBs[[i]]
	sp4 = read_rds(Files[[i]])
	Mother_infant_formatting(sp4 ) ->  subset_mother_infant_pairs
	#print(subset_mother_infant_pairs)
	for (Timepoint in unique(subset_mother_infant_pairs$infant_timepoint_categorical)) {
		t_subset_mother_infant_pairs = filter(subset_mother_infant_pairs, infant_timepoint_categorical==Timepoint)
		#if (dim(t_subset_mother_infant_pairs)[1] < 5  ){  Avr_sharing_table = tibble( SGB = SGB_name, Timepoint = Timepoint, Avg_sharing = NA, N_pairs = dim(t_subset_mother_infant_pairs)[1] ) %>% rbind(Avr_sharing_table, .) ; next  }
		Avr_sharing = 100* dim(filter(t_subset_mother_infant_pairs, Strain_sharing == "yes"))[1]  / dim(t_subset_mother_infant_pairs)[1]
		cat(SGB_name, Timepoint, Avr_sharing, "\n")
		Avr_sharing_table = tibble( SGB = SGB_name, Timepoint = Timepoint, Avg_sharing = Avr_sharing, N_pairs = dim(t_subset_mother_infant_pairs)[1] ) %>% rbind(Avr_sharing_table, .)		
	}
}

print(Avr_sharing_table)

Info = read_tsv("metadata/link_SGB_species_names.txt")
Avr_sharing_table %>% left_join(Info, by="SGB") -> Avr_sharing_table

write_tsv(Avr_sharing_table, "Results/SGB_perc_sharing.tsv")



} else { read_tsv("Results/SGB_perc_sharing.tsv") -> Avr_sharing_table }
#SGB               Timepoint Avg_sharing N_pairs
timepoint_order <- c("W2", "M1", "M2", "M3", "M6", "M9", "M12")
Avr_sharing_table %>% mutate(Timepoint = factor( Timepoint, levels = timepoint_order) ) -> Avr_sharing_table
#Avr_sharing_table %>% ggplot(aes(x = Timepoint, y = Avg_sharing) ) + geom_boxplot(outlier.shape = NA) + geom_sina(alpha=0.5)  + theme_bw() -> Plot1
Avr_sharing_table %>% group_by(SGB) %>% summarise( MA = max(Avg_sharing), MI = min(Avg_sharing), Sample= mean(N_pairs) ) %>% mutate(Change = MA -MI)  %>% arrange(desc(Change))
Avr_sharing_table %>% mutate(Trend = ifelse(SGB %in% c("t__SGB1855_group", "t__SGB1832" ), "high_sharing", ifelse(SGB %in% c("t__SGB8007_group", "t__SGB8002"), "low_sharing", ifelse(SGB %in% c("t__SGB2303", "t__SGB17248"), "high_sharing_change", ifelse(SGB %in% c("t__SGB10068", "t__SGB17247", "t__SGB1815" ), "other_sharing",  "other") ) ) ) ) -> Avr_sharing_table
Avr_sharing_table %>% arrange(desc(SGB), desc(N_pairs)) -> Avr_sharing_table
Avr_sharing_table %>% filter(N_pairs >= 7) -> Avr_sharing_table
Avr_sharing_table %>% mutate(`Number of mother infant pairs` = N_pairs ) %>% ggplot(aes(x=Timepoint, y = Avg_sharing )) + geom_boxplot(outlier.shape = NA) + geom_point(aes(size=`Number of mother infant pairs`, color=Trend ), alpha = ifelse(Avr_sharing_table$Trend == "other", 0.4, 1 ) ) + theme_bw() + 
  geom_line(data = Avr_sharing_table, aes(x = Timepoint, y = Avg_sharing, group = SGB, color= Trend ), alpha = 0.4, linetype= 2 ) + scale_color_manual(values = c("high_sharing" = "#136F63", "low_sharing" = "#032B43", "high_sharing_change" ="#D00000","other_sharing" = "#CC9200" ,  "other" = "grey" ) ) +
  geom_text_repel(data = Avr_sharing_table %>% filter(Timepoint != "M12") %>% distinct(SGB, Trend, .keep_all = TRUE) %>%  filter(Trend != "other"), aes(label = paste0( str_replace(SGB, "t__", "") ,"\n", Species) , color = Trend), hjust = -0.07, vjust = 0.5) + 
  theme(legend.position = "bottom") + ylab( "SGB strain sharing rate (%)"  ) + xlab("Timepoint") -> Plot1


ggsave("plots_associations/SGB_averagesahring_perTime.pdf", Plot1, width=5.5, height=4.5)


#Do the same but taking into account mother timepoints 

if (! file.exists("Results/SGB_perc_sharing_mothertimepoint.tsv") ) {
Avr_sharing_table2 = tibble()
for (i in seq(1, length(SGBs))){
        SGB_name = SGBs[[i]]
        sp4 = read_rds(Files[[i]])
        Mother_infant_formatting(sp4, select_mother = F) ->  subset_mother_infant_pairs
        #print(subset_mother_infant_pairs)
        for (Timepoint in unique(subset_mother_infant_pairs$infant_timepoint_categorical)) {
		t_subset_mother_infant_pairs = filter(subset_mother_infant_pairs, infant_timepoint_categorical==Timepoint)
		for (MotherTimepoint in unique(t_subset_mother_infant_pairs$mother_timepoint_categorical) ){
			t_subset_mother_infant_pairs2 = filter(t_subset_mother_infant_pairs, mother_timepoint_categorical==MotherTimepoint)
                	#if (dim(t_subset_mother_infant_pairs2)[1] < 5  ){ next }
                	Avr_sharing = 100* dim(filter(t_subset_mother_infant_pairs2, Strain_sharing == "yes"))[1]  / dim(t_subset_mother_infant_pairs2)[1]
                	cat(SGB_name, Timepoint,MotherTimepoint, Avr_sharing, "\n")
                	Avr_sharing_table2 = tibble( SGB = SGB_name, Timepoint = Timepoint, Mother_timepoint = MotherTimepoint , Avg_sharing = Avr_sharing, N_pairs = dim(t_subset_mother_infant_pairs2)[1] ) %>% rbind(Avr_sharing_table2, .)        
        	}   
	}
}
Info = read_tsv("metadata/link_SGB_species_names.txt")
Avr_sharing_table2 %>% left_join(Info, by="SGB") -> Avr_sharing_table2
write_tsv(Avr_sharing_table2, "Results/SGB_perc_sharing_mothertimepoint.tsv")
} else {Avr_sharing_table2 = read_tsv("Results/SGB_perc_sharing_mothertimepoint.tsv") }

timepoint_order <- c("W2", "M1", "M2", "M3", "M6", "M9", "M12")
Avr_sharing_table2 %>% mutate(Timepoint = factor( Timepoint, levels = timepoint_order) ) -> Avr_sharing_table2
Avr_sharing_table2 %>% mutate(Mother_timepoint= factor( Mother_timepoint, levels = c("P28", "B" ,"M3") ) ) -> Avr_sharing_table2

Avr_sharing_table2 %>% mutate(Timepoint = factor( Timepoint, levels = timepoint_order) ) -> Avr_sharing_table2

Avr_sharing_table2 %>% filter(! is.na(Mother_timepoint) ) %>% ggplot(aes(x = Timepoint, y = Avg_sharing) ) + facet_wrap(~Mother_timepoint) + geom_boxplot(outlier.shape = NA) + geom_sina(alpha=0.5)  + theme_bw() -> Plot2
ggsave("plots_associations/SGB_averagesahring_perTime_perMotherTime.pdf", Plot2)


###Mother-infant sharing rate distribution
read_rds("tmp/Complete_long_sharing.rds") -> Info
Info %>% group_by(next_id_infant, infant_timepoint_categorical) %>% summarise(Sharing_rate = 100*mean(Strain_sharing), Species_n = n() ) %>% ungroup() -> Info
Info %>% filter(Species_n>=3) %>% group_by(infant_timepoint_categorical) %>% mutate(M = median(Sharing_rate) ) %>% ungroup() %>%  ggplot(aes(x=infant_timepoint_categorical, y=Sharing_rate)) + geom_boxplot(outlier.shape = NA) +
  ggforce::geom_sina(aes(size=Species_n), alpha=0.3)  + theme_bw() +
  labs(size = "Number of species used",x = "Timepoint", y = "Family sharing rate (%)") +   theme(legend.position = "bottom") -> Plot3

ggsave("plots_associations/MotherInfant_sharing_perTime.pdf", Plot3, width=5.5, height=4.5)










