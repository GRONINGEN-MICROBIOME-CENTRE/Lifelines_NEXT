library(tidyverse)
library(lmerTest)
library(ggforce)

Long_sharing_file = "tmp/Complete_long_sharing.rds"

Read_all = function(){

        file_list <- list.files("Distance_tables", pattern = "*.rds", full.names = TRUE)
        SGB_list = c()
        for (Entry in file_list){
                SGB_list = c(SGB_list, str_remove(basename(Entry), ".rds") )
        }
        return(list(file_list, SGB_list))
}
source("Common_function.R") #Mother_infant_formatting

Merge_all = function(){
	Read_all() -> Info
	Merged_table = tibble()
	file_list = Info[[1]]
	SGB_list = Info[[2]]
	for(i in 1:length(file_list)){
		SGB_name = SGB_list[i]  
		sp4 = read_rds(file_list[i])
		Mother_infant_formatting(sp4) ->  subset_mother_infant_pairs
		if (dim(subset_mother_infant_pairs)[1] == 0){ next }	
		subset_mother_infant_pairs %>% mutate(SGB =  SGB_name) %>% select(SGB,  next_id_infant, Age_days_infant, infant_timepoint_categorical, mother_timepoint_categorical, Strain_sharing ) %>% mutate(Strain_sharing = ifelse(Strain_sharing=="yes", 1, ifelse(Strain_sharing == "no", 0, NA))) %>% rbind(Merged_table, . ) -> Merged_table
	}
	write_rds(Merged_table,Long_sharing_file)
	return(Merged_table)
}
if (!file.exists(Long_sharing_file) ){
	Merge_all() -> df
} else {
	read_rds(Long_sharing_file) -> df
}


df %>% drop_na()  %>% group_by(next_id_infant, infant_timepoint_categorical, mother_timepoint_categorical) %>% summarise( N = n(), N_share = sum(Strain_sharing), Age_days_infant  ) %>% mutate(Percentage_sharing = N_share/N)  -> Counts
Counts$Percentage_sharing %>% summary()
Counts %>% filter(Percentage_sharing == 0) %>% group_by(infant_timepoint_categorical) %>% summarise(n()) %>% print()
Counts %>% filter(Percentage_sharing == 1) %>%  group_by(infant_timepoint_categorical) %>% summarise(n()) %>% print()

#Applying a filter of at least 5 SGBs per infant
Counts %>% filter(N>=5) %>% ungroup() %>% mutate(sample_id_infant = paste0(next_id_infant,"_", infant_timepoint_categorical) )%>%
  left_join(. , dynamic_prioritized_phenos) %>% left_join(. , prioritized_phenos) -> Counts


Run_linear_model = function(Counts){

	Results_association_strainSharing = tibble()
	Pheno_list = c("Age_days_infant", "mother_birthcard_duration_water_broken_min" , "birth_birthcard_total_duration_delivery_min", "mother_birthcardself_gestational_age_weeks","birth_deliverybirthcard_mode_binary", "birth_deliverybirthcard_place_delivery_simple", "birth_delivery_mode_simple",  "birth_delivery_mode_complex", "mother_birthcard_parity",  "family_pets_any","infant_birthcard_feeding_mode_after_birth", "infant_ffq_feeding_mode_simple",  "infant_ffq_feeding_mode_complex")     
	for ( Pheno in Pheno_list){
		subset_mother_infant_pairs_it = Counts
		if (Pheno %in%  c("mother_birthcard_duration_water_broken_min" , "birth_birthcard_total_duration_delivery_min", "birth_deliverybirthcard_place_delivery_simple" ) ){ subset_mother_infant_pairs_it %>% filter(birth_deliverybirthcard_mode_binary == "VG") -> subset_mother_infant_pairs_it }
		subset_mother_infant_pairs_it %>% filter(!is.na(!!sym(Pheno))) -> subset_mother_infant_pairs_it
		print(paste0("Iterating through pheno: ", Pheno))
		#2. build the model
		if (! Pheno %in% c("infant_timepoint_categorical", "Age_days_infant") ){
			Model0 =  as.formula(paste0( "Percentage_sharing ~ (1|next_id_infant) + infant_timepoint_categorical + mother_timepoint_categorical" ))
			Model = as.formula(paste0( "Percentage_sharing ~ (1|next_id_infant) + infant_timepoint_categorical + mother_timepoint_categorical +  ", Pheno ))
		} else { 
			Model0 = as.formula(paste0( "Percentage_sharing ~ mother_timepoint_categorical + (1|next_id_infant)"))
			Model = as.formula(paste0( "Percentage_sharing ~ mother_timepoint_categorical + (1|next_id_infant) +", Pheno)) 
		}
		Result_model <- tryCatch({
			lmer(formula = Model, data = subset_mother_infant_pairs_it)
			}, error = function(e) {
				print(paste0("Model fitting failed for phenotype: ", Pheno, " with error: ", e$message))
				return(NA)
			}) 
		if (class(Result_model) == "logical") { next }
     		#Check if convergence
     		Convergence = summary(Result_model)$optinfo$conv$lme4
     		if (length(Convergence) == 0 ){ Convergence = "yes"} else { Convergence = "no" }
     
     		as.data.frame(summary(Result_model)$coefficients) -> Result_coefficients
     		rownames(Result_coefficients)[grepl(Pheno, rownames(Result_coefficients))] -> Rows_to_take
   		P_anova = NA
     		if (length(Rows_to_take) > 2){ 
       			Result_model0 <- tryCatch({
         			lmer(formula = Model0,data = subset_mother_infant_pairs_it)
       				}, error = function(e) {
         			print(paste0("Model fitting failed for phenotype: ", Pheno, " with error: ", e$message))
         			return(NA)
       				}) 
       			if (class(Result_model0) == "logical") { next}
       			anova(Result_model,Result_model0) -> P_anova
       			P_anova$`Pr(>Chisq)`[2] -> P_anova
     		}
     		Result_coefficients[Rows_to_take,] %>% rownames_to_column("Phenotype") %>% as_tibble() %>% mutate(P_ANOVA = P_anova, Overall_name = Pheno, Model_converge = Convergence) %>% rbind(Results_association_strainSharing, . ) -> Results_association_strainSharing
 	}
	write_tsv(Results_association_strainSharing, "Results/OverallSharing_summaryLMERTEST_COV_infanttime+mothertime+nextid.tsv")
	return(Results_association_strainSharing)
}
Output = "Results/OverallSharing_summaryLMERTEST_COV_infanttime+mothertime+nextid.tsv"
if (!file.exists(Output) ){
	Results = Run_linear_model(Counts)
} else{
	Results = read_tsv(Output)
}
Results  %>% mutate(FDR_all = p.adjust(`Pr(>|t|)`, "fdr")) -> Summary_statistics
Summary_statistics %>% arrange(`Pr(>|t|)`) -> Summary_statistics
print(Summary_statistics)
write_tsv(Summary_statistics, Output) 
Counts %>% mutate(Percentage_sharing = 100*Percentage_sharing) %>% mutate(Time = ifelse(infant_timepoint_categorical %in% c("W2", "M1", "M2", "M3"), "early", "late" )) -> Counts
print(Counts)
##Make some plots
print("Making plots")
#Top association: Age
print("Age")
Counts %>% mutate(Number_SGB_in_both = N ) %>% ggplot(aes( x = infant_timepoint_categorical, y = Percentage_sharing  )  ) + geom_boxplot(outlier.shape = NA) + geom_sina(alpha=0.2, aes(size=Number_SGB_in_both)) + theme_bw()  -> Plot_timepoint
ggsave("plots_associations/Sharingrates_vs_timepoint.pdf", Plot_timepoint)
#Second top association (not significant) delivery mdoe
print("Delivery")
Counts %>% mutate(Number_SGB_in_both = N ) %>% filter(!is.na(birth_deliverybirthcard_mode_binary)) %>%group_by(birth_deliverybirthcard_mode_binary, Time) %>% mutate(M=median(Percentage_sharing)) %>% ungroup()  %>%  ggplot(aes( x= birth_deliverybirthcard_mode_binary, y = Percentage_sharing) ) + geom_boxplot(outlier.shape = NA) + geom_sina(alpha=0.2, aes(size=Number_SGB_in_both) ) + theme_bw()+geom_point(aes(x=birth_deliverybirthcard_mode_binary, y=M), col="#f45b69", shape=18, size=4 ) + labs(size = "Number of species used",x = "infant timepoint categorical", y = "Mother-infant sharing rate (%)") + theme(legend.position = "bottom") + facet_wrap(~Time)  -> Plot_deliverymode
ggsave("plots_associations/Sharingrates_vs_deliverymode.pdf" ,Plot_deliverymode, width=4, height=2 )
#Association feeding mode
print("Fedding")
Counts %>% mutate(Number_SGB_in_both = N ) %>% filter(!is.na(infant_ffq_feeding_mode_complex)) %>% group_by(infant_ffq_feeding_mode_complex, Time) %>% mutate(M=median(Percentage_sharing)) %>% ungroup()  %>% ggplot(aes( x= infant_ffq_feeding_mode_complex, y = Percentage_sharing) ) + geom_boxplot(outlier.shape = NA) + geom_sina(alpha=0.2, aes(size=Number_SGB_in_both) ) + theme_bw() +geom_point(aes(x=infant_ffq_feeding_mode_complex, y=M), col="#f45b69", shape=18, size=4 ) + labs(size = "Number of species used",x = "Feeding Mode", y = "Mother-infant sharing rate (%)") + theme(legend.position = "bottom") + facet_wrap(~Time)  -> Plot_feedingmode
ggsave("plots_associations/Sharingrates_vs_feedingmode.pdf" ,Plot_feedingmode, width=4, height=2 )



print("Sharing rate in early timepoints")
Counts %>% filter(infant_timepoint_categorical %in% c("W2", "M1")) %>%  summarise(median(Percentage_sharing)) %>% print()
Counts %>% filter(infant_timepoint_categorical %in% c("M12")) %>%  summarise(median(Percentage_sharing)) %>% print()










