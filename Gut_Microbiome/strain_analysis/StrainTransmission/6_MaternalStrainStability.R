library(tidyverse)


#mother consistency during pregnancy

#Get RDS file from command line
args <- commandArgs(trailingOnly = TRUE)
file = args[1]
sp4 = read_rds(file)
#Get SGB name
SGB_name = str_remove(basename(file), "\\.rds") 
#Keep only longitudinal comparison in mothers
sp4 %>% filter(related == "related" & (Type_1 == "mother" & Type_2 == "mother") ) %>% filter(NEXT_ID_short_1 == NEXT_ID_short_2 ) -> motherinfo
motherinfo %>% filter(! Timepoint_categorical_1 == Timepoint_categorical_2 ) -> motherinfo
#Calculate exact time differences between timepoints
motherinfo %>% mutate(TimeDif = abs(exact_age_months_at_collection_1 - exact_age_months_at_collection_2) ) -> motherinfo
#Preare df
motherinfo %>% select(NEXT_ID_short_1, NEXT_ID_short_2, Timepoint_categorical_1, Timepoint_categorical_2, Strain_sharing, nGD, TimeDif)  -> motherinfo
motherinfo %>% mutate(Combined_Timepoint = paste0(pmin(Timepoint_categorical_1, Timepoint_categorical_2), "-", pmax(Timepoint_categorical_1, Timepoint_categorical_2)))  -> motherinfo
Order_timepoints = c("P12" ,"P28", "B", "M1", "M2", "M3")
motherinfo %>% mutate(Timepoint_categorical_1 = factor(Timepoint_categorical_1, levels=Order_timepoints ), Timepoint_categorical_2 = factor(Timepoint_categorical_2, levels=Order_timepoints ) ) -> motherinfo
#Create two groups comparisons between pregnancy time poitns, and after birth timpeoints
Pre = c("P12", "P28", "B")
Post = c("B", "M1", "M2", "M3")
motherinfo %>% mutate( Time_period = ifelse( Timepoint_categorical_1 %in% Pre & Timepoint_categorical_2 %in% Pre, "Pregnancy", ifelse(Timepoint_categorical_1 %in% Post & Timepoint_categorical_2 %in% Post, "After", "Mixed") ) ) -> motherinfo
#Compare if the probability of strain sharing is different before birth than after birht. Control for natural strain replacement due to time
motherinfo %>% mutate(Strain_sharing = factor(Strain_sharing, levels=c("no", "yes") )) %>% filter(!  Time_period == "Mixed") -> motherinfo
motherinfo %>% glm( Strain_sharing ~ TimeDif + Time_period , ., family=binomial() ) %>% summary() -> Model
#Save results
Model$coefficients %>% as.data.frame() %>% rownames_to_column("Coefficient") %>% filter(Coefficient == "Time_periodPregnancy") %>% as_tibble() %>% 
  mutate(SGB = SGB_name, Total_n = dim(motherinfo)[1], Mother_n = length(unique(motherinfo$NEXT_ID_short_1)), Pregnancy_group=dim(filter(motherinfo, Time_period == "Pregnancy"))[1], PostBirth_group =dim(filter(motherinfo, Time_period == "After"))[1],
         Pregnancy_group_mothers = length(unique(filter(motherinfo, Time_period == "Pregnancy")$NEXT_ID_short_1 )), PostBirth_group_mothers = length(unique(filter(motherinfo, Time_period == "After")$NEXT_ID_short_1 ))   ) -> Result

write_tsv(Result, paste0("Results/MotherConsistency/", SGB_name, ".tsv" ) )




