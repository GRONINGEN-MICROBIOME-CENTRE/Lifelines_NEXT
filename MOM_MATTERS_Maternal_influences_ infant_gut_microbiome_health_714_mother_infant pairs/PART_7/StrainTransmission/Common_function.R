metadata<-read.delim("metadata/LLNEXT_metadata_03_01_2024.txt")
cross_phenos<-read.delim("metadata/masterfile_cross_sectional_2023_11_15.txt")
prioritized_phenos <- cross_phenos %>%
  select(next_id_infant,
        mother_birthcardself_gestational_age_weeks, 
         birth_deliverybirthcard_mode_binary,
         birth_deliverybirthcard_place_delivery_simple,
         birth_delivery_mode_simple,
         birth_delivery_mode_complex,
         mother_birthcard_parity,
         family_pets_any, 
         infant_birthcard_feeding_mode_after_birth,mother_birthcard_duration_water_broken_min, birth_birthcard_total_duration_delivery_min) %>% as_tibble() %>% mutate(next_id_infant= as.character(next_id_infant))
dynamic_phenos<-read.delim("metadata/masterfile_longitudinal_2023_09_29.txt")
prioritized_phenos %>% mutate(mother_birthcard_parity = ifelse(mother_birthcard_parity == 0, 0, ifelse(mother_birthcard_parity==1, 1, ifelse(mother_birthcard_parity>1, 2, NA)))) -> prioritized_phenos

dynamic_prioritized_phenos<-dynamic_phenos %>%
  select(next_id_infant,
         SAMPLE_ID,
         infant_ffq_feeding_mode_simple,
         infant_ffq_feeding_mode_complex) %>% as_tibble() %>% mutate(next_id_infant= as.character(next_id_infant))

names (dynamic_prioritized_phenos)[2]<-"sample_id_infant"



Mother_infant_formatting = function(sp4, select_mother = T ){
  sp5<-sp4 %>% filter(related == "related")
  mother_infant_pairs <- sp5 %>%  filter((Type_1 == "mother" & Type_2 == "infant") | (Type_1 == "infant" & Type_2 == "mother"))
  #Just continue if there are enough mother_infant pairs
  if (dim(mother_infant_pairs)[1] <= 2) { return(tibble()) }
  
  ######################
  ####Formatting
  ######################
  # Define whether mother's timepoint is pre or post-birth
  mother_infant_pairs <- mother_infant_pairs %>% mutate(mother_timeframe = case_when( (Type_1 == "mother" & (Timepoint_categorical_1 %in% c("P12", "P28", "B"))) | 
                                                                                        (Type_2 == "mother" & (Timepoint_categorical_2 %in% c("P12", "P28", "B"))) ~ "pre_birth",
                                                                                      (Timepoint_categorical_1 %in% c("M1", "M2", "M3")) |   (Timepoint_categorical_2 %in% c("M1", "M2", "M3")) ~ "post_birth", TRUE ~ NA_character_  ))
  
  # Define exact timepoint categorical of infant 
  mother_infant_pairs %>% mutate(infant_timepoint_categorical = ifelse( Type_1 == "infant",  Timepoint_categorical_1,  ifelse( Type_2 == "infant",  Timepoint_categorical_2, NA ) ) ) -> mother_infant_pairs
  # Define exact timepoint categorical of mother 
  mother_infant_pairs %>% mutate(mother_timepoint_categorical = ifelse( Type_1 == "mother", Timepoint_categorical_1, ifelse( Type_2 == "mother", mother_infant_pairs$Timepoint_categorical_2, NA  ) ) ) -> mother_infant_pairs
  
  
  # Define early or late timepoint  of infant 
  mother_infant_pairs %>% mutate( infant_early_late_timepoint = ifelse( infant_timepoint_categorical %in% c("W2", "M1", "M2", "M3"), "early", ifelse( infant_timepoint_categorical %in% c("M6", "M9", "M12"), "late", NA ) ) ) -> mother_infant_pairs 
  
  timepoint_order <- c("W2", "M1", "M2", "M3", "M6", "M9", "M12")
  mother_infant_pairs %>% mutate(infant_timepoint_categorical = factor( infant_timepoint_categorical, levels = timepoint_order) ) -> mother_infant_pairs
  timepoint_order_mothers<-c("P12", "P28", "B", "M1", "M2", "M3")
  mother_infant_pairs %>% mutate(mother_timepoint_categorical  = factor(mother_timepoint_categorical, levels = timepoint_order_mothers ) ) -> mother_infant_pairs
  
  # For strain sharing with mother and different infant timepoints take into account only the birth timepoint and if not available the P28 or the M3 timepoint
  if (select_mother == T) {
  subset_mother_infant_pairs =  mother_infant_pairs %>%  group_by(FAMILY_1) %>%
    mutate(  selected_timepoint = ifelse("B" %in% mother_timepoint_categorical, "B", ifelse("P28" %in% mother_timepoint_categorical, "P28", NA_character_)) ) %>%
    #mutate(  selected_timepoint = ifelse("B" %in% mother_timepoint_categorical, "B", ifelse("P28" %in% mother_timepoint_categorical, "P28", ifelse("M3" %in% mother_timepoint_categorical, "M3", NA_character_))) ) %>%
    ungroup() %>% filter(!is.na(selected_timepoint)) %>% filter( as.character(selected_timepoint) == mother_timepoint_categorical)
 } else {
   subset_mother_infant_pairs =  mother_infant_pairs
 }  

  #############################################################
  # Merging with phenotypes (first pulling out the infant ID)
  ##############################################################
  subset_mother_infant_pairs %>% mutate(next_id_infant = ifelse( Type_1 == "infant", NEXT_ID_short_1, ifelse( Type_2 == "infant", NEXT_ID_short_2, NA)) ) %>% mutate(next_id_infant= as.character(next_id_infant))  -> subset_mother_infant_pairs
  subset_mother_infant_pairs<-left_join(subset_mother_infant_pairs, prioritized_phenos, by = "next_id_infant")
  
  subset_mother_infant_pairs %>% mutate(sample_id_infant = paste0(next_id_infant,"_", infant_timepoint_categorical) )  %>% mutate(next_id_infant= as.character(next_id_infant))  -> subset_mother_infant_pairs
  subset_mother_infant_pairs<-left_join(subset_mother_infant_pairs, dynamic_prioritized_phenos)
  
  #releveling
  subset_mother_infant_pairs %>% mutate (infant_birthcard_feeding_mode_after_birth = factor( infant_birthcard_feeding_mode_after_birth, levels = c("BF", "MF", "FF") ) ) -> subset_mother_infant_pairs
  subset_mother_infant_pairs %>% mutate(infant_ffq_feeding_mode_complex = factor(infant_ffq_feeding_mode_complex, levels = c("BF", "MF", "FF")) ) -> subset_mother_infant_pairs
  
  subset_mother_infant_pairs %>% mutate( Age_days_infant = ifelse(Type_1 == "infant", exact_age_months_at_collection_1,  ifelse(Type_2 == "infant", exact_age_months_at_collection_2, NA ) ) ) ->  subset_mother_infant_pairs
  
  return(subset_mother_infant_pairs)
  
} 


Plot_stuff = function(subset_mother_infant_pairs, Phenotype ,out_name, Plot_type="Bar",Do_facet = F ){
        if (Plot_type =="Box"){
                subset_mother_infant_pairs %>% ggplot(aes_string(x="Strain_sharing", y=Phenotype)) +  geom_boxplot(outlier.shape = NA) + ggforce::geom_sina(alpha=0.5) + theme_bw() -> Plot
        }else if( Plot_type == "Bar"){
                if ( class( as_vector(subset_mother_infant_pairs[,Phenotype])) == "integer") { subset_mother_infant_pairs[,Phenotype] = as.factor(as_vector(subset_mother_infant_pairs[,Phenotype])) }


		if (Do_facet == F){
                counts <- subset_mother_infant_pairs %>%
                        count(!!sym(Phenotype), Strain_sharing) %>%
                        group_by(!!sym(Phenotype)) %>%
                        mutate(proportion = n / sum(n)) %>% drop_na()
		} else{
			counts <- subset_mother_infant_pairs %>%
   				 count(!!sym(Phenotype), Strain_sharing, mother_timepoint_categorical) %>%
  				  group_by(!!sym(Phenotype), mother_timepoint_categorical ) %>%
   				 mutate(proportion = n / sum(n)) %>% drop_na()
		}
                ggplot(counts, aes_string(x = Phenotype, y = "proportion", fill = "Strain_sharing" )) +
                geom_bar(stat = "identity", position = "stack") +
                geom_text(aes(label = n, y = proportion), col="white", vjust = -0.5, size = 5) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values= c("#243b74", "#cb343f")) + theme_bw() -> Plot

		if (Do_facet == T){ Plot + facet_wrap( ~mother_timepoint_categorical) -> Plot }

        }
        ggsave(out_name, Plot, width=4, height=2)
}


Fit_glmer = function( Model, For_analysis){
	Problem_fit = F
	glmer(Model, For_analysis, family=binomial() ) -> Res 
	summary(Res)$coefficients %>% as.data.frame() -> Res_info
	
	N_param = length(Res_info$`Pr(>|z|)`)
	if (N_param > 2){
		threshold = N_param -2
	} else { threshold=N_param }
	if ( length( unique( Res_info$`Pr(>|z|)` ) ) < threshold  ) { 
		glmer(Model, For_analysis, family=binomial(),control= glmerControl(optimizer="nlminbwrap") ) -> Res 
        	summary(Res)$coefficients %>% as.data.frame() -> Res_info 
        	if ( length( unique( Res_info$`Pr(>|z|)` ) ) <= threshold ) { 
                	glmer(Model, For_analysis, family=binomial(),control= glmerControl(optimizer="bobyqa") ) -> Res 
                	summary(Res)$coefficients %>% as.data.frame() -> Res_info
                	if ( length( unique( Res_info$`Pr(>|z|)` ) ) < threshold  ) { print(Res_info);  print("Error in Pvalues...") ; Problem_fit = T }
        	}   
	} 
	return(list(Res, Problem_fit))
} 
