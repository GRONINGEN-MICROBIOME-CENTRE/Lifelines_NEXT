################ Cazyme preparation ################

library(tidyverse)


#### Load functions ####
Prepare_cazymes = function(File = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes/NEXT_cayman_infant.txt", Prevalence_min = 0.1, Transform = "log" ){
  caz = read_tsv(File)
  caz %>% filter(! feature %in% c("total_reads", "filtered_reads")) -> caz
  #Get unique counts + ambigous counts normalized
  caz %>% select(ID, feature, combined_rpkm ) %>% spread(key = feature, value = combined_rpkm) -> caz
  #Change ID
  caz$ID = caz$ID 
  #Get summary
  summary_as_df <- function(x) {
    result <- summary(x)
    data.frame(stat = names(result), value = as.vector(result))
  }
  summary_df <- caz %>%
    select(-ID) %>%
    purrr::map(summary_as_df) %>%
    bind_rows(.id = "Cazyme") %>% spread(stat, value) %>% as_tibble()
  summary_df = summary_df %>% mutate( `NA's` = ifelse( is.na(`NA's`), 0,  `NA's` ), N = dim(caz)[1], Present= N - `NA's`)
  summary_df = summary_df %>% mutate(Present_perc = Present/N ) %>% mutate(Keep = ifelse(Present_perc>=Prevalence_min, T, F ) )
  #Make NA into 0
  caz %>% mutate_all(~ replace(., is.na(.), 0)) -> caz
  #Remove Cazymes with high NA prop
  caz %>% select(c("ID", summary_df$Cazyme[summary_df$Keep==T] )) -> caz_anal
  if (Transform == "log"){
    Log_transform(caz_anal)  ->caz_anal_tranf
  } else if (Transform == "clr"){
    CLR_transformation(caz_anal)  -> caz_anal_tranf
  } else { caz_anal_tranf = NULL}
  
  return(list("caz" = caz, "caz_filtered" = caz_anal, "caz_tranf"=caz_anal_tranf  ,"Summary_caz" = summary_df ) )
}



# Data loading and transformations ----------------------------------------
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes")
metadata<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"

caz_tables = Prepare_cazymes(File = "NEXT_cayman_infant.txt", Prevalence_min = 0.3, Transform = "NULL" )
caz_filtered = as.data.frame (caz_tables[["caz_filtered"]])
row.names(caz_filtered)<-caz_filtered$ID
caz_filtered$ID<-NULL


## Loading microbiome data 
metaphlan = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
metaphlan_clr = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
metaphlan_infants = metaphlan[match(as.character(metadata_infants$NG_ID),rownames(metaphlan)),grep("t__",colnames(metaphlan))]
colnames(metaphlan_infants) = sub(".*s__","",colnames(metaphlan_infants))
metaphlan_infants = metaphlan_infants[,colnames(metaphlan_clr)] # Filtering to only the columns in the filtered file 
metaphlan_infants <- metaphlan_infants[rownames(metaphlan_infants) %in% rownames(caz_filtered), ]

# Data transformations 
cazyme_alr = do_clr_externalWeighting(caz_filtered,metaphlan_infants)
cazyme_alr<-as.data.frame(cazyme_alr)

# Write files 

write.table(cazyme_alr, "LLNEXT_cazyme_ALR_transformed_prev_30_07_07_2025.txt", sep = "\t")
write.table(caz_filtered, "LLNEXT_cazyme_raw_filtered_prev_30_07_07_2025.txt", sep = "\t")


