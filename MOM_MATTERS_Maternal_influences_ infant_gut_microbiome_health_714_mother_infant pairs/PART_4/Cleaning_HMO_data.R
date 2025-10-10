############ CLEANING HMO DATA #################################


# Load metadata and phenotypes
hmo<-read.delim("230810_HMO_data_Molgenis_milk_groups_ugml_cleaned_n1542.txt")
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
linkage<-read.delim("~/Desktop/LLNEXT/Analysis/linkage/20230414_Linkingfile_fam_longID_formatted_TS_SB.txt")
names (hmo)[1]<-"next_id_mother"
names (hmo)[5]<-"sample_id_mother"

# Modify HMO data to get infant SAMPLE_ID's 
hmo_linkage<-left_join(hmo, linkage)
hmo_linkage$SAMPLE_ID<-paste0(hmo_linkage$next_id_infant,"_", hmo_linkage$time_point )

# Merging files and selecting infant relevant data
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"
metadata_infants$BATCH_NUMBER<-as.factor(metadata_infants$BATCH_NUMBER)

# For this specific data remove duplicates for categorical timepoints 
metadata_infants <- metadata_infants %>%
  distinct(SAMPLE_ID, .keep_all = TRUE) 

# Merging HMO with microbiome metadata 
hmo_met<-left_join(metadata_infants, hmo_linkage)


names (hmo_met)[2]<-"NEXT_ID"
hmo_met$NEXT_ID=as.factor(hmo_met$NEXT_ID)
row.names(hmo_met)<-hmo_met$NG_ID

hmo_met$twin_pair=NULL
hmo_met$infant_relations=NULL
hmo_met$Type=NULL
hmo_met$Sequenced=NULL
hmo_met$sequence_control=NULL
hmo_met$isolation_control=NULL
hmo_met$metaphlan4_unclassified_high_contaminants_factor_75=NULL
hmo_met$NEXT_ID_mother_simple=NULL
hmo_met$pregnancy=NULL
hmo_met$sample_id_mother=NULL
hmo_met$time_point=NULL
hmo_met$mother_milk_HMO_UPLC_ID=NULL
hmo_met <- hmo_met[!(hmo_met$Timepoint_categorical %in% c("M6", "M9", "M12")), ]

hmo_met$mother_milk_HMO_milk_group <- factor(hmo_met$mother_milk_HMO_milk_group, levels = c("Le+Se+", "Le+Se-", "Le-Se+", "Le-Se-"))

# Remove mixed feeding and only do breastfeeding 
dynamic_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_longitudinal_2023_09_29.txt")
feeding <-dynamic_phenotypes[, c("SAMPLE_ID", "infant_ffq_feeding_mode_simple")]

hmo_met_feed<-left_join(hmo_met, feeding)
hmo_met_feed <- hmo_met_feed[(hmo_met_feed$infant_ffq_feeding_mode_simple %in% c("excl_BF")), ]
hmo_met_feed <- hmo_met_feed[!is.na(hmo_met_feed$mother_milk_HMO_milk_group), ]
row.names(hmo_met_feed)<-hmo_met_feed$NG_ID
num_unique_infants <- length(unique(hmo_met_feed$NEXT_ID))
# 277

write.table(hmo_met_feed, "HMOs_matched_infant_MGS_breastfeeding_only_early_timepoints.txt", sep = "\t", row.names = F)

