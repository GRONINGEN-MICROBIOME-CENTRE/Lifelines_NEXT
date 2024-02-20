setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/02.SCRIPTS/01.METADATA_HARMONIZATION/')

#############################################################
# Here we prepapre metadats for further analysis
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################
library(readxl)
library(reshape2)

##############################
# Input data & data transformation
##############################
Chiliadal_sequencing <- read_xlsx('~/Desktop/Projects_2021/NEXT_virome/06.Sequencing/Chiliadal_virome_sequencing_batch_I.xlsx')
Chiliadal_sequencing2 <- read_xlsx('~/Desktop/Projects_2021/NEXT_virome/06.Sequencing/Chiliadal_virome_sequencing_batch_II.xlsx')
Chiliadal_sequencing2 <- Chiliadal_sequencing2[!is.na(Chiliadal_sequencing2$NEXT_ID),]

Chiliadal_sequencing <- rbind(Chiliadal_sequencing, Chiliadal_sequencing2)
Chiliadal_sequencing[is.na(Chiliadal_sequencing$Sequencing_comments),]$Sequencing_comments <- 'Sequenced'
Chiliadal_sequenced_samples <- Chiliadal_sequencing[Chiliadal_sequencing$NEXT_ID!='BLANK' & Chiliadal_sequencing$Sequencing_comments=='Sequenced',]
# removing unnecessary columns
Chiliadal_sequenced_samples$ALIQUOTE_NR <- NULL
Chiliadal_sequenced_samples$...15 <- NULL
Chiliadal_sequenced_samples$ADDITION <- NULL
Chiliadal_sequenced_samples$Sequencing_ID...20 <- NULL
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$DNA_CONC=='too low',]$DNA_CONC <- '0.001'
Chiliadal_sequenced_samples$DNA_CONC <- as.numeric(Chiliadal_sequenced_samples$DNA_CONC)
# renaming columns for convenience
colnames(Chiliadal_sequenced_samples)[16] <- 'Sequencing_ID'
colnames(Chiliadal_sequenced_samples)[c(1:3,7,16,17)] <- paste0(colnames(Chiliadal_sequenced_samples)[c(1:3,7,16,17)], '_VLP')

# steps to harmonize VLP metadata before merging with the metadata of Big Gut NEXT paper
Chiliadal_sequenced_samples$NEXT_ID_VLP <- paste0('LLNEXT', Chiliadal_sequenced_samples$NEXT_ID_VLP)

# changing the NEXT ID because it was incorrectly assigned during collection (father ID was used instead child ID)
# this info is coming from Trishla's script: https://umcgonline.sharepoint.com/:u:/r/sites/LifelinesNEXTGutMicrobiome/Shared%20Documents/NEXT_metagenomics/QC/QC_script_MGS_TS.R?csf=1&web=1&e=tgeDh7
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$NEXT_ID_VLP=='LLNEXT011062',]$NEXT_ID_VLP <- 'LLNEXT011053'
Chiliadal_sequenced_samples$Universal_ID <- paste0(Chiliadal_sequenced_samples$NEXT_ID_VLP, '_', Chiliadal_sequenced_samples$Timepoint_VLP)

# metadata from Big Gut Next Paper (last date of editing: October 25, 2022)
bgnp_metadata <- read.table('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/01.METADATA/final_metadata_NEXT_25_10_2022.txt', sep='\t', header=T)

# due to the reassignment of Timepoint in the Big Gut Next paper based on exact age at the time of collection,
# we first need to reconstruct the original Timepoints to merge with VLP metadata
bgnp_metadata$Timepoint_regular <- substr(bgnp_metadata$NG_ID, 2,3)
bgnp_metadata[bgnp_metadata$Timepoint_regular=='01',]$Timepoint_regular <- 'M1'
bgnp_metadata[bgnp_metadata$Timepoint_regular=='02',]$Timepoint_regular <- 'M2'
bgnp_metadata[bgnp_metadata$Timepoint_regular=='03',]$Timepoint_regular <- 'M3'
bgnp_metadata[bgnp_metadata$Timepoint_regular=='06',]$Timepoint_regular <- 'M6'
bgnp_metadata[bgnp_metadata$Timepoint_regular=='09',]$Timepoint_regular <- 'M9'
bgnp_metadata[bgnp_metadata$Timepoint_regular=='12',]$Timepoint_regular <- 'M12'
bgnp_metadata[bgnp_metadata$Timepoint_regular=='MB',]$Timepoint_regular <- 'B'
bgnp_metadata[bgnp_metadata$Timepoint_regular=='P3',]$Timepoint_regular <- 'P12'
bgnp_metadata[bgnp_metadata$Timepoint_regular=='P7',]$Timepoint_regular <- 'P28'

# due to the reassignment of Type in the Big Gut Next paper based on the sample mix up during collection,
# we need to make a list for VLP samples, where the same reassignment might be needed
bgnp_metadata$type_id_derived <- substr(bgnp_metadata$NG_ID, 1,1)
bgnp_metadata[bgnp_metadata$type_id_derived%in% c('A', 'B'),]$type_id_derived <- 'mother'
bgnp_metadata[bgnp_metadata$type_id_derived%in% c('C', 'D', 'Y'),]$type_id_derived <- 'infant'

bgnp_metadata$check_VLP_for_mixup <- FALSE

for (i in bgnp_metadata[paste0(bgnp_metadata$Type, '_', bgnp_metadata$Timepoint)!=paste0(bgnp_metadata$type_id_derived, '_', bgnp_metadata$Timepoint_regular),]$NG_ID) {
  if (bgnp_metadata[bgnp_metadata$NG_ID==i,]$Type!=bgnp_metadata[bgnp_metadata$NG_ID==i,]$type_id_derived) {
    bgnp_metadata[bgnp_metadata$NG_ID==i,]$Timepoint_regular <- bgnp_metadata[bgnp_metadata$NG_ID==i,]$Timepoint_categorical
    bgnp_metadata[bgnp_metadata$NG_ID==i,]$check_VLP_for_mixup <- TRUE
  }
}
bgnp_metadata$type_id_derived <- NULL

# renaming the column names
colnames(bgnp_metadata)[5] <- 'STATUS'
bgnp_metadata$STATUS <- 'Sequenced'

# some samples from Big Gut Next paper either failed to be sequenced, or were excluded based on the low
# read depth, or being an outlier at the PCA
excluded_bgnp <- read_xlsx('~/Desktop/Projects_2021/NEXT_virome/08.Big_gut_NEXT_paper_data/Removed_MGS_samples/QC_MGS_TS.xlsx', sheet = 3)
excluded_bgnp$Timepoint_categorical <- excluded_bgnp$Timepoint
excluded_bgnp$Timepoint_regular <- excluded_bgnp$Timepoint
excluded_bgnp$Timepoint <- NULL
excluded_bgnp$SAMPLE_ID <- paste0(excluded_bgnp$NEXT_ID, '_', excluded_bgnp$Timepoint_regular)
excluded_bgnp$check_VLP_for_mixup <- ifelse(excluded_bgnp$STATUS=='FAILED' | excluded_bgnp$STATUS=='LOW_RD', F, T)

# extending metadata from Big Gut Next paper with samples that were excluded during QC:
bgnp_with_excluded <- rbind(bgnp_metadata, excluded_bgnp)

# harmonization:
bgnp_with_excluded[bgnp_with_excluded$Type=='mother',]$Type <- 'M'
bgnp_with_excluded[bgnp_with_excluded$Type=='infant',]$Type <- 'K'
bgnp_with_excluded$Universal_ID <- paste0(bgnp_with_excluded$NEXT_ID, '_', bgnp_with_excluded$Timepoint_regular)

# adding '_MGS' to some column names for convenience:
colnames(bgnp_with_excluded)[c(2,3,5:11,13:17,22)] <- paste0( colnames(bgnp_with_excluded)[c(2,3,5:11,13:17,22)], '_', 'MGS' )

# the full basic metadata for sequenced VLP samples:
metadata_chiliadal <- merge(Chiliadal_sequenced_samples, bgnp_with_excluded, by='Universal_ID', all.x = T)

metadata_chiliadal_w_imputation_all <- merge(Chiliadal_sequenced_samples, bgnp_with_excluded, by='Universal_ID', all = T)

# future steps: add aliquot id (if necessary? check if this could potentially affect the exact age)
# add the weight of the fecal aliquot for VLP
# add the aliquoting comments and isolation comments if any 

# creating a metadata for families that have at least one VLP sample and MGS samples
families_to_keep <- list()
for (i in unique(metadata_chiliadal_w_imputation_all$FAMILY)) {
  if ( sum(!is.na(metadata_chiliadal_w_imputation_all[metadata_chiliadal_w_imputation_all$FAMILY==i,]$Sequencing_ID_VLP)) != 0) {
    families_to_keep[i] <- i
  }
}

metadata_chiliadal_w_imputation <- metadata_chiliadal_w_imputation_all[metadata_chiliadal_w_imputation_all$FAMILY %in% unname(unlist(families_to_keep)),]

# spliting the metadata containing VLP samples and MGS samples:

metadata_vertical_vlp <- metadata_chiliadal[,c("Universal_ID", "NEXT_ID_VLP","Type_VLP", "Timepoint_VLP", "FAMILY", "Sequencing_ID_VLP")]
colnames(metadata_vertical_vlp)[c(2:6)] <- c('NEXT_ID', 'Type', 'Timepoint', 'FAMILY', 'Sequencing_ID')
metadata_vertical_vlp$Timepoint <- paste0(metadata_vertical_vlp$Timepoint, '_VLP')

# also excluding failed samples:
metadata_vertical_mgs <- bgnp_with_excluded[bgnp_with_excluded$STATUS_MGS!='FAILED' & bgnp_with_excluded$FAMILY %in% metadata_vertical_vlp$FAMILY,c("Universal_ID", "NEXT_ID_MGS", "Type_MGS", "Timepoint_regular", "FAMILY", "NG_ID")]
colnames(metadata_vertical_mgs)[c(2:6)] <- c('NEXT_ID', 'Type', 'Timepoint', 'FAMILY', 'Sequencing_ID')
metadata_vertical_mgs$Timepoint <- paste0(metadata_vertical_mgs$Timepoint, '_MGS')

metadata_vertical <- rbind(metadata_vertical_vlp, metadata_vertical_mgs)
metadata_vertical$Timepoint <- factor(metadata_vertical$Timepoint, levels=c("P12_VLP", "P28_VLP", "B_VLP", "M1_VLP", "M3_VLP", "M6_VLP", "M12_VLP",
                                                                            "P12_MGS", "P28_MGS", "B_MGS", "W2_MGS", "M1_MGS", "M2_MGS", "M3_MGS", "M6_MGS", "M9_MGS", "M12_MGS"),
                                                                    ordered=T)


##### OUTPUT #####
write.table(metadata_chiliadal, '../../01.METADATA/Chiliadal_metadata_ver_01_06042023.txt', sep='\t', quote=F, row.names=F)

