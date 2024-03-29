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
library(vegan)

# library(lme4)
# library(RLRsim)
# library(lmerTest)
# 
# library(ggplot2)
# library(ggforce)
##############################
# Input data & data transformation
##############################
### VLP datasets' metadata:
## Chiliadal_virome_sequencing_batch_I.xlsx: list of Chiliadal VLP extractions sent for sequencing in August 2022
## Chiliadal_virome_sequencing_batch_II.xlsx: list of Chiliadal VLP extractions sent for sequencing in December 2022
## Samples_BaseClear_pilot_sequencing.xlsx: list of NEXT samples included in the Big Gut NEXT paper that were extracted and sequenced to test the new VLP protocol

### MGS datasets' metadata
## LLNEXT_metadata_03_01_2024.txt: metadata for NEXT samples included in the Big Gut NEXT paper
## LLNEXT_metadata_03_03_2023.txt: metadata for NEXT samples included in the Big Gut NEXT paper that includes samples that were excluded based on Sphingomonas contamination
## QC_MGS_TS.xlsx (sheet 3): all previously excluded from the Big Gut NEXT paper samples
## LLNEXT_metaphlan_4_complete_10_02_2023.txt: merged MetaPhlAn4 output file for NEXT samples included in the Big Gut NEXT paper
## DNA_CONC_ALL_WITH_DUPLICATES_MERGED_15_07_2022_TS.txt: DNA concentration measurements by Qubit received from Novogene sample QC reports


 Chiliadal_sequencing <- read_xlsx('~/Desktop/Projects_2021/NEXT_virome/06.Sequencing/Chiliadal_virome_sequencing_batch_I.xlsx')
 Chiliadal_sequencing2 <- read_xlsx('~/Desktop/Projects_2021/NEXT_virome/06.Sequencing/Chiliadal_virome_sequencing_batch_II.xlsx')
 Chiliadal_sequencing2 <- Chiliadal_sequencing2[!is.na(Chiliadal_sequencing2$NEXT_ID),]
 
 Prechiliadal_sequencing <- read_xlsx('../../01.METADATA/Samples_BaseClear_pilot_sequencing.xlsx', sheet = 2)
 
 Chiliadal_sequencing <- rbind(Chiliadal_sequencing, Chiliadal_sequencing2, Prechiliadal_sequencing)
 Chiliadal_sequencing[is.na(Chiliadal_sequencing$Sequencing_comments),]$Sequencing_comments <- 'Sequenced'
 Chiliadal_sequenced_samples <- Chiliadal_sequencing[Chiliadal_sequencing$Sequencing_comments=='Sequenced',]
 # removing unnecessary columns
 Chiliadal_sequenced_samples$ALIQUOTE_NR <- NULL
 Chiliadal_sequenced_samples$...17 <- NULL
 Chiliadal_sequenced_samples$ADDITION <- NULL
 Chiliadal_sequenced_samples$Sequencing_ID...21 <- NULL
 Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$DNA_CONC=='too low',]$DNA_CONC <- '0.001'
 Chiliadal_sequenced_samples$DNA_CONC <- as.numeric(Chiliadal_sequenced_samples$DNA_CONC)
 # renaming columns for convenience
 colnames(Chiliadal_sequenced_samples)[grep("Sequencing_ID...22", colnames(Chiliadal_sequenced_samples))] <- 'Sequencing_ID'
 # adding "_VLP" to the column names to match with BGNP metadata
 colnames(Chiliadal_sequenced_samples)[c(1:3,7,18,19)] <- paste0(colnames(Chiliadal_sequenced_samples)[c(1:3,7,18,19)], '_VLP')

 # steps to harmonize VLP metadata before merging with the metadata of Big Gut NEXT paper
 Chiliadal_sequenced_samples$NEXT_ID_VLP <- paste0('LLNEXT', Chiliadal_sequenced_samples$NEXT_ID_VLP)

 # changing the NEXT ID because it was incorrectly assigned during collection (father ID was used instead child ID)
 # this info is coming from Trishla's script: https://umcgonline.sharepoint.com/:u:/r/sites/LifelinesNEXTGutMicrobiome/Shared%20Documents/NEXT_metagenomics/QC/QC_script_MGS_TS.R?csf=1&web=1&e=tgeDh7
 Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$NEXT_ID_VLP=='LLNEXT011062',]$NEXT_ID_VLP <- 'LLNEXT011053'
 Chiliadal_sequenced_samples$Universal_ID <- paste0(Chiliadal_sequenced_samples$NEXT_ID_VLP, '_', Chiliadal_sequenced_samples$Timepoint_VLP)

# metadata from Big Gut Next Paper (last date of editing: January 3, 2024)
bgnp_metadata <- read.table('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/01.METADATA/LLNEXT_metadata_03_01_2024.txt', sep='\t', header=T)

# removing Shannon index calculated at s__ + t__ taxa levels
bgnp_metadata$shannon <- NULL

# have to also load previous version of metadta due to the undocumented exclusions from the current
# version of BGNP metadata
early_bgnp <- read.table("../../01.METADATA/LLNEXT_metadata_03_03_2023.txt", sep='\t', header=T)

# adding a column that will encompass inclusion/exclusion from BGNP
early_bgnp$BGNP <- "Included"
early_bgnp[! early_bgnp$NG_ID %in% bgnp_metadata$NG_ID,]$BGNP <- "Excluded: > 75% contamination"
early_bgnp[(! early_bgnp$NG_ID %in% bgnp_metadata$NG_ID) & early_bgnp$metaphlan4_unclassified_with_contaminants < 75,]$BGNP <- "Excluded: Unknown reason"

# adding this comprehensive column allows us to remove non-informative "Sequenced" column:
early_bgnp$Sequenced <- NULL

# some samples from Big Gut Next paper either failed to be sequenced, or were excluded based on the low
# read depth, or being an outlier at the PCA
excluded_bgnp <- read_xlsx('~/Desktop/Projects_2021/NEXT_virome/08.Big_gut_NEXT_paper_data/Removed_MGS_samples/QC_MGS_TS.xlsx', sheet = 3)

excluded_bgnp$BGNP <- paste0("Excluded: ", excluded_bgnp$STATUS)
excluded_bgnp$STATUS <- NULL

# removing and renaming columns that do not allow rbinding
excluded_bgnp$Timepoint_categorical <- excluded_bgnp$Timepoint
excluded_bgnp$Timepoint <- NULL
excluded_bgnp$Timepoint_regular <- NULL
excluded_bgnp$check_VLP_for_mixup <- NULL
colnames(excluded_bgnp)[grep('metaphlan_unknown_perc', colnames(excluded_bgnp))] <- "metaphlan4_unclassified"

# to add columns that are present in early_bgnp, we need metaphlan4 data:
metaphlan <- read.table("../../../08.Big_gut_NEXT_paper_data/LLNEXT_metaphlan_4_complete_10_02_2023.txt", sep='\t', header=T)

metaphlan$NCBI_tax_id <- NULL
colnames(metaphlan) <- substr(colnames(metaphlan), 0, 12) 
row.names(metaphlan) <- metaphlan$clade_name
metaphlan$clade_name <- NULL

# calculating shannong index at SGB level (we will need it later)
strains <- metaphlan[grep("t__",row.names(metaphlan)),]
bac_diversity_shannon <- data.frame(diversity(t(strains), "shannon"))
colnames(bac_diversity_shannon) <- "diversity"

# it is essential to pick either only SGB or species level:
contaminants <- as.data.frame(t(metaphlan[  c( grep('UNCLASSIFIED', row.names(metaphlan)),
                                                grep('Sphingomonas_sp_FARSPH', row.names(metaphlan))[1],
                                                grep('Phyllobacterium_myrsinacearum', row.names(metaphlan))[1]), ] ))

contaminants$sum <- rowSums(contaminants)

# adding columns that are present in early_bgnp
excluded_bgnp$metaphlan4_unclassified <- contaminants$UNCLASSIFIED[match(excluded_bgnp$NG_ID, row.names(contaminants))]
excluded_bgnp$contaminant_1_Sphingomonas_sp_FARSPH <- contaminants$`k__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Sphingomonadales|f__Sphingomonadaceae|g__Sphingomonas|s__Sphingomonas_sp_FARSPH`[match(excluded_bgnp$NG_ID, row.names(contaminants))]
excluded_bgnp$contaminant_2_Phyllobacterium_myrsinacearum <- contaminants$`k__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Rhizobiales|f__Phyllobacteriaceae|g__Phyllobacterium|s__Phyllobacterium_myrsinacearum`[match(excluded_bgnp$NG_ID, row.names(contaminants))]
excluded_bgnp$metaphlan4_unclassified_with_contaminants <- contaminants$sum[match(excluded_bgnp$NG_ID, row.names(contaminants))]
excluded_bgnp$metaphlan4_unclassified_high_contaminants_factor_75 <- ifelse(excluded_bgnp$metaphlan4_unclassified_with_contaminants>=75, "high", "low")

# adding missing data to samples that were excluded early in the pipeline
excluded_bgnp$reads_lost_QC <- 1 - excluded_bgnp$clean_reads_FQ_1/excluded_bgnp$raw_reads_FQ_1

excluded_bgnp[grep('X', excluded_bgnp$NG_ID), "isolation_method"] <- "fsk_vortex_gro" 
excluded_bgnp[grep('X', excluded_bgnp$NG_ID, invert = T), "isolation_method"] <- "fsk_bb_kiel" 

# DNA concentration as measured by Novogene:
dna_conc <- read.table("../../../08.Big_gut_NEXT_paper_data/DNA_CONC_ALL_WITH_DUPLICATES_MERGED_15_07_2022_TS.txt", sep='\t', header = T)
dna_conc <- dna_conc[dna_conc$ID!="Sample Name",]
dna_conc$DNA_CONC <- as.numeric(dna_conc$DNA_CONC)
excluded_bgnp[is.na(excluded_bgnp$dna_conc), "dna_conc"] <- dna_conc$DNA_CONC[match(excluded_bgnp[is.na(excluded_bgnp$dna_conc), ]$NG_ID, dna_conc$ID)]

excluded_bgnp[is.na(excluded_bgnp$NG_ID_short),]$NG_ID_short <- substr(excluded_bgnp[is.na(excluded_bgnp$NG_ID_short),]$NG_ID, 0, 8)

excluded_bgnp$SAMPLE_ID <- paste0(excluded_bgnp$NEXT_ID, '_', excluded_bgnp$Timepoint_categorical)

excluded_bgnp <- excluded_bgnp[,colnames(early_bgnp)]

# ready to combine excluded samples togther with the currently included samples:
bgnp_metadata_all <- rbind(early_bgnp, excluded_bgnp)

# due to the reassignment of Timepoint in the Big Gut Next paper based on exact age at the time of collection,
# we first need to reconstruct the original Timepoints to merge with VLP metadata
bgnp_metadata_all$Timepoint_original <- substr(bgnp_metadata_all$NG_ID, 2,3)
bgnp_metadata_all[bgnp_metadata_all$Timepoint_original=='01',]$Timepoint_original <- 'M1'
bgnp_metadata_all[bgnp_metadata_all$Timepoint_original=='02',]$Timepoint_original <- 'M2'
bgnp_metadata_all[bgnp_metadata_all$Timepoint_original=='03',]$Timepoint_original <- 'M3'
bgnp_metadata_all[bgnp_metadata_all$Timepoint_original=='06',]$Timepoint_original <- 'M6'
bgnp_metadata_all[bgnp_metadata_all$Timepoint_original=='09',]$Timepoint_original <- 'M9'
bgnp_metadata_all[bgnp_metadata_all$Timepoint_original=='12',]$Timepoint_original <- 'M12'
bgnp_metadata_all[bgnp_metadata_all$Timepoint_original=='MB',]$Timepoint_original <- 'B'
bgnp_metadata_all[bgnp_metadata_all$Timepoint_original=='P3',]$Timepoint_original <- 'P12'
bgnp_metadata_all[bgnp_metadata_all$Timepoint_original=='P7',]$Timepoint_original <- 'P28'
bgnp_metadata_all[bgnp_metadata_all$Timepoint_original=='BB',]$Timepoint_original <- 'B'

# creating Universal_ID that includes NEXT_ID and ORIGINAL TIMEPOINT
bgnp_metadata_all$Universal_ID <- paste0(bgnp_metadata_all$NEXT_ID, '_', bgnp_metadata_all$Timepoint_original)

# annotating the isolation & sequencing controls:
bgnp_metadata_all[bgnp_metadata_all$Universal_ID=="LLNEXT005816_M6","isolation_control"] <- "yes"
bgnp_metadata_all[bgnp_metadata_all$Universal_ID=="LLNEXT008093_M9","sequence_control"] <- "yes"

# excluding the duplicated isolation and sequencing controls (failed controls of above):
bgnp_metadata_all <- bgnp_metadata_all[! bgnp_metadata_all$NG_ID %in% c("C06X009601G1", "C09F047214A9"),]

# filling some missing data:
bgnp_metadata_all[is.na(bgnp_metadata_all$isolation_control), "isolation_control"] <- "no"
bgnp_metadata_all[is.na(bgnp_metadata_all$sequence_control), "sequence_control"] <- "no"

# from parsed NEXT pilot kneaddata log:
bgnp_metadata_all[bgnp_metadata_all$NG_ID=="C01X01943A03",c("human_reads_FQ_1", "human_reads_FQ_2")] <- 30548/2

# CW2F002016H7 was presumably corrupted during BBDUK for the current metaphlan4 run, 
# making AdapTr_1 != AdapTr_2, so make sure to use properly cleaned data for any of your analyses
bgnp_metadata_all[bgnp_metadata_all$NG_ID=="CW2F002016H7",]$BGNP <- "EXCLUDE OR RERUN"

# Exact ages at time of collection:

# WAIT FOR THE FINAL FILE FROM TRISHLA
# exact_ages <- read_xlsx("../../../08.Big_gut_NEXT_paper_data/2022_07_11_NEXT_metagenomics_septemeber_freeze_2021_batch_2_update.xlsx")
# exact_ages <- exact_ages[!is.na(exact_ages$exact_age_days_at_collection),]
# exact_ages$SAMPLE_ID <- paste0('LLNEXT', exact_ages$next_id, '_', exact_ages$time_point)
# 
# exact_ages$exact_age_days_at_collection[match(excluded_bgnp$SAMPLE_ID, exact_ages$SAMPLE_ID)]

# Removing some redundancy from columns:
bgnp_metadata_all$raw_reads <- bgnp_metadata_all$raw_reads_FQ_1 + bgnp_metadata_all$raw_reads_FQ_2
bgnp_metadata_all$human_reads <- bgnp_metadata_all$human_reads_FQ_1 + bgnp_metadata_all$human_reads_FQ_2
bgnp_metadata_all$clean_reads <- bgnp_metadata_all$clean_reads_FQ_1 + bgnp_metadata_all$clean_reads_FQ_2

bgnp_metadata_all[,grep('FQ_', colnames(bgnp_metadata_all))] <- NULL

##### PROCEED WITH HARMONIZATION:
# 3. COMPARE TO CURRENT VERSION OF BGNP METADATA ADD SHANNON INDEX BASED ON t__

# arguably, samples that failed sequencing or were outliers in MGS analysis, 
# should also be checked thoroughly in Chiliadal virome
bgnp_metadata_all$check_VLP_for_mixup <- ifelse(bgnp_metadata_all$BGNP=='Included', F, T)

# due to the reassignment of Type in the Big Gut Next paper based on the sample mix up during collection,
# we need to make a list for VLP samples, where the same reassignment might be needed
bgnp_metadata_all$type_id_derived <- substr(bgnp_metadata_all$NG_ID, 1,1)
bgnp_metadata_all[bgnp_metadata_all$type_id_derived%in% c('A', 'B'),]$type_id_derived <- 'mother'
bgnp_metadata_all[bgnp_metadata_all$type_id_derived%in% c('C', 'D', 'Y'),]$type_id_derived <- 'infant'

# NG_IDs of samples were either type or timepoint differ from the one assigned to them during DNA extraction
reassigned_type_or_timepoint <- bgnp_metadata_all[paste0(bgnp_metadata_all$Type, '_', bgnp_metadata_all$Timepoint_categorical)!=paste0(bgnp_metadata_all$type_id_derived, '_', bgnp_metadata_all$Timepoint_original),]$NG_ID

for (i in reassigned_type_or_timepoint) {
  
  # since some timepoints were reassigned due to the exact ages of infants at the moment of stool collection, 
  # we focus only on those were the type of the sample differs
  if (bgnp_metadata_all[bgnp_metadata_all$NG_ID==i,]$Type!=bgnp_metadata_all[bgnp_metadata_all$NG_ID==i,]$type_id_derived) {
    
    bgnp_metadata_all[bgnp_metadata_all$NG_ID==i,]$Timepoint_original <- bgnp_metadata_all[bgnp_metadata_all$NG_ID==i,]$Timepoint_categorical
    bgnp_metadata_all[bgnp_metadata_all$NG_ID==i,]$check_VLP_for_mixup <- TRUE
    
  }
}

bgnp_metadata_all$type_id_derived <- NULL

# reordering the columns:
bgnp_metadata_all <- bgnp_metadata_all[,c("NG_ID", "NEXT_ID", "Type", "FAMILY", 
                                          "BGNP", "raw_reads", "human_reads", "clean_reads",
                                          "reads_lost_QC", "sequence_control", "isolation_control", "dna_conc",
                                          "isolation_method", "NG_ID_short", "exact_age_days_at_collection", "exact_age_months_at_collection",
                                          "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID", "metaphlan4_unclassified",
                                          "contaminant_1_Sphingomonas_sp_FARSPH", "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants", "metaphlan4_unclassified_high_contaminants_factor_75",
                                          "Timepoint_original", "Universal_ID", "check_VLP_for_mixup")]
# ADD MISSING BATCH NUMBERS:
# STILL SOME NA, WAIT FOR THE FILE FROM TRISHLA
bgnp_metadata_all$BATCH_NUMBER <- bgnp_metadata$BATCH_NUMBER[match(bgnp_metadata_all$NG_ID, bgnp_metadata$NG_ID)]

# ADD MISSING INFANT RELATIONS:
bgnp_metadata_all$infant_relations <- bgnp_metadata$infant_relations[match(bgnp_metadata_all$NEXT_ID, bgnp_metadata$NEXT_ID)]
# FAM0875 and FAM0259 both had singleton infants
bgnp_metadata_all[is.na(bgnp_metadata_all$infant_relations) & bgnp_metadata_all$Type=='infant',"infant_relations"] <- "singleton"

# ADD Modified_NEXT_ID_without_preg_number
bgnp_metadata_all$Modified_NEXT_ID_without_preg_number <- sub('_.', '', bgnp_metadata_all$NEXT_ID)

# ADD days_from_first_collection
# STILL SOME NA, WAIT FOR THE FILE FROM TRISHLA
bgnp_metadata_all$days_from_first_collection <- bgnp_metadata$days_from_first_collection[match(bgnp_metadata_all$NG_ID, bgnp_metadata$NG_ID)]

# ADD Prokaryotic Shannon diversity index
bgnp_metadata_all$shannon <- bac_diversity_shannon$diversity[match(bgnp_metadata_all$NG_ID, row.names(bac_diversity_shannon))]

# harmonization:
bgnp_metadata_all[bgnp_metadata_all$Type=='mother',]$Type <- 'M'
bgnp_metadata_all[bgnp_metadata_all$Type=='infant',]$Type <- 'K'

# adding '_MGS' to some column names for convenience:
colnames(bgnp_metadata_all)[c(2,3,6:13, 18:19, 25, 28, 32)] <- paste0( colnames(bgnp_metadata_all)[c(2,3,6:13, 18:19, 25, 28, 32)], '_', 'MGS' )

row.names(bgnp_metadata_all) <- NULL #to harmonize row names, somewhere these were NG_IDs

##### MAIN OUTPUT: the full basic metadata for sequenced VLP samples: #####
metadata_chiliadal <- merge(Chiliadal_sequenced_samples, bgnp_metadata_all, by='Universal_ID', all.x=T)

# future steps: add aliquot id (if necessary? check if this could potentially affect the exact age)
# add the weight of the fecal aliquot for VLP
# add the aliquoting comments and isolation comments if any 
# add samples from pilot sequencing in baseclear


##### MAIN OUTPUT: the full basic metadata for all MGS samples that have at least 1 VLP sample from the family
metadata_chiliadal_w_imputation_all <- merge(Chiliadal_sequenced_samples, bgnp_metadata_all, by='Universal_ID', all = T)

# creating a metadata for families that have at least one VLP sample and MGS samples
families_to_keep <- list()
for (i in unique(metadata_chiliadal_w_imputation_all$FAMILY)) {
  if ( sum(!is.na(metadata_chiliadal_w_imputation_all[metadata_chiliadal_w_imputation_all$FAMILY==i,]$Sequencing_ID_VLP)) != 0) {
    families_to_keep[i] <- i
  }
}

metadata_chiliadal_w_imputation <- metadata_chiliadal_w_imputation_all[metadata_chiliadal_w_imputation_all$FAMILY %in% unname(unlist(families_to_keep)),]

# check how many orphan samples (only one member of the family is present among VLPs):
N_orphans <- as.data.frame(t(data.frame(families_to_keep)))
row.names(N_orphans) <- NULL
colnames(N_orphans) <- "FAMILY"
N_orphans$N_samples <- NA
N_orphans$N_fam_members <- NA

twin_families <- unique(metadata_chiliadal[metadata_chiliadal$Type_MGS=="K" & !is.na(metadata_chiliadal$Type_MGS) & metadata_chiliadal$infant_relations=='twins',"FAMILY"])
N_orphans$Twin_pregnancy <- ifelse(N_orphans$FAMILY %in% twin_families, "Yes", "No")

# second pregnancy families:
P2_families <- unique(bgnp_metadata_all[grep('_2', bgnp_metadata_all$NEXT_ID_MGS),"FAMILY"])
N_orphans$P2_present <- ifelse(N_orphans$FAMILY %in% P2_families, "Yes", "No")

# Preset mother and infant variables for all families
N_orphans$Mother_P1 <- "Yes" 
N_orphans$Mother_P2 <- "No" # the majority of the families do not have 2nd pregnancies
N_orphans$Infant_P1 <- "Yes"
N_orphans$Infant_P2 <- "No" # the majority of the families do not have 2nd pregnancies

# Loop through each family
for (i in N_orphans$FAMILY) {
  
  # NUMBER OF SAMPLES PER FAMILY
  N_orphans[N_orphans$FAMILY==i,"N_samples"] <- sum(metadata_chiliadal$FAMILY==i, na.rm = T)
  
  # NUMBER OF FAMILY MEMBERS
  N_orphans[N_orphans$FAMILY==i,"N_fam_members"] <- length(unique(metadata_chiliadal[metadata_chiliadal$FAMILY==i & !is.na(metadata_chiliadal$FAMILY),]$NEXT_ID_VLP))
  
  # NEXT IDs OF MATERNAL SAMPLES (for infering N pregnancies present in Chiliadal)
  maternal_samples <- unique(metadata_chiliadal[metadata_chiliadal$FAMILY==i & metadata_chiliadal$Type_VLP=="M","NEXT_ID_VLP"])
  
  # CHECK IF MATERNAL SAMPLES ARE PRESENT AT ALL
  if (  length(maternal_samples) < 1 ) {
    N_orphans[N_orphans$FAMILY==i, "Mother_P1"] <- "No"
  } else {
    # IF 2ND PREGNANCY IS PRESENT, CHECK IF MATERNAL SAMPLE FOR 1ST PREGNANCY IS PRESENT
    if ( (i %in% P2_families) & length(grep('_2', maternal_samples))!=0 & length(maternal_samples) < 2 ) {
      N_orphans[N_orphans$FAMILY==i, "Mother_P1"] <- "No"
      N_orphans[N_orphans$FAMILY==i, "Mother_P2"] <- "Yes"
    } else {
      # IF BOTH 1ST AND 2ND PREGNANCY SAMPLES ARE PRESENT
      if (length(maternal_samples) == 2) {
        N_orphans[N_orphans$FAMILY==i, "Mother_P2"] <- "Yes"
      }
    }
    
  }
  
  # NEXT IDs OF INFANT SAMPLES (for infering N pregnancies present in Chiliadal)
  infant_samples <- unique(metadata_chiliadal[metadata_chiliadal$FAMILY==i & metadata_chiliadal$Type_VLP=="K","NEXT_ID_VLP"])
  
  # DUE TO THE PRESENCE OF TWINS, WE ALSO NEED SAMPLE NAMES
  infant_samples_IDs <- metadata_chiliadal[metadata_chiliadal$FAMILY==i & metadata_chiliadal$Type_VLP=="K","NG_ID_short"]
  
  # CHECK IF INFANT SAMPLES ARE PRESENT AT ALL
  if (  length(infant_samples) < 1 ) {
    N_orphans[N_orphans$FAMILY==i, "Infant_P1"] <- "No"
  } else {
    # IF 2ND PREGNANCY IS PRESENT, CHECK IF INFANT SAMPLE FOR 1ST PREGNANCY IS PRESENT
    if ( (i %in% P2_families) & length(grep('C', infant_samples_IDs))<1 & length(infant_samples) < 2 ) {
      N_orphans[N_orphans$FAMILY==i, "Infant_P1"] <- "No"
      N_orphans[N_orphans$FAMILY==i, "Infant_P2"] <- "Yes"
    } else {
      # If there are two or more infant samples but the family is not marked as having twins,
      # assume a second pregnancy exists
      if (length(infant_samples) >=2 & !(i %in% twin_families)) {
        N_orphans[N_orphans$FAMILY==i, "Infant_P2"] <- "Yes"
      }
    }

  }
  
}

# NUMBER OF MOTHER-INFANT DYADS:
length(N_orphans[N_orphans$Twin_pregnancy=="Yes" & N_orphans$Mother_P1=="Yes", "FAMILY"]) #10 twin families
length(N_orphans[N_orphans$Twin_pregnancy=="No" & N_orphans$Mother_P1=="Yes" & N_orphans$Infant_P1=="Yes", "FAMILY"]) # 258
length(N_orphans[N_orphans$Twin_pregnancy=="No" & N_orphans$Mother_P2=="Yes" & N_orphans$Infant_P2=="Yes", "FAMILY"]) # 9
# Total number of mother-infant dyads with at least 1 maternal and 1 infant sample: 10*2 + 258 + 9 = 287

N_orphans$Mother_P1_MGS <- ifelse(N_orphans$FAMILY %in% unique(bgnp_metadata_all[bgnp_metadata_all$Type_MGS=="M" & bgnp_metadata_all$Modified_NEXT_ID_without_preg_number==bgnp_metadata_all$NEXT_ID_MGS,]$FAMILY), "Yes", "No")
N_orphans$Infant_P1_MGS <- ifelse(N_orphans$FAMILY %in% unique(bgnp_metadata_all[bgnp_metadata_all$Type_MGS=="K",]$FAMILY), "Yes", "No")
N_orphans$Mother_P2_MGS <- ifelse(N_orphans$FAMILY %in% unique(bgnp_metadata_all[bgnp_metadata_all$Type_MGS=="M" & bgnp_metadata_all$Modified_NEXT_ID_without_preg_number!=bgnp_metadata_all$NEXT_ID_MGS,]$FAMILY), "Yes", "No")
N_orphans$Infant_P2_MGS <- ifelse(N_orphans$FAMILY %in% unique(bgnp_metadata_all[bgnp_metadata_all$Type_MGS=="K" & substr(bgnp_metadata_all$NG_ID, 1, 1)=="Y",]$FAMILY), "Yes", "No")

# Total number of mother-infant dyads with at least 1 VLP sample per pair and MGS-imputable rest: 332 = 12*2 (twins) + 290 (no P2) + 8 (P1 at least mom or infant) + 10 (P2 at least mom or infant) 


# # All these calculations of the 2nd families just confuse me, I am thinking of artificially change their family IDs
# metadata_chiliadal$FAMILY_ed <- metadata_chiliadal$FAMILY
# metadata_chiliadal[grep('_2', metadata_chiliadal$NEXT_ID_MGS),]$FAMILY_ed <- sub('FAM0', 'FAM2', metadata_chiliadal[grep('_2', metadata_chiliadal$NEXT_ID_MGS),]$FAMILY_ed)
# metadata_chiliadal[grep('Y', metadata_chiliadal$NG_ID),]$FAMILY_ed <- sub('FAM0', 'FAM2', metadata_chiliadal[grep('Y', metadata_chiliadal$NG_ID),]$FAMILY_ed)
# 
# tmp <- data.frame(unique(metadata_chiliadal[!is.na(metadata_chiliadal$FAMILY_ed),]$FAMILY_ed))
# row.names(tmp) <- NULL
# colnames(tmp) <- "FAMILY"
# tmp$N_samples <- NA
# tmp$N_fam_members <- NA
# tmp$Twin_pregnancy <- ifelse(tmp$FAMILY %in% twin_families, "Yes", "No")
# tmp$Mother <- "Yes" 
# tmp$Infant <- "Yes" 
# 
# for (i in tmp$FAMILY) {
#   
#   # NUMBER OF SAMPLES PER FAMILY
#   tmp[tmp$FAMILY==i,"N_samples"] <- sum(metadata_chiliadal$FAMILY_ed==i, na.rm = T)
#   
#   # NUMBER OF FAMILY MEMBERS
#   tmp[tmp$FAMILY==i,"N_fam_members"] <- length(unique(metadata_chiliadal[metadata_chiliadal$FAMILY_ed==i & !is.na(metadata_chiliadal$FAMILY_ed),]$NEXT_ID_VLP))
#   
#   # NEXT IDs OF MATERNAL SAMPLES (for infering N pregnancies present in Chiliadal)
#   maternal_samples <- unique(metadata_chiliadal[metadata_chiliadal$FAMILY_ed==i & metadata_chiliadal$Type_VLP=="M","NEXT_ID_VLP"])
#   
#   # CHECK IF MATERNAL SAMPLES ARE PRESENT AT ALL
#   if (  length(maternal_samples) < 1 ) {
#     tmp[tmp$FAMILY==i, "Mother"] <- "No"
#   } 
#   
#   # NEXT IDs OF INFANT SAMPLES (for infering N pregnancies present in Chiliadal)
#   infant_samples <- unique(metadata_chiliadal[metadata_chiliadal$FAMILY_ed==i & metadata_chiliadal$Type_VLP=="K","NEXT_ID_VLP"])
#   
#    # CHECK IF INFANT SAMPLES ARE PRESENT AT ALL
#   if (  length(infant_samples) < 1 ) {
#     tmp[tmp$FAMILY==i, "Infant"] <- "No"
#   } 
# }

# # spliting the metadata containing VLP samples and MGS samples:
# 
# metadata_vertical_vlp <- metadata_chiliadal[,c("Universal_ID", "NEXT_ID_VLP","Type_VLP", "Timepoint_VLP", "FAMILY", "Sequencing_ID_VLP")]
# colnames(metadata_vertical_vlp)[c(2:6)] <- c('NEXT_ID', 'Type', 'Timepoint', 'FAMILY', 'Sequencing_ID')
# metadata_vertical_vlp$Timepoint <- paste0(metadata_vertical_vlp$Timepoint, '_VLP')
# 
# # also excluding failed samples:
# metadata_vertical_mgs <- bgnp_with_excluded[bgnp_with_excluded$STATUS_MGS!='FAILED' & bgnp_with_excluded$FAMILY %in% metadata_vertical_vlp$FAMILY,c("Universal_ID", "NEXT_ID_MGS", "Type_MGS", "Timepoint_regular", "FAMILY", "NG_ID")]
# colnames(metadata_vertical_mgs)[c(2:6)] <- c('NEXT_ID', 'Type', 'Timepoint', 'FAMILY', 'Sequencing_ID')
# metadata_vertical_mgs$Timepoint <- paste0(metadata_vertical_mgs$Timepoint, '_MGS')
# 
# metadata_vertical <- rbind(metadata_vertical_vlp, metadata_vertical_mgs)
# metadata_vertical$Timepoint <- factor(metadata_vertical$Timepoint, levels=c("P12_VLP", "P28_VLP", "B_VLP", "M1_VLP", "M3_VLP", "M6_VLP", "M12_VLP",
#                                                                             "P12_MGS", "P28_MGS", "B_MGS", "W2_MGS", "M1_MGS", "M2_MGS", "M3_MGS", "M6_MGS", "M9_MGS", "M12_MGS"),
#                                       ordered=T)


### CHECKING DISCERPANCIES IN BGNP METADATA ###
# species <- metaphlan[grep("s__",row.names(metaphlan)),]
# species <- row.names(species)[grep("t__", row.names(species), invert = T)]
# species_only_diversity <- data.frame(diversity(t(species), index = "shannon"))
# colnames(species_only_diversity) <- "Shannon_species"

# bgnp_metadata$shannon_species <- species_only_diversity$Shannon_species[match(bgnp_metadata$NG_ID, row.names(species_only_diversity))]
# 
# bgnp_metadata$delta_shannon <- bgnp_metadata$shannon - bgnp_metadata$shannon_species

# type_mod0 <- lm(delta_shannon ~ Type + dna_conc + clean_reads_FQ_1, data = bgnp_metadata)
# type_mod1  = lmer(delta_shannon ~ Type + dna_conc + clean_reads_FQ_1 + (1|NEXT_ID), REML = F, data = bgnp_metadata)
# BIC(type_mod0, type_mod1)
# exactLRT(type_mod1,type_mod0)
# summary(type_mod1)
# coef(summary(type_mod1))

# pdf('../../05.PLOTS/01.METADATA_HARMONIZATION/Shannon_diff_level_difference.pdf', width=6/2.54, height=10/2.54)
# ggplot(bgnp_metadata, aes(Type, delta_shannon, fill=Type)) + 
#   geom_sina(aes(color=Type)) +
#   geom_boxplot(alpha=0.5, width=0.5, outlier.shape = NA) +
#   labs(x="", y="Delta between shannon and shannon_species") +
#   theme_bw() +
#   theme(legend.position = "none", 
#         axis.title = element_text(face="bold"))
# dev.off()

# there are some samples that got double SAMPLE_IDs due to timepoint reassignment
View(bgnp_metadata[bgnp_metadata$SAMPLE_ID %in% bgnp_metadata[which(duplicated(bgnp_metadata$SAMPLE_ID)),"SAMPLE_ID"],])

##### OUTPUT #####
#### RESOLVED: due to the undocumented exclusion of some BGNP samples, had to track them down
#### write.table(metadata_chiliadal[is.na(metadata_chiliadal$sequence_control_MGS),], "../../03.OUTPUT/01.METADATA_HARMONIZATION/01.IN_PROGRESS/Missing_data_from_BGNP_03_01_2024.txt", sep='\t', row.names=F, quote=F)
#### RESOLVED: due to the undocumented exclusion of some BGNP samples, previous version of metadata was not complete
#### write.table(metadata_chiliadal, '../../01.METADATA/Chiliadal_metadata_ver_01_06042023.txt', sep='\t', quote=F, row.names=F)

write.table(bgnp_metadata_all, '../../03.OUTPUT/01.METADATA_HARMONIZATION/01.IN_PROGRESS/LLNEXT_metadata_220324_with_excluded.txt', sep='\t', quote=F, row.names = F)
write.table(metadata_chiliadal, '../../01.METADATA/Chiliadal_metadata_ver_03_29032024.txt', sep='\t', quote=F, row.names=F)


