setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Preparation of EGA metadata
#############################################################

##############################
# Loading libraries
##############################
library(tidyr)
library(dplyr)


# To harmonize MGS metadata: 
Total_MGS <- read.table('~/Resilio Sync/NEXT Pilot/rebuttal_new_analysis_microbiome/EGA_template/NEXT_PILOT_MGS_EGA.txt', sep='\t', header=T)

MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)

MGS_metadata$gender <- MGS_metadata$infant_sex
MGS_metadata[MGS_metadata$Type=='Mother',]$gender <- 'female'
MGS_metadata$Sequencing_type <- 'total_MGS'
MGS_fqs <- read.table('02.CLEAN_DATA/final_clean_total_MGS_fqs_to_EGA_07082023')
MGS_fqs$V1 <- gsub('.*/','', MGS_fqs$V1)
MGS_fqs$NG_ID <- substr(MGS_fqs$V1, 1, 12)

MGS_fqs_wide <- MGS_fqs %>%
  group_by(NG_ID) %>%
  mutate(Fastq_num = paste0("fastq_R", row_number())) %>%
  pivot_wider(names_from = Fastq_num, values_from = V1) %>%
  select(fastq_R1, fastq_R2, NG_ID)

MGS_metadata <- merge(MGS_metadata, MGS_fqs_wide, by='NG_ID')

# To harmonize VLP metadata:
VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)
VLP_metadata$Sequencing_type <- 'MVS'
VLP_metadata$gender <- VLP_metadata$infant_sex
VLP_metadata[VLP_metadata$Type=='Mother',]$gender <- 'female'
# adding fastqs names:
VLP_fqs <- read.table('02.CLEAN_DATA/final_clean_MVS_fqs_to_EGA_07082023')
VLP_fqs$NG_ID <- substr(VLP_fqs$V1, 1, 14)

VLP_fqs_wide <- VLP_fqs %>%
  group_by(NG_ID) %>%
  mutate(Fastq_num = paste0("fastq_R", row_number())) %>%
  pivot_wider(names_from = Fastq_num, values_from = V1) %>%
  select(fastq_R1, fastq_R2, NG_ID)

VLP_metadata <- merge(VLP_metadata, VLP_fqs_wide, by='NG_ID', all=T)

columns_to_EGA <- c('Old_ID', "Universal_fecal_ID", "NEXT_ID", "Individual_ID",
                    "NG_ID", "gender", "Sequencing_type",
                    "Type", "Timepoint", "FAM_ID",
                    "Raw_reads", "Clean_reads", "DNA_CONC", 
                    "infant_type_pregnancy", "infant_gestational_age", "infant_birthweight",
                    "infant_mode_delivery", "mother_age_years", "infant_ever_never_breastfed",
                    "fastq_R1", "fastq_R2")

EGA_metadata_combined <- rbind(VLP_metadata[, columns_to_EGA], MGS_metadata[, columns_to_EGA])

# check if there are inconsistencies between Universal_fecal_ID and NEXT_ID:

unique_combinations <- unique(EGA_metadata_combined[, c("Individual_ID", "NEXT_ID")])
duplicated_ids <- duplicated(unique_combinations$Individual_ID)
if (any(duplicated_ids)) {
  cat("Not all IDs from column1 have unique corresponding IDs in column2.\n")
} else {
  cat("Each ID from column1 has a unique corresponding ID in column2.\n")
  EGA_metadata_combined$NEXT_ID <- NULL # removing NEXT_IDs
}

# Further harmonizing (changing columnnames to self-explanatory names where needed):
colnames(EGA_metadata_combined)[grep('Old_ID', colnames(EGA_metadata_combined))] <- 'alias'
colnames(EGA_metadata_combined)[grep('Individual_ID', colnames(EGA_metadata_combined))] <- 'subjectId'
colnames(EGA_metadata_combined)[grep('NG_ID', colnames(EGA_metadata_combined))] <- 'bioSampleId'
colnames(EGA_metadata_combined)[grep('Type', colnames(EGA_metadata_combined))] <- 'Family_member'
colnames(EGA_metadata_combined)[grep('FAM_ID', colnames(EGA_metadata_combined))] <- 'Family_ID'

# adding missing columns:
EGA_metadata_combined$Title <- 'LLNEXT_pilot'
EGA_metadata_combined$description <- 'dna_stool_Lifelines_NEXT_pilot'
EGA_metadata_combined[EGA_metadata_combined$bioSampleId=='LN_7C08_VL_405',]$description <- 'dna_negative_control_Lifelines_NEXT_pilot'
EGA_metadata_combined$caseOrControl <- NA
EGA_metadata_combined$organismPart <- NA
EGA_metadata_combined$cellLine <- NA
EGA_metadata_combined$region <- NA

columns_to_EGA_order <- c("Title", "alias", "description",
                          "subjectId", "bioSampleId", "caseOrControl",
                          "organismPart", "cellLine",
                          "region", "Universal_fecal_ID","Sequencing_type", "gender", 
                          "Family_member", "Family_ID", "Timepoint", 
                          "DNA_CONC", "Raw_reads", "Clean_reads",
                          "infant_type_pregnancy", "infant_gestational_age", "infant_birthweight",
                          "infant_mode_delivery", "mother_age_years", "infant_ever_never_breastfed",
                          "fastq_R1", "fastq_R2")

EGA_metadata_combined <- EGA_metadata_combined[,columns_to_EGA_order]
# adding missing data for NGCTRL:
EGA_metadata_combined[EGA_metadata_combined$bioSampleId=='LN_7C08_VL_405',]$alias <- 'TVL405F1NOTR'
EGA_metadata_combined[EGA_metadata_combined$bioSampleId=='LN_7C08_VL_405',]$subjectId <- 'NGCTRL'
EGA_metadata_combined[EGA_metadata_combined$bioSampleId=='LN_7C08_VL_405',]$Sequencing_type <- 'MVS'
EGA_metadata_combined[EGA_metadata_combined$bioSampleId=='LN_7C08_VL_405',]$DNA_CONC <- 'too_low'
EGA_metadata_combined[EGA_metadata_combined$bioSampleId=='LN_7C08_VL_405',]$Raw_reads <- 24403158
EGA_metadata_combined[EGA_metadata_combined$bioSampleId=='LN_7C08_VL_405',]$Clean_reads <- 20305936

# changing DNA_concnetration back to Qubit output with too_low:
EGA_metadata_combined[EGA_metadata_combined$DNA_CONC==0.01 & EGA_metadata_combined$Sequencing_type=='MVS',]$DNA_CONC <- 'too_low'

write.table(EGA_metadata_combined,'05.MANUSCRIPT/NEXT_PILOT_EGA_combined.txt', sep='\t', quote=F, row.names = F)
write.table(EGA_metadata_combined[EGA_metadata_combined$Sequencing_type=="MVS",],'05.MANUSCRIPT/NEXT_PILOT_EGA_metaviromes.txt', sep='\t', quote=F, row.names = F)
write.table(EGA_metadata_combined[EGA_metadata_combined$Sequencing_type=="total_MGS",],'05.MANUSCRIPT/NEXT_PILOT_EGA_metagenomes.txt', sep='\t', quote=F, row.names = F)


