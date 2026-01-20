setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/02.SCRIPTS/01.METADATA_HARMONIZATION/')

#############################################################
# Here we prepare and organize LONG metadata for further 
# analyses
#############################################################

# There are many tables loaded and used only once in the script,
# delete them after use to declutter the global environment

#############################################################
# 0. Used file source
#############################################################

## VLP datasets' metadata:
## Chiliadal_virome_sequencing_batch_I.xlsx: list of Chiliadal VLP extractions sent for sequencing in August 2022
## Chiliadal_virome_sequencing_batch_II.xlsx: list of Chiliadal VLP extractions sent for sequencing in December 2022
## Samples_BaseClear_pilot_sequencing.xlsx: list of NEXT samples included in the Big Gut NEXT paper that were extracted and sequenced to test the new VLP protocol
## Controls_BaseClear_sequencing.xlsx: list of all sequencing control samples from BaseClear
## genomics:
## N_raw_reads_all: concatenated raw reads counts derived from per-sample bbduk.log
## N_human_reads_all_27052025: concatenated human reads counts derived from per-sample kneaddata.log
## viromeqc_all_CHILIADAL.txt: concatenated viromeQC output for all Chiliadal samples
## Contigs_stat_all: concatenated Quast output for all Chiliadal sample assemblies

## MGS datasets' metadata
## LLNEXT_metadata_03_01_2024.txt: metadata for NEXT samples included in the Big Gut NEXT paper
## LLNEXT_metadata_03_03_2023.txt: metadata for NEXT samples included in the Big Gut NEXT paper that includes samples that were excluded based on Sphingomonas contamination
## QC_MGS_TS.xlsx (sheet 3): all previously excluded from the Big Gut NEXT paper samples
## LLNEXT_metaphlan_4_complete_10_02_2023.txt: merged MetaPhlAn4 output file for NEXT samples included in the Big Gut NEXT paper
## DNA_CONC_ALL_WITH_DUPLICATES_MERGED_15_07_2022_TS.txt: DNA concentration measurements by Qubit received from Novogene sample QC reports
## genomics:
## Contigs_stat_1110_mgs: concatenated Quast output for all Chiliadal-connected MGS sample assemblies
## viromeqc_1110mgs: concatenated viromeQC output for all Chiliadal-connected MGS samples
#############################################################
# Functions
#############################################################

#############################################################
# 1. Loading libraries
#############################################################
library(readxl)
library(reshape2)
library(vegan)
library(dplyr)
library(tidyverse)

#############################################################
# 2. Load Input Data
#############################################################

# Chiliadal sequencing batches
Chiliadal_sequencing <- read_xlsx('~/Desktop/Projects_2021/NEXT_virome/06.Sequencing/Chiliadal_virome_sequencing_batch_I.xlsx')
Chiliadal_sequencing2 <- read_xlsx('~/Desktop/Projects_2021/NEXT_virome/06.Sequencing/Chiliadal_virome_sequencing_batch_II.xlsx') %>%
  filter(!is.na(NEXT_ID))

Prechiliadal_sequencing <- read_xlsx('../../01.METADATA/Samples_BaseClear_pilot_sequencing.xlsx', sheet = 2)

BaseClear_controls <- read_xlsx('../../01.METADATA/Controls_BaseClear_sequencing.xlsx')

# BGNP metadata
bgnp_metadata <- read.table('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/01.METADATA/LLNEXT_metadata_03_01_2024.txt', sep='\t', header=T)
early_bgnp <- read.table("../../01.METADATA/LLNEXT_metadata_03_03_2023.txt", sep='\t', header=T)

# Excluded BGNP samples
excluded_bgnp <- read_xlsx('~/Desktop/Projects_2021/NEXT_virome/08.Big_gut_NEXT_paper_data/Removed_MGS_samples/QC_MGS_TS.xlsx', sheet = 3)
excluded_bgnp[excluded_bgnp$NEXT_ID=="LLNEXT005780", "Type"] <- "Mother"
excluded_bgnp[excluded_bgnp$NEXT_ID=="LLNEXT005780", "FAMILY"] <- "FAM0236"

# Metaphlan BGNP data
metaphlan <- read.table("../../../08.Big_gut_NEXT_paper_data/LLNEXT_metaphlan_4_complete_10_02_2023.txt", sep='\t', header=T)
# Merge metaphlan data for contamination info
metaphlan <- metaphlan %>%
  select(-NCBI_tax_id) %>% 
  column_to_rownames(var = "clade_name") %>%
  set_names(substr(names(.), 1, 12))

# BGNP DNA concentration
dna_conc <- read.table("../../../08.Big_gut_NEXT_paper_data/DNA_CONC_ALL_WITH_DUPLICATES_MERGED_15_07_2022_TS.txt", sep='\t', header = T)
dna_conc <- dna_conc[dna_conc$ID!="Sample Name",]
dna_conc$DNA_CONC <- as.numeric(dna_conc$DNA_CONC)

# chiliadal sample genomic characteristics:
vlp_raw_reads <- read.table('../../04.RAW_DATA/02.Read_stats/N_raw_reads_all', sep='\t', header=F, col.names = c('Sequencing_ID', 'raw_reads'))
vlp_human_reads <- read.table('../../04.RAW_DATA/02.Read_stats/N_human_reads_all_27052025', sep='\t', header=F, col.names = c('Sequencing_ID', 'human_reads'))

vlp_clean_reads <- read.table('../../04.RAW_DATA/02.Read_stats/viromeqc_all_CHILIADAL.txt', sep = '\t', header=T)
colnames(vlp_clean_reads)[grep('Sample', colnames(vlp_clean_reads))] <- 'Sequencing_ID'
colnames(vlp_clean_reads)[grep('Reads', colnames(vlp_clean_reads))] <- 'clean_reads'

reads_list <- list(vlp_raw_reads, vlp_human_reads, 
                   vlp_clean_reads[,c("Sequencing_ID", "clean_reads")])

VLP_reads <- reads_list %>% reduce(full_join, by='Sequencing_ID')
VLP_reads$reads_lost_QC <- 1 - VLP_reads$clean_reads/VLP_reads$raw_reads

contig_stat_vlp <- read.delim('../../04.RAW_DATA/02.Read_stats/Contigs_stat_all', sep = '\t', header=T)
colnames(contig_stat_vlp) <- janitor::make_clean_names(colnames(contig_stat_vlp))
colnames(contig_stat_vlp)[grep('assembly', colnames(contig_stat_vlp))] <- 'Sequencing_ID'

vqc_vlp <- vlp_clean_reads
colnames(vqc_vlp)[-1] <- janitor::make_clean_names(colnames(vqc_vlp[,-1]))

genomics_list <- list(VLP_reads, 
                      contig_stat_vlp[,c("Sequencing_ID", "contigs_0_bp", "contigs_1000_bp")],
                      vqc_vlp[,c("Sequencing_ID", "bacterial_markers_alignment_rate")])

chili_genomics <- genomics_list %>% reduce(full_join, by = "Sequencing_ID")

BC_test_key <- read_xlsx('../../01.METADATA/Samples_BaseClear_pilot_sequencing.xlsx', sheet=1)

chili_genomics[grep('VNP', chili_genomics$Sequencing_ID),]$Sequencing_ID <- BC_test_key$Sequencing_ID[match(BC_test_key$SAMPLE, chili_genomics[grep('VNP', chili_genomics$Sequencing_ID),]$Sequencing_ID)]

rm(list = c('reads_list', 'vlp_clean_reads', 'vlp_human_reads', 
            'vlp_raw_reads', 'BC_test_key', 'contig_stat_vlp', 
            'genomics_list', 'VLP_reads', 'vqc_vlp'))

# pentachiliadal sample genomic characteristics:
contig_stat_mgs <- read.delim('../../04.RAW_DATA/02.Read_stats/Contigs_stat_1110_mgs', sep='\t', header=T)
colnames(contig_stat_mgs) <- janitor::make_clean_names(colnames(contig_stat_mgs))
colnames(contig_stat_mgs)[grep('assembly', colnames(contig_stat_mgs))] <- 'Sequencing_ID'

vqc_mgs <- read.table('../../04.RAW_DATA/02.Read_stats/viromeqc_1110mgs', sep='\t', header=T)
colnames(vqc_mgs) <- janitor::make_clean_names(colnames(vqc_mgs))
colnames(vqc_mgs)[grep('sample', colnames(vqc_mgs))] <- 'Sequencing_ID'

penta_genomics <- full_join(contig_stat_mgs[,c("Sequencing_ID", "contigs_0_bp", "contigs_1000_bp")], 
                            vqc_mgs[,c("Sequencing_ID", "bacterial_markers_alignment_rate")], by="Sequencing_ID")
rm(contig_stat_mgs, vqc_mgs)
#############################################################
# 3. Combine and clean Chiliadal metadata
#############################################################
Chiliadal_sequenced_samples <- rbind(Chiliadal_sequencing, Chiliadal_sequencing2, 
                              Prechiliadal_sequencing,  BaseClear_controls) %>%
  mutate(Sequencing_comments = ifelse(is.na(Sequencing_comments), 'Sequenced', Sequencing_comments),
         NEXT_ID = paste0('LLNEXT', NEXT_ID),
         NEXT_ID = ifelse(NEXT_ID == 'LLNEXT011062', 'LLNEXT011053', NEXT_ID), # changing the NEXT ID because it was incorrectly assigned during collection (father ID was used instead child ID); also in BGNP: QC_script_MGS_TS.R 
         DNA_CONC = ifelse(DNA_CONC == 'too low' & !is.na(DNA_CONC), '0.001', DNA_CONC), #assigning the lowest concentration possible to the samples where Qubit did not manage to measure the concentration
         DNA_CONC = as.numeric(DNA_CONC)
         ) %>%
  select(
    -ALIQUOTE_NR,
    -`...17`,
    -ADDITION,
    -`Sequencing_ID...21`
  ) %>%
  rename(Sequencing_ID = `Sequencing_ID...22`) %>%
  mutate(Universal_ID = paste0(NEXT_ID, '_', Timepoint))

rm(list = c('BaseClear_controls','Chiliadal_sequencing',
            'Chiliadal_sequencing2', 'Prechiliadal_sequencing'))
#############################################################
# 3.1 Fixing sample swaps
#############################################################
#known sample swap (confirmed by the downstream VLP analysis)
Chiliadal_sequenced_samples$NEXT_ID[Chiliadal_sequenced_samples$Sequencing_ID == "CHV024503E02"] <- "LLNEXT206301"
Chiliadal_sequenced_samples$Type[Chiliadal_sequenced_samples$Sequencing_ID == "CHV024503E02"] <- "K"

Chiliadal_sequenced_samples$NEXT_ID[Chiliadal_sequenced_samples$Sequencing_ID == "CHV020702E09"] <- "LLNEXT206295"
Chiliadal_sequenced_samples$Type[Chiliadal_sequenced_samples$Sequencing_ID == "CHV020702E09"] <- "M"
#############################################################
# 4. Process BGNP Metadata
#############################################################
# Add BGNP inclusion/exclusion status
early_bgnp <- early_bgnp %>%
  mutate(BGNP = ifelse(NG_ID %in% bgnp_metadata$NG_ID, "Included",
                       ifelse(metaphlan4_unclassified_with_contaminants < 75, "Excluded: Unknown reason", "Excluded: > 75% contamination"))) %>%
  select(-Sequenced)
  
# Some samples from Big Gut Next paper either failed to be sequenced, or were excluded based on the low
# read depth, or being an outlier at the PCA
excluded_bgnp <- excluded_bgnp %>%
  mutate(
         BGNP = paste0("Excluded: ", STATUS),
         reads_lost_QC = 1 - clean_reads_FQ_1/raw_reads_FQ_1,
         isolation_method = ifelse(grepl('X', NG_ID), 'fsk_vortex_gro', 'fsk_bb_kiel'),
         NG_ID_short = ifelse(is.na(NG_ID_short), substr(NG_ID, 1, 8), NG_ID_short),
         SAMPLE_ID = paste0(NEXT_ID, '_', Timepoint_categorical)
                                   ) %>%
  select(-STATUS, -Timepoint_regular, -Timepoint_categorical, -check_VLP_for_mixup) %>%
  rename('Timepoint_categorical' = 'Timepoint')

# Bacterial Shannon index (subspecies level)
strains <- metaphlan[grep("t__", rownames(metaphlan)),]
bac_diversity_shannon <- data.frame(diversity(t(strains), "shannon")) %>%
  set_names("diversity")

# Contamination estimation using metaphlan:
contaminants <- as.data.frame(t(metaphlan[  c( grep('UNCLASSIFIED', row.names(metaphlan)),
                                               grep('Sphingomonas_sp_FARSPH', row.names(metaphlan))[1],
                                               grep('Phyllobacterium_myrsinacearum', row.names(metaphlan))[1]), ] ))

contaminants$sum <- rowSums(contaminants)

# Match and integrate excluded BGNP with early BGNP
excluded_bgnp$metaphlan4_unclassified <- contaminants$UNCLASSIFIED[match(excluded_bgnp$NG_ID, rownames(contaminants))]
excluded_bgnp$contaminant_1_Sphingomonas_sp_FARSPH <- contaminants[[1]][match(excluded_bgnp$NG_ID, rownames(contaminants))]
excluded_bgnp$contaminant_2_Phyllobacterium_myrsinacearum <- contaminants[[2]][match(excluded_bgnp$NG_ID, rownames(contaminants))]
excluded_bgnp$metaphlan4_unclassified_with_contaminants <- contaminants$sum[match(excluded_bgnp$NG_ID, rownames(contaminants))]
excluded_bgnp$metaphlan4_unclassified_high_contaminants_factor_75 <- ifelse(excluded_bgnp$metaphlan4_unclassified_with_contaminants>=75, "high", "low")

# DNA concentration updates
excluded_bgnp$dna_conc <- dna_conc$DNA_CONC[match(excluded_bgnp$NG_ID, dna_conc$ID)]

# Sample ID and other harmonizations
excluded_bgnp$SAMPLE_ID <- paste0(excluded_bgnp$NEXT_ID, '_', excluded_bgnp$Timepoint_categorical)


# Combine early and excluded metadata
excluded_bgnp <- excluded_bgnp[,colnames(early_bgnp)]
bgnp_metadata_all <- bind_rows(early_bgnp, excluded_bgnp)

# Assign Universal_ID
bgnp_metadata_all <- bgnp_metadata_all %>%
  mutate(Timepoint_original = case_when(
    substr(NG_ID, 2, 3) == '01' ~ 'M1',
    substr(NG_ID, 2, 3) == '02' ~ 'M2',
    substr(NG_ID, 2, 3) == '03' ~ 'M3',
    substr(NG_ID, 2, 3) == '06' ~ 'M6',
    substr(NG_ID, 2, 3) == '09' ~ 'M9',
    substr(NG_ID, 2, 3) == '12' ~ 'M12',
    substr(NG_ID, 2, 3) == 'MB' ~ 'B',
    substr(NG_ID, 2, 3) == 'P3' ~ 'P12',
    substr(NG_ID, 2, 3) == 'P7' ~ 'P28',
    substr(NG_ID, 2, 3) == 'BB' ~ 'B',
    TRUE ~ substr(NG_ID, 2, 3)
  ),
  Universal_ID = paste0(NEXT_ID, '_', Timepoint_original))


# adding some missing info
bgnp_metadata_all[bgnp_metadata_all$NG_ID == "C01X01943A03", c("human_reads_FQ_1", "human_reads_FQ_2")] <- 30548 / 2
bgnp_metadata_all[bgnp_metadata_all$NG_ID == "CW2F002016H7", "BGNP"] <- "EXCLUDE OR RERUN"

# Remove known failed controls
bgnp_metadata_all <- bgnp_metadata_all[! bgnp_metadata_all$NG_ID %in% c("C06X009601G1", "C09F047214A9"),]

# Add sequencing controls info
bgnp_metadata_all <- bgnp_metadata_all %>%
  mutate(isolation_control = ifelse(Universal_ID == "LLNEXT005816_M6", "yes", "no"),
         sequence_control = ifelse(Universal_ID == "LLNEXT008093_M9", "yes", "no"),
         isolation_control = ifelse(is.na(isolation_control), "no", isolation_control),
         sequence_control = ifelse(is.na(sequence_control), "no", sequence_control),
         
         # Remove redundant read columns by aggregating then dropping raw/clean/human reads
         raw_reads = raw_reads_FQ_1 + raw_reads_FQ_2,
         human_reads = human_reads_FQ_1 + human_reads_FQ_2,
         clean_reads = clean_reads_FQ_1 + clean_reads_FQ_2,
         check_VLP_for_mixup = ifelse(BGNP == "Included", FALSE, TRUE),
         type_id_derived = case_when(
           substr(NG_ID, 1, 1) %in% c("A", "B") ~ "mother",
           substr(NG_ID, 1, 1) %in% c("C", "D", "Y") ~ "infant",
           TRUE ~ NA_character_
         )
  ) %>%
  select(-matches("FQ_"))

# Handle reassigned types or timepoints
reassigned_type_or_timepoint <- bgnp_metadata_all %>%
  filter(paste0(Type, '_', Timepoint_categorical) != paste0(type_id_derived, '_', Timepoint_original)) %>%
  pull(NG_ID)

for (i in reassigned_type_or_timepoint) {
  if (bgnp_metadata_all[bgnp_metadata_all$NG_ID == i, "Type"] != bgnp_metadata_all[bgnp_metadata_all$NG_ID == i, "type_id_derived"]) {
    bgnp_metadata_all[bgnp_metadata_all$NG_ID == i, "Timepoint_original"] <- bgnp_metadata_all[bgnp_metadata_all$NG_ID == i, "Timepoint_categorical"]
    bgnp_metadata_all[bgnp_metadata_all$NG_ID == i, "check_VLP_for_mixup"] <- TRUE
  }
}

# Drop helper column
bgnp_metadata_all <- bgnp_metadata_all %>% select(-type_id_derived)

# Add metadata fields from original BGNP
bgnp_metadata_all$BATCH_NUMBER <- bgnp_metadata$BATCH_NUMBER[match(bgnp_metadata_all$NG_ID, bgnp_metadata$NG_ID)]
bgnp_metadata_all$infant_relations <- bgnp_metadata$infant_relations[match(bgnp_metadata_all$NEXT_ID, bgnp_metadata$NEXT_ID)]
bgnp_metadata_all$infant_relations[is.na(bgnp_metadata_all$infant_relations) & bgnp_metadata_all$Type == 'infant'] <- 'singleton'
bgnp_metadata_all$Modified_NEXT_ID_without_preg_number <- sub('_.', '', bgnp_metadata_all$NEXT_ID)
bgnp_metadata_all$days_from_first_collection <- bgnp_metadata$days_from_first_collection[match(bgnp_metadata_all$NG_ID, bgnp_metadata$NG_ID)]


# Harmonize Type column
bgnp_metadata_all <- bgnp_metadata_all %>%
  mutate(Type = case_when(
    Type == 'mother' ~ 'M',
    Type == 'infant' ~ 'K',
    TRUE ~ Type
  ))

# Reorder columns
bgnp_metadata_all <- bgnp_metadata_all %>%
  select(NG_ID, NEXT_ID, Type, FAMILY, BGNP, raw_reads, human_reads, clean_reads,
         reads_lost_QC, sequence_control, isolation_control, dna_conc, isolation_method, 
         NG_ID_short, exact_age_days_at_collection, exact_age_months_at_collection, BATCH_NUMBER,
         exact_age_years_at_collection, Timepoint_categorical, SAMPLE_ID, metaphlan4_unclassified,
         contaminant_1_Sphingomonas_sp_FARSPH, contaminant_2_Phyllobacterium_myrsinacearum,
         metaphlan4_unclassified_with_contaminants, metaphlan4_unclassified_high_contaminants_factor_75,
         Timepoint_original, Universal_ID, check_VLP_for_mixup, Modified_NEXT_ID_without_preg_number, 
         infant_relations) %>%
  add_row(NG_ID = "FSKXNCTR20B3",
          NEXT_ID = "LLNEXTBLANK",
          Type = "BLANK", 
          FAMILY = "Controls",
          BGNP = "Excluded: LOW_RD",
          raw_reads = 4387732*2,
          clean_reads = 3746094*2,
          human_reads = 371328*2,
          reads_lost_QC = (4387732 - 3746094)/4387732,
          sequence_control = "no",
          isolation_control = "no", 
          dna_conc = 0.01, 
          isolation_method = "fsk_vortex_gro",
          NG_ID_short = "FSKXNCTR",
          Timepoint_categorical = "BLANK", 
          SAMPLE_ID ="FSKXNCTR_BLANK", 
          Timepoint_original = "BLANK",
          Universal_ID = "LLNEXTBLANK_BLANK",
          check_VLP_for_mixup = F,
          Modified_NEXT_ID_without_preg_number = "LLNEXTBLANK"
          )
  
# Clean up intermediate variables
rm(list = c("bgnp_metadata", "early_bgnp", "excluded_bgnp", 
            "contaminants", "strains",  
            "dna_conc"))
row.names(bgnp_metadata_all) <- NULL

bac_diversity_shannon$Universal_ID <- bgnp_metadata_all$Universal_ID[match(row.names(bac_diversity_shannon), bgnp_metadata_all$NG_ID)]


##############################
# 6. Finalize and Output Cleaned Metadata
##############################
# Long metadata creation:

Chiliadal_sequenced_samples$FAMILY <- bgnp_metadata_all$FAMILY[match(Chiliadal_sequenced_samples$Universal_ID,
                                                                     bgnp_metadata_all$Universal_ID)]

Chiliadal_sequenced_samples[is.na(Chiliadal_sequenced_samples$FAMILY), ]$FAMILY <- "Controls"

# preparing columns for long metadata merging:

# renaming the column to make it uniform
colnames(Chiliadal_sequenced_samples)[grep('DNA_CONC', colnames(Chiliadal_sequenced_samples))] <- 'dna_conc'

# adding the column with sequencing type
bgnp_metadata_all$seq_type <- 'MGS'
Chiliadal_sequenced_samples$seq_type <- 'VLP'

# splitting the bgnp column into status and reason
bgnp_metadata_all$bgnp_status <- ifelse(bgnp_metadata_all$BGNP=="Included", "Included", "Excluded")
bgnp_metadata_all$excl_bgnp_reason <- bgnp_metadata_all$BGNP
bgnp_metadata_all[bgnp_metadata_all$bgnp_status=="Included",]$excl_bgnp_reason <- NA
bgnp_metadata_all$BGNP <- NULL

Chiliadal_sequenced_samples$bgnp_status <- 'VLP'
Chiliadal_sequenced_samples$excl_bgnp_reason <- NA

# adding a provisional column about chiliadal samples to both:
bgnp_metadata_all$chiliadal_status <- 'MGS'
Chiliadal_sequenced_samples$chiliadal_status <- 'Included'
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$Sequencing_comments!="Sequenced",]$chiliadal_status <- 'Excluded'

# assigning "excluded" to the biological replicates:
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$Sequencing_ID=="CHV300415A05", ]$chiliadal_status <- 'Excluded'
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$Sequencing_ID=="CHV300715A08", ]$chiliadal_status <- 'Excluded'

# adding a column with the reasons for exclusions
Chiliadal_sequenced_samples$excl_chiliadal_reason <- NA
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$Sequencing_comments!="Sequenced",]$excl_chiliadal_reason <- Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$Sequencing_comments!="Sequenced",]$Sequencing_comments
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$Sequencing_ID=="CHV300415A05", ]$excl_chiliadal_reason <- 'Isolation replicate w lower read depth'
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$Sequencing_ID=="CHV300715A08", ]$excl_chiliadal_reason <- 'Isolation replicate w lower read depth'
bgnp_metadata_all$excl_chiliadal_reason <- NA

# adding a column on replicate_type
bgnp_metadata_all$replicate_type <- "none"
bgnp_metadata_all[bgnp_metadata_all$isolation_control=="yes",]$replicate_type <- "isolation_rep"
bgnp_metadata_all[bgnp_metadata_all$sequence_control=="yes",]$replicate_type <- "sequencing_rep"
bgnp_metadata_all$isolation_control <- NULL
bgnp_metadata_all$sequence_control <- NULL

Chiliadal_sequenced_samples$replicate_type <- "none"
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$Universal_ID=="LLNEXT123468_M12",]$replicate_type <- "isolation_rep"
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$Universal_ID=="LLNEXT003818_B",]$replicate_type <- "isolation_rep"
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$Universal_ID=="LLNEXT007878_M3",]$replicate_type <- "isolation_rep"

# renaming the column w isolation method:
colnames(Chiliadal_sequenced_samples)[grep('Tube', colnames(Chiliadal_sequenced_samples))] <- 'isolation_method'
Chiliadal_sequenced_samples[is.na(Chiliadal_sequenced_samples$isolation_method),"isolation_method"] <- 'BC_control'
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$isolation_method!='BC_control',]$isolation_method <- paste0('vlp_', Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$isolation_method!='BC_control',]$isolation_method)

# renaming the column w batch:
colnames(bgnp_metadata_all)[grep('BATCH_NUMBER', colnames(bgnp_metadata_all))] <- 'sequencing_batch'
bgnp_metadata_all[!is.na(bgnp_metadata_all$sequencing_batch),"sequencing_batch"] <- paste0('mgs_', stringr::str_pad(bgnp_metadata_all[!is.na(bgnp_metadata_all$sequencing_batch),"sequencing_batch"], 2, pad=0))
bgnp_metadata_all$sequencing_batch[bgnp_metadata_all$NEXT_ID == "LLNEXTBLANK"] <- "mgs_17"

colnames(Chiliadal_sequenced_samples)[grep('BATCH_SEQ', colnames(Chiliadal_sequenced_samples))] <- 'sequencing_batch'
Chiliadal_sequenced_samples[is.na(Chiliadal_sequenced_samples$sequencing_batch),"sequencing_batch"] <- 'BC_control'
Chiliadal_sequenced_samples$sequencing_batch <- paste0('vlp_', stringr::str_pad(gsub('SB', '', Chiliadal_sequenced_samples$sequencing_batch), 2, pad=0))

# renaming the column w position in batch & adding to bgnp:
bgnp_metadata_all$isolation_position <- NA
colnames(Chiliadal_sequenced_samples)[grep('ISOLATION', colnames(Chiliadal_sequenced_samples))] <- 'isolation_position'
Chiliadal_sequenced_samples[is.na(Chiliadal_sequenced_samples$isolation_position),"isolation_position"] <- 'BC_control'

# removing non-informative column BOX w physical aliquote location & aliquote NRs:
Chiliadal_sequenced_samples$BOX <- NULL
Chiliadal_sequenced_samples$Location <- NULL
Chiliadal_sequenced_samples$Aliquote_NR <- NULL
Chiliadal_sequenced_samples$DNA_NR <- NULL
bgnp_metadata_all$NG_ID_short <- NULL

# rename isolation batch:
bgnp_metadata_all$isolation_batch <- NA

colnames(Chiliadal_sequenced_samples)[grep('BATCH_ISO', colnames(Chiliadal_sequenced_samples))] <- 'isolation_batch'
Chiliadal_sequenced_samples[is.na(Chiliadal_sequenced_samples$isolation_batch),"isolation_batch"] <- 'BC_control'

# rename RT-step plate column:
bgnp_metadata_all$plate_RT <- NA
colnames(Chiliadal_sequenced_samples)[grep('Plate_N', colnames(Chiliadal_sequenced_samples))] <- 'plate_RT'
Chiliadal_sequenced_samples[is.na(Chiliadal_sequenced_samples$plate_RT),"plate_RT"] <- 'BC_control'

# rename RT-step plate position:
bgnp_metadata_all$plate_RT_pos <- NA
colnames(Chiliadal_sequenced_samples)[grep('Plate_position', colnames(Chiliadal_sequenced_samples))] <- 'plate_RT_pos'
Chiliadal_sequenced_samples[is.na(Chiliadal_sequenced_samples$plate_RT_pos),"plate_RT_pos"] <- 'BC_control'

# adding other missing columns to BGNP & renaming:
bgnp_metadata_all$COMMENTS <- NA
colnames(bgnp_metadata_all)[grep('NG_ID', colnames(bgnp_metadata_all))] <- 'Sequencing_ID'
Chiliadal_sequenced_samples$SAMPLE_ID <- NA

# adding other missing columns to BGNP & renaming:
bgnp_metadata_all$DATE_VLP <- NA
bgnp_metadata_all$DATE_DNA <- NA

# adding other missing columns to BGNP & renaming:
colnames(Chiliadal_sequenced_samples)[grep('Sequencing_comments', colnames(Chiliadal_sequenced_samples))] <- 'sequencing_status'
Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$sequencing_status!='Sequenced',]$sequencing_status <- 'Failed'

bgnp_metadata_all$sequencing_status <- 'Sequenced'
bgnp_metadata_all[grepl("Excluded: FAILED", bgnp_metadata_all$excl_bgnp_reason),]$sequencing_status <- 'Failed'

# taking into account 2nd pregnancies:
Chiliadal_sequenced_samples$Modified_NEXT_ID_without_preg_number <- gsub('_2', '', Chiliadal_sequenced_samples$NEXT_ID)

# adding columns lacking in Chiliadal:
Chiliadal_sequenced_samples$exact_age_days_at_collection <- bgnp_metadata_all$exact_age_days_at_collection[match(Chiliadal_sequenced_samples$Universal_ID, bgnp_metadata_all$Universal_ID)]
Chiliadal_sequenced_samples$exact_age_months_at_collection <- bgnp_metadata_all$exact_age_months_at_collection[match(Chiliadal_sequenced_samples$Universal_ID, bgnp_metadata_all$Universal_ID)]
Chiliadal_sequenced_samples$exact_age_years_at_collection <- bgnp_metadata_all$exact_age_years_at_collection[match(Chiliadal_sequenced_samples$Universal_ID, bgnp_metadata_all$Universal_ID)]

Chiliadal_sequenced_samples$Timepoint_categorical <- NA #will be edited upon reassigning

colnames(Chiliadal_sequenced_samples)[grep('Timepoint$', colnames(Chiliadal_sequenced_samples))] <- 'Timepoint_original'

Chiliadal_sequenced_samples$days_from_first_collection <- bgnp_metadata_all$days_from_first_collection[match(Chiliadal_sequenced_samples$Universal_ID, bgnp_metadata_all$Universal_ID)]

Chiliadal_sequenced_samples$metaphlan4_unclassified <- NA
Chiliadal_sequenced_samples$contaminant_1_Sphingomonas_sp_FARSPH <- NA
Chiliadal_sequenced_samples$contaminant_2_Phyllobacterium_myrsinacearum <- NA
Chiliadal_sequenced_samples$metaphlan4_unclassified_with_contaminants <- NA
Chiliadal_sequenced_samples$metaphlan4_unclassified_high_contaminants_factor_75 <- NA

Chiliadal_sequenced_samples$check_VLP_for_mixup <- bgnp_metadata_all$check_VLP_for_mixup[match(Chiliadal_sequenced_samples$Universal_ID, bgnp_metadata_all$Universal_ID)]


# changing the infant_relations column to Family structure to better reflect it:
bgnp_metadata_all$Family_structure <- NA
bgnp_metadata_all[bgnp_metadata_all$Type=='M',]$Family_structure <- 'MP1'
bgnp_metadata_all[bgnp_metadata_all$Type=='M' & grepl('_2', bgnp_metadata_all$NEXT_ID),]$Family_structure <- 'MP2'

bgnp_metadata_all[bgnp_metadata_all$Type=='K',]$Family_structure <- 'K1P1'

# twins
twin_fams <- unique(na.omit(bgnp_metadata_all[bgnp_metadata_all$infant_relations=="twins",]$FAMILY))
for (fam in twin_fams) {
  
  twins <- unique(bgnp_metadata_all[bgnp_metadata_all$FAMILY==fam & bgnp_metadata_all$Type=="K",]$NEXT_ID)
  twins <- as.numeric(gsub('LLNEXT', '', twins))
  
  T1 <- paste0('LLNEXT', stringr::str_pad(min(twins), 6, pad=0))
  T2 <- paste0('LLNEXT', stringr::str_pad(max(twins), 6, pad=0))
  
  bgnp_metadata_all[bgnp_metadata_all$NEXT_ID==T1,]$Family_structure <- 'K1P1'
  bgnp_metadata_all[bgnp_metadata_all$NEXT_ID==T2,]$Family_structure <- 'K2P1'
  
  
}


# multiple pregnancies

mult_preg_fams <- unique(na.omit(bgnp_metadata_all[bgnp_metadata_all$Family_structure=="MP2",]$FAMILY))
for (fam in  mult_preg_fams) {
  
  siblings <- unique(bgnp_metadata_all[bgnp_metadata_all$FAMILY==fam & bgnp_metadata_all$Type=="K",]$NEXT_ID)
  siblings <- as.numeric(gsub('LLNEXT', '', siblings))
  
  S1 <- paste0('LLNEXT', stringr::str_pad(min(siblings), 6, pad=0))
  S2 <- paste0('LLNEXT', stringr::str_pad(max(siblings), 6, pad=0))
  
  bgnp_metadata_all[bgnp_metadata_all$NEXT_ID==S1,]$Family_structure <- 'K1P1'
  bgnp_metadata_all[bgnp_metadata_all$NEXT_ID==S2,]$Family_structure <- 'K1P2'
  
}

bgnp_metadata_all$infant_relations <- NULL

Chiliadal_sequenced_samples$Family_structure <- bgnp_metadata_all$Family_structure[match(Chiliadal_sequenced_samples$Universal_ID, bgnp_metadata_all$Universal_ID)]

# adding reads info etc:
Chiliadal_sequenced_samples <- full_join(Chiliadal_sequenced_samples, chili_genomics, by='Sequencing_ID')

rm(chili_genomics)

# keeping only MGS samples that are connected to the VLP samples
mgs_smeta <- bgnp_metadata_all[bgnp_metadata_all$FAMILY %in% Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$sequencing_status=="Sequenced",]$FAMILY,]

# adding contig and vqc data to full overlap data:
mgs_smeta <- full_join(mgs_smeta, penta_genomics, by="Sequencing_ID")
rm(penta_genomics)

# cleaning up:
rm(bgnp_metadata_all)

# there are no VLP samples that are not present in BGNP except for NCs
Chiliadal_sequenced_samples[!Chiliadal_sequenced_samples$Universal_ID %in% mgs_smeta$Universal_ID, 'Universal_ID']

# full overlap MGS and VLP samples:
# sequenced & included in Chiliadal:


A <- Chiliadal_sequenced_samples[Chiliadal_sequenced_samples$sequencing_status=="Sequenced" & 
                                     Chiliadal_sequenced_samples$chiliadal_status=="Included" &
                                     Chiliadal_sequenced_samples$Type %in% c('K', 'M'),]$Universal_ID

# sequenced in BGNP:
B <- mgs_smeta[mgs_smeta$sequencing_status=="Sequenced" & 
                   #mgs_smeta$bgnp_status=="Included" & # some BGNP exclusion criteria are irrelevant for Chiliadal study
                   mgs_smeta$Type %in% c('K', 'M'),]$Universal_ID

UID_overlap <- intersect(A, B) # 1,110 samples; 1,097 otherwise (if only BGNP_included samples are used)

# merging into long metadata:
long_metadata <- rbind(Chiliadal_sequenced_samples, mgs_smeta)
long_metadata$full_overlap <- F

# adding a full_overlap column:
long_metadata[long_metadata$Universal_ID %in% UID_overlap,]$full_overlap <- T

# adding a non_dyad column: find orphan fams

nondyad_candidates <- c()
dummy <- c('K1P1', 'MP1')

for (fam in unique(na.omit(long_metadata$FAMILY))) {
  
 N_members <- length(unique(long_metadata[long_metadata$FAMILY==fam & long_metadata$sequencing_status=="Sequenced",]$NEXT_ID))
 
 members <- sort(unique(long_metadata[long_metadata$FAMILY==fam & long_metadata$sequencing_status=="Sequenced",]$Family_structure))

if (N_members!= 2 | all(members!=dummy)) {
  nondyad_candidates[fam] <- fam 
}
  
}

names(nondyad_candidates) <- NULL

nondyad <- setdiff(nondyad_candidates, c(twin_fams, mult_preg_fams, 'Controls'))

# adding a column reflecting the diff fam structures for an easy selection:
long_metadata$non_dyads <- "Dyad"
long_metadata[long_metadata$FAMILY %in% twin_fams,"non_dyads"] <- "Twin family"
long_metadata[long_metadata$FAMILY %in% mult_preg_fams,"non_dyads"] <- "Multiple pregnancy family"
long_metadata[long_metadata$FAMILY %in% nondyad,"non_dyads"] <- "Non-dyad"
long_metadata[long_metadata$FAMILY == 'Controls',"non_dyads"] <- "Controls"

# Sorting the columns:

column_order <- c('NEXT_ID', 'Type', 'Timepoint_original',
                  'FAMILY', 'Universal_ID', 'Sequencing_ID',
                  'seq_type', 'dna_conc', 'isolation_method',
                  'sequencing_status', 'replicate_type', 'raw_reads', 
                  'human_reads', 'clean_reads', 'reads_lost_QC',
                  'bacterial_markers_alignment_rate', 'contigs_0_bp', 'contigs_1000_bp',
                  'full_overlap', 'chiliadal_status',
                  'excl_chiliadal_reason', 'check_VLP_for_mixup',
                  'isolation_batch', 'sequencing_batch', 'DATE_VLP', 'DATE_DNA',
                  'isolation_position', 'COMMENTS', 'plate_RT', 'plate_RT_pos',
                  'bgnp_status', 'excl_bgnp_reason', 'SAMPLE_ID', 'metaphlan4_unclassified',
                  'contaminant_1_Sphingomonas_sp_FARSPH', 'contaminant_2_Phyllobacterium_myrsinacearum',
                  'metaphlan4_unclassified_with_contaminants', 'metaphlan4_unclassified_high_contaminants_factor_75',
                  'Modified_NEXT_ID_without_preg_number', 'Family_structure', 'non_dyads',
                  'exact_age_days_at_collection', 'exact_age_months_at_collection', 
                  'exact_age_years_at_collection','Timepoint_categorical')

long_metadata <- long_metadata[,column_order]
long_metadata$bacShannon <- bac_diversity_shannon$diversity[match(long_metadata$Universal_ID, bac_diversity_shannon$Universal_ID)]

rm(bac_diversity_shannon)
                  
# Files to be exported & shared:

# Long metadata for all fecal samples from Chiliadal-included individuals:
extended_long_metadata <- long_metadata

# Long metadata for sequenced-only and Chiliadal included samples:
use_long_metadata <- long_metadata[long_metadata$sequencing_status=="Sequenced" &
                                     long_metadata$chiliadal_status %in% c('Included', 'MGS'),]

use_long_metadata$sequencing_status <- NULL
use_long_metadata$chiliadal_status <- NULL
use_long_metadata$excl_chiliadal_reason <- NULL

# Basic metadata for sequenced-only Chiliadal samples:
use_VLP_metadata <- long_metadata[long_metadata$sequencing_status=="Sequenced" &
                                     long_metadata$chiliadal_status=='Included' &
                                    long_metadata$Type %in% c('M', 'K') &
                                    long_metadata$full_overlap==T,] # excludes 7 samples

remove_use_VLP <- c('seq_type', 'sequencing_status', 'replicate_type', 'full_overlap', 
                    'chiliadal_status', 'excl_chiliadal_reason', 'bgnp_status', 'excl_bgnp_reason',
                    'SAMPLE_ID', 'metaphlan4_unclassified', 'contaminant_1_Sphingomonas_sp_FARSPH',
                    'contaminant_2_Phyllobacterium_myrsinacearum', 'metaphlan4_unclassified_with_contaminants',
                    'metaphlan4_unclassified_high_contaminants_factor_75',
                    'Timepoint_categorical')

use_VLP_metadata <- use_VLP_metadata %>%
  select(-all_of(remove_use_VLP))

# Basic long metadata for VLP-MGS connection comparison:
long_VLP_MGS_overlap_metadata <- long_metadata[long_metadata$full_overlap == T & 
                                                 long_metadata$chiliadal_status %in% c('Included', 'MGS'),]

remove_use_VLP_MGS <- c('sequencing_status', 'replicate_type', 'chiliadal_status', 
                        'excl_chiliadal_reason')

long_VLP_MGS_overlap_metadata <- long_VLP_MGS_overlap_metadata %>%
  select(-all_of(remove_use_VLP_MGS))
  
# calculating VLP enrichment efficacy using own bacterial marker alignment:

# running calculation for VLPs and MGS sepratately (MGS is relative to own median, which is similar to the median in the original paper, 0.7)
wider_df <- long_VLP_MGS_overlap_metadata %>%
  pivot_wider(
    id_cols = Universal_ID,
    names_from = seq_type,
    values_from = bacterial_markers_alignment_rate
  ) %>%
  mutate(
    own_total_VLP_enrich = pmin(MGS / VLP, 100),
    own_total_MGS_enrich = pmin(median(MGS, na.rm = TRUE) / MGS, 100)
  )

# adding the new metrics to the table:
long_VLP_MGS_overlap_metadata <- long_VLP_MGS_overlap_metadata %>%
  left_join(
    wider_df %>%
      select(Universal_ID, own_total_VLP_enrich, own_total_MGS_enrich),
    by = "Universal_ID"
  ) %>%
  mutate(
    own_vir_enrich = case_when(
      seq_type == "VLP" ~ own_total_VLP_enrich,
      seq_type == "MGS" ~ own_total_MGS_enrich,
      TRUE ~ NA_real_
    )
  ) %>%
  select(-c(own_total_VLP_enrich, own_total_MGS_enrich))

rm(wider_df)
# metaphlan table for concurrent MGS:
metaphlan_overlap <- metaphlan[,colnames(metaphlan) %in% long_metadata[long_metadata$full_overlap==T & long_metadata$seq_type=="MGS",]$Sequencing_ID]
metaphlan_overlap <- metaphlan_overlap[rowSums(metaphlan_overlap) > 0,]
rm(metaphlan)

# for later: figure out where to reassign
reass <- long_metadata[long_metadata$seq_type=="MGS" & long_metadata$SAMPLE_ID!=long_metadata$Universal_ID, ]$Universal_ID
View(use_VLP_metadata[use_VLP_metadata$Universal_ID %in% reass, ])


# check missing:
skimr::skim(long_metadata)

#############################################################
# X. OUTPUT
#############################################################
write.table(mgs_smeta[mgs_smeta$sequencing_status=="Sequenced", "Sequencing_ID"], '../../06.CLEAN_DATA/Intermediate/MGS_NG_IDs_Chiliadal.txt', sep='\t', row.names = F, col.names = F, quote=F)

# Long metadata for all fecal samples from Chiliadal-included individuals:
write.table(extended_long_metadata, '../../06.CLEAN_DATA/02.FINAL/Chiliadal_meta_Ext_v05.txt', sep='\t', row.names = F, quote=F)

# Long metadata for all successfully sequenced unique fecal samples from Chiliadal-included individuals (both available VLP and MGS samples)
write.table(use_long_metadata, '../../06.CLEAN_DATA/02.FINAL/Chiliadal_meta_ExtFiltered_v05.txt', sep='\t', row.names = F, quote=F)

# Short VLP metadata for all unique VLP samples from Chiliadal that have a paired MGS sample:
write.table(use_VLP_metadata, '../../06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLPmatched_v04.txt', sep='\t', row.names = F, quote=F)

# Long metadata for all unique fecal samples from Chiliadal-included individuals (VLP-MGS overlap, 2220 samples)
write.table(long_VLP_MGS_overlap_metadata, '../../06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v04.txt', sep='\t', row.names = F, quote=F)

# metaphlan unfiltered full taxonomy for 1110 MGS samples with full overlap to Chiliadal:
write.table(metaphlan_overlap, '../../06.CLEAN_DATA/02.FINAL/MGS_Chiliadal_metaphlan_full_taxonomy_ver_01_07102025.txt', sep='\t', quote=F)



