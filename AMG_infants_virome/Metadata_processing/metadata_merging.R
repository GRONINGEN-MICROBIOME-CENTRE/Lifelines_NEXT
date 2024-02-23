library(readr)
library(dplyr)
library(openxlsx)

setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\metadata")

# Loading metadata from the papers and linking files (either from the paper or from ENA)
meta_garmaeva <- as.data.frame(read_tsv("VLP_metadata_final_10_05_2023.txt"))
meta_walters1 <- as.data.frame(read_tsv("metadata_mapping_stork_viromics.txt"))
meta_maqsood1 <- read.xlsx("40168_2019_766_MOESM7_ESM.xlsx")
meta_liang1 <- read.xlsx("41586_2020_2192_MOESM2_ESM.xlsx", sheet = 1, rows=c(3:352), colNames=T)
meta_liang2 <- read.xlsx("41586_2020_2192_MOESM2_ESM.xlsx", sheet = 8, rows=c(3:679), colNames=T)
meta_maqsood2 <- as.data.frame(read_tsv("linking_file_maqsood_edited.txt"))
meta_maqsood3 <- as.data.frame(read_tsv("filereport_read_run_PRJEB33578_tsv.txt"))
meta_walters2 <- as.data.frame(read_tsv("filereport_read_run_PRJNA916952_tsv.txt"))

# Loading sample lists that I was running
sample_list_garmaeva <- readLines("garmaeva_sample_list.txt")
sample_list_shah <- readLines("shah_sample_list.txt")
sample_list_liang <- readLines("liang_sample_list.txt")
sample_list_maqsood <- readLines("maqsood_sample_list.txt")
sample_list_walters <- readLines("walters_sample_list.txt")

#Processing Liang metadata
meta_liang2 <- subset(meta_liang2, Accession %in% sample_list_liang)
liang_lost <- sample_list_liang[!(sample_list_liang %in% meta_liang2$Accession)]

meta_liang <- merge(x=meta_liang1, y=meta_liang2, by="Sample_id")
liang_to_drop <- c("Bioproject_accession", "Biosample_accession", "Library_ID", "Library_type")
meta_liang <- meta_liang[, !(names(meta_liang) %in% liang_to_drop)]
meta_liang <- meta_liang %>%
  mutate(Cohort = "liang",
         Sample_id = paste0(Sample_id, "_", substr(Note, 1, 1)))

col_names_liang <- c("Sample_ID", "Subject_ID", "Type", "Timepoint", "infant_feeding_mode", "infant_formula_type", 
                     "infant_mode_delivery", "infant_sex", "Race", "mother_bmi_category", "timepoint_continuous_hours",
                     "mother_pregnancy_weight_gain", "mother_age_years", "infant_gestational_age", "infant_birthweight", 
                     "household_number", "household_underage", "mother_pregnancy_induce_HTN_diabetes", "mother_chorioamnionitis",
                     "Cohort", "Sample_name", "COMMENTS")
colnames(meta_liang) <- col_names_liang

# Processing Walters metadata
meta_walters2 <- meta_walters2[, names(meta_walters2) %in% c("run_accession", "library_name")]
meta_walters2$library_name <- sub("NA$", "", meta_walters2$library_name)
meta_walters2$library_name <- sub("_", "", meta_walters2$library_name, fixed = TRUE)
meta_walters2$library_name <- sub("-", ".", meta_walters2$library_name, fixed = TRUE)
names(meta_walters1)[names(meta_walters1) == "#SampleID"] <- "library_name"

walters_rename <- meta_walters2$library_name[!(meta_walters2$library_name %in% meta_walters1$library_name)]
walters_rename <- walters_rename[walters_rename != "M2201B3_D"]
meta_walters2$library_name[meta_walters2$library_name %in% walters_rename] <- sub("_D", "_R", walters_rename)
meta_walters <- merge(x=meta_walters1, y=meta_walters2, by="library_name")
walters_lost <- meta_walters2$run_accession[!(meta_walters2$library_name %in% meta_walters1$library_name)]
meta_walters <- meta_walters %>%
  mutate(Cohort = "walters")
walters_to_drop <- c("VisitGroup", "VisitNum", "CombinedExtractionTypeVisitNum", "VisitIDQIIME", "CombinedExtractionTypeVisitID")
meta_walters <- meta_walters[, !(names(meta_walters) %in% walters_to_drop)]
col_names_walters <- c("Sample_ID", "FAM_ID", "Date", "Visit_ID", "timepoint_continuous_days", "Type", "Extraction_Type", 
                       "Sample_name", "Cohort")
colnames(meta_walters) <- col_names_walters
meta_walters <- meta_walters %>%
  mutate(Subject_ID = paste0(FAM_ID, substr(Type, 1, 1)))

# Processing Maqsood metadata
meta_maqsood2$new_name <- sub("_R1.fastq.gz", "", meta_maqsood2$new_name, fixed = TRUE)
meta_maqsood2$new_name <- sub("_R2.fastq.gz", "", meta_maqsood2$new_name, fixed = TRUE)
meta_maqsood2$new_name <- sub("_", "", meta_maqsood2$new_name, fixed = TRUE)
meta_maqsood2$new_name <- sub("_", "", meta_maqsood2$new_name, fixed = TRUE)
meta_maqsood2 <- meta_maqsood2[meta_maqsood2$new_name %in% sample_list_maqsood, ]

meta_maqsood3 <- meta_maqsood3[, names(meta_maqsood3) %in% c("submitted_ftp", "sample_alias")]
meta_maqsood3$submitted_ftp <- sub(".*/", "", meta_maqsood3$submitted_ftp)
names(meta_maqsood2)[names(meta_maqsood2) == "old_name"] <- "submitted_ftp"

meta_maqsood22 <- merge(x=meta_maqsood2, y=meta_maqsood3, by="submitted_ftp")
meta_maqsood22$sample_alias <- sub("_viral", "", meta_maqsood22$sample_alias)
names(meta_maqsood22)[names(meta_maqsood22) == "sample_alias"] <- "Subject_ID"

meta_maqsood <- merge(x=meta_maqsood1, y=meta_maqsood22, by="Subject_ID")
meta_maqsood <- meta_maqsood %>%
  mutate(Cohort = "maqsood",
         infant_type_pregnancy = case_when(InfantMother == "infant" ~ "twin", 
                                           InfantMother == "mother" ~ NA))
maqsood_to_drop <- c("Inclusion.status.for.bacterial.analysis", "Inclusion.status.for.viral.analysis", "submitted_ftp" )
meta_maqsood <- meta_maqsood[, !(names(meta_maqsood) %in% maqsood_to_drop)]
maqsood_lost <- sample_list_maqsood[!(sample_list_maqsood %in% meta_maqsood$new_name)]

# Editing the colnames for the further merging
col_names_maqsood <- c("Subject_ID", "FAM_ID", "Type", "infant_mode_delivery", 
                       "infant_feeding_mode", "infant_zygosity", "timepoint_continuous_days", 
                       "timepoint_continuous_hours", "infant_site_delivery", "mother_age_years",
                       "mother_bmi_category", "Race", "Sample_name", "Cohort","infant_type_pregnancy")
colnames(meta_maqsood) <- col_names_maqsood

control_samples <- data.frame(matrix(NA, ncol = ncol(meta_maqsood), nrow = 8))
control_samples$X13 <- maqsood_lost
control_samples$X3 <- "Neg_ctlr"
control_samples$X14 <- "maqsood"
control_samples$X1 <- c("buffer_control_10", "buffer_control_11", "buffer_control_12", "buffer_control_13", 
                        "orsay_control_4633", "orsay_control_4654", "orsay_control_4676", "orsay_control_4699")
colnames(control_samples) <- col_names_maqsood
meta_maqsood <- rbind(meta_maqsood, control_samples)
meta_maqsood$Sample_ID <- meta_maqsood$Subject_ID

# Processing Garmaeva metadata
meta_garmaeva <- meta_garmaeva %>%
  mutate(Cohort = "garmaeva")

garmaeva_to_drop <- c("Universal_fecal_ID", "Old_ID", "Raw_reads", "Clean_reads", "Human_reads", "reads_lost_QC", 
                      "Short_sample_ID", "Short_sample_ID_bact", "Individual_ID", "N_viral_contigs", "Length_all_viral_contigs",
                      "contigs...0bp.", "contigs...1000bp.", "Total_length...1000bp.", "N50", "L50", 
                      "bacterial_contamination_perc_reads", "perc_reads_aligned", "perc_viral_contigs", "proportion_viral_length", 
                      "temperate_richness", "temperate_RA", "phatyp_temperate_RA", "phatyp_temperate_ricnhess", "phatyp_temperate_RA_all", 
                      "phatyp_temperate_RA_most", "phatyp_temperate_ricnhess_most", "phatyp_temperate_RA_genomad", 
                      "phatyp_temperate_ricnhess_genomad", "phatyp_temperate_ricnhess_all")

meta_garmaeva <- meta_garmaeva[, !(names(meta_garmaeva) %in% garmaeva_to_drop)]
col_names_garmaeva <- c("Sample_ID", "Sample_name", "Subject_ID", "Type", "Timepoint", "FAM_ID", "DNA_CONC", "infant_type_pregnancy", 
                        "infant_sex", "infant_gestational_age", "infant_birthweight", "infant_place_delivery", 
                        "infant_mode_delivery", "mother_age_years", "infant_feeding_mode_imputed_W2", 
                        "infant_feeding_mode_imputed_M1", "infant_feeding_mode_imputed_M2", "infant_feeding_mode_imputed_M3", 
                        "infant_ffq_feeding_mode_derived_M6", "infant_ffq_feeding_mode_derived_M9", 
                        "infant_ffq_feeding_mode_derived_M12", "infant_ever_never_breastfed", 
                        "infant_feeding_mode_imputed_B_to_M3", "infant_food_solid_intro_M6", "Age_days", "Age_months", 
                        "Age_years", "Timepoint_continuous", "N_timepoints", "viral_richness", "viral_alpha_diversity", 
                        "bacterial_alpha_diversity", "NUMBER", "COMMENTS", "Isolation_batch", 
                        "infant_ffq_feeding_mode_simple", "infant_ffq_feeding_mode_complex", "Cohort" )
colnames(meta_garmaeva) <- col_names_garmaeva

garmaeva_lost <- sample_list_garmaeva[!(sample_list_garmaeva %in% meta_garmaeva$Sample_name)]
control_samples <- data.frame(matrix(NA, ncol = ncol(meta_garmaeva), nrow = 1))
control_samples$X2 <- garmaeva_lost
control_samples$X4 <- "Neg_ctlr"
control_samples$X38 <- "garmaeva"
colnames(control_samples) <- col_names_garmaeva
meta_garmaeva <- rbind(meta_garmaeva, control_samples)

# Merging the tables
meta_garmaeva_liang <- merge(x=meta_garmaeva, y=meta_liang, all=TRUE)
meta_walters_maqsood <- merge(x=meta_walters, y=meta_maqsood, all=TRUE)
meta_all <- merge(x=meta_garmaeva_liang, y=meta_walters_maqsood, all=TRUE)

# Fixing the columns one by one
meta_all_test <- meta_all %>%
  mutate(Sample_ID = ifelse(grepl("Mitomycin", COMMENTS), paste0(Sample_ID, "IPM"), Sample_ID),
         Sample_ID = ifelse(grepl("no ind", COMMENTS), paste0(Sample_ID, "IPN"), Sample_ID),
         Sample_ID = ifelse(grepl("with Car", COMMENTS), paste0(Sample_ID, "IPCa"), Sample_ID),
         Sample_ID = ifelse(grepl("with Cip", COMMENTS), paste0(Sample_ID, "IPCi"), Sample_ID),
         Sample_ID = ifelse(grepl("in LB", COMMENTS), paste0(Sample_ID, "LB"), Sample_ID),
         Sample_ID = ifelse(grepl("in BSM", COMMENTS), paste0(Sample_ID, "BSM"), Sample_ID),
         Type = ifelse(grepl("virome of infant stool", COMMENTS), "Infant", Type),
         Type = ifelse(Timepoint == "Negative Control" & complete.cases(Timepoint), "Neg_ctlr", Type),
         Type = ifelse(Timepoint == "Positive Control" & complete.cases(Timepoint), "Pos_ctlr", Type),
         Type = ifelse(Timepoint == "Germ-free Mice" & complete.cases(Timepoint), "GFMice", Type),
         Type = ifelse(Type == "mother", "Mother", Type),
         Type = ifelse(Type == "infant", "Infant", Type),
         Type = ifelse(grepl("induced phages", COMMENTS), "Induced_phages", Type),
         infant_mode_delivery = ifelse(grepl("aginal", infant_mode_delivery), "VG", infant_mode_delivery),
         infant_mode_delivery = ifelse(grepl("C", infant_mode_delivery), "CS", infant_mode_delivery)
         )

## Columns done
# Sample_ID
# Sample_name
# Subject_ID
# Type
# infant_mode_delivery

## Questions to discuss:
# mother_age_years -> should I round it? Maqsood has it rounded, and us and liang don't
# 


