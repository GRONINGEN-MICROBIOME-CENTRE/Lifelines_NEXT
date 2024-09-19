library(readr)
library(dplyr)
library(openxlsx)
library(stringr)
library(reshape2)
library(ggplot2)
library(tidyr)

setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\metadata")

qc_assembly_vd_ra_stats <- as.data.frame(read_tsv("qc_assembly_vd_ra_stats.tsv"))
names(qc_assembly_vd_ra_stats)[names(qc_assembly_vd_ra_stats) == 'sample_id'] <- 'Sample_name'


# Loading metadata from the papers and linking files (either from the paper or from ENA)
meta_garmaeva <- as.data.frame(read_tsv("VLP_metadata_final_10_05_2023.txt"))
meta_garmaeva1 <- as.data.frame(read_tsv("C:\\Users\\Natal\\Documents\\UMCG\\TCAM_BGNP\\masterfile_cross_sectional_2023_11_15.txt"))
meta_maqsood1 <- read.xlsx("40168_2019_766_MOESM7_ESM.xlsx")
meta_liang1 <- read.xlsx("41586_2020_2192_MOESM2_ESM.xlsx", sheet = 1, rows=c(3:352), colNames=T)
meta_liang2 <- read.xlsx("41586_2020_2192_MOESM2_ESM.xlsx", sheet = 8, rows=c(3:679), colNames=T)
meta_maqsood2 <- as.data.frame(read_tsv("linking_file_maqsood_edited.txt"))
meta_maqsood3 <- as.data.frame(read_tsv("filereport_read_run_PRJEB33578_tsv.txt"))
meta_walters <- as.data.frame(read_tsv("SRA_PRJNA916952_added_metadata.txt"))

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
  mutate(Subcohort = Cohort,
         Cohort = "liang",
         Sample_id = paste0("L", Sample_id, "_", substr(Note, 1, 1)),
         Subject_id = paste0("L", Subject_id),
         FAM_ID = Subject_id)

col_names_liang <- c("Sample_ID", "Subject_ID", "Type", "Timepoint", "infant_feeding_mode", "infant_formula_type", 
                     "infant_mode_delivery", "infant_sex", "Race", "mother_bmi_category", "timepoint_continuous_hours",
                     "mother_pregnancy_weight_gain", "mother_age_years", "infant_gestational_age", "infant_birthweight", 
                     "household_number", "household_underage", "mother_pregnancy_induce_HTN_diabetes", "mother_chorioamnionitis",
                     "Cohort", "Sample_name", "COMMENTS", "Subcohort", "FAM_ID")
colnames(meta_liang) <- col_names_liang

# Processing Walters metadata

meta_walters <- meta_walters %>%
  mutate(study_accession = NULL,
         secondary_study_accession = NULL,
         sample_accession = NULL,
         secondary_sample_accession = NULL,
         experiment_accession = NULL,
         submission_accession = NULL, 
         tax_id = NULL,
         scientific_name = NULL,
         instrument_platform = NULL,
         instrument_model = NULL,
         nominal_length = NULL,
         library_name = NULL,
         library_layout = NULL,
         library_strategy = NULL,
         library_source = NULL,
         library_selection = NULL,
         library_source = NULL,            
         library_selection = NULL,
         read_count = NULL,
         base_count = NULL,
         center_name = NULL,
         first_public = NULL,
         last_updated = NULL,
         experiment_title = NULL,
         study_title = NULL,
         study_alias = NULL,
         run_alias = NULL,
         fastq_bytes = NULL,
         fastq_md5 = NULL,
         fastq_ftp = NULL,
         fastq_aspera = NULL,
         fastq_galaxy = NULL,
         submitted_bytes = NULL,
         submitted_md5 = NULL,
         submitted_ftp = NULL,
         submitted_aspera = NULL,
         submitted_galaxy = NULL,
         submitted_format = NULL,
         sra_bytes = NULL,
         sra_md5 = NULL,
         sra_ftp = NULL,
         sra_aspera = NULL,
         sra_galaxy = NULL,
         broker_name = NULL,
         nominal_sdev = NULL,
         first_created = NULL,
         experiment_alias= NULL,
         sample_title = NULL,
         FathersWeightkg = NULL,
         FathersHeightcm = NULL,
         Language = NULL,
         Ethnicity = NULL,
         Date = NULL,
         BirthLengthcm = NULL,
         BabyWeightkg = NULL,
         Type = ifelse(grepl("M", ID), "Mother", "Infant"),
         FAM_ID = paste0("W", str_replace(ID, "M", "")),
         mother_age_years = floor((AgeMotherDays - AgeBabyDays)/ 365),
         ID = NULL,
         AgeMotherDays = NULL,
         Subject_ID = paste0(FAM_ID, substr(Type, 1, 1)),
         GestationWks = round(GestationWks, 0),
         infant_feeding_mode_complex = ifelse(Breast_milk == 1, "BF", NA),
         infant_feeding_mode_complex = ifelse(Formula == 1, "FF", infant_feeding_mode_complex),
         infant_feeding_mode_complex = ifelse(Breast_milk + Formula == 2, "MF", infant_feeding_mode_complex),
         infant_bithweight = BirthWeightkg * 1000,
         BirthWeightkg = NULL,
         Breast_milk = NULL,
         Formula = NULL,
         Timepoint = ifelse(VisitID == "M1", "Mtrim2", NA),
         Timepoint = ifelse(VisitID == "M2", "Mtrim3", Timepoint),
         Timepoint = ifelse(!is.na(AgeBabyDays), paste0("M", round(AgeBabyDays/30, 0)), Timepoint),
         Cohort = "walters"
  )
  
family_id_motherage <- unique(meta_walters[!is.na(meta_walters$mother_age_years), c("FAM_ID", "mother_age_years")])

family_id_visitid_timepoint <- unique(meta_walters[!is.na(meta_walters$Timepoint), c("FAM_ID", "VisitID", "Timepoint")])

meta_walters <- merge(meta_walters, family_id_motherage, by="FAM_ID", all.x = TRUE)
meta_walters <- merge(meta_walters, family_id_visitid_timepoint, by=c("FAM_ID", "VisitID"), all.x = TRUE)

meta_walters <- meta_walters %>%
  mutate(mother_age_years = ifelse(Type == "Mother", mother_age_years.y, NA),
         mother_age_years.x = NULL,
         mother_age_years.y = NULL,
         VisitID = NULL,
         Timepoint.x = NULL)


names(meta_walters) <- c("FAM_ID", "Sample_name", "Sample_ID", "Age_days", "infant_sex", "Race", "infant_place_delivery", 
                       "infant_gestational_age", "infant_mode_delivery", "Cows_milk", "Type", "Subject_ID", 
                       "infant_feeding_mode_complex", "infant_birthweight", "Cohort", "Timepoint", "mother_age_years")

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
                                           InfantMother == "mother" ~ NA),
         Subject_ID = paste0("M", Subject_ID),
         Family_Group = paste0("M", Family_Group),
         )
maqsood_to_drop <- c("Inclusion.status.for.bacterial.analysis", "Inclusion.status.for.viral.analysis", "submitted_ftp" )
meta_maqsood <- meta_maqsood[, !(names(meta_maqsood) %in% maqsood_to_drop)]
maqsood_lost <- sample_list_maqsood[!(sample_list_maqsood %in% meta_maqsood$new_name)]

# Editing the colnames for the further merging
col_names_maqsood <- c("Subject_ID", "FAM_ID", "Type", "infant_mode_delivery", 
                       "infant_feeding_mode", "infant_zygosity", "Age_days", 
                       "timepoint_continuous_hours", "infant_site_delivery", "mother_age_years",
                       "mother_bmi_category", "Race", "Sample_name", "Cohort","infant_type_pregnancy")
colnames(meta_maqsood) <- col_names_maqsood

control_samples <- data.frame(matrix(NA, ncol = ncol(meta_maqsood), nrow = 8))
control_samples$X13 <- maqsood_lost
control_samples$X3 <- "Neg_ctrl"
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
                      "phatyp_temperate_ricnhess_genomad", "phatyp_temperate_ricnhess_all", "infant_feeding_mode_imputed_W2", 
                      "infant_feeding_mode_imputed_M1", "infant_feeding_mode_imputed_M2", "infant_feeding_mode_imputed_M3", 
                      "infant_ffq_feeding_mode_derived_M6", "infant_ffq_feeding_mode_derived_M9", 
                      "infant_ffq_feeding_mode_derived_M12", "infant_ever_never_breastfed", "infant_feeding_mode_imputed_B_to_M3",
                      "viral_richness", "viral_alpha_diversity", "bacterial_alpha_diversity", "NUMBER", "Timepoint_continuous", "N_timepoints",
                      "DNA_CONC", "Isolation_batch")
meta_garmaeva <- meta_garmaeva[, !(names(meta_garmaeva) %in% garmaeva_to_drop)]

# Deriving maternal BMI from llnext data
llnext_pheno <- cross_phenotypes <- read.delim("C:\\Users\\Natal\\Documents\\UMCG\\TCAM_BGNP\\masterfile_cross_sectional_2023_11_15.txt")
llnext_pheno$NEXT_ID <- llnext_pheno$next_id_mother

meta_garmaeva <-merge(x=meta_garmaeva, y=unique(llnext_pheno[c("NEXT_ID", "mother_lifelinesbirthcards_prepreg_bmi_kg_m2")]), all.x = TRUE)
meta_garmaeva <- meta_garmaeva %>%
  mutate(mother_lifelinesbirthcards_prepreg_bmi_kg_m2 = ifelse((mother_lifelinesbirthcards_prepreg_bmi_kg_m2 >= 18 & mother_lifelinesbirthcards_prepreg_bmi_kg_m2 < 25), "normal", mother_lifelinesbirthcards_prepreg_bmi_kg_m2),
         mother_lifelinesbirthcards_prepreg_bmi_kg_m2 = ifelse((mother_lifelinesbirthcards_prepreg_bmi_kg_m2 >= 25 & mother_lifelinesbirthcards_prepreg_bmi_kg_m2 < 30), "overweight", mother_lifelinesbirthcards_prepreg_bmi_kg_m2),
         mother_lifelinesbirthcards_prepreg_bmi_kg_m2 = ifelse((mother_lifelinesbirthcards_prepreg_bmi_kg_m2 >= 30 & mother_lifelinesbirthcards_prepreg_bmi_kg_m2 < 35), "obese", mother_lifelinesbirthcards_prepreg_bmi_kg_m2))

col_names_garmaeva <- c("Subject_ID", "Sample_ID", "Sample_name", "Type", "Timepoint", "FAM_ID", "infant_type_pregnancy", 
                        "infant_sex", "infant_gestational_age", "infant_birthweight", "infant_place_delivery", 
                        "infant_mode_delivery", "mother_age_years", "infant_food_solid_intro_M6", "Age_days", "Age_months", 
                        "Age_years", "COMMENTS", "infant_feeding_mode_simple", "infant_feeding_mode_complex", "Cohort", "mother_bmi_category")
colnames(meta_garmaeva) <- col_names_garmaeva

garmaeva_lost <- sample_list_garmaeva[!(sample_list_garmaeva %in% meta_garmaeva$Sample_name)]
control_samples <- data.frame(matrix(NA, ncol = ncol(meta_garmaeva), nrow = 1))
control_samples$X3 <- garmaeva_lost
control_samples$X4 <- "Neg_ctrl"
control_samples$X21 <- "garmaeva"
colnames(control_samples) <- col_names_garmaeva
meta_garmaeva <- rbind(meta_garmaeva, control_samples)
meta_garmaeva$infant_type_pregnancy[meta_garmaeva$Sample_name %in% c("LN_7B03_VL_391", "LN_7B04_VL_392")] <- "single"

# Merging the tables
meta_garmaeva_liang <- merge(x=meta_garmaeva, y=meta_liang, all=TRUE)
meta_walters_maqsood <- merge(x=meta_walters, y=meta_maqsood, all=TRUE)
meta_all <- merge(x=meta_garmaeva_liang, y=meta_walters_maqsood, all=TRUE)

meta_all_with_qc <- merge(x=meta_all, y=qc_assembly_vd_ra_stats, by="Sample_name", all = T)
# samples_exclude <- c("SRR22944375", "SRR8653201", "SRR8800143", "SRR8800149")
# meta_all_with_qc <- meta_all_with_qc[!(meta_all_with_qc$Sample_name %in% samples_exclude), ]

# Fixing the columns one by one
meta_all_with_qc_curated <- meta_all_with_qc %>%
  mutate(Sample_ID = ifelse(grepl("Mitomycin", COMMENTS), paste0(Sample_ID, "IPM"), Sample_ID),
         Sample_ID = ifelse(grepl("no ind", COMMENTS), paste0(Sample_ID, "IPN"), Sample_ID),
         Sample_ID = ifelse(grepl("with Car", COMMENTS), paste0(Sample_ID, "IPCa"), Sample_ID),
         Sample_ID = ifelse(grepl("with Cip", COMMENTS), paste0(Sample_ID, "IPCi"), Sample_ID),
         Sample_ID = ifelse(grepl("in LB", COMMENTS), paste0(Sample_ID, "LB"), Sample_ID),
         Sample_ID = ifelse(grepl("in BSM", COMMENTS), paste0(Sample_ID, "BSM"), Sample_ID),
         Type = ifelse(grepl("virome of infant stool", COMMENTS), "Infant", Type),
         Type = ifelse(Timepoint == "Negative Control" & complete.cases(Timepoint), "Neg_ctrl", Type),
         Type = ifelse(Timepoint == "Positive Control" & complete.cases(Timepoint), "Pos_ctlr", Type),
         Type = ifelse(Timepoint == "Germ-free Mice" & complete.cases(Timepoint), "GFMice", Type),
         Type = ifelse(Type == "mother", "Mother", Type),
         Type = ifelse(Type == "infant", "Infant", Type),
         Type = ifelse(grepl("induced phages", COMMENTS), "Induced_phages", Type),
         Type = ifelse(grepl("kid", Sample_name), "Infant", Type),
         Type = ifelse(grepl("ctl", Sample_name), "Neg_ctrl", Type),
         Type = ifelse(grepl("formula", COMMENTS), "Formula_virome", Type),
         Cohort = ifelse(grepl("kid", Sample_name), "shah", Cohort),
         Cohort = ifelse(grepl("ctl", Sample_name), "shah", Cohort),
         infant_mode_delivery = ifelse(grepl("aginal", infant_mode_delivery), "VG", infant_mode_delivery),
         infant_mode_delivery = ifelse(grepl("C", infant_mode_delivery), "CS", infant_mode_delivery),
         mother_age_years = floor(mother_age_years),
         Subject_ID = ifelse(Type %in% c("Infant", "Mother"), Subject_ID, NA),
         FAM_ID = ifelse(Type %in% c("Infant", "Mother"), FAM_ID, NA),
         infant_type_pregnancy = ifelse(Cohort == "liang" & Type == "Infant", "single", infant_type_pregnancy),
         infant_type_pregnancy = ifelse(Cohort == "walters" & Type == "Infant", "single", infant_type_pregnancy),
         Race = ifelse(Race %in% c("African American", "AfricanAmerican"), "Black", Race),
         Race = ifelse(Cohort == "garmaeva" & Type != "Neg_ctrl", "White", Race),
         Race = ifelse(Race == "Native American", "Native_American", Race),
         Race = ifelse(Race == "Mix", "Other", Race),
         mother_bmi_category = ifelse(mother_bmi_category == "Lean", "normal", mother_bmi_category),
         mother_bmi_category = ifelse(mother_bmi_category == "Obese", "obese", mother_bmi_category),
         Timepoint = ifelse(Timepoint %in% c("Germ-free Mice", "Positive Control", "Negative Control"), NA, Timepoint),
         Timepoint = ifelse(Type == "Induced_phages", NA, Timepoint),
         Timepoint = str_replace(Timepoint, "Month ", "M"),
         Timepoint = str_replace(Timepoint, "Year 2~5", "Y2-5"),
         Timepoint = str_replace(Timepoint, "B", "M0"),
         Timepoint = ifelse(grepl("kid", Sample_name), "M12", Timepoint),
         Timepoint = ifelse(Timepoint == "P7", "Mtrim3", Timepoint),
         Age_months = ifelse(Cohort != "garmaeva" & !(is.na(Age_days)) & Type %in% c("Infant", "Mother"), Age_days/30, Age_months),
         Age_years = ifelse(Cohort != "garmaeva" & !(is.na(Age_months)) & Type %in% c("Infant", "Mother"), Age_months/12, Age_years),
         Timepoint = ifelse(is.na(Timepoint) & Type %in% c("Infant", "Mother") & !(is.na(Age_months)), paste0("M", round(Age_months, 0)), Timepoint),
         Timepoint = ifelse(is.na(Timepoint) & Type == "Mother" & !(is.na(Age_days)) & Age_days < 0, "M0", Timepoint),
         infant_place_delivery = ifelse(Cohort == "liang" & Type == "Infant", "hospital", infant_place_delivery),
         infant_place_delivery = ifelse(Cohort == "maqsood" & Type == "Infant", "hospital", infant_place_delivery),
         infant_sex = ifelse(infant_sex == "male" & (!is.na(infant_sex)), "Male", infant_sex),
         infant_sex = ifelse(infant_sex == "female" & (!is.na(infant_sex)), "Female", infant_sex),
         infant_gestational_age = round(infant_gestational_age, 0),
         infant_feeding_mode_complex = ifelse(infant_feeding_mode %in% c("Mixed", "Mix", "mix"), "MF", infant_feeding_mode_complex),
         infant_feeding_mode_complex = ifelse(infant_feeding_mode %in% c("Breastmilk", "breast_milk"), "BF", infant_feeding_mode_complex),
         infant_feeding_mode_complex = ifelse(infant_feeding_mode %in% c("Formula", "formula"), "FF", infant_feeding_mode_complex),
         infant_feeding_mode_simple = ifelse(infant_feeding_mode_complex %in% c("BF", "MF"), "ever_BF", infant_feeding_mode_simple),
         infant_feeding_mode_simple = ifelse(infant_feeding_mode_complex == "FF" & cohort == "maqsood", "never_BF", infant_feeding_mode_simple),
         infant_feeding_mode_simple = ifelse(infant_feeding_mode_simple == "excl_FF", "never_BF", infant_feeding_mode_simple),
         infant_zygosity = ifelse(Subject_ID %in% c("LLNEXT002053", "LLNEXT002062"), "Dizygotic", infant_zygosity),  # based on the cross_phenotypes column "mother_birthcard_type_twin"
         infant_zygosity = ifelse(Subject_ID %in% c("LLNEXT003908", "LLNEXT003917"), "Monozygotic", infant_zygosity),  # based on the cross_phenotypes column "mother_birthcard_type_twin"
         mother_chorioamnionitis = NULL,
         infant_feeding_mode = NULL,
         timepoint_continuous_hours = NULL,
         infant_food_solid_intro_M6 = NULL,
         Age_months = NULL,
         Age_years = NULL,
         Cohort = NULL,
         infant_site_delivery = NULL,
         mother_pregnancy_weight_gain = NULL,
         household_number = NULL,
         household_underage = NULL,
         mother_pregnancy_induce_HTN_diabetes = NULL,
         infant_formula_type = NULL
         )

breastfed <- unique(meta_all_with_qc_curated$Subject_ID[(meta_all_with_qc_curated$infant_feeding_mode_simple == "ever_BF") &
                                                          !is.na(meta_all_with_qc_curated$infant_feeding_mode_simple) &
                                                          !is.na(meta_all_with_qc_curated$Subject_ID)])

meta_all_with_qc_curated <- meta_all_with_qc_curated %>%
  mutate(infant_feeding_mode_simple = ifelse(Subject_ID %in% breastfed, "ever_BF", infant_feeding_mode_simple),
         infant_feeding_mode_simple = ifelse((is.na(infant_feeding_mode_simple) & infant_feeding_mode_complex == "FF"), "never_BF", infant_feeding_mode_simple))

for_infant_count <- unique(meta_all_with_qc_curated[meta_all_with_qc_curated$Type =="Infant" & !(is.na(meta_all_with_qc_curated$Type)), 
                                      c("FAM_ID", "Subject_ID")])
for_infant_count_duplicated <- for_infant_count[duplicated(for_infant_count$FAM_ID), ]
meta_all_with_qc_curated$infant_number_infamily <- NA
meta_all_with_qc_curated$infant_number_infamily[meta_all_with_qc_curated$Type =="Infant" & !(is.na(meta_all_with_qc_curated$Type))] <- "B1"
meta_all_with_qc_curated$infant_number_infamily[meta_all_with_qc_curated$Subject_ID %in% 
                                                  for_infant_count_duplicated$Subject_ID] <- "B2"


stats <- summary(meta_all_with_qc_curated)
na_count <- as.data.frame(colSums(data.frame(sapply(meta_all_with_qc_curated, is.na))))


setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\data")
nmds_coords <- read.delim("data.scores.txt")

Q_ndms1 <- quantile(nmds_coords$NMDS1, probs=c(.25, .75), na.rm = FALSE)
low_ndms1 <- Q_ndms1[1]-1.5*IQR(nmds_coords$NMDS1)
up_ndms1 <-  Q_ndms1[2]+1.5*IQR(nmds_coords$NMDS1)

Q_ndms2 <- quantile(nmds_coords$NMDS2, probs=c(.25, .75), na.rm = FALSE)
low_ndms2 <- Q_ndms2[1]-1.5*IQR(nmds_coords$NMDS2)
up_ndms2 <-  Q_ndms2[2]+1.5*IQR(nmds_coords$NMDS2)

clean_nmds <- subset(nmds_coords, nmds_coords$NMDS1 > low_ndms1 & nmds_coords$NMDS1 < up_ndms1 &
                       nmds_coords$NMDS2 > low_ndms2 & nmds_coords$NMDS2 < up_ndms2)

outliers <- row.names(nmds_coords)[!(row.names(nmds_coords) %in% row.names(clean_nmds))]

nmds_coords$Sample_name <- row.names(nmds_coords)
meta_all_with_qc_curated <- merge(meta_all_with_qc_curated, nmds_coords, by="Sample_name", all = TRUE)

meta_all_with_qc_curated$outliers <- 0
meta_all_with_qc_curated$outliers[meta_all_with_qc_curated$Sample_name %in% outliers] <- 1

meta_all_with_qc_curated$inclusion_status <- 1
meta_all_with_qc_curated$inclusion_status[meta_all_with_qc_curated$Type == "GFMice"] <- 0
meta_all_with_qc_curated$inclusion_status[meta_all_with_qc_curated$Type == "Pos_ctlr"] <- 0
meta_all_with_qc_curated$inclusion_status[meta_all_with_qc_curated$clean_reads_comb == 0] <- 0
meta_all_with_qc_curated$inclusion_status[is.na(meta_all_with_qc_curated$NMDS1)] <- 0

meta_all_with_qc_curated$inclusion_status_nmds <- 1
meta_all_with_qc_curated$inclusion_status_nmds[meta_all_with_qc_curated$Type == "GFMice"] <- 0
meta_all_with_qc_curated$inclusion_status_nmds[meta_all_with_qc_curated$Type == "Induced_phages"] <- 0
meta_all_with_qc_curated$inclusion_status_nmds[meta_all_with_qc_curated$Type == "Pos_ctlr"] <- 0
meta_all_with_qc_curated$inclusion_status_nmds[meta_all_with_qc_curated$clean_reads_comb == 0] <- 0
meta_all_with_qc_curated$inclusion_status_nmds[meta_all_with_qc_curated$Sample_name %in% outliers] <- 0
meta_all_with_qc_curated$inclusion_status_nmds[is.na(meta_all_with_qc_curated$NMDS1)] <- 0
meta_all_with_qc_curated$inclusion_status_nmds[meta_all_with_qc_curated$Type == "Formula_virome"] <- 0


summary_neg_ctls <- meta_all_with_qc_curated[meta_all_with_qc_curated$Type == "Neg_ctrl" & !(is.na(meta_all_with_qc_curated$Type)), 
                                            c("Sample_name", "COMMENTS", "Type", "cohort", "contigs_total", "total_viruses_discovered", 
                                              "outliers", "NMDS1", "NMDS2", "inclusion_status", "inclusion_status_nmds")]


write.table(meta_all_with_qc_curated, "C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\metadata\\metadata_with_qc_v3.tsv", sep='\t', row.names=F, col.names=T, quote=F)


meta_all_with_qc_curated_NCP <- meta_all_with_qc_curated[meta_all_with_qc_curated$Type %in% c("Neg_ctrl", "Infant", "Mother") & 
                                                           !(is.na(meta_all_with_qc_curated$Type)) & meta_all_with_qc_curated$cohort != "walters", ]

meta_all_with_qc_curated_NCP <- meta_all_with_qc_curated_NCP[, !colnames(meta_all_with_qc_curated_NCP) %in%
                                                               c("infant_sex", "infant_gestational_age", "infant_birthweight", "infant_mode_delivery", "mother_age_years",
                                                                 "Age_days", "mother_bmi_category", "infant_type_pregnancy", "infant_place_delivery",
                                                                 "infant_feeding_mode_complex", "Race", "infant_feeding_mode_simple", "Cows_milk", 
                                                                 "infant_zygosity", "NMDS1", "NMDS2", "outliers", "inclusion_status", "inclusion_status_nmds",
                                                                 "infant_number_infamily")]

meta_all_with_qc_curated_NCP$reads_mapped <- 1
meta_all_with_qc_curated_NCP$reads_mapped[meta_all_with_qc_curated_NCP$Sample_name %in% 
                                            c("SRR9660140", "SRR9660156", "SRR9660164", "SRR9660165",
                                              "SRR9660192", "SRR9660199", "SRR9660205", "SRR9660207")] <- 0

write.table(meta_all_with_qc_curated_NCP, "C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\metadata\\metadata_with_qc_NCPv2.tsv", sep='\t', row.names=F, col.names=T, quote=F)