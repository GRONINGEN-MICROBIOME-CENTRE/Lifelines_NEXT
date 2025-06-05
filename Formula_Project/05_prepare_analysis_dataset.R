# Load libraries
library(here)
library(tidyverse)
library(forcats)
library(dplyr)

# Load data
data <- readRDS(here("scripts/GITHUB/output/df_NG_ID_continuous_metadata.rds")) %>% as_tibble()

# Separate FOS and GOS
fos_gos_data <- data %>%
  filter(variable %in% c("FOS (gr/100mL)", "GOS (gr/100mL)"))

other_data <- data %>%
  filter(!variable %in% c("FOS (gr/100mL)", "GOS (gr/100mL)"))

# Pivot both datasets
pivot_vars <- c("NG_ID", "timepoint", "formula", "NEXT_ID", "Type", "FAMILY", "Sequenced",
                "raw_reads_FQ_1", "raw_reads_FQ_2", "human_reads_FQ_1", "human_reads_FQ_2",
                "clean_reads_FQ_1", "clean_reads_FQ_2", "reads_lost_QC", "sequence_control",
                "isolation_control", "dna_conc", "isolation_method", "NG_ID_short",
                "exact_age_days_at_collection", "exact_age_months_at_collection", "exact_age_years_at_collection",
                "Timepoint_categorical", "SAMPLE_ID", "metaphlan4_unclassified",
                "contaminant_1_Sphingomonas_sp_FARSPH", "contaminant_2_Phyllobacterium_myrsinacearum",
                "metaphlan4_unclassified_with_contaminants", "metaphlan4_unclassified_high_contaminants_factor_75",
                "BATCH_NUMBER", "infant_relations", "Modified_NEXT_ID_without_preg_number",
                "days_from_first_collection", "shannon")

fos_gos_data2 <- fos_gos_data %>%
  mutate(Timepoint_categorical = fct_relevel(Timepoint_categorical, "W2", "M1", "M2", "M3")) %>%
  pivot_wider(id_cols = all_of(pivot_vars), names_from = variable, values_from = c("Continuous", "categorical", "median group"))

other_data2 <- other_data %>%
  mutate(Timepoint_categorical = fct_relevel(Timepoint_categorical, "W2", "M1", "M2", "M3")) %>%
  pivot_wider(id_cols = all_of(pivot_vars), names_from = variable, values_from = c("Continuous", "median group"))

# Combine datasets
data2 <- full_join(fos_gos_data2, other_data2, by = pivot_vars)

# Check for NAs
cat("Total NAs in data2:", sum(is.na(data2)), "\n")
na_columns_data2 <- colSums(is.na(data2))
cat("Columns with NAs:\n")
print(na_columns_data2[na_columns_data2 > 0])

# Define original long names to rename (columns 35 to 180)
names_components <- names(data2)[35:180]

# Remove low-quality or irrelevant columns
exclude_columns <- c(
  "Continuous_Trans Fatty acids (gr/100mL)", "median group_Trans Fatty acids (gr/100mL)",
  "Continuous_Cholesterol (mg/100mL)", "median group_Cholesterol (mg/100mL)",
  "Continuous_Beta palmitate (mg/100mL)", "median group_Beta palmitate (mg/100mL)",
  "Continuous_Added sugars (gr/100mL)", "median group_Added sugars (gr/100mL)",
  "Continuous_Glucose (gr/100mL)", "median group_Glucose (gr/100mL)",
  "Continuous_Maltose (gr/100mL)", "median group_Maltose (gr/100mL)",
  "Continuous_Polysaccharides (gr/100mL)", "median group_Polysaccharides (gr/100mL)",
  "Continuous_Other (gr/mL)", "median group_Other (gr/mL)",
  "Continuous_Starch (gr/100mL)", "median group_Starch (gr/100mL)",
  "Continuous_X2-FL Oligosaccharide (mg/100mL)", "median group_X2-FL Oligosaccharide (mg/100mL)",
  "Continuous_Molybdenum (µg/100ml)", "median group_Molybdenum (µg/100ml)",
  "Continuous_Chromium (µg/100ml)", "median group_Chromium (µg/100ml)",
  "Continuous_Ash.content..gr.", "median group_Ash.content..gr.",
  "Continuous_B..infantis", "median group_B..infantis",
  "Continuous_B..longum", "median group_B..longum",
  "Continuous_B.lactis", "median group_B.lactis",
  "Continuous_Bifidus..BL..culture", "median group_Bifidus..BL..culture",
  "Continuous_Lactobacillus.rhamnosus.GG.culture", "median group_Lactobacillus.rhamnosus.GG.culture",
  "Continuous_Limosilactobacillus.fermentum.hereditum", "median group_Limosilactobacillus.fermentum.hereditum"
)

exclude_columns <- exclude_columns[exclude_columns %in% names(data2)]
data2 <- data2 %>% dplyr::select(-all_of(exclude_columns))
names_components <- setdiff(names_components, exclude_columns)

# Define the new names for the columns
names_components_abbr <- c(
  "continuous_gos_g", "continuous_fos_g", "categorical_gos_g", "categorical_fos_g", "median_group_gos_g", "median_group_fos_g",
  "continuous_energy_kcal", "continuous_fats_g", "continuous_sat_fats_g", "continuous_mono_unsat_g", "continuous_poly_unsat_g",
  "continuous_linoleic_omega6_g","continuous_linolenic_omega3_g", "continuous_AA_omega6_mg", "continuous_DHA_omega3_mg","continuous_carbs_g", "continuous_sugars_g",
  "continuous_lactose_g","continuous_inositol_mg", "continuous_fiber_g","continuous_X3_GL_galactosyllactose_mg", "continuous_proteins_g", "continuous_casein_g", "continuous_whey_g",
  "continuous_whey_casein_ratio","continuous_salt_g", "continuous_vitamin_a_ug", "continuous_vitamin_d_ug",
  "continuous_vitamin_e_mg", "continuous_vitamin_k_ug", "continuous_vitamin_c_mg", "continuous_vitamin_b1_mg",
  "continuous_vitamin_b2_mg", "continuous_niacin_mg", "continuous_vitamin_b6_mg", "continuous_folic_acid_ug",
  "continuous_vitamin_b12_ug", "continuous_biotin_ug", "continuous_pantothenic_acid_mg", "continuous_folate_ug",
  "continuous_sodium_mg", "continuous_potassium_mg", "continuous_chloride_mg", "continuous_calcium_mg",
  "continuous_phosphorus_mg", "continuous_magnesium_mg", "continuous_iron_mg", "continuous_zinc_mg",
  "continuous_copper_mg", "continuous_manganese_mg", "continuous_fluoride_mg", "continuous_selenium_ug",
  "continuous_iodine_ug", "continuous_choline_mg", "continuous_taurine_mg", "continuous_carnitine_mg",
  "continuous_nucleotides_mg", "median_group_energy_kcal", "median_group_fats_g", "median_group_sat_fats_g",
  "median_group_mono_unsat_g", "median_group_poly_unsat_g", "median_group_linoleic_omega6_g", "median_group_linolenic_omega3_g",
  "median_group_AA_omega6_mg", "median_group_DHA_omega3_mg","median_group_carbs_g", "median_group_sugars_g",
  "median_group_lactose_g","median_group_inositol_mg", "median_group_fiber_g", "median_group_X3_GL_galactosyllactose_mg",
  "median_group_proteins_g", "median_group_casein_g", "median_group_whey_g", "median_group_whey_casein_ratio",
  "median_group_salt_g","median_group_vitamin_a_ug", "median_group_vitamin_d_ug", "median_group_vitamin_e_mg",
  "median_group_vitamin_k_ug", "median_group_vitamin_c_mg", "median_group_vitamin_b1_mg", "median_group_vitamin_b2_mg",
  "median_group_niacin_mg", "median_group_vitamin_b6_mg", "median_group_folic_acid_ug", "median_group_vitamin_b12_ug",
  "median_group_biotin_ug", "median_group_pantothenic_acid_mg", "median_group_folate_ug", "median_group_sodium_mg",
  "median_group_potassium_mg", "median_group_chloride_mg", "median_group_calcium_mg", "median_group_phosphorus_mg",
  "median_group_magnesium_mg", "median_group_iron_mg", "median_group_zinc_mg", "median_group_copper_mg",
  "median_group_manganese_mg", "median_group_fluoride_mg", "median_group_selenium_ug", "median_group_iodine_ug",
 "median_group_choline_mg", "median_group_taurine_mg", "median_group_carnitine_mg", "median_group_nucleotides_mg"
)

# Rename columns
if (length(names_components_abbr) == length(names_components)) {
  name_mapping <- setNames(names_components_abbr, names_components)
  names(data2)[35:(34 + length(names_components_abbr))] <- name_mapping[names(data2)[35:(34 + length(names_components_abbr))]]
} else {
  stop("The lengths of names_components_abbr and names_components do not match.")
}

# Remove any leftover placeholder columns
cols_to_remove <- c("Continuous_X", "Continuous_X.1", "median group_X", "median group_X.1")
data2 <- data2 %>% dplyr::select(-any_of(cols_to_remove))

# Save cleaned dataset
saveRDS(data2, file = "data2.rds")
