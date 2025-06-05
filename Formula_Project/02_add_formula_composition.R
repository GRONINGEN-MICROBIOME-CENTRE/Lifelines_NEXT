# Load libraries
library(dplyr)
library(tidyverse)
library(stringr)
library(tidyr)
library(tibble)

# Load data
data_combined <- readRDS("output/df_data_combined.rds")

# Filter to only used formulas (value == 1 and not "other")

filtered_data_combined <- data_combined %>%
  filter(value == 1, formula != "other") %>%
  dplyr::select(-value)


# Load and clean formula composition data
formula_compositions <- read.csv2("data/Formula Ingredients .csv", header = TRUE, fileEncoding = "latin1") %>%
  rename(formula = ï..Formula)

# Standardize formula names
clean_formula_names <- function(df) {
  df$formula[df$formula == "Nutrilon_Omneo.comfort_1"] <- "Nutrilon_Omneo_comfort_1"
  df$formula[df$formula == "Nutrilon_Omneo.comfort_ 1"] <- "Nutrilon_Omneo_comfort_1"
  df$formula[df$formula == "Nutramigen_1_LG"] <- "Nutramigen_1_LGG"
  df
}

filtered_data_combined <- clean_formula_names(filtered_data_combined)
formula_compositions <- clean_formula_names(formula_compositions)

# Filter formulas used between W2–M3
formulas_W2_M3 <- c("Hero_standaard_1", "HiPP_1_Bio", "Nutrilon_Forte_1", "Nutrilon_Nenatal_1", "Nutrilon_Omneo_comfort_1", 
                    "Nutrilon_Pepti_1", "Nutrilon_Standaard_1", "Nutrilon_Hypo_Allergeen_1", "AH_standaard_1", 
                    "Jumbo_standaard_1", "Kruidvat_1", "HIPP_PRE_HA", "Nutramigen_1_LGG", "AH_bio_1", "Hero_HA", 
                    "Enfamil_AR_1", "Friso_standaard_1")

# Convert nutrient values to numeric
convert_nutrient_values <- function(df) {
  df %>%
    filter(formula %in% formulas_W2_M3) %>%
    mutate(across(-formula, ~ as.numeric(gsub(",", ".", 
                                              ifelse(grepl("^<", .), 
                                                     as.numeric(gsub(",", ".", sub("^<", "", .))) - 0.001, 
                                                     ifelse(. == "Yes", NA, .))))))
}

filtered_formula_compositions <- convert_nutrient_values(formula_compositions)

# Define numeric columns
numeric_columns <- c(
  "Energy.kcal.100ml", "Fats..gr.", "Saturated.fatty.acids..gr.", "Single.unsaturated..gr.",
  "Polyunsaturated..gr.", "Linoleic.acid.Omega6..gr.", "Linolenic.acid.Omega3.gr.", "Arachidonic.acid..AA..Omega6.mg.",
  "Docosahexaenoic.acid..DHA..Omega3.mg.", "Carbohydrates..gr.", "Sugars..gr.", "Lactose..gr.", "Glucose..gr.",
  "Maltose..gr.", "Inositol..mg.", "Polysaccharides..gr.", "Fiber..gr.", "GOS..gr.", "FOS..gr.",
  "X2..FL..Oligosaccharide.2..FL..mg.", "X3.GL..3..Galactosyllactose..mg.", "Proteins..gr.", "Casein..gr.",
  "Whey..gr.", "Whey.Casein", "Salt..gr.", "Vitamin.A..Î.g.", "Vitamin.D..Î.g.", "Vitamin.E..mg.",
  "Vitamin.K..Î.g.", "Vitamin.C..mg.", "Vitamin.B1..mg.", "Vitamin.B2..mg.", "Niacin..mg.", "Vitamin.B6..mg.",
  "Folic.acid..Î.g.", "Vitamin.B12..Î.g.", "Biotin..Î.g.", "Pantothenic.acid..mg.", "Folate..Î.g.", "Natrium..mg.",
  "Potassium..mg.", "Chloride..mg.", "Calcium..mg.", "Phosphorus..mg.", "Magnesium..mg.", "Iron..mg.", "Zinc..mg.",
  "Copper..mg.", "Manganese..mg.", "Fluoride..mg.", "Selenium..Î.g.", "Iodine..Î.g.", "Molybdenum..Î.g.",
  "Chromium..Î.g.", "Choline..mg.", "Taurine..mg.", "L.carnitine.mg.", "Nucleotides..mg."
)

# Create binary columns based on median
create_median_indicators <- function(df, columns) {
  for (col in columns) {
    median_val <- median(df[[col]], na.rm = TRUE)
    df[[paste0(col, "_median_group")]] <- ifelse(df[[col]] < median_val, 0, 1)
  }
  df
}

filtered_formula_compositions <- create_median_indicators(filtered_formula_compositions, numeric_columns)

# Select presence/absence columns
presence_cols <- c("GOS..gr.", "FOS..gr.", "X2..FL..Oligosaccharide.2..FL..mg.", "X3.GL..3..Galactosyllactose..mg.", "Whey.Casein")
selected_columns <- filtered_formula_compositions %>%
  dplyr::select(formula, dplyr::all_of(presence_cols))

# Merge with filtered data
merged_data_combined <- filtered_data_combined %>%
  left_join(selected_columns, by = "formula") %>%
  filter(formula != "Nutrilon_Nutriton")

# Convert presence columns to binary
merged_data_combined <- merged_data_combined %>%
  mutate(across(all_of(presence_cols), ~ ifelse(!is.na(.), 1, 0)))

# Add median group columns
median_columns <- filtered_formula_compositions %>%
  dplyr::select(formula, ends_with("_median_group"))

merged_data_combined2 <- merged_data_combined %>%
  left_join(median_columns, by = "formula")

# Save final dataframe
saveRDS(merged_data_combined2, file = "output/df_NG_ID_nutrient_categorical_median_wide.rds")
