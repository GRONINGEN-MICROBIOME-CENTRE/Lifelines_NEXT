# Load libraries
library(dplyr)
library(tidyverse)
library(stringr)
library(tidyr)

# Load data
merged_data_combined2 <- readRDS("output/df_NG_ID_nutrient_categorical_median_wide.rds")
formula_compositions <- read.csv2("data/Formula Ingredients .csv", header = TRUE, fileEncoding = "latin1") %>%
  rename(formula = ï..Formula)

# Rename the columns
names(merged_data_combined2)[names(merged_data_combined2) == 'GOS..gr.'] <- 'GOS (gr/100mL)_categorical'
names(merged_data_combined2)[names(merged_data_combined2) == 'FOS..gr.'] <- 'FOS (gr/100mL)_categorical'
names(merged_data_combined2)[names(merged_data_combined2) == 'X2..FL..Oligosaccharide.2..FL..mg.'] <- 'X2-FL Oligosaccharide (mg/100mL)_categorical'
names(merged_data_combined2)[names(merged_data_combined2) == 'X3.GL..3..Galactosyllactose..mg.'] <- 'X3-GL Galactosyllactosyllactose (mg/100mL)_categorical'

# Delete the column Whey.Casein
merged_data_combined2 <- merged_data_combined2[, !names(merged_data_combined2) %in% 'Whey.Casein']


median_names <- c(
  "Energy.kcal.100ml_median_group" = "Energy (kcal/100ml)_median group", "Trans.fatty.acids..gr._median_group" = "Trans Fatty acids (gr/100mL)_median group","Cholesterol..mg._median_group" = "Cholesterol (mg/100mL)_median group",
  "Beta.palmitate..mg._median_group" = "Beta palmitate (mg/100mL)_median group","Added.sugars..gr._median_group" = "Added sugars (gr/100mL)_median group","Other..gr._median_group" = "Other (gr/mL)_median group",
  "Starch..gr._median_group" = "Starch (gr/100mL)_median group", "Fats..gr._median_group" = "Fats (gr/100ml)_median group","Saturated.fatty.acids..gr._median_group" = "Saturated Fatty Acids (gr/100mL)_median group", "Single.unsaturated..gr._median_group" = "Single Unsaturated (gr/100mL)_median group", "Polyunsaturated..gr._median_group" = "Polyunsaturated (gr/100mL)_median group", "Linoleic.acid.Omega6..gr._median_group" = "Linoleic Acid (Omega-6) (gr/100mL)_median group", "Linolenic.acid.Omega3.gr._median_group" = "Linolenic Acid (Omega-3) (gr/100mL)_median group",
  "Arachidonic.acid..AA..Omega6.mg._median_group" = "Arachidonic Acid (AA) (Omega-6) (mg/100mL)_median group", "Docosahexaenoic.acid..DHA..Omega3.mg._median_group" = "Docosahexaenoic Acid (DHA) (Omega-3) (mg/100mL)_median group",
  "Carbohydrates..gr._median_group" = "Carbohydrates (gr/100mL)_median group", "Sugars..gr._median_group" = "Sugars (gr/100mL)_median group", "Lactose..gr._median_group" = "Lactose (gr/100mL)_median group",
  "Glucose..gr._median_group" = "Glucose (gr/100mL)_median group", "Maltose..gr._median_group" = "Maltose (gr/100mL)_median group", "Inositol..mg._median_group" = "Inositol (mg/100mL)_median group",
  "Polysaccharides..gr._median_group" = "Polysaccharides (gr/100mL)_median group", "Fiber..gr._median_group" = "Fiber (gr/100mL)_median group","GOS..gr._median_group" = "GOS (gr/100mL)_median group", "FOS..gr._median_group" = "FOS (gr/100mL)_median group",
  "X2..FL..Oligosaccharide.2..FL..mg._median_group" = "X2-FL Oligosaccharide (mg/100mL)_median group","X3.GL..3..Galactosyllactose..mg._median_group" = "X3-GL Galactosyllactose (mg/100mL)_median group",
  "Proteins..gr._median_group" = "Proteins (gr/100mL)_median group","Casein..gr._median_group" = "Casein (gr/100mL)_median group","Whey..gr._median_group" = "Whey (gr/100mL)_median group",
  "Whey.Casein_median_group" = "Whey/Casein Ratio_median group","Salt..gr._median_group" = "Salt (gr/100mL)_median group", "Vitamin.A..Î.g._median_group" = "Vitamin A (µg/100mL)_median group", "Vitamin.D..Î.g._median_group" = "Vitamin D (µg/100mL)_median group","Vitamin.E..mg._median_group" = "Vitamin E (mg/100mL)_median group","Vitamin.K..Î.g._median_group" = "Vitamin K (µg/100mL)_median group", "Vitamin.C..mg._median_group" = "Vitamin C (mg/100mL)_median group","Vitamin.B1..mg._median_group" = "Vitamin B1 (mg/100mL)_median group","Vitamin.B2..mg._median_group" = "Vitamin B2 (mg/100mL)_median group","Niacin..mg._median_group" = "Niacin (mg/100mL)_median group","Vitamin.B6..mg._median_group" = "Vitamin B6 (mg/100mL)_median group","Folic.acid..Î.g._median_group" = "Folic Acid (µg/100mL)_median group", "Vitamin.B12..Î.g._median_group" = "Vitamin B12 (µg/100mL)_median group", "Biotin..Î.g._median_group" = "Biotin (µg/100mL)_median group","Pantothenic.acid..mg._median_group" = "Pantothenic Acid (mg/100mL)_median group","Folate..Î.g._median_group" = "Folate (µg/100mL)_median group","Natrium..mg._median_group" = "Natrium (mg/100mL)_median group", "Potassium..mg._median_group" = "Potassium (mg/100mL)_median group", "Chloride..mg._median_group" = "Chloride (mg/100ml)_median group", "Calcium..mg._median_group" = "Calcium (mg/100ml)_median group", "Phosphorus..mg._median_group" = "Phosphorus (mg/100ml)_median group","Magnesium..mg._median_group" = "Magnesium (mg/100ml)_median group","Iron..mg._median_group" = "Iron (mg/100ml)_median group","Zinc..mg._median_group" = "Zinc (mg/100ml)_median group", "Copper..mg._median_group" = "Copper (mg/100ml)_median group","Manganese..mg._median_group" = "Manganese (mg/100ml)_median group", "Fluoride..mg._median_group" = "Fluoride (mg/100ml)_median group","Selenium..Î.g._median_group" = "Selenium (µg/100ml)_median group","Iodine..Î.g._median_group" = "Iodine (µg/100ml)_median group", "Molybdenum..Î.g._median_group" = "Molybdenum (µg/100ml)_median group", "Chromium..Î.g._median_group" = "Chromium (µg/100ml)_median group","Choline..mg._median_group" = "Choline (mg/100ml)_median group", "Taurine..mg._median_group" = "Taurine (mg/100ml)_median group","L.carnitine.mg._median_group" = "L-carnitine (mg/100ml)_median group","Nucleotides..mg._median_group" = "Nucleotides (mg/100ml)_median group")


#Changing the column names for the categorical median variables 

colnames(merged_data_combined2) <- ifelse(colnames(merged_data_combined2) %in% names(median_names), median_names[colnames(merged_data_combined2)], colnames(merged_data_combined2))


#Load formula composition data continuous 

#Load the formula composition data 
formula_compositions <- read.csv2("C:/Users/User/Desktop/Formula_project/data/Formula Ingredients .csv",header = TRUE, fileEncoding = "latin1") 
formula_compositions <-
  formula_compositions %>%
  rename(formula = ï..Formula )

# Changing the name of the formula "Nutrilon_Omneo.comfort_1" to "Nutrilon_Omneo_comfort_1"

formula_compositions$formula[formula_compositions$formula == "Nutrilon_Omneo.comfort_ 1"] <- "Nutrilon_Omneo_comfort_1"


#selecting only the formulas given from W2 to M3 (no Nutrilon_Nutriton)

formulas_W2_M3 <- c("Hero_standaard_1","HiPP_1_Bio","Nutrilon_Forte_1","Nutrilon_Nenatal_1","Nutrilon_Omneo_comfort_1","Nutrilon_Pepti_1","Nutrilon_Standaard_1","Nutrilon_Hypo_Allergeen_1", "AH_standaard_1", "Jumbo_standaard_1", "Kruidvat_1", "HIPP_PRE_HA", "Nutramigen_1_LGG", "AH_bio_1","Hero_HA", "Enfamil_AR_1", "Friso_standaard_1")

filtered_formula_compositions <- formula_compositions[formula_compositions$formula %in% formulas_W2_M3, ] %>%
  mutate(across(-formula, ~ as.numeric(gsub(",", ".", ifelse(grepl("^<", .), as.numeric(gsub(",", ".", sub("^<", "", .))) - 0.001, ifelse(. == "Yes", NA, .))))))


# Define the display names for the columns
display_names <- c(
  "Energy.kcal.100ml" = "Energy (kcal/100ml)", "Trans.fatty.acids..gr." = "Trans Fatty acids (gr/100mL)", "Cholesterol..mg." = "Cholesterol (mg/100mL)", "Beta.palmitate..mg." = "Beta palmitate (mg/100mL)", "Added.sugars..gr." = "Added sugars (gr/100mL)", "Other..gr." = "Other (gr/mL)", "Starch..gr." = "Starch (gr/100mL)", "GOS..gr..x" = "GOS (gr/100mL)", "FOS..gr..x (gr/100mL)", "X2..FL..Oligosaccharide.2..FL..mg..x" = "X2'FL Oligosaccharide (mg/100ml)", "X3.GL..3..Galactosyllactose..mg..x" = "X3' GL Galactosyllactose (mg/100mL) ", "GOS..gr..y" = "GOS_categorical", "FOS..gr..y" = "FOS_categorical", "X2..FL..Oligosaccharide.2..FL..mg..y" = "X2_FL_Oligosaccharide_categorical", "X3.GL..3..Galactosyllactose..mg..y" = "X3_GL_Galactosyllactose_categorical", "Whey.Casein.y" = "Whey/Cazein ratio", "Fats..gr." = "Fats (gr/100ml)", "Saturated.fatty.acids..gr." = "Saturated Fatty Acids (gr/100mL)",
  "Single.unsaturated..gr." = "Single Unsaturated (gr/100mL)", "Polyunsaturated..gr." = "Polyunsaturated (gr/100mL)", "Linoleic.acid.Omega6..gr." = "Linoleic Acid (Omega-6) (gr/100mL)",
  "Linolenic.acid.Omega3.gr." = "Linolenic Acid (Omega-3) (gr/100mL)", "Arachidonic.acid..AA..Omega6.mg." = "Arachidonic Acid (AA) (Omega-6) (mg/100mL)",
  "Docosahexaenoic.acid..DHA..Omega3.mg." = "Docosahexaenoic Acid (DHA) (Omega-3) (mg/100mL)", "Carbohydrates..gr." = "Carbohydrates (gr/100mL)", "Sugars..gr." = "Sugars (gr/100mL)",
  "Lactose..gr." = "Lactose (gr/100mL)", "Glucose..gr." = "Glucose (gr/100mL)", "Maltose..gr." = "Maltose (gr/100mL)", "Inositol..mg." = "Inositol (mg/100mL)",
  "Polysaccharides..gr." = "Polysaccharides (gr/100mL)", "Fiber..gr." = "Fiber (gr/100mL)", "GOS..gr." = "GOS (gr/100mL)", "FOS..gr." = "FOS (gr/100mL)",
  "X2..FL..Oligosaccharide.2..FL..mg." = "X2-FL Oligosaccharide (mg/100mL)", "X3.GL..3..Galactosyllactose..mg." = "X3-GL Galactosyllactose (mg/100mL)",
  "Proteins..gr." = "Proteins (gr/100mL)", "Casein..gr." = "Casein (gr/100mL)", "Whey..gr." = "Whey (gr/100mL)", "Whey.Casein" = "Whey/Casein Ratio",
  "Salt..gr." = "Salt (gr/100mL)", "Vitamin.A..Î.g." = "Vitamin A (µg/100mL)", "Vitamin.D..Î.g." = "Vitamin D (µg/100mL)", "Vitamin.E..mg." = "Vitamin E (mg/100mL)",
  "Vitamin.K..Î.g." = "Vitamin K (µg/100mL)", "Vitamin.C..mg." = "Vitamin C (mg/100mL)", "Vitamin.B1..mg." = "Vitamin B1 (mg/100mL)", "Vitamin.B2..mg." = "Vitamin B2 (mg/100mL)",
  "Niacin..mg." = "Niacin (mg/100mL)", "Vitamin.B6..mg." = "Vitamin B6 (mg/100mL)", "Folic.acid..Î.g." = "Folic Acid (µg/100mL)", "Vitamin.B12..Î.g." = "Vitamin B12 (µg/100mL)",
  "Biotin..Î.g." = "Biotin (µg/100mL)", "Pantothenic.acid..mg." = "Pantothenic Acid (mg/100mL)", "Folate..Î.g." = "Folate (µg/100mL)", "Natrium..mg." = "Natrium (mg/100mL)",
  "Potassium..mg." = "Potassium (mg/100mL)", "Chloride..mg." = "Chloride (mg/100ml)", "Calcium..mg." = "Calcium (mg/100ml)", "Phosphorus..mg." = "Phosphorus (mg/100ml)",
  "Magnesium..mg." = "Magnesium (mg/100ml)", "Iron..mg." = "Iron (mg/100ml)", "Zinc..mg." = "Zinc (mg/100ml)", "Copper..mg." = "Copper (mg/100ml)",
  "Manganese..mg." = "Manganese (mg/100ml)", "Fluoride..mg." = "Fluoride (mg/100ml)", "Selenium..Î.g." = "Selenium (µg/100ml)", "Iodine..Î.g." = "Iodine (µg/100ml)",
  "Molybdenum..Î.g." = "Molybdenum (µg/100ml)", "Chromium..Î.g." = "Chromium (µg/100ml)", "Choline..mg." = "Choline (mg/100ml)", "Taurine..mg." = "Taurine (mg/100ml)",
  "L.carnitine.mg." = "L-carnitine (mg/100ml)", "Nucleotides..mg." = "Nucleotides (mg/100ml)"
)

#Changing the column names for the continuous variables 

colnames(filtered_formula_compositions) <- ifelse(colnames(filtered_formula_compositions) %in% names(display_names),display_names[colnames(filtered_formula_compositions)], colnames(filtered_formula_compositions))


#Merge the two data frames based on the formula column

merged_final <- merge(filtered_formula_compositions, merged_data_combined2, by = "formula")

#Changing the order of the columns 

merged_final <- merged_final %>% dplyr::select(NG_ID, timepoint, formula, everything())


# Reshape the dataframe into long format

merged_final_long <- pivot_longer(merged_final, cols = -c(NG_ID, timepoint, formula), names_to = "variable", values_to = "value")


# Separate the median group columns
merged_final_long <- merged_final_long %>%
  mutate(median_group = grepl("_median group", variable)) %>%
  mutate(variable = gsub("_median group", "", variable))


# Identify categorical columns
categorical_columns <- merged_final_long %>%
  filter(grepl("_categorical$", variable)) %>%
  mutate(variable = gsub("_categorical$", "", variable))

# Add categorical values to the main dataframe
merged_final_long <- merged_final_long %>%
  left_join(categorical_columns, by = c("NG_ID", "timepoint", "formula", "variable"), suffix = c("", "_categorical")) %>%
  mutate(categorical = ifelse(!is.na(value_categorical), value_categorical, NA)) %>%
  dplyr::select(-value_categorical)

# Filter out the original categorical rows
merged_final_long <- merged_final_long %>%
  filter(!grepl("_categorical$", variable))

# Pivot the dataframe to get the desired long format
merged_final_pivot <- merged_final_long %>%
  pivot_wider(names_from = median_group, values_from = value, names_prefix = "group_") %>%
  rename(Continuous = group_FALSE, `median group` = group_TRUE)


#Delete median_group_categorical column

merged_final_pivot <- merged_final_pivot %>% dplyr::select(-median_group_categorical)


#Change the order of the columns 

merged_final_pivot <- merged_final_pivot %>% dplyr::select(-categorical, categorical)



# Save final long-format dataframe
saveRDS(merged_final_pivot, file = "output/df_NG_ID_continuous_median_long.rds")












