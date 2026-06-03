# 📊 The signature of formula components on the gut microbiome of exclusively formula-fed infants in the Lifelines NEXT cohort


## Overview

The repository includes the full pipeline used to explore how infant formula composition 
is associated with: 

-Gut microbiome diversity (shannon index, richness)
-Bacterial species abundance 


This workflow includes:
  
- **Data preparation** from questionnaires and formula composition tables 
- **Integration** of phenotype, metagenomic, and formula data  
- **Computation** of diversity metrics 
- **Linear mixed-effects modeling** to assess associations between:
  - Nutrient concentrations (continuous and median-grouped)
  - Microbiome diversity (Shannon, richness)
  - Species abundance
  

---

## 📜 Analysis Scripts

### 🧹 Data Preparation

| Script | Description |
|--------|-------------|
| `01_prepare_formula_data.R` | Loads and prepares formula data from questionnaires |
| `02_add_formula_composition.R` | Adds nutrient composition to formula data |
| `03_clean_and_reshape.R` | Cleans and reshapes data for merging |
| `04_merge_metadata.R` | Merges questionnaire and metadata |
| `05_prepare_analysis_dataset.R` | Finalizes dataset for modeling |
| `06_merge_species_filtered.R` | Merges filtered species abundance data |
| `07_add_richness.R` | Adds richness estimates to the dataset|

### 📈 Modeling and Analysis

| Script | Description |
|--------|-------------|
| `01_clean_and_merge_for_modeling.R` | Prepares final dataset with phenotype and diversity data |
| `02_run_mixed_models_shannon_continuous.R` | Runs mixed models: Shannon ~ Continuous nutrient concentrations |
| `03_run_mixed_models_shannon_median_group.R` | Runs mixed models: Shannon ~ Median-grouped nutrient values |
| `04_run_mixed_models_richness_continuous.R` | Runs mixed models: Richness ~ Continuous nutrient concentrations |
| `05_run_mixed_models_richness_median_group.R` | Runs mixed models: Richness ~ Median-grouped nutrient values |
| `06_run_mixed_models_species_continuous.R` | Runs mixed models: Species abundance ~ Continuous nutrient concentrations |
| `07_run_mixed_models_species_median_group.R` | Runs mixed models: Species abundance ~ Median-grouped nutrient values |


Each script performs: 

-Removal of outliers
-Z-score data transformation 
-Mixed-effects model fitting 
-Likelihood ratio testing 
-BH FDR correction 
-Saving results for downstream generation of plots 


---

## ▶️ How to Run

1. Start with the data preparation scripts (`01_prepare_formula_data.R` to `07_add_richness.R`)
2. Then run the modeling scripts depending on your analysis goal (e.g., Shannon, richness, species)
3. Results will be saved in the `output/` directory.

---

## 📦 Dependencies

This project uses the following R packages:

- `tidyverse`
- `here`
- `broom` / `broom.mixed`
- `lmerTest`
- `MuMIn`
- `parallel`
- `forcats`

Install them using:

```r
install.packages(c(
  "tidyverse", "here", "broom", "lmerTest", "MuMIn", "parallel", "forcats"
))

