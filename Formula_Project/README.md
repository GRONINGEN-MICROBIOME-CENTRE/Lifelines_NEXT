# 📊 The effect of infant formula on the gut microbiome

## Overview
This project includes:
  
- **Data preparation** from questionnaire and formula composition files  
- **Integration** of phenotype data and metagenomic (species abundance) data  
- **Calculation** of diversity metrics (Shannon, richness)  
- **Linear mixed-effects modeling** to assess associations between:
  - Nutrient concentrations (continuous and median-grouped)
  - Microbiome diversity (Shannon, richness)
  - Species-level abundance

---

## 📜 Analysis Scripts

### 🧹 Data Preparation

| Script | Description |
|--------|-------------|
| `01_prepare_formula_data.R` | Loads and prepares formula data from questionnaires |
| `02_add_formula_composition.R` | Adds nutrient composition to formula data |
| `03_clean_and_reshape.R` | Cleans and reshapes data for merging |
| `04_merge_metadata.R` | Merges questionnaire and metadata |
| `05_prepare_analysis_dataset.R` | Finalizes data set for analysis |
| `06_merge_species_filtered.R` | Merges filtered species abundance data |
| `07_add_richness.R` | Merges richness to the finalized data|

### 📈 Modeling and Analysis

| Script | Description |
|--------|-------------|
| `01_clean_and_merge_for_modeling.R` | Prepares final dataset with phenotype and diversity data |
| `02_run_mixed_models.R` | Runs mixed models: Shannon ~ Continuous nutrient values |
| `03_run_mixed_models_median_group.R` | Runs mixed models: Shannon ~ Median-grouped nutrient values |
| `04_run_mixed_models_richness_continuous.R` | Runs mixed models: Richness ~ Continuous nutrient values |
| `05_run_mixed_models_richness_median_group.R` | Runs mixed models: Richness ~ Median-grouped nutrient values |
| `06_run_mixed_models_species_continuous.R` | Runs mixed models: Species abundance ~ Continuous nutrient values |
| `07_run_mixed_models_species_median_group.R` | Runs mixed models: Species abundance ~ Median-grouped nutrient values |

---

## ▶️ How to Run

1. Start with the data preparation scripts (`01_prepare_formula_data.R` to `07_add_richness.R`)
2. Then run the modeling scripts depending on your analysis goal (e.g., Shannon, richness, species)
3. Results will be saved in the `output/` folder.

---

## 📦 Dependencies

This project uses the following R packages:

- `tidyverse`
- `here`
- `broom` / `broom.mixed`
- `lmerTest`
- `MuMIn`
- `parallel`

Install them using:

```r
install.packages(c("tidyverse", "here", "broom", "lmerTest", "MuMIn", "parallel"))

