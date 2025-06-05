# 06_run_mixed_models_species_continuous.R

# Load libraries
library(here)
library(tidyverse)
library(broom.mixed)
library(lmerTest)
library(MuMIn)
library(dplyr)

# Load data
data2 <- readRDS(here("scripts/GITHUB/output/NG_ID_shannon_richness_continuous_median_species_birth_mode.rds"))

# Define variables
names_components <- names(data2)[grepl("^continuous_", names(data2))]
bacteria_columns <- names(data2)[grepl("^\\.s__", names(data2))]

# Function to run mixed models
run_mixed_mod <- function(df, component, bacteria_col) {
  message("Processing: ", component, " ~ ", bacteria_col)
  
  data_to_model <- df %>%
    dplyr::select(all_of(component), Timepoint_categorical, dna_conc, NEXT_ID, clean_reads_FQ_1, birth_deliverybirthcard_mode_binary, BATCH_NUMBER, all_of(bacteria_col)) %>%
    mutate(across(c(Timepoint_categorical, NEXT_ID), as.factor)) %>%
    filter(!if_any(everything(), is.na)) %>%
    mutate(across(all_of(component), scale)) %>%
    mutate_if(is.factor, fct_drop)
  
  m1 <- tryCatch({
    lmer(as.formula(paste(bacteria_col, "~ Timepoint_categorical + dna_conc + BATCH_NUMBER + clean_reads_FQ_1 + birth_deliverybirthcard_mode_binary + (1|NEXT_ID)")), data = data_to_model)
  }, error = function(e) NULL)
  
  m2 <- tryCatch({
    lmer(as.formula(paste(bacteria_col, "~", component, "* Timepoint_categorical + dna_conc + BATCH_NUMBER + clean_reads_FQ_1 + birth_deliverybirthcard_mode_binary + (1|NEXT_ID)")), data = data_to_model)
  }, error = function(e) NULL)
  
  if (!is.null(m1) && !is.null(m2)) {
    comparison <- tryCatch(anova(m1, m2, test = "LRT"), error = function(e) NULL)
    r2 <- tryCatch(MuMIn::r.squaredGLMM(m2), error = function(e) NULL)
    
    list(
      component = component,
      bacteria_col = bacteria_col,
      model = m2,
      comparison = comparison,
      r2_marginal = if (!is.null(r2)) r2[1, "R2m"] else NA,
      r2_conditional = if (!is.null(r2)) r2[1, "R2c"] else NA,
      p_value = if (!is.null(comparison)) broom.mixed::tidy(comparison)$p.value[2] else NA
    )
  } else {
    NULL
  }
}

# Run models
model_results <- list()
for (component in names_components) {
  for (bacteria_col in bacteria_columns) {
    result <- tryCatch(run_mixed_mod(data2, component, bacteria_col), error = function(e) NULL)
    if (!is.null(result)) {
      model_results[[paste(component, bacteria_col, sep = "_")]] <- result
    }
  }
}

# Save raw model results
saveRDS(model_results, file = here("scripts/GITHUB/output/model_results_continuous_species_corrected_birth_mode_standardized_final.rds"))

# Extract results
res_continuous <- map_dfr(model_results, function(res) {
  if (!is.null(res$model)) {
    broom.mixed::tidy(res$model) %>%
      mutate(
        component = res$component,
        bacteria_col = res$bacteria_col,
        AIC = AIC(res$model),
        npar = attr(logLik(res$model), "df"),
        p_value = res$p_value,
        r2_marginal = res$r2_marginal,
        r2_conditional = res$r2_conditional
      )
  }
})

# Save full results
dir.create(here("scripts/GITHUB/output/association_analysis_results_species_corrected_birth_mode"), showWarnings = FALSE)
saveRDS(res_continuous, file = here("scripts/GITHUB/output/association_analysis_results_species_corrected_birth_mode", "df_association_analysis_continuous_species_corrected_birth_mode_standardized_final.rds"))

# BH correction
res_continuous_sig_BH <- res_continuous %>%
  group_by(component, bacteria_col) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  filter(p_adj < 0.05)

# Save significant results
saveRDS(res_continuous_sig_BH, file = here("scripts/GITHUB/output/df_res_continuous_species_sig_BH_corrected_birth_mode_standardized_final.rds"))

# Add BH results back to full results for plotting
res_combined <- res_continuous %>%
  left_join(res_continuous_sig_BH %>% dplyr::select(component, bacteria_col, p_adj), by = c("component", "bacteria_col")) %>%
  filter(p_adj < 0.05)

saveRDS(res_combined, file = here("scripts/GITHUB/output/df_res_continuous_species_sig_BH_corrected_birth_mode_for_plots_standardized_final.rds"))
