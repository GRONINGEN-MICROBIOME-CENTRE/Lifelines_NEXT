# 07_run_mixed_models_species_median_group.R

# Load libraries
library(here)
library(tidyverse)
library(broom.mixed)
library(lmerTest)
library(MuMIn)

# Load data
data2 <- readRDS(here("output/NG_ID_shannon_richness_continuous_median_species_birth_mode.rds"))

# Define variables
exclude_columns <- c("median_group_X3_GL_galactosyllactosyllactose_mg", "median_group_casein_g")
names_components <- names(data2)[grepl("^median_", names(data2)) & !names(data2) %in% exclude_columns]
bacteria_columns <- names(data2)[grepl("^\\.s__", names(data2))]

# Function to run mixed models
run_mixed_mod <- function(df, component, bacteria_col) {
  message("Processing: ", component, " ~ ", bacteria_col)
  
  data_to_model <- df %>%
    select(all_of(component), Timepoint_categorical, dna_conc, NEXT_ID, clean_reads_FQ_1, birth_deliverybirthcard_mode_binary, BATCH_NUMBER, all_of(bacteria_col)) %>%
    mutate(across(c(Timepoint_categorical, NEXT_ID), as.factor)) %>%
    filter(!if_any(everything(), is.na)) %>%
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
saveRDS(model_results, file = here("output", "model_results_median_species_corrected_birth_mode_final.rds"))

# Extract results
res_median <- map_dfr(model_results, function(res) {
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
dir.create(here("output/association_analysis_results_species_corrected_birth_mode"), showWarnings = FALSE)
saveRDS(res_median, file = here("output/association_analysis_results_species_corrected_birth_mode", "df_association_analysis_median_species_corrected_birth_mode_final.rds"))

# BH correction
res_median_sig_BH <- res_median %>%
  group_by(component, bacteria_col) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  filter(p_adj < 0.05)

# Save significant results
saveRDS(res_median_sig_BH, file = here("output", "df_res_median_species_sig_BH_corrected_birth_mode_final.rds"))

# Add BH results back to full results for plotting
res_median_plots <- res_median %>%
  left_join(res_median_sig_BH %>% select(component, bacteria_col, p_adj), by = c("component", "bacteria_col")) %>%
  filter(p_adj < 0.05)

saveRDS(res_median_plots, file = here("output", "df_res_median_species_sig_BH_corrected_birth_mode_for_plots_final.rds"))
