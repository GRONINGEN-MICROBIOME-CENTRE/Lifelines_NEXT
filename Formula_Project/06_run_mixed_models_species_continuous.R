# 06_run_mixed_models_species_continuous.R

# Load libraries
library(here)
library(tidyverse)
library(broom)
library(broom.mixed)
library(lmerTest)
library(MuMIn)
library(dplyr)
library(forcats)


# Load data
data2 <- readRDS(here("scripts/GITHUB/output/NG_ID_shannon_richness_continuous_median_species_birth_mode.rds"))


# Identify continuous nutrient components
names_components <- names(data2)[grepl("^continuous_", names(data2))]

# Identify species columns
bacteria_columns <- data2 %>% select(starts_with(".s__")) %>% names()

# ------------------------------------------------------------
# Mixed model function
# ------------------------------------------------------------

run_mixed_mod <- function(df, component, bacteria_col) {
  message("Processing: ", component, " with bacteria column: ", bacteria_col)
  
  # Prep data
  data_to_model <- df %>%
    select(
      all_of(component),
      Timepoint_categorical,
      dna_conc,
      NEXT_ID,
      clean_reads_FQ_1,
      birth_deliverybirthcard_mode_binary,
      BATCH_NUMBER,
      all_of(bacteria_col)
    ) %>%
    mutate(
      Timepoint_categorical = as.factor(Timepoint_categorical),
      NEXT_ID = as.factor(NEXT_ID)
    ) %>%
    filter(!if_any(everything(), is.na)) %>%
    mutate(
      z_score_component = scale(.data[[component]]) %>% as.numeric(),
      z_score_species = scale(.data[[bacteria_col]]) %>% as.numeric()
    ) %>%
    filter(
      abs(z_score_component) <= 3.5,
      abs(z_score_species) <= 3.5
    ) %>%
    mutate(across(where(is.factor), fct_drop))
  
  # Null model
  m1 <- tryCatch({
    lmer(
      as.formula(
        paste0(
          bacteria_col,
          " ~ Timepoint_categorical + dna_conc + BATCH_NUMBER + ",
          "clean_reads_FQ_1 + birth_deliverybirthcard_mode_binary + (1|NEXT_ID)"
        )
      ),
      data = data_to_model,
      REML = FALSE
    )
  }, error = function(e) NULL)
  
  # Full model
  m2 <- tryCatch({
    lmer(
      as.formula(
        paste0(
          bacteria_col,
          " ~ z_score_component + Timepoint_categorical + dna_conc + ",
          "BATCH_NUMBER + clean_reads_FQ_1 + birth_deliverybirthcard_mode_binary + (1|NEXT_ID)"
        )
      ),
      data = data_to_model,
      REML = FALSE
    )
  }, error = function(e) NULL)
  
  # Model comparison
  model_comparison <- if (!is.null(m1) && !is.null(m2)) {
    tryCatch(anova(m1, m2, test = "LRT"), error = function(e) NULL)
  } else NULL
  
  model_comparison_tidied <- if (!is.null(model_comparison)) {
    broom.mixed::tidy(model_comparison) %>%
      mutate(component = component, bacteria_col = bacteria_col)
  } else NULL
  
  # R2 values
  r2_values <- if (!is.null(m2)) MuMIn::r.squaredGLMM(m2) else NA
  r2_marginal <- if (!is.null(r2_values)) r2_values[1, "R2m"] else NA
  r2_conditional <- if (!is.null(r2_values)) r2_values[1, "R2c"] else NA
  
  list(
    component = component,
    bacteria_col = bacteria_col,
    m1 = m1,
    m2 = m2,
    model_comparison_tidied = model_comparison_tidied,
    p_value = if (!is.null(model_comparison_tidied)) model_comparison_tidied$p.value[2] else NA,
    r2_marginal = r2_marginal,
    r2_conditional = r2_conditional
  )
}

# ------------------------------------------------------------
# Run models
# ------------------------------------------------------------

model_results <- list()

for (component in names_components) {
  for (bacteria_col in bacteria_columns) {
    tryCatch({
      model_results[[paste(component, bacteria_col, sep = "_")]] <-
        run_mixed_mod(data2, component, bacteria_col)
    }, error = function(e) {
      message("Error in iteration: ", e$message)
      NULL
    })
  }
}

saveRDS(
  model_results,
  file = "model_results_continuous_species_corrected_birth_mode_standardized_final.rds"
)

# ------------------------------------------------------------
# Extract results
# ------------------------------------------------------------

extract_m2_results <- function(results_list) {
  extracted <- list()
  
  for (res in results_list) {
    if (!is.null(res$m2) && res$m2$success) {
      model <- res$m2$model
      
      coef_df <- broom.mixed::tidy(model) %>%
        mutate(
          component = res$component,
          bacteria_col = res$bacteria_col
        )
      
      coef_df <- coef_df %>%
        mutate(
          AIC = AIC(model),
          npar = attr(logLik(model), "df"),
          p_value = res$p_value
        )
      
      extracted[[length(extracted) + 1]] <- coef_df
    }
  }
  
  bind_rows(extracted)
}

res_continuous <- extract_m2_results(model_results)

saveRDS(
  res_continuous,
  file = "df_association_analysis_continuous_species_corrected_birth_mode_standardized_final.rds"
)

# ------------------------------------------------------------
# FDR correction
# ------------------------------------------------------------

res_continuous_sig_BH <- res_continuous %>%
  group_by(component, bacteria_col) %>%
  slice(1) %>%   # LRT row
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  filter(p_adj < 0.05)

saveRDS(
  res_continuous_sig_BH,
  file = "df_res_continuous_species_sig_BH_corrected_birth_mode_standardized_final.rds"
)

# ------------------------------------------------------------
# Add FDR results back to full results for plotting
# ------------------------------------------------------------

res_continuous_sig_BH2 <- res_continuous %>%
  group_by(component, bacteria_col) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  select(component, bacteria_col, p_adj)

res_continuous_combined <- res_continuous %>%
  left_join(res_continuous_sig_BH2, by = c("component", "bacteria_col")) %>%
  filter(p_adj < 0.05)

saveRDS(
  res_continuous_combined,
  file = "df_res_continuous_species_sig_BH_corrected_birth_mode_for_plots_standardized_final.rds"
)
