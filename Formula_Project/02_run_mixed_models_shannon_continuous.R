# 02_run_mixed_models_shannon_continuous.R


# Load libraries

library(here)
library(tidyverse)
library(broom)
library(broom.mixed)
library(lmerTest)
library(MuMIn)
library(parallel)

# ------------------------------------------------------------
# Load cleaned dataset
# ------------------------------------------------------------

data2 <- readRDS(here("scripts/GITHUB/output/NG_ID_shannon_richness_continuous_median_species_birth_mode.rds"))

# Identify nutrient components
names_components <- names(data2) %>% keep(~ str_detect(.x, "continuous"))
names_components <- setdiff(names_components, "continuous_X3_galactosyllactose_mg")

# ------------------------------------------------------------
# Mixed model function (updated)
# ------------------------------------------------------------

run_mixed_mod <- function(df, component) {
  message("Processing: ", component)
  
  # Select variables
  data_selected <- df %>%
    select(
      all_of(component),
      Timepoint_categorical,
      shannon,
      dna_conc,
      NEXT_ID,
      clean_reads_FQ_1,
      birth_deliverybirthcard_mode_binary,
      BATCH_NUMBER
    ) %>%
    filter(!if_any(everything(), is.na)) %>%
    mutate(across(where(is.factor), fct_drop))
  
  # Z-score + outlier removal
  data_to_model <- data_selected %>%
    mutate(
      Timepoint_categorical = as.factor(Timepoint_categorical),
      NEXT_ID = as.factor(NEXT_ID),
      z_score_shannon = scale(shannon) %>% as.numeric(),
      !!sym(paste0("z_score_", component)) := scale(.data[[component]]) %>% as.numeric()
    ) %>%
    filter(
      abs(z_score_shannon) <= 3,
      abs(.data[[paste0("z_score_", component)]]) <= 3
    )
  
  # Null model
  m1 <- tryCatch({
    lmer(
      shannon ~ Timepoint_categorical + dna_conc + BATCH_NUMBER +
        clean_reads_FQ_1 + birth_deliverybirthcard_mode_binary +
        (1 | NEXT_ID),
      data = data_to_model,
      REML = FALSE
    )
  }, error = function(e) NULL)
  
  # Full model
  m2 <- tryCatch({
    lmer(
      as.formula(
        paste0(
          "shannon ~ z_score_", component,
          " + Timepoint_categorical + dna_conc + BATCH_NUMBER + ",
          "clean_reads_FQ_1 + birth_deliverybirthcard_mode_binary + (1|NEXT_ID)"
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
    broom.mixed::tidy(model_comparison) %>% mutate(component = component)
  } else NULL
  
  # R2 values
  r2_values <- if (!is.null(m2)) MuMIn::r.squaredGLMM(m2) else NA
  r2_marginal <- if (!is.null(r2_values)) r2_values[1, "R2m"] else NA
  r2_conditional <- if (!is.null(r2_values)) r2_values[1, "R2c"] else NA
  
  list(
    component = component,
    m1 = m1,
    m2 = m2,
    model_comparison = model_comparison,
    model_comparison_tidied = model_comparison_tidied,
    r2_marginal = r2_marginal,
    r2_conditional = r2_conditional
  )
}

# ------------------------------------------------------------
# Run models
# ------------------------------------------------------------

model_results <- list()

for (component in names_components) {
  tryCatch({
    model_results[[component]] <- run_mixed_mod(data2, component)
  }, error = function(e) {
    message("Error in component ", component, ": ", e$message)
  })
}

# ------------------------------------------------------------
# Extract model results
# ------------------------------------------------------------

extract_m2_results <- function(results_list) {
  extracted <- list()
  
  for (res in results_list) {
    if (!is.null(res$m2)) {
      model <- res$m2
      
      coef_df <- broom.mixed::tidy(model) %>%
        mutate(component = res$component)
      
      model_aic <- AIC(model)
      model_npar <- attr(logLik(model), "df")
      
      model_pvalue <- if (!is.null(res$model_comparison_tidied)) {
        res$model_comparison_tidied$p.value[2]
      } else NA_real_
      
      coef_df <- coef_df %>%
        mutate(
          AIC = model_aic,
          npar = model_npar,
          p_value = model_pvalue
        )
      
      extracted[[length(extracted) + 1]] <- coef_df
    }
  }
  
  bind_rows(extracted)
}

res_continuous_shannon <- extract_m2_results(model_results)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

output_dir <- here("scripts/GITHUB/output/association_analysis_results_shannon_corrected_birth_mode")
dir.create(output_dir, showWarnings = FALSE)

saveRDS(
  res_continuous_shannon,
  file = file.path(output_dir, "df_association_analysis_continuous_shannon_corrected_birth_mode_standardized_final.rds")
)

# BH correction
res_continuous_shannon_sig_BH <- res_continuous_shannon %>%
  group_by(component) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  filter(p_adj < 0.05)

saveRDS(
  res_continuous_shannon_sig_BH,
  file = file.path(output_dir, "df_significant_results_BH_corrected.rds")
)
