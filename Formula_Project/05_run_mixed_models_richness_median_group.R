# 05_run_mixed_models_richness_median_group.R

# Load libraries
library(here)
library(tidyverse)
library(broom)
library(broom.mixed)
library(lmerTest)
library(MuMIn)
library(parallel)
library(forcats)

# Load data
data2 <- readRDS(here("scripts/GITHUB/output/NG_ID_shannon_richness_continuous_median_species_birth_mode.rds"))


# Identify median nutrient components
names_components <- names(data2)[grepl("^median_", names(data2))]

# Exclude specific columns
exclude_columns <- c(
  "median_group_X3_galactosyllactose_mg",
  "median_group_casein_g"
)

names_components <- setdiff(names_components, exclude_columns)

# ------------------------------------------------------------
# Mixed model function
# ------------------------------------------------------------

run_mixed_mod <- function(df, component) {
  message("Processing: ", component)
  
  # Select variables
  data_selected <- df %>%
    select(
      all_of(component),
      Timepoint_categorical,
      richness,
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
      z_score_richness = scale(richness) %>% as.numeric()
    ) %>%
    filter(abs(z_score_richness) <= 3.5)
  
  # Null model
  m1 <- tryCatch({
    lmer(
      richness ~ Timepoint_categorical + dna_conc + BATCH_NUMBER +
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
          "richness ~ ", component,
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

res_median_richness <- extract_m2_results(model_results)
res_median_richness_sig <- res_median_richness %>% filter(p_value < 0.05)

saveRDS(res_median_richness, file = "df_association_analysis_median_richness_removal_richness_outliers.rds")

# Reload from final location
res_median_richness <- readRDS(
  here("dataframes/association_analysis_results_richness_birth_mode/df_association_analysis_median_richness_removal_richness_outliers.rds")
)

# ------------------------------------------------------------
# FDR correction
# ------------------------------------------------------------

res_median_richness_sig_BH <- res_median_richness %>%
  group_by(component) %>%
  slice(2) %>%   # second row = nutrient effect
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Extract only FDR-significant components
res_median_richness_sig_BH2 <- res_median_richness_sig_BH %>%
  select(component, p_adj)

# Merge back into full results
res_median_plots <- res_median_richness %>%
  left_join(res_median_richness_sig_BH2, by = "component") %>%
  filter(p_adj < 0.05)

# Save complete dataframe
saveRDS(
  res_median_plots,
  file = "df_association_analysis_median_richness_removal_richness_outliers_FDR_complete.rds"
)

# ------------------------------------------------------------
# Save Excel output
# ------------------------------------------------------------

library(writexl)

data <- readRDS(
  here("dataframes/association_analysis_results_richness_birth_mode/df_association_analysis_median_richness_removal_richness_outliers.rds")
)

write_xlsx(data, "df_association_analysis_median_richness_removal_richness_outliers.xlsx")

# Print sig results
res_median_richness_sig <- res_median_richness %>% filter(p_value < 0.05)
print(res_median_richness_sig, n = nrow(res_median_richness_sig))
