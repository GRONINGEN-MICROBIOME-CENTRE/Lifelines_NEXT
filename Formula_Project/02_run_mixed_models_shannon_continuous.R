# 02_run_mixed_models_UPDATED.R

# Load libraries
library(here)
library(tidyverse)
library(broom.mixed)
library(lmerTest)
library(MuMIn)

# Load cleaned dataset
data2 <- readRDS(here("scripts/GITHUB/output/NG_ID_shannon_richness_continuous_median_species_birth_mode.rds"))

# Define nutrient variables
names_components <- names(data2) %>% keep(~ str_detect(.x, "continuous"))

# Function to run and save mixed model results
run_mixed_mod <- function(df, component) {
  print(paste("Processing:", component))
  
  data_to_model <- df %>%
    dplyr::select(all_of(component), Timepoint_categorical, shannon, dna_conc, NEXT_ID, clean_reads_FQ_1, birth_deliverybirthcard_mode_binary, BATCH_NUMBER) %>%
    mutate(
      Timepoint_categorical = as.factor(Timepoint_categorical),
      NEXT_ID = as.factor(NEXT_ID)
    ) %>%
    filter(!if_any(everything(), is.na)) %>%
    mutate(across(all_of(component), scale)) %>%
    mutate_if(is.factor, fct_drop)
  
  m1 <- tryCatch({
    model <- lmerTest::lmer(
      as.formula(paste("shannon ~ Timepoint_categorical + dna_conc + BATCH_NUMBER + clean_reads_FQ_1 + birth_deliverybirthcard_mode_binary + (1|NEXT_ID)")),
      data = data_to_model
    )
    list(success = TRUE, model = model)
  }, error = function(e) list(success = FALSE, error_message = e$message, model = NULL))
  
  m2 <- tryCatch({
    model <- lmerTest::lmer(
      as.formula(paste("shannon ~", component, "* Timepoint_categorical + dna_conc + BATCH_NUMBER + clean_reads_FQ_1 + birth_deliverybirthcard_mode_binary + (1|NEXT_ID)")),
      data = data_to_model
    )
    list(success = TRUE, model = model)
  }, error = function(e) list(success = FALSE, error_message = e$message, model = NULL))
  
  model_comparison <- NULL
  model_comparison_tidied <- NULL
  
  if (m1$success & m2$success) {
    model_comparison <- tryCatch({
      anova(m1$model, m2$model, test = "LRT")
    }, error = function(e) NULL)
    
    if (!is.null(model_comparison)) {
      model_comparison_tidied <- broom.mixed::tidy(model_comparison) %>%
        mutate(component = component)
    }
  }
  
  r2_marginal <- NA
  r2_conditional <- NA
  
  if (m2$success) {
    r2_values <- MuMIn::r.squaredGLMM(m2$model)
    r2_marginal <- r2_values[1, "R2m"]
    r2_conditional <- r2_values[1, "R2c"]
  }
  
  model_results <- list(
    component = component,
    m1 = m1,
    m2 = m2,
    model_comparison = model_comparison,
    model_comparison_tidied = model_comparison_tidied,
    r2_marginal = r2_marginal,
    r2_conditional = r2_conditional
  )
}

# Run models
model_results <- list()

for (i in names_components) {
  component <- i
  df <- data2
  
  tryCatch({
    model_results[[component]] <- run_mixed_mod(df, component)
  }, error = function(e) {
    message(paste("Error in iteration:", e$message))
    return(NULL)
  })
}

# Extract results
extract_m2_results <- function(results_list) {
  extracted_results <- list()
  
  for (i in seq_along(results_list)) {
    result <- results_list[[i]]
    
    if (!is.null(result) && result$m2$success) {
      model <- result$m2$model
      
      coef_df <- broom.mixed::tidy(model) %>%
        mutate(component = as.character(attr(model@frame, "names")[2]))
      
      model_aic <- AIC(model)
      model_npar <- attr(logLik(model), "df")
      
      model_pvalue <- if (!is.null(result$model_comparison_tidied) &&
                          is.data.frame(result$model_comparison_tidied)) {
        result$model_comparison_tidied$p.value[2]
      } else {
        NA_real_
      }
      
      coef_df <- coef_df %>%
        mutate(AIC = model_aic,
               npar = model_npar,
               p_value = model_pvalue)
      
      extracted_results[[i]] <- coef_df
    }
  }
  
  final_df <- dplyr::bind_rows(extracted_results)
  return(final_df)
}

res_continuous_shannon <- extract_m2_results(model_results)


# Save full results
output_dir <- here("scripts/GITHUB/output/association_analysis_results_shannon_corrected_birth_mode")
dir.create(output_dir, showWarnings = FALSE)
saveRDS(res_continuous_shannon, file = file.path(output_dir, "df_association_analysis_continuous_shannon_corrected_birth_mode_standardized_final.rds"))


# Filter significant results
res_continuous_shannon_sig_BH <- res_continuous_shannon %>%
  group_by(component) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  filter(p_adj < 0.05)

# Save significant results
saveRDS(res_continuous_shannon_sig_BH, file = file.path(output_dir, "df_significant_results_BH_corrected.rds"))
