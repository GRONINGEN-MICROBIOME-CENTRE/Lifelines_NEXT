# here glmmPen will choose the lambda grid.

# Changes to this script include scaling the model matrix to make predictors on the same scale.
library(tidyverse)
library(MOFA2)
library(here)

# track time
start_time <- Sys.time()
cat("Script started at:", format(start_time), "\n")

#data_pheno <- read_delim(file=here("all_data/data_phenotypes/infants_MEFISTO_cross_sectional_phenotypes_max_NA_40_5_variance.txt"))
data_pheno <- read_delim(file=here("all_data/data_phenotypes/infants_MEFISTO_cross_sectional_phenotypes_max_NA_15_5_2428_samples.txt"))

data_pheno <- data_pheno %>% rename(NEXT_ID = next_id_infant)


k3 <- load_model(here("all_data/MEFISTO_K3_FULL.hdf5"), load_interpol_Z = FALSE)

Zmat <- as.data.frame(k3@expectations$Z$single_group)

Zmat$NEXT_ID <- rownames(Zmat)

data_Z <- as_tibble(Zmat) %>% 
  mutate(Timepoint_categorical = str_extract(NEXT_ID, "_.*") %>% str_remove("^_")) %>%
  mutate(NEXT_ID = str_remove(NEXT_ID, "_.*")) %>%
  mutate(
    Timepoint_categorical = fct_relevel(
      Timepoint_categorical,
      "W2", "M1", "M2", "M3", "M6", "M9", "M12"
    )
  )

data_full <- 
  data_Z %>%
  left_join(., data_pheno, by = c("NEXT_ID", "Timepoint_categorical")) %>%
  filter(!is.na(NG_ID)) # filter out those samples that do not have metadata because they were missing (and never sequenced).

# across all phenotypes there are no complete cases. 
phenotype_names <- names(data_full)[11:length(names(data_full))]
covariates <- names(data_full)[c(7,8,9,10)] 

# scale time ----
# Below is the old scaling of time
# data_complete_cases <-
#   data_full %>%
#   filter(complete.cases(.)) %>%
#   mutate(
#     age_c  = scale(exact_age_months_at_collection, center = TRUE, scale = FALSE)[,1],
#     age_c2 = age_c^2,
#     age_c3 = age_c^3
#   )

data_complete_cases <-
data_full %>%
  filter(complete.cases(.)) %>%
  mutate(
    age_c = poly(data_full$exact_age_months_at_collection,3)[,1], 
    age_c2 = poly(data_full$exact_age_months_at_collection,3)[,2], 
    age_c3 = poly(data_full$exact_age_months_at_collection,3)[,3]) %>%
  mutate(
    across(c(age_c, age_c2, age_c3), scale)
  )

############################################################
## Complete glmmPen pipeline 
############################################################

## -----------------------------
## 0. Libraries
## -----------------------------
library(glmmPen)
library(parallel)
library(dplyr)

## -----------------------------
## 1. Responses
## -----------------------------
responses <- c("Factor1", "Factor2", "Factor3")

## -----------------------------
## 2. Predictors (fixed effects) [as input for scaling function]
## -----------------------------
model_variables <- setdiff(
  c(phenotype_names, covariates, "age_c", "age_c2", "age_c3"),
  "exact_age_months_at_collection"
)

## -----------------------------
## 3. Standardize variables so (penalized) effects sizes of glmmPen are comparable
## -----------------------------
source(here("all_data/scripts/helper_functions/prep_and_scale_model_predictors.R")) # time data is already going into the function scaled.

head(data_scaled)

## -----------------------------
## 4. Model Variables (fixed effects), as input for glmmPen
## -----------------------------

model_variables_glmmPen <- setdiff(
  colnames(data_scaled),
  c("NEXT_ID", "NG_ID", responses)   # drop grouping + responses
)

# (optional) sanity check: ensure time columns are present
stopifnot(all(c("age_c","age_c2","age_c3") %in% model_variables))


## -----------------------------
## 5. Alpha grid (you still choose alpha)
## -----------------------------
# If you were previously using names(lambdas), just define the alpha values directly:
alpha_values <- c(1, 0.5)   # <-- put whatever set you want


## -----------------------------
## 6. glmmPen wrapper (NO lambda0_seq provided)
## -----------------------------
run_glmmPen <- function(response, alpha_value) {
  
  # build formula using SCALED predictors
  formula <- as.formula(
    paste(
      response, "~",
      paste(model_variables_glmmPen, collapse = " + "),
      "+ (1 | NEXT_ID)"
    )
  )
  
  glmmPen(
    formula = formula,
    data    = data_scaled,     # <--- USE SCALED DATA
    family  = "gaussian",
    
    penalty = "lasso",
    alpha   = alpha_value,
    
    optim_options = optimControl(
      var_restriction = "fixef",
      nMC_start = 250,
      nMC_max   = 1000,
      maxitEM   = 50
    ),
    
    tuning_options = selectControl(
      # lambda0_seq omitted on purpose -> glmmPen chooses grid
      lambda1_seq = 0,
      search = "abbrev",
      BIC_option = "BIC",
      logLik_calc = TRUE
    ),
    
    progress = FALSE
  )
}

#run_glmmPen("Factor1", "alpha_1")


## -----------------------------
## 7. Run in parallel
## -----------------------------
n_cores <- 4

tasks <- expand.grid(
  response   = responses,
  alpha_name = alpha_values,
  stringsAsFactors = FALSE
)

fits_list <- mclapply(
  seq_len(nrow(tasks)),
  function(i) run_glmmPen(tasks$response[i], tasks$alpha_name[i]),
  mc.cores = n_cores
)

names(fits_list) <- paste(tasks$response, tasks$alpha_name, sep = "__")

# split by alpha_name first
fits <- split(fits_list, tasks$alpha_name)

# within each alpha, split by response
fits <- lapply(
  fits,
  function(alpha_block) {
    # names look like "Factor1__alpha_0.5" -> extract response part
    resp_names <- sub("__.*$", "", names(alpha_block))
    setNames(split(alpha_block, resp_names), unique(resp_names))
  }
)

saveRDS(
  fits,
  file = here(
    paste0(
      "results_glmmPen/results_glmmPen_chooses_lambda",
      format(Sys.Date(), "%d-%m-%Y"),
      ".rds"
    )
  )
)

# end time
end_time <- Sys.time()
cat("Script finished at:", format(end_time), "\n")
cat("Total runtime:", difftime(end_time, start_time, units = "mins"), "\n")

fits <- readRDS(file = here("results_glmmPen/results_glmmPen_chooses_lambda31-01-2026.rds")) # time IS scaled with more MCMC draws



## -----------------------------
## 8. Convergence check
## -----------------------------
check_pglmm_convergence <- function(fits, conv_tol = 0.0015) {
  out <- list()
  k <- 0L

  for (fit_i in seq_along(fits)) {
    fit_set <- fits[[fit_i]]

    # e.g. Factor1, Factor2, ...
    factor_names <- names(fit_set)
    if (is.null(factor_names)) next

    for (fac in factor_names) {
      fac_obj <- fit_set[[fac]]
      if (is.null(fac_obj) || is.null(names(fac_obj))) next

      # e.g. Factor1__alpha_0.5, Factor1__alpha_0.1, ...
      model_names <- names(fac_obj)

      for (m in model_names) {
        model <- fac_obj[[m]]

        opt <- tryCatch(model$optinfo, error = function(e) NULL)

        iter <- if (!is.null(opt$iter)) opt$iter else NA_integer_
        conv <- if (!is.null(opt$conv)) opt$conv else NA_real_

        warn <- NA_character_
        if (!is.null(opt$warnings)) {
          if (length(opt$warnings) == 0) {
            warn <- NA_character_
          } else {
            warn <- paste(opt$warnings, collapse = " | ")
          }
        }

        # Simple "converged?" rule:
        # - conv exists AND conv <= tol
        # - and no warnings (if warnings field exists)
        converged <- !is.na(conv) && conv <= conv_tol && (is.na(warn) || warn == "")

        k <- k + 1L
        out[[k]] <- data.frame(
          fits_index = fit_i,
          factor = fac,
          model = m,
          iter = iter,
          conv = conv,
          conv_tol = conv_tol,
          warnings = warn,
          converged = converged,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  res <- do.call(rbind, out)
  if (is.null(res)) {
    res <- data.frame(
      fits_index = integer(),
      factor = character(),
      model = character(),
      iter = integer(),
      conv = numeric(),
      conv_tol = numeric(),
      warnings = character(),
      converged = logical(),
      stringsAsFactors = FALSE
    )
  }
  res[order(res$fits_index, res$factor, res$model), ]
}

# usage:
conv_tbl <- check_pglmm_convergence(fits, conv_tol = 0.0015)
conv_tbl

# quick checks:
all(conv_tbl$converged)
subset(conv_tbl, !converged)



print(converged)

lapply(fits, function(f) f$optinfo$warnings)


# ## -----------------------------
# ## 9. Extract selected variables
# ## -----------------------------
stability_table <- dplyr::bind_rows(
  lapply(names(fits), function(alpha_name) {
    dplyr::bind_rows(
      lapply(names(fits[[alpha_name]]), function(response) {

        fit <- fits[[alpha_name]][[response]][[1]]  # pglmmObj

        fixef <- fit$fixef
        n_nonzero <- sum(fixef[names(fixef) != "(Intercept)"] != 0)

        # lambda selected (robust to structure)
        lambda0_selected <- NA_real_
        if (!is.null(fit$results_optim) && "lambda0" %in% colnames(fit$results_optim)) {
          lambda0_selected <- as.numeric(fit$results_optim[1, "lambda0"])
        }

        data.frame(
          alpha = alpha_name,
          response = response,
          n_nonzero_fixef = n_nonzero,
          lambda0_selected = lambda0_selected,
          stringsAsFactors = FALSE
        )
      })
    )
  })
)

print(stability_table)

# ## -----------------------------
# ## 10. Stability summaries
# ## -----------------------------

stability_table <- dplyr::bind_rows(
  lapply(names(fits), function(alpha_name) {
    dplyr::bind_rows(
      lapply(names(fits[[alpha_name]]), function(response) {

        fit <- fits[[alpha_name]][[response]][[1]]  # pglmmObj

        fixef <- fit$fixef
        n_nonzero <- sum(fixef[names(fixef) != "(Intercept)"] != 0)

        # lambda selected 
        lambda0_selected <- NA_real_
        if (!is.null(fit$results_optim) && "lambda0" %in% colnames(fit$results_optim)) {
          lambda0_selected <- as.numeric(fit$results_optim[1, "lambda0"])
        }

        data.frame(
          alpha = alpha_name,
          response = response,
          n_nonzero_fixef = n_nonzero,
          lambda0_selected = lambda0_selected,
          stringsAsFactors = FALSE
        )
      })
    )
  })
)

print(stability_table)
#
## -----------------------------
## 11. plots (fixed)
## -----------------------------
library(ggplot2)
library(tidytext)
library(dplyr)

# renaming variables
new_names <- c(
  "birth_deliverybirthcard_mode_binaryVG" = "Delivery mode: vaginal",
  "birth_deliverybirthcard_place_delivery_simplehospital" = "Delivery place: hospital",
  "infant_ffq_ever_never_breastfednever_BF" = "Feeding Mode: never BF",
  "infant_growth_birth_weight_kg" = "Birth weight",
  "infant_growth_standardized_weight_slope_kg" = "Infant weight gain slope",
  "mother_birthcardself_gestational_age_weeks" = "Gestational age",
  "clean_reads_FQ_1" = "Sequencing depth",
  "dna_conc" = "DNA concentration",
  "age_c" = "Infant age",
  "age_c2" = "Infant age\u00B2",
  "age_c3" = "Infant age\u00B3",
  "infant_misc_sexmale" = "Infant sex: male",
  "mother_deliverybirthcard_preg_comp_hypertensionyes" = "Pregnancy hypertension: yes",
  "BATCH_NUMBER" = "Batch"
)

plot_data <- dplyr::bind_rows(
  lapply(names(fits), function(alpha_name) {
    dplyr::bind_rows(
      lapply(names(fits[[alpha_name]]), function(response) {
        
        fit <- fits[[alpha_name]][[response]][[1]]
        
        lambda0 <- NA_real_
        if (!is.null(fit$results_optim) && "lambda0" %in% colnames(fit$results_optim)) {
          lambda0 <- as.numeric(fit$results_optim[1, "lambda0"])
        }
        
        data.frame(
          alpha = alpha_name,
          response = response,
          covariate = names(fit$fixef),
          effect_size = as.numeric(fit$fixef),
          lambda0 = lambda0,
          stringsAsFactors = FALSE
        )
      })
    )
  })
) %>%
  dplyr::filter(
    effect_size != 0,
    covariate != "(Intercept)"
  ) %>%
  dplyr::mutate(
    facet_label = paste0(
      sub("Factor(\\d+)", "Factor \\1", response),
      " | \u03B1 = ", alpha, "\n",
      "\u03BB = ", signif(lambda0, 3)),
    covariate_ordered = tidytext::reorder_within(
      covariate,
      effect_size,
      facet_label
    )
  )

plot_data <- plot_data %>%
  dplyr::mutate(
    covariate_nice = dplyr::recode(covariate, !!!new_names, .default = covariate),
    covariate_ordered = tidytext::reorder_within(covariate_nice, effect_size, facet_label)
  )

plot_data <- plot_data %>%
  dplyr::mutate(
    factor_num = as.integer(stringr::str_extract(response, "\\d+")),
    alpha_num  = as.numeric(gsub("[^0-9.]", "", alpha))  # works if alpha is like "0.5" or "alpha=0.5"
  ) %>%
  dplyr::mutate(
    facet_label = paste0(
      "Factor ", factor_num,
      " | \u03B1 = ", alpha_num, "\n",
      "\u03BB = ", signif(lambda0, 3)
    )
  ) %>%
  as_tibble() 


plot_data2 <- plot_data %>%
  mutate(
    factor_num = as.integer(str_extract(response, "\\d+")),
    factor_col = paste0("Factor ", factor_num),
    
    alpha_num  = as.numeric(gsub("[^0-9.]", "", alpha)),
    alpha_row  = factor(alpha_num, levels = c(0.5, 1.0),
                        labels = c("\u03B1 = 0.5", "\u03B1 = 1.0")),
    
    lambda_lab = paste0("\u03BB = ", signif(lambda0, 3)),
    
    # panel label
    panel_lab = paste0(factor_col, "\n", as.character(alpha_row), "\n", lambda_lab)
  ) %>%
  mutate(
    panel_lab = factor(
      panel_lab,
      levels = panel_lab[order(alpha_num, factor_num)] |> unique()
    ),
    
    covariate_nice = recode(covariate, !!!new_names, .default = covariate),
    
    # reorder within each panel
    covariate_ordered = tidytext::reorder_within(covariate_nice, effect_size, panel_lab)
  )



plt <- ggplot(plot_data2, aes(x = effect_size, y = covariate_ordered, fill = effect_size > 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_col(width = 0.7) +
  facet_wrap(~ panel_lab, ncol = 3, scales = "free_y") +
  tidytext::scale_y_reordered() +
  scale_fill_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"), guide = "none") +
  theme_bw(base_size = 8) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "effect size (penalized coefficient)", y = "predictor")

plt



####
date_tag <- Sys.Date()

ggsave(
  filename = here(
    paste0("all_data/figure/Supplmental_Figure_A_", date_tag, ".png")
  ),
  plot   = plt,
  width  = 180,
  height = 107.4,
  units  = "mm",
  dpi    = 600
)

ggsave(
  filename = here(
    paste0("all_data/figure/Supplmental_Figure_A_", date_tag, ".pdf")
  ),
  plot   = plt,
  width  = 180,
  height = 107.4,
  units  = "mm",
  dpi    = 600,
  device = cairo_pdf
)



# best fit
best_per_factor <- do.call(rbind, lapply(responses, function(resp) {

  best_alpha <- NA_character_
  best_bic   <- Inf
  best_lam   <- NA_real_

  for (a in names(fits)) {
    fit <- fits[[a]][[resp]][[1]]
    ro  <- fit$results_optim

    if (!is.null(ro) && "BIC" %in% colnames(ro)) {
      bic <- as.numeric(ro[1, "BIC"])
      lam <- if ("lambda0" %in% colnames(ro)) as.numeric(ro[1, "lambda0"]) else NA_real_

      if (!is.na(bic) && bic < best_bic) {
        best_bic   <- bic
        best_alpha <- a
        best_lam   <- lam
      }
    }
  }

  data.frame(
    response = resp,
    best_alpha = best_alpha,
    best_lambda0 = best_lam,
    best_BIC = best_bic,
    stringsAsFactors = FALSE
  )
}))

print(best_per_factor)

# print best model for each Factor (without heavy copying )
best_models <- setNames(vector("list", length(responses)), responses)

for (i in seq_len(nrow(best_per_factor))) {
  resp <- best_per_factor$response[i]
  a    <- best_per_factor$best_alpha[i]

  cat("\n=============================\n")
  cat("Best model for", resp, "\n")
  cat("alpha:", a, "\n")
  cat("lambda0:", best_per_factor$best_lambda0[i], "\n")
  cat("BIC:", best_per_factor$best_BIC[i], "\n")
  cat("=============================\n")

  best_models[[resp]] <- fits[[a]][[resp]][[1]]

  # If printing the full object is heavy, comment this out and print only fixef
  print(best_models[[resp]])
}

# Optional: save just the best models
# saveRDS(best_models, here("results_glmmPen/best_models_per_factor.rds"))

# Make Plots
best_per_factor <- do.call(rbind, lapply(responses, function(resp) {

  best_alpha <- NA_character_
  best_bic   <- Inf
  best_lam   <- NA_real_

  for (a in names(fits)) {
    fit <- fits[[a]][[resp]][[1]]
    ro  <- fit$results_optim

    if (!is.null(ro) && "BIC" %in% colnames(ro)) {
      bic <- as.numeric(ro[1, "BIC"])
      lam <- if ("lambda0" %in% colnames(ro)) as.numeric(ro[1, "lambda0"]) else NA_real_

      if (!is.na(bic) && bic < best_bic) {
        best_bic   <- bic
        best_alpha <- a
        best_lam   <- lam
      }
    }
  }

  data.frame(
    response = resp,
    best_alpha = best_alpha,
    best_lambda0 = best_lam,
    best_BIC = best_bic,
    stringsAsFactors = FALSE
  )
}))

print(best_per_factor)

library(dplyr)
library(tidytext)
library(ggplot2)

plot_data_best <- dplyr::bind_rows(
  lapply(seq_len(nrow(best_per_factor)), function(i) {

    resp <- best_per_factor$response[i]
    a    <- best_per_factor$best_alpha[i]

    fit <- fits[[a]][[resp]][[1]]

    lambda0 <- NA_real_
    if (!is.null(fit$results_optim) && "lambda0" %in% colnames(fit$results_optim)) {
      lambda0 <- as.numeric(fit$results_optim[1, "lambda0"])
    }

    data.frame(
      alpha = a,
      response = resp,
      covariate = names(fit$fixef),
      effect_size = as.numeric(fit$fixef),
      lambda0 = lambda0,
      stringsAsFactors = FALSE
    )
  })
) %>%
  dplyr::filter(
    effect_size != 0,
    covariate != "(Intercept)"
  ) %>%
  dplyr::mutate(
    facet_label = paste0(
      response, " | BEST: ", alpha, "\n",
      "lambda = ", signif(lambda0, 3)
    ),
    covariate_ordered = tidytext::reorder_within(
      covariate,
      effect_size,
      facet_label
    )
  )


ggplot(
  plot_data_best,
  aes(
    x = effect_size,
    y = covariate_ordered,
    fill = effect_size > 0
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_col(width = 0.7) +
  facet_wrap(~ facet_label, scales = "free_y") +
  tidytext::scale_y_reordered() +
  scale_fill_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"), guide = "none") +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Best glmmPen model per factor (selected by BIC across alphas)",
    x = "Effect size (penalized coefficient)",
    y = "Predictor"
  )




# make tables----
library(dplyr)
library(tidyr)
library(stringr)

# 1) Extract the "best" row 
extract_results_optim <- function(fits, alpha, factor_name) {
  alpha <- as.character(alpha)
  
  obj <- fits[[alpha]][[factor_name]]
  if (is.null(obj)) stop("Can't find factor: ", factor_name, " under alpha: ", alpha)
  
  # find the subkey like "Factor1__0.5"
  subkey <- grep(paste0("^", factor_name, "__", alpha, "$"), names(obj), value = TRUE)
  if (length(subkey) != 1) {
    stop("Expected exactly one subkey like '", factor_name, "__", alpha, "'. Found: ",
         paste(names(obj), collapse = ", "))
  }
  
  res <- obj[[subkey]]$results_optim
  if (is.null(res)) stop("No $results_optim found at that location.")
  
  as.data.frame(res) # ensure data.frame
}

make_penalized_coef_table_from_results <- function(res_optim,
                                                   zero_tol = 1e-12,
                                                   include_intercept = TRUE,
                                                   include_gamma0 = FALSE,
                                                   digits = 4) {
  meta_cols <- c(
    "lambda0", "lambda1",
    "BICh", "BIC", "BICq", "BICNgrp",
    "LogLik",
    "Non_0_fef", "Non_0_ref", "Non_0_coef",
    "EM_iter", "Converged"
  )
  
  coef_cols <- setdiff(names(res_optim), meta_cols)
  
  if (!include_intercept) {
    coef_cols <- setdiff(coef_cols, c("(Intercept)", "Intercept"))
  }
  if (!include_gamma0) {
    coef_cols <- setdiff(coef_cols, "Gamma0")
  }
  
  # model-level info (one row)
  model_info <- res_optim[1, intersect(meta_cols, names(res_optim)), drop = FALSE]
  
  new_names <- c(
    "birth_deliverybirthcard_mode_binaryVG" = "Delivery mode: vaginal",
    "birth_deliverybirthcard_place_delivery_simplehospital" = "Delivery place: hospital",
    "infant_ffq_ever_never_breastfednever_BF" = "Breastfed: never",
    "infant_growth_birth_weight_kg" = "Birth weight",
    "infant_growth_standardized_weight_slope_kg" = "Infant weight gain slope",
    "mother_birthcardself_gestational_age_weeks" = "Gestational age",
    "clean_reads_FQ_1" = "Sequencing depth",
    "dna_conc" = "DNA concentration",
    "age_c" = "Infant age",
    "age_c2" = "Infant age\u00B2",
    "age_c3" = "Infant age\u00B3",
    "infant_misc_sexmale" = "Infant sex: male",
    "mother_deliverybirthcard_preg_comp_hypertensionyes" = "Pregnancy hypertension: yes",
    "BATCH_NUMBER" = "Batch"
  )
  
  new_names["(Intercept)"] <- "Intercept"
  
  
  tab <- res_optim %>%
    dplyr::select(dplyr::all_of(coef_cols)) %>%
    tidyr::pivot_longer(
      cols = everything(),
      names_to = "Predictor",
      values_to = "Penalized_coefficient"
    ) %>%
    mutate(
      Predictor = dplyr::recode(Predictor, !!!new_names, .default = Predictor),
      Selected_nonzero = abs(Penalized_coefficient) > zero_tol,
      Direction = ifelse(!Selected_nonzero, "\u2014",
                         ifelse(Penalized_coefficient > 0, "+", "\u2212")),
      Penalized_coefficient = round(Penalized_coefficient, digits),
      Selected_nonzero = ifelse(Selected_nonzero, "Yes", "No")
    )
  
  tab <- dplyr::bind_cols(
    model_info[rep(1, nrow(tab)), c("lambda0", "BIC"), drop = FALSE],
    tab
  )
  
  list(
    coef_table = tab,
    model_info = model_info
  )
}


# 3) wrapper: give alpha + factor, get the table
make_table_for_fit <- function(fits, alpha, factor_name,
                               zero_tol = 1e-12,
                               include_intercept = TRUE,
                               include_gamma0 = FALSE,
                               digits = 4) {
  res <- extract_results_optim(fits, alpha = alpha, factor_name = factor_name)
  out <- make_penalized_coef_table_from_results(
    res_optim = res,
    zero_tol = zero_tol,
    include_intercept = include_intercept,
    include_gamma0 = include_gamma0,
    digits = digits
  )
  
  # add identifiers (alpha / factor) as cols
  out$coef_table <- out$coef_table %>%
    mutate(alpha = as.character(alpha), factor = factor_name) %>%
    select(factor, alpha, lambda0, BIC, Predictor, Penalized_coefficient, Selected_nonzero)
  
  out$model_info <- out$model_info %>%
    mutate(alpha = as.character(alpha), factor = factor_name)
  
  out
}

make_all_tables <- function(fits,
                            zero_tol = 1e-12,
                            include_intercept = TRUE,
                            include_gamma0 = FALSE,
                            digits = 4) {
  alphas <- names(fits)
  factors <- unique(unlist(lapply(fits, names)))  # e.g., Factor1/2/3
  
  all_coef <- list()
  all_info <- list()
  
  for (a in alphas) {
    for (f in factors) {
      if (!is.null(fits[[a]][[f]])) {
        tmp <- make_table_for_fit(
          fits, alpha = a, factor_name = f,
          zero_tol = zero_tol,
          include_intercept = include_intercept,
          include_gamma0 = include_gamma0,
          digits = digits
        )
        all_coef[[paste(a, f, sep = "_")]] <- tmp$coef_table
        all_info[[paste(a, f, sep = "_")]] <- tmp$model_info
      }
    }
  }
  
  list(
    coef_table = bind_rows(all_coef),
    model_info = bind_rows(all_info)
  )
}

all_out <- make_all_tables(fits)

clean_model_info_table <- function(model_info) {
  model_info %>%
    dplyr::select(
      factor,
      alpha,
      lambda  = lambda0,
      BIC,
      Number_selected_fef = Non_0_fef
    ) %>%
    dplyr::arrange(factor, alpha)
}

clean_model_info <- clean_model_info_table(all_out$model_info) 
print(clean_model_info)

all_out$coef_table

table_final <- all_out$coef_table 

rownames(table_final) <- NULL

# only one lambda and BIC per model
table_final2 <- table_final %>%
  dplyr::arrange(factor, alpha, Predictor) %>%
  dplyr::group_by(factor, alpha) %>%
  dplyr::mutate(
    lambda0 = if_else(
      dplyr::row_number() == 1L,
      format(signif(lambda0, 3), scientific = FALSE),
      ""
    ),
    BIC = if_else(
      dplyr::row_number() == 1L,
      format(signif(BIC, 3), scientific = FALSE),
      ""
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    lambda = lambda0,
    variable_retained = Selected_nonzero
  )



write_csv(
  table_final2,
  file = here(
    paste0(
      "all_data/final_table/Supplemental_Table_B_",
      Sys.Date(),
      ".csv"
    )
  )
)


library(kableExtra)

table_final2 %>%
  kbl(format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      caption = "Penalized model results. $\\lambda$ and BIC shown once per model.") %>%
  kable_styling(latex_options = c("hold_position"))



# two seperate tables
model_summary <- table_final2 %>%
  group_by(factor, alpha) %>%
  summarise(
    lambda = first(na.omit(lambda)),
    BIC    = first(na.omit(BIC)),
    n_retained = sum(Penalized_coefficient != 0),
    .groups = "drop"
  )

model_summary

write_csv(
  model_summary,
  file = here(
    paste0(
      "all_data/final_table/Supplemental_Table_B_model_summary",
      Sys.Date(),
      ".csv"
    )
  )
)


coef_table <- table_final2 %>%
  select(factor, alpha, Predictor, Penalized_coefficient, variable_retained)

coef_table

write_csv(
  coef_table,
  file = here(
    paste0(
      "all_data/final_table/Supplemental_Table_C_coefficient_table",
      Sys.Date(),
      ".csv"
    )
  )
)

