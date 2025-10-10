##########################
### AUTHOR: Cyrus A. Mallon
### Date: 4 Dec 2024
### Info: Script to replicate factor analysis used in LL-NEXT study, using demo data.
### N.B. For an exact replication of code used in study, see MEFISTO folder located in LL-NEXT github repo
### LL-NEXT/MEFISTO repo: <https://github.com/GRONINGEN-MICROBIOME-CENTRE/Lifelines_NEXT/tree/main/MEFISTO>
##########################


###Extended Info### 
#The purpose of this script is to showcase the analysis strategy used to understand how microbial community variation changes over time in the LL-NEXT cohort,
#and whether this variation is associated with mother-infant phenotypes. In order to understand how the trajectories of microbial community variation (i.e., beta diversity) 
#changes over time, we employed the MEFISTO package. MEFISTO is a multimodal factor analysis that incorporates time as a covariate into the model, allowing one to understand how variation in the dataset of interest changes over time.
#In order to associate phenotypes to factors, we employed a step wise regression strategy. More info is found in the methods of the manuscript.

#The script is split into two sections. (i) Upstream and (ii) downstream analysis. Upstream analysis involves creating the MEFISTO model. 
#Downstream analysis involves vizualizing the results of the MEFISTO model and the stepwise regression of phenotypes to factors (i.e., latent variables estimated via the MEFISTO model).


###############
####UPSTREAM###
###############

# In the downstream part of this script we will:
# (i) tidy the meta and sequencing data in preparation for MEFISTO
# (ii) create a MEFISTO object, setting MEFISTO model and training options
# (iii) run/train the MEFISTO model

# load libraries
library(MOFA2)
library(here)
library(tidyverse)
library(stringi)


# Part (i): load and tidy data----
## N.B. simulated_metadata.txt contains both meta data and phenotypic data

## load meta data
meta_data <- read.delim(here("demo/demo_data/simulated_metadata.txt"), header = TRUE)

## select infant data
meta_data <-  
  meta_data %>%
  rownames_to_column(var = "sample") %>%
  relocate(sample) %>%
  as_tibble() %>%
  filter(Type == "Infant")

n_infants <- nrow(meta_data)

## load abundance data 
## N.B. Real analysis was performed on all taxonomic levels and centered log-ratio transformation was applied to the data. 

data <- read.delim(here("demo/demo_data/simulated_metaphlan.txt"), header = TRUE)

data <- 
  data %>%
  rownames_to_column(var = "sample") %>%
  relocate(sample) %>%
  as_tibble() %>%
  slice(1:n_infants) 

## create feature matrix from data
feature_matrix <-
  as_tibble(colnames(data)[-1]) %>% 
  rename(taxon = value) %>%
  mutate(taxonID = NA) %>%
  mutate(taxonID = purrr::map_chr(taxonID, ~ stri_rand_strings(1, length = 7, pattern = "[A-Za-z0-9]"))) %>% #assign unique ID to each taxon 
  relocate(taxonID)

## check unique codes to make sure there are no duplicates
potential_duplicates <- feature_matrix %>%
  filter(duplicated(taxonID) | duplicated(taxonID, fromLast = TRUE))

if (any(duplicated(potential_duplicates$taxonID))) {
  stop("Duplicated taxonID values found. Script stopped.")
}

## create input data for MOFA/MEFISTO
data_bacteria <- 
  data %>%
  pivot_longer(-sample) %>%
  rename(taxon = name) %>% #nrow = # nrow = 30195680
  full_join(., feature_matrix, by = "taxon") %>% # nrow = 30195680
  full_join(.,meta_data, by = "sample") %>% # nrow = 30195680
  rename(sample = sample,
         feature = taxonID,
         value = value) %>% 
  mutate(time = if_else(Timepoint == "W2", 0.5, as.numeric(str_extract(Timepoint, "\\d+")))) %>% 
  mutate(view = "microbiome",
         feature = as.character(feature)) %>%
  relocate(view, sample, time, feature, value)

# test:  do all samples contain the same number of features?
tib <-
  data_bacteria %>%
  group_by(sample) %>%
  reframe(n_features_per_sample = n()) 

are_all_values_the_same <- all(length(unique(tib$n_features_per_sample)) == 1)

if (are_all_values_the_same == FALSE) {
  stop("ERROR: The number of features is not equal across all samples.")
}

# put all data views together 
data_for_MOFA <- data_bacteria

# Part (ii) ----
# create MOFA/MEFISTO object
untrained_object <- create_mofa(data = data_for_MOFA)

# set covariates
untrained_object <- set_covariates(untrained_object, covariates = "time")

# set data options
data_opts <- get_default_data_options(untrained_object)

#set model options
model_opts <- get_default_model_options(untrained_object)
model_opts$num_factors <- 10

#set MEFISTO optioms
mefisto_opts <- get_default_mefisto_options(untrained_object)
mefisto_opts$n_grid <- 10
mefisto_opts$start_opt <- 50
mefisto_opts$opt_freq <- 50
mefisto_opts$model_groups <- FALSE # fast option (no optimization of group relationships) 

# set training options
train_opts <- get_default_training_options(untrained_object)
train_opts$seed <- 2020
train_opts$convergence_mode <- "slow"
train_opts$maxiter <- 1000
train_opts$verbose <- TRUE

# pass options to object
MOFAobject_untrained <- prepare_mofa(
  object = untrained_object,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts,
  mefisto_options = mefisto_opts
) 


# Part (iii): train the model

# Uncommenting and running the commands below will train a MEFISTO model on the demo data.
# This will take some time. A pretrained model is provided in the demo folder, should one not want to train the model and
# proceed directly to the downstream analysis.


# outfile <- paste0(here("LLNEXT_microbiome_model_DEMO.hdf5"))

# run_mofa(MOFAobject_untrained, outfile = outfile, use_basilisk = TRUE) # this will become the MOFAobject 


#########################
###DOWNSTREAM ANALYSIS###
#########################

# clear environment
rm(list = ls())

# load libraries
library(kableExtra)
library(lme4)
library(lmerTest)
library(broom)
library(broom.mixed)
library(tidyverse)
library(MOFA2)
library(here)
library(openxlsx)


MOFAobject <- load_model(here("LLNEXT_microbiome_model_DEMO.hdf5"), load_interpol_Z = FALSE)

# Explore model output----

## plot variance explained 
res_variance <- calculate_variance_explained(MOFAobject)

## how many factors do I want to plot?
start_factor <- 1
end_factor <- 10

## create a vector of factor names
factor_names <- paste0("factor_", seq(start_factor, end_factor))

## plot as bar charts
p_variance_bar_charts <-
  as_tibble(res_variance$r2_per_factor$single_group) %>% 
  mutate(factor = all_of(factor_names)) %>%
  mutate(factor = factor(factor, levels = factor_names )) %>%
  ggplot() +
  geom_bar(aes(x = factor, y = microbiome), stat = "identity") +
  ylab("variance explained (%)") +
  xlab("factor") +
  ggtitle("Variance Explained by Factor") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

p_variance_bar_charts

## examine variance total
p <- plot_variance_explained(MOFAobject, plot_total = T)

## variance decomposition by factor
p_variance_by_factor <-
  p[[1]] + theme(axis.text.x = element_text(angle = 90))

## total variance explained per view 
p_variance_total <-
  p[[2]] + theme(axis.text.x = element_text(angle = 90))

## plot projections
p_dimensions_reduced <-
  plot_factors(MOFAobject, color_by = "time")

## show plots
p
p_variance_by_factor
p_variance_total
p_dimensions_reduced

data_to_plot <-
  data_to_test %>%
  pivot_wider(names_from = factor, values_from = value)

# Plot Factor 1 Trajectory ----
#colors
col_vec2 <- c("#FFC20A", "#0C7BDC","#004D40")

fac1_trajectory_sp <-
  data_to_plot %>%
  ggplot(aes(x = Timepoint, y = Factor1, group = Factor_phenotype_static, color = Factor_phenotype_static)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = col_vec2) +
  stat_summary(
    fun = mean,   # calc mean
    geom = "line",  
    aes(group = Factor_phenotype_static),  # draw lines for each group
    size = 1  
  ) +
  theme_classic() +
  xlab("time (month)") +
  ylab("factor 1 score") +
  labs(color = "Factor_phenotype_static")

fac1_trajectory_sp


# Associate phenotypes to factors ----
## N.B. simulated_metadata.txt contains both meta data and phenotypic data

meta_data <- read.delim(here("demo/demo_data/simulated_metadata.txt"), header = TRUE)

## select infant data
meta_data <-  
  meta_data %>%
  rownames_to_column(var = "sample") %>%
  relocate(sample) %>%
  as_tibble() %>%
  filter(Type == "Infant")

data_pheno <-
  meta_data %>% 
  select(sample, Family, Timepoint,Binary_phenotype_static, Factor_phenotype_static,Quant_phenotype_static, Covariate1, Covariate2)

phenotypes_to_test <- c("Binary_phenotype_static", "Factor_phenotype_static", "Quant_phenotype_static")

data_factors <-
  get_factors(MOFAobject, as.data.frame = TRUE, factors = c(1:10)) 

data_to_test <-
  data_pheno %>% 
  left_join(., data_factors, by = "sample") 

## function to run mixed models with stepwise correction 
run_and_correct_MM <- function(df, phenotypes, factor_name) {
  initial_results <- list()
  corrected_results <- list()
  
  # run initial mixed models
  for (i in seq_along(phenotypes)) {
    phenotype <- phenotypes[i]
    
    # remove NAs from data and select needed data 
    df_filtered <- df %>%
      ungroup() %>%
      filter(factor == factor_name) %>%
      mutate(Timepoint = factor(Timepoint, levels = c("W2", "M1","M3","M6","M12"))) %>%
      filter(if_all(everything(), ~ !is.na(.)))
    
    # fit model with phenotype
    fit1 <- tryCatch({
      model <- lmerTest::lmer(as.formula(paste("value ~", "Covariate1 + Covariate2 +", phenotype, "* Timepoint + (1|Family)")), data = df_filtered, REML = FALSE)
      list(success=TRUE, model=model, summary=broom::glance(model), tidied=broom::tidy(model) %>% mutate(phenotype = phenotype, factor = factor_name))
    }, error=function(e) {
      list(success=FALSE, error_message=conditionMessage(e), model=NULL, summary=NULL, tidied=NULL)
    }, warning = function(w) {
      list(success=FALSE, warning_message=conditionMessage(w), model=NULL, summary=NULL, tidied=NULL)
    })
    
    # fit model without phenotype
    fit0 <- tryCatch({
      model <- lmerTest::lmer(as.formula(paste("value ~", "Covariate1 + Covariate2 + Timepoint + (1|Family)")), data = df_filtered, REML = FALSE)
      list(success=TRUE, model=model, summary=broom::glance(model), tidied=broom::tidy(model) %>% mutate(phenotype = phenotype, factor = factor_name))
    }, error=function(e) {
      list(success=FALSE, error_message=conditionMessage(e), model=NULL, summary=NULL, tidied=NULL)
    }, warning = function(w) {
      list(success=FALSE, warning_message=conditionMessage(w), model=NULL, summary=NULL, tidied=NULL)
    })
    
    # compare models if both were successfully fitted
    if (fit1$success && fit0$success) {
      model_comparison <- anova(fit1$model, fit0$model, test = "LRT")
      model_comparison_tidied <- broom::tidy(model_comparison) %>% mutate(phenotype = phenotype, factor = factor_name)
    } else {
      model_comparison <- NULL
      model_comparison_tidied <- NULL
    }
    
    # store all results 
    initial_results[[i]] <- list(fit1=fit1, fit0=fit0, model_comparison=model_comparison, model_comparison_tidied=model_comparison_tidied)
  }
  
  # summarize the initial model comparisons 
  summary_df <- data.frame()
  for (i in seq_along(initial_results)) {
    model_comparison_tidied <- initial_results[[i]]$model_comparison_tidied
    
    if (!is.null(model_comparison_tidied)) {
      summary_df <- bind_rows(summary_df, model_comparison_tidied) %>%  filter(!is.na(p.value)) %>%  mutate(p.adj = p.adjust(p.value, method = "bonferroni")) 
    }
  }
  
  # identify the phenotype with the lowest p-value
  lowest_p_phenotype <- summary_df %>% filter(!is.na(p.value)) %>% arrange(p.adj) %>% slice(1) %>% pull(phenotype)
  
  # corrected Associations: correct other phenotypes using the phenotype with the lowest p-value
  for (phenotype in phenotypes) {
    if (phenotype != lowest_p_phenotype) {
      # Remove NAs from data
      df_filtered <- df %>%
        ungroup() %>%
        filter(factor == factor_name) %>%
        mutate(Timepoint = factor(Timepoint, levels = c("W2", "M1","M3","M6","M12"))) %>%
        filter(if_all(everything(), ~ !is.na(.)))
      
      # fit model with phenotype
      fit1 <- tryCatch({
        model <- lmerTest::lmer(as.formula(paste("value ~", "Covariate1 + Covariate2 +", lowest_p_phenotype, "+", phenotype, "* Timepoint + (1|Family)")), data = df_filtered, REML = FALSE)
        list(success=TRUE, model=model, summary=broom::glance(model), tidied=broom::tidy(model) %>% mutate(corrected_phenotype = phenotype, factor = factor_name))
      }, error=function(e) {
        list(success=FALSE, error_message=conditionMessage(e), model=NULL, summary=NULL, tidied=NULL)
      }, warning = function(w) {
        list(success=FALSE, warning_message=conditionMessage(w), model=NULL, summary=NULL, tidied=NULL)
      })
      
      # fit model without phenotype
      fit0 <- tryCatch({
        model <- lmerTest::lmer(as.formula(paste("value ~", "Covariate1 + Covariate2 +", lowest_p_phenotype, "+ Timepoint + (1|Family)")), data = df_filtered, REML = FALSE)
        list(success=TRUE, model=model, summary=broom::glance(model), tidied=broom::tidy(model) %>% mutate(corrected_phenotype = phenotype, factor = factor_name))
      }, error=function(e) {
        list(success=FALSE, error_message=conditionMessage(e), model=NULL, summary=NULL, tidied=NULL)
      }, warning = function(w) {
        list(success=FALSE, warning_message=conditionMessage(w), model=NULL, summary=NULL, tidied=NULL)
      })
      
      # compare models if both were successfully fitted
      if (fit1$success && fit0$success) {
        model_comparison <- anova(fit1$model, fit0$model, test = "LRT")
        model_comparison_tidied <- broom::tidy(model_comparison) %>% mutate(corrected_phenotype = phenotype, factor = factor_name)
      } else {
        model_comparison <- NULL
        model_comparison_tidied <- NULL
      }
      
      # store the result
      corrected_results[[phenotype]] <- list(fit1=fit1, fit0=fit0, model_comparison=model_comparison, model_comparison_tidied=model_comparison_tidied)
    }
  }
  
  # return results
  return(list(initial_results=initial_results, corrected_results=corrected_results, lowest_p_phenotype=lowest_p_phenotype))
}

## helper function to summarize model comparisons
summarize_model_comparisons <- function(results) {
  summary_df <- data.frame()
  
  for (i in seq_along(results)) {
    model_comparison_tidied <- results[[i]]$model_comparison_tidied
    
    if (!is.null(model_comparison_tidied)) {
      summary_df <- bind_rows(summary_df, model_comparison_tidied) %>% filter(!is.na(p.value)) %>%  mutate(p.adj = p.adjust(p.value, method = "bonferroni")) %>%
        select(npar,df,p.value,matches("pheno*"),factor,p.adj) %>% arrange(p.adj)
    }
  }
  
  return(summary_df)
}

# Perform associations by factor----

## Factor 1 
Fac1_models <-run_and_correct_MM(df = data_to_test, phenotypes = phenotypes_to_test, factor_name = "Factor1")

## summarize results
Fac1_initial_summary <- summarize_model_comparisons(Fac1_models$initial_results) %>% mutate(source = "initial_model")
Fac1_corrected_summary <- summarize_model_comparisons(Fac1_models$corrected_results) %>%   mutate(source = "corrected_model")
Fac1_combined_summary <- bind_rows(Fac1_initial_summary, Fac1_corrected_summary) %>%  mutate_all(~ ifelse(is.na(.), "", .)) %>% select(factor, source, npar, df, phenotype, corrected_phenotype, p.value, p.adj)

## table of results
kable(Fac1_combined_summary, format = "html", table.attr = "class='table table-bordered'")

## Factor 2 
Fac2_models <-run_and_correct_MM(df = data_to_test, phenotypes = phenotypes_to_test, factor_name = "Factor2")

## summarize results
Fac2_initial_summary <- summarize_model_comparisons(Fac2_models$initial_results) %>% mutate(source = "initial_model")
Fac2_corrected_summary <- summarize_model_comparisons(Fac2_models$corrected_results) %>%   mutate(source = "corrected_model")
Fac2_combined_summary <- bind_rows(Fac2_initial_summary, Fac2_corrected_summary) %>%  mutate_all(~ ifelse(is.na(.), "", .)) %>% select(factor, source, npar, df, phenotype, corrected_phenotype, p.value, p.adj)

## table of results
kable(Fac2_combined_summary, format = "html", table.attr = "class='table table-bordered'")


## Factor 3 
Fac3_models <-run_and_correct_MM(df = data_to_test, phenotypes = phenotypes_to_test, factor_name = "Factor3")

## summarize results
Fac3_initial_summary <- summarize_model_comparisons(Fac3_models$initial_results) %>% mutate(source = "initial_model")
Fac3_corrected_summary <- summarize_model_comparisons(Fac3_models$corrected_results) %>%   mutate(source = "corrected_model")
Fac3_combined_summary <- bind_rows(Fac3_initial_summary, Fac3_corrected_summary) %>%  mutate_all(~ ifelse(is.na(.), "", .)) %>% select(factor, source, npar, df, phenotype, corrected_phenotype, p.value, p.adj)

## table of results
kable(Fac3_combined_summary, format = "html", table.attr = "class='table table-bordered'")

## Factor 4 
Fac4_models <-run_and_correct_MM(df = data_to_test, phenotypes = phenotypes_to_test, factor_name = "Factor4")

## summarize results
Fac4_initial_summary <- summarize_model_comparisons(Fac4_models$initial_results) %>% mutate(source = "initial_model")
Fac4_corrected_summary <- summarize_model_comparisons(Fac4_models$corrected_results) %>%   mutate(source = "corrected_model")
Fac4_combined_summary <- bind_rows(Fac4_initial_summary, Fac4_corrected_summary) %>%  mutate_all(~ ifelse(is.na(.), "", .)) %>% select(factor, source, npar, df, phenotype, corrected_phenotype, p.value, p.adj)

## table of results
kable(Fac4_combined_summary, format = "html", table.attr = "class='table table-bordered'")

## Factor 5 
Fac5_models <-run_and_correct_MM(df = data_to_test, phenotypes = phenotypes_to_test, factor_name = "Factor5")

## summarize results
Fac5_initial_summary <- summarize_model_comparisons(Fac5_models$initial_results) %>% mutate(source = "initial_model")
Fac5_corrected_summary <- summarize_model_comparisons(Fac5_models$corrected_results) %>%   mutate(source = "corrected_model")
Fac5_combined_summary <- bind_rows(Fac5_initial_summary, Fac5_corrected_summary) %>%  mutate_all(~ ifelse(is.na(.), "", .)) %>% select(factor, source, npar, df, phenotype, corrected_phenotype, p.value, p.adj)

## table of results
kable(Fac5_combined_summary, format = "html", table.attr = "class='table table-bordered'")

## Factor 6
Fac6_models <-run_and_correct_MM(df = data_to_test, phenotypes = phenotypes_to_test, factor_name = "Factor6")

## summarize results
Fac6_initial_summary <- summarize_model_comparisons(Fac6_models$initial_results) %>% mutate(source = "initial_model")
Fac6_corrected_summary <- summarize_model_comparisons(Fac6_models$corrected_results) %>%   mutate(source = "corrected_model")
Fac6_combined_summary <- bind_rows(Fac6_initial_summary, Fac6_corrected_summary) %>%  mutate_all(~ ifelse(is.na(.), "", .)) %>% select(factor, source, npar, df, phenotype, corrected_phenotype, p.value, p.adj)

## table of results
kable(Fac6_combined_summary, format = "html", table.attr = "class='table table-bordered'")

## Factor 7 
Fac7_models <-run_and_correct_MM(df = data_to_test, phenotypes = phenotypes_to_test, factor_name = "Factor7")

## summarize results
Fac7_initial_summary <- summarize_model_comparisons(Fac7_models$initial_results) %>% mutate(source = "initial_model")
Fac7_corrected_summary <- summarize_model_comparisons(Fac7_models$corrected_results) %>%   mutate(source = "corrected_model")
Fac7_combined_summary <- bind_rows(Fac7_initial_summary, Fac7_corrected_summary) %>%  mutate_all(~ ifelse(is.na(.), "", .)) %>% select(factor, source, npar, df, phenotype, corrected_phenotype, p.value, p.adj)

## table of results
kable(Fac7_combined_summary, format = "html", table.attr = "class='table table-bordered'")

## Factor 8 
Fac8_models <-run_and_correct_MM(df = data_to_test, phenotypes = phenotypes_to_test, factor_name = "Factor8")

## summarize results
Fac8_initial_summary <- summarize_model_comparisons(Fac8_models$initial_results) %>% mutate(source = "initial_model")
Fac8_corrected_summary <- summarize_model_comparisons(Fac8_models$corrected_results) %>%   mutate(source = "corrected_model")
Fac8_combined_summary <- bind_rows(Fac8_initial_summary, Fac8_corrected_summary) %>%  mutate_all(~ ifelse(is.na(.), "", .)) %>% select(factor, source, npar, df, phenotype, corrected_phenotype, p.value, p.adj)

## table of results
kable(Fac8_combined_summary, format = "html", table.attr = "class='table table-bordered'")

## Factor 9 
Fac9_models <-run_and_correct_MM(df = data_to_test, phenotypes = phenotypes_to_test, factor_name = "Factor9")

## summarize results
Fac9_initial_summary <- summarize_model_comparisons(Fac9_models$initial_results) %>% mutate(source = "initial_model")
Fac9_corrected_summary <- summarize_model_comparisons(Fac9_models$corrected_results) %>%   mutate(source = "corrected_model")
Fac9_combined_summary <- bind_rows(Fac9_initial_summary, Fac9_corrected_summary) %>%  mutate_all(~ ifelse(is.na(.), "", .)) %>% select(factor, source, npar, df, phenotype, corrected_phenotype, p.value, p.adj)

## table of results
kable(Fac9_combined_summary, format = "html", table.attr = "class='table table-bordered'")

## Factor 10
Fac10_models <-run_and_correct_MM(df = data_to_test, phenotypes = phenotypes_to_test, factor_name = "Factor10")

## summarize results
Fac10_initial_summary <- summarize_model_comparisons(Fac10_models$initial_results) %>% mutate(source = "initial_model")
Fac10_corrected_summary <- summarize_model_comparisons(Fac10_models$corrected_results) %>%   mutate(source = "corrected_model")
Fac10_combined_summary <- bind_rows(Fac10_initial_summary, Fac10_corrected_summary) %>%  mutate_all(~ ifelse(is.na(.), "", .)) %>% select(factor, source, npar, df, phenotype, corrected_phenotype, p.value, p.adj)

## table of results
kable(Fac10_combined_summary, format = "html", table.attr = "class='table table-bordered'")


###########
# put it all together
Fac_summaries_vec <- paste0("Fac", 1:10, "_combined_summary")

# Access each variable dynamically and store in a list
summaries_list <- lapply(Fac_summaries_vec, get)

# final table of results
final_table <-
  bind_rows(summaries_list) %>% 
  group_by(factor, source) %>%
  filter(p.adj < 0.05)

# create latex table
kable(final_table, "latex", booktabs = TRUE) 

###########################





