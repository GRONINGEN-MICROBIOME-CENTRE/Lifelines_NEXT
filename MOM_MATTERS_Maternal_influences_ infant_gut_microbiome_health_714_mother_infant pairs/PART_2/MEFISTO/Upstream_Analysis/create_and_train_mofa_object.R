# Purpose: Set MEFISTO and model option. Launch training.

# opening line ----
print("creating MOFA object...")

# create untrained MOFA object ----
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

# make file name
file_name <- paste0(here("MOFAobject_untrained_"), current_date, ".rda")


save(MOFAobject_untrained, file = file_name)

# train the model
outfile <- paste0(here("LLNEXT_microbiome_model_"), current_date, ".hdf5")

run_mofa(MOFAobject_untrained, outfile = outfile, use_basilisk = TRUE) #This will become the MOFAobject

# closing line ----
print("creating MOFA object COMPLETED.")
