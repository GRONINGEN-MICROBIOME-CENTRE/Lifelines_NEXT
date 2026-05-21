#!/usr/bin/env Rscript
library(MOFA2)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
K <- as.integer(args[1])

data_for_MEFISTO <- readRDS("data_for_MEFISTO.rds")

subdata <- data_for_MEFISTO

obj <- create_mofa(subdata)
obj <- set_covariates(obj, covariates="time")

data_opts <- get_default_data_options(obj)
model_opts <- get_default_model_options(obj)
model_opts$num_factors <- K

mefisto_opts <- get_default_mefisto_options(obj)
mefisto_opts$n_grid <- 10
mefisto_opts$start_opt <- 50
mefisto_opts$opt_freq <- 50

train_opts <- get_default_training_options(obj)
train_opts$seed <- 2020 + K
train_opts$convergence_mode <- "slow"
train_opts$maxiter <- 1000
train_opts$verbose <- FALSE

untrained <- prepare_mofa(
  object=obj,
  data_options=data_opts,
  model_options=model_opts,
  training_options=train_opts,
  mefisto_options=mefisto_opts
)

outfile <- sprintf("MEFISTO_K%d_FULL.hdf5", K)

trained <- run_mofa(untrained, outfile = outfile, use_basilisk = TRUE)


