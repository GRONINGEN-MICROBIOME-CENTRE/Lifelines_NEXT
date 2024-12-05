#load libraries
library(here)
library(tidyverse)
library(MOFA2)
library(reshape2)
library(cowplot)
library(magrittr)
library(stringi)

# print run info
print("Only bacteria. No Groups. Factor variance threshold cutoff is 0.025. Using default number of factors.")

#create run directory----

factors<-"DEFAULTFAC"
views<- "BAC"
groups<- "NOGRP"
current_date <- format(Sys.Date(), format = "%d_%m_%Y")
current_time <- format(Sys.time(), "%H:%M:%S")

directory_path <- paste0(getwd(), "/run_", 
                         factors, "_", 
                         views, "_",
                         groups,"_",
                         current_time, "_",
                         current_date)

dir.create(directory_path)

# load and merge meta data and cross-sectional phenotype data ----
source(here("scripts/add_cross_sectional_phenotypes_to_metadata.R"))

# load bacterial sequencing data and make feature matrix ----
source(here("scripts/import_bacterial_relative_abundance_data.R"))

# create data that will be fed to MOFA ----
source(here("scripts/create_data_for_mofa.R"))

# create and train MOFA object ----
source(here("scripts/create_and_train_mofa_object.R"))

