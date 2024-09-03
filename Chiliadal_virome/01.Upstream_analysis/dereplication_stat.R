##########################################
# UNDER DEVELOPMENT
# Parsing VC data and VC size after
# dereplication
##########################################
args = commandArgs(trailingOnly=TRUE)
##############################
# Loading libraries
##############################
library(tidyverse)
##############################
# Functions
##############################

##############################
# Input data
##############################
viral_clusters <- read_delim(args[1], col_names = F)
##############################
# ANALYSIS
##############################

# Transforming into long format to merge with Extended_TOF
long_df <- viral_clusters %>%
  separate_rows(X2, sep = ",\\s*") %>%
  arrange(X2, X1) %>%
  rename(Representative = X1, Cluster_member = X2)

# Getting the size of dereplication cluster
cluster_counts <- long_df %>%
  count(Representative) %>%
  rename(Cluster_size = n)

write.table(long_df, paste0(sub('\\.tsv', '', args[1]), '_long_format.txt'), sep='\t', row.names = F, quote=F)
write.table(cluster_counts, paste0(sub('\\.tsv', '', args[1]), '_size.txt'), sep='\t', row.names = F, quote=F)
