.libPaths("/home1/p309176/R_libraries")
library(readr)
library(vegan)
library(dplyr)
library(stringr)
library(reshape2)
library(tidyr)

meta_all_with_qc_curated <- as.data.frame(read_tsv('../../metadata_with_qc_NCPv2.tsv'))
negative_controls_all <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$Type == "Neg_ctrl"]

viral_clusters_initial <- as.data.frame(read_tsv('/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/initial_dereplication/NCP_viral_clusters.tsv',
                                                 col_names = FALSE))

viral_clusters_98 <- as.data.frame(read_tsv('/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/initial_dereplication/NCP_viral_clusters_98.tsv',
                                                 col_names = FALSE))

viral_clusters_99 <- as.data.frame(read_tsv('/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/initial_dereplication/NCP_viral_clusters_99.tsv',
                                            col_names = FALSE))

contains_negative_control <- function(text, negative_controls_all) {
  any(sapply(negative_controls_all, grepl, text))
}


nc_viral_clusters_98 <- viral_clusters_98 %>%
  filter(sapply(X2, contains_negative_control, negative_controls_all))

nc_viral_clusters_99 <- viral_clusters_99 %>%
  filter(sapply(X2, contains_negative_control, negative_controls_all))

all_values_98 <- unlist(strsplit(as.character(nc_viral_clusters_98$X2), split = ","))
all_values_99 <- unlist(strsplit(as.character(nc_viral_clusters_99$X2), split = ","))

writeLines(all_values_98, "viruses_remove_98.txt")
writeLines(all_values_99, "viruses_remove_99.txt")
