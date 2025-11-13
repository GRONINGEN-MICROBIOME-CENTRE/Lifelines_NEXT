library(tidyverse)

# Get the KO annotattions first via wget: version 115.1 (released on 01.08.2025, downloaded 13.08.2025)
# download.file("https://rest.kegg.jp/link/pathway/ko", destfile = "../data/AMGP_v1/AMGP_EMP_PRG/KEGG_pathway_ko.txt", method = "curl")
# download.file("https://rest.kegg.jp/list/pathway/", destfile = "../data/AMGP_v1/AMGP_EMP_PRG/KEGG_pathway.txt", method = "curl")
# download.file("https://rest.kegg.jp/link/ko/module", destfile = "../data/AMGP_v1/AMGP_EMP_PRG/KEGG_module_ko.txt", method = "curl")
# download.file("https://rest.kegg.jp/list/module", destfile = "../data/AMGP_v1/AMGP_EMP_PRG/KEGG_module.txt", method = "curl")
# The manual's in here: https://www.kegg.jp/kegg/rest/keggapi.html

# Load the data
KEGG <- read.delim("../data/AMGP_v2/ko-annotations.tsv", header=T)
KEGG_pathway_ko <- read.delim("../data/AMGP_v1/AMGP_EMP_PRG/KEGG_pathway_ko.txt", header=F)
KEGG_pathway <- read.delim("../data/AMGP_v1/AMGP_EMP_PRG/KEGG_pathway.txt", header=F)

PRG <- read.delim("../data/AMGP_v2/HMMER_processed_table.tsv", header=T)
PRG_annot <- read.delim("../data/AMGP_v1/AMGP_EMP_PRG/phrog_annot_v4.tsv", header=T)

# FILTERING AND PREPROCESSING KEGG ANNOTATION - PART 1
KEGG <- KEGG[KEGG$thrshld != "-------" & KEGG$thrshld != "",]
KEGG <- KEGG %>%
  mutate(score = as.numeric(score),
         thrshld = as.numeric(thrshld)) %>%
  filter(score >= thrshld) %>%
  group_by(gene.name) %>%
  slice_min(order_by = E.value, with_ties = TRUE) %>%
  slice_max(order_by = score, with_ties = TRUE) %>%  # Check this code if running on other dfs; it works here because both annotations for duplicate proteins are OK, so it coild be a random choice as well
  ungroup() %>%
  select(-X.) %>%
  rename(gene_ID = gene.name,
         KO_ID = KO,
         KO_function = KO.definition)

# FILTERING AND PREPROCESSING PHROG ANNOTATION
PRG_annot$phrog <- paste0("phrog_", PRG_annot$phrog)
PRG_annot$color <- NULL
names(PRG)[names(PRG) == 'Match'] <- 'phrog'
PRG <- merge(PRG, PRG_annot, by='phrog', all.x=T)
PRG <- PRG[PRG$category != "unknown function", ]

# For unified filtering
set.seed(13)

PRG <- PRG %>%
  group_by(Protein_ID) %>%
  filter(E_value < 1e-5) %>%
  slice_max(order_by = Score, with_ties = TRUE) %>%
  slice_min(order_by = E_value, with_ties = TRUE) %>%
  slice_sample(n = 1) %>%  # there were 10 duplicates left with the same category (only 1 had head/tail distinct assigned categories)
  ungroup()

# FILTERING AND PREPROCESSING KEGG ANNOTATION - PART 2
KEGG <- KEGG[!is.na(KEGG$KO_ID), ]

colnames(KEGG_pathway_ko) <- c("ko", "path")
colnames(KEGG_pathway) <- c("path", "name")
KEGG_pathway_ko <- KEGG_pathway_ko %>%
  mutate(ko = str_replace(ko, "ko:", ""),
         path = str_replace(path, "path:", ""))

# Selecting only the metabolic pathways: checked manually
# The order correspond to the order in here: https://www.kegg.jp/kegg/pathway.html
KEGG_pathway_metabolism <- KEGG_pathway[1:192, ]
KOs_metabolism <- unique(KEGG_pathway_ko$ko[KEGG_pathway_ko$path %in% KEGG_pathway_metabolism$path])

# Getting Phrog and KEGG assignments ready for merging
KEGG_processed <- KEGG %>%
  mutate(KO_MG = ifelse(KO_ID %in% KOs_metabolism, "yes", "no"),
         KO_VIR = "no") %>%
  rename(Protein_ID = gene_ID)

keywords_vir <- c('portal', 'terminase', 'spike', 'capsid', 'sheath', 'tail', 'virion', 'holin', 'base plate', 'baseplate', 'lysozyme', 'head', 'structural', 'phage', 'vir')
for (i in keywords_vir){
  KEGG_processed$KO_VIR[grep(i, KEGG_processed$KO_function, ignore.case = TRUE)] <- 'yes'
}
unique(KEGG_processed$KO_function[KEGG_processed$KO_VIR == 'yes']) 

KEGG_processed$KO_VIR[KEGG_processed$KO_ID %in% c("K07271", "K18234", "K02168", "K03810")] <- 'no'

PRG_processed <- PRG %>%
  mutate(PHROG_VIR = ifelse(category %in% c("head and packaging", "tail", "lysis", "connector", "integration and excision"), "yes", "no"),
         E_value = NULL,
         Score = NULL,
         Description = NULL,
  ) %>%
  rename(PHROG_ID = phrog,
         PHROG_function = annot,
         PHROG_category = category)

PRG_KEGG <- merge(PRG_processed, KEGG_processed, by="Protein_ID", all = T)

PRG_KEGG <- PRG_KEGG %>%
  mutate(VIR = case_when(
    !is.na(PHROG_VIR) &  is.na(KO_VIR) ~ PHROG_VIR,
    is.na(PHROG_VIR) & !is.na(KO_VIR) ~ KO_VIR,
    PHROG_VIR == KO_VIR ~ PHROG_VIR,
    PHROG_VIR == "yes" & KO_VIR == "no" ~ PHROG_VIR,
    PHROG_VIR == "no" & KO_VIR == "yes" ~ KO_VIR,  # mostly nucleic acid metabolism here, that was assigned as viral with the keyword search
    TRUE ~ NA_character_
  ),
  MG = ifelse(KO_MG == "yes" & !(is.na(KO_MG)), "yes", "no"))

write.table(PRG_KEGG, "../results/AMGs_v2.tsv", sep='\t', row.names=F, col.names=T, quote=F)  # final table output