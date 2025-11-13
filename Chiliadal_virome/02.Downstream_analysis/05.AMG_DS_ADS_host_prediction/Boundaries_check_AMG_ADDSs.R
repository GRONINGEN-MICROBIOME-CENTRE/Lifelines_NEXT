library(tidyverse)

## Load the data
DF_PDC_DBA_ini <- read.delim("../results/Defense_antidefense_syss_combined_v2.tsv", header=T)  # out file of ADDSs_prediction_DF_PDC_DBA.R
PRG_KEGG_ini <- read.delim("../results/AMGs_v2.tsv", header=T)  # out file of AMG_prediction_PRG_KEGG.R

linking <- read.delim("../data/AMGP_v2/linking_proteins.txt", header=F)

KEGG_pathway_ko <- read.delim("../data/AMGP_v1/AMGP_EMP_PRG/KEGG_pathway_ko.txt", header=F)
KEGG_pathway <- read.delim("../data/AMGP_v1/AMGP_EMP_PRG/KEGG_pathway.txt", header=F)

# Preprocess the linking file
linking <- linking %>%
  rename(Protein_ID = V1,
         New_protein_ID = V2) %>%
  mutate(
    Protein_ID = sub(" .*", "", Protein_ID),
    virus = sub("_CDS.*$", "", Protein_ID),
    protein_start = str_extract(Protein_ID, "\\d+(?=\\.\\.)"),
    protein_end = str_extract(Protein_ID, "(?<=\\.\\.)\\d+"),
    protein_start = as.integer(protein_start),
    protein_end = as.integer(protein_end)
  )

# Preprocess the data
antidefense_proteins_adds <- unlist(strsplit(DF_PDC_DBA_ini$protein_in_syst[DF_PDC_DBA_ini$activity == "Antidefense"], split = ","))  # retrieve all proteins assigned as antidefense system 
defense_proteins_adds <- unlist(strsplit(DF_PDC_DBA_ini$protein_in_syst[DF_PDC_DBA_ini$activity == "Defense"], split = ","))  # retrieve all proteins assigned as defense system 

PRG_KEGG <- merge(PRG_KEGG_ini, linking, by="Protein_ID", all=T)

PRG_KEGG$VIR[is.na(PRG_KEGG$VIR)] <- "no"
PRG_KEGG$VIR[PRG_KEGG$New_protein_ID %in% antidefense_proteins_adds] <- "yes"
PRG_KEGG$VIR[PRG_KEGG$New_protein_ID %in% defense_proteins_adds] <- "no"
PRG_KEGG$MG[is.na(PRG_KEGG$MG)] <- "no"

DF_PDC_DBA_long <- DF_PDC_DBA_ini %>%
  separate_rows(protein_in_syst, sep = ",") %>%
  rename(New_protein_ID = protein_in_syst)

DF_PDC_DBA_long <- merge(DF_PDC_DBA_long, linking[c("New_protein_ID", "Protein_ID")], by="New_protein_ID", all.x=T)
DF_PDC_DBA_long <- DF_PDC_DBA_long[c("New_protein_ID", "Protein_ID", "sys_id", "activity", "type", "subtype", "genes_count")]

DF_PDC_DBA_PRG_KEGG <- merge(PRG_KEGG, DF_PDC_DBA_long, all=T, by="Protein_ID")

## Preparing for boundary check
DF_PDC_DBA_PRG_KEGG <- DF_PDC_DBA_PRG_KEGG %>%
  mutate(VIR = ifelse(is.na(VIR), "no", VIR),
         MG = ifelse(is.na(MG), "no", MG),
         ADS = ifelse(activity == "Antidefense" & !is.na(activity), "yes", "no"),
         DS = ifelse(activity == "Defense" & !is.na(activity), "yes", "no"))%>%
  group_by(virus) %>%
  filter(
    sum(VIR == "yes", na.rm = TRUE) >= 2,  # at least 2 viral proteins
  ) %>%
  ungroup()

## Boundary check
DF_PDC_DBA_PRG_KEGG <- DF_PDC_DBA_PRG_KEGG %>%
  group_by(virus) %>%
  mutate(
    locus_start = min(protein_start[VIR == "yes"], na.rm = TRUE),
    locus_end = max(protein_start[VIR == "yes"], na.rm = TRUE),
    location = ifelse(
      protein_start >= locus_start & protein_start <= locus_end,
      "inside",
      "outside"
    )
  ) %>%
  select(-locus_start, -locus_end) %>% 
  ungroup() 

inside_systems <- unique(DF_PDC_DBA_PRG_KEGG$sys_id[!is.na(DF_PDC_DBA_PRG_KEGG$sys_id) & DF_PDC_DBA_PRG_KEGG$location == "inside"])
inside_amgs <- unique(DF_PDC_DBA_PRG_KEGG$Protein_ID[DF_PDC_DBA_PRG_KEGG$MG == "yes" & DF_PDC_DBA_PRG_KEGG$location == "inside"])
DF_PDC_DBA_fin <- DF_PDC_DBA_ini %>%
  mutate(location = ifelse(sys_id %in% inside_systems, "inside", "outside"))

DF_PDC_DBA_fin <- DF_PDC_DBA_fin[DF_PDC_DBA_fin$location == "inside", ]
DF_PDC_DBA_fin$location <- NULL


## MG refinement
MG_simple_filtered <- DF_PDC_DBA_PRG_KEGG[DF_PDC_DBA_PRG_KEGG$MG == "yes" & DF_PDC_DBA_PRG_KEGG$location == "inside", ]
MG_simple_filtered <- MG_simple_filtered[c("Protein_ID", "PHROG_ID", "PHROG_function", "PHROG_category", "KO_ID", "KO_function", "VIR",
  "virus",  "protein_start", "protein_end", "ADS", "DS")]

colnames(KEGG_pathway_ko) <- c("ko", "path")
colnames(KEGG_pathway) <- c("path", "name")
KEGG_pathway_ko <- KEGG_pathway_ko %>%
  mutate(ko = str_replace(ko, "ko:", ""),
         path = str_replace(path, "path:", ""))

# Selecting only the metabolic pathways: checked manually
# The order correspond to the order in here: https://www.kegg.jp/kegg/pathway.html
KEGG_pathway_metabolism <- KEGG_pathway[1:192, ]
KEGG_pathway_metabolism$group <- NA
KEGG_pathway_metabolism$group[1:13] <- "Global and overview maps"
KEGG_pathway_metabolism$group[14:27] <- "Carbohydrate metabolism"
KEGG_pathway_metabolism$group[28:35] <- "Energy metabolism"
KEGG_pathway_metabolism$group[36:52] <- "Lipid metabolism"
KEGG_pathway_metabolism$group[53:54] <- "Nucleotide metabolism"
KEGG_pathway_metabolism$group[55:68] <- "Amino acid metabolism"
KEGG_pathway_metabolism$group[69:75] <- "Metabolism of other amino acids"
KEGG_pathway_metabolism$group[76:98] <- "Glycan biosynthesis and metabolism"
KEGG_pathway_metabolism$group[99:110] <- "Metabolism of cofactors and vitamins"
KEGG_pathway_metabolism$group[111:131] <- "Metabolism of terpenoids and polyketides"
KEGG_pathway_metabolism$group[132:162] <- "Biosynthesis of other secondary metabolites"
KEGG_pathway_metabolism$group[163:183] <- "Xenobiotics biodegradation and metabolism"
KEGG_pathway_metabolism$group[184:192] <- "Chemical structure transformation maps"

KEGG_pathway_KO_description <- merge(KEGG_pathway_ko, KEGG_pathway_metabolism, by="path")
colnames(KEGG_pathway_KO_description) <- c("PATH_ID", "KO_ID", "PATH_NAME", "PATH_GROUP")
KOs_metabolism <- unique(KEGG_pathway_ko$ko[KEGG_pathway_ko$path %in% KEGG_pathway_metabolism$path])

KOs_nucleotide_metabolism <- unique(KEGG_pathway_ko$ko[KEGG_pathway_ko$path %in% c("map01232", "map00230", "map00240")])

MG_simple_filtered <- MG_simple_filtered %>%
  mutate(MG_category = ifelse(KO_ID %in% KOs_nucleotide_metabolism, "VDBG", "ABG"),  # VDBG (virus direct benefit gene) ABG (ambiguous benefit gene)
         MG_category = ifelse(VIR == "yes", "VDBG", MG_category),
         MG_category = ifelse(ADS == "yes", "VDBG", MG_category))

KOs_VBG <- unique(MG_simple_filtered$KO_ID[MG_simple_filtered$MG_category == "VDBG"])
MG_simple_filtered <- MG_simple_filtered %>%
  mutate(MG_category = ifelse(KO_ID %in% KOs_VBG, "VDBG", MG_category))

MG_simple_filtered <- merge(MG_simple_filtered, KEGG_pathway_KO_description[KEGG_pathway_KO_description$PATH_GROUP != "Global and overview maps", ], by="KO_ID", all.x = T)

write.table(DF_PDC_DBA_fin, "../results/Defense_antidefense_syss_combined_after_boundaries_check_v2.tsv", sep='\t', row.names=F, col.names=T, quote=F)  # final table output
write.table(MG_simple_filtered, "../results/MGs_after_boundaries_check_v2.tsv", sep='\t', row.names=F, col.names=T, quote=F)

write.table(DF_PDC_DBA_PRG_KEGG, "../results/Defense_antidefense_syss_MG_combined_full_v2.tsv", sep='\t', row.names=F, col.names=T, quote=F)
