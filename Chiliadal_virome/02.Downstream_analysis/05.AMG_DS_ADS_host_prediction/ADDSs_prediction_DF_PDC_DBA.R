library(tidyverse)

# Merging outs from combined Padloc & DfenseFinder and  and dbAPIS
DBA_hmm_full <- read.delim("../data/DSP_v2/hmmscan.out.parsed.tsv", header=TRUE)
DBA_dmnd_full <- read.delim("../data/DSP_v2/diamond.out.parsed.tsv", header=TRUE)

DF_PDC <- read.delim("../data/DSP_v2/deduplicated_defense_systems.tsv", header=TRUE)
DF_ads <- read.delim("../data/DSP_v2/High_quality_vOTUs_PNT_ed_defense_finder_systems_antidefense.tsv", header=TRUE)

# For unified filtering
set.seed(13)

# Processing dbAPIS output
# Filter HMM out
DBA_hmm <- DBA_hmm_full %>%
  filter(Domain.c.evalue < 1e-5, Domain.score >= 50) %>%
  group_by(Query) %>%
  slice_max(order_by = Domain.score, with_ties = TRUE) %>%
  slice_min(order_by = Domain.c.evalue, with_ties = TRUE) %>%
  ungroup()

DBA_hmm <- DBA_hmm %>%
  group_by(Query, Hit.CLAN) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  mutate(tool_original = "DBA_hmm")

# Filter DIAMOND out
DBA_dmnd <- DBA_dmnd_full %>%
  select(-seqid) %>%
  filter(evalue < 1e-5, bitscore >= 50) %>%
  distinct() %>%
  group_by(qseqid) %>%
  slice_max(order_by = bitscore, with_ties = TRUE) %>%
  slice_min(order_by = evalue, with_ties = TRUE) %>%
  ungroup()

DBA_dmnd <- DBA_dmnd %>%
  group_by(qseqid, famid) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  mutate(tool_original = "DBA_dmnd")

# Merge the two
# Select by the e-value
colnames(DBA_dmnd) <- c("Query", "Hit.family", "Defense.type", "Hit.CLAN", "Hit.CLAN.defense.type", "pident", 
                        "align.length", "Domain.c.evalue", "Domain.score", "qcov", "scov", "tool_original")

DBA_dmnd <- DBA_dmnd[c("Query", "Hit.family", "Defense.type", "Hit.CLAN", "Hit.CLAN.defense.type",
                       "Domain.c.evalue", "Domain.score", "tool_original")]

DBA_hmm <- DBA_hmm[c("Query", "Hit.family", "Defense.type", "Hit.CLAN", "Hit.CLAN.defense.type",
                       "Domain.c.evalue", "Domain.score", "tool_original")]

DBA <- rbind(DBA_dmnd, DBA_hmm)

DBA <- DBA %>%
  group_by(Query) %>%
  slice_min(order_by = Domain.c.evalue, with_ties = TRUE) %>%
  ungroup()

DBA <- DBA %>%
  group_by(Query, Hit.CLAN) %>%
  slice_sample(n = 1) %>%
  ungroup()
# Check at this point - some duplicates might be left!

# Removing the duplicates manually after cross-check with the defensefinder out
DBA <- DBA[!(DBA$Query == "GPD_uvig_152846_N0_L38737_K0.0_E0_P0_F0_39" & DBA$Hit.family == "APIS029"), ]
DBA <- DBA[!(DBA$Query == "IMGVR_UViG_3300029549_000118_N0_L45016_K0.0_E0_P1_F1_20" & DBA$Hit.family == "APIS029"), ]
DBA <- DBA[!(DBA$Query == "NEXT_MA03F028153A3_N325_L43919_K33.4_E1_P0_F0_25" & DBA$Hit.family == "APIS029"), ]
DBA <- DBA[!(DBA$Query == "NEXT_V0249_N223_L36234_K6.4_E0_P0_F0_45" & DBA$Hit.family == "APIS029"), ]
DBA <- DBA[!(DBA$Query == "MGV-GENOME-0260989_N0_L38626_K0.0_E0_P0_F0_13" & DBA$Hit.family == "APIS052"), ]

# Getting dbAPIS output unified with the DefenseFinder/Padloc
DBA_preped <- DBA %>%
  filter(Hit.family != "APIS028") %>%  # This family was excluded according to the dbAPIS readme
  mutate(seqid = str_replace(Query, "_[^_]*$", ""),
         Domain.score = NULL,
         Hit.CLAN.defense.type = ifelse(Hit.CLAN.defense.type == "", "CRISPR-Cas evasion", Hit.CLAN.defense.type),
         activity = "Antidefense",
         subtype = "",
         tool_final = "DBA",
         Domain.c.evalue = NULL,
         tool_original = "DBA"
         ) %>%
  rename(type = Hit.CLAN.defense.type)

DBA_preped <- DBA_preped %>%
  mutate(type = case_when(
    type == "CRISPR-Cas evasion" ~ "Anti_CRISPR",
    type == "Thoeris" ~ "Anti_Thoeris",
    type == "NAD+ reconstitution pathway (NARP)" ~ "NADP",
    type == "Gabija" ~ "Anti_Gabija",
    type == "toxin-antitoxin (TA)" ~ "Anti_TA",
    type == "Retron" ~ "Anti_Retron",
    type == "CRISPR-Cas evasion by DNA repair" ~ "Anti_CRISPR",
    type == "restriction-modification (RM)" ~ "Anti_RM",
    type == "RecBCD" ~ "Anti_RecBCD",
    type == "cyclic oligonucleotide-based antiphage signaling system (CBASS)" ~ "Anti_CBASS",
    type == "Dnd" ~ "Anti_Dnd",
    type == "pyrimidine cyclase system for antiphage resistance (Pycsar)" ~ "Anti_Pycsar",
    type == "superinfection exclusion" ~ "Superinfection_exclusion",
    type == "SOS response" ~ "Anti_SOS_response",
    type == "restriction-modification (RM);bacteriophage exclusion (BREX)" ~ "Anti_RM",
    type == "Gasdermin" ~ "Anti_Gasdermin",
    type == "TIR-STING" ~ "Anti_TIR-STING",
    type == "O-antigen-based barrier" ~ "Anti_O-antigen-based_barrier",
    type == "broad-spectrum counter-defense" ~ "Anti_defense_broad",
    type == "CBASS, Pycsar and CRISPR–Cas (type III)" ~ "Anti_defense_broad",
    TRUE ~ type
  ))

DBA_preped <- DBA_preped %>%
  mutate(
    prot_num = as.integer(str_extract(Query, "(?<=_)\\d+$"))
  )

DBA_preped <- DBA_preped %>%
  group_by(seqid, type) %>%
  arrange(prot_num, .by_group = TRUE) %>%
  mutate(
    chain_break = prot_num - lag(prot_num, default = first(prot_num)) != 1,
    chain_id_local = cumsum(replace_na(chain_break, TRUE))
  ) %>%
  ungroup() %>%
  mutate(
    sys_id = paste0(seqid, "_", type, "_", chain_id_local)
  ) %>%
  group_by(sys_id) %>%
  mutate(
    sequential = if_else(n() > 1, "yes", "no"),
    protein_in_syst = paste(Query, collapse = ","),
    name_of_profiles_in_sys = paste(Hit.family, collapse = ", ")
  ) %>%
  ungroup()

DBA_preped <- DBA_preped %>%
  mutate(
    sys_beg = str_extract(protein_in_syst, "^[^,]+"),
    sys_end = str_extract(protein_in_syst, "[^,]+$"),
    genes_count = str_count(protein_in_syst, ",") + 1
  )

DBA_preped <- DBA_preped[colnames(DF_PDC)]
DBA_preped <- DBA_preped[!duplicated(DBA_preped), ]

# Processing DF_PDC
# Merging DF_PDC out with the antidefense systems that were found by defense finder
DF_ads <- DF_ads %>%
  mutate(seqid = str_replace(sys_beg, "_[^_]*$", ""),
         tool_original = "defensefinder",
         tool_final = "defensefinder")
DF_ads <- DF_ads[!duplicated(DF_ads), ]

DF_ads <- DF_ads[, colnames(DF_PDC)]
DF_PDC <- rbind(DF_PDC, DF_ads)

DF_PDC <- DF_PDC %>%
  filter(subtype != "VSPR",
         subtype != "DMS_other") %>%
  mutate(subtype = case_when(
    subtype == "cas_type_IV-B" ~ "CAS_Class1-Type-IV-B",
    subtype == "ietAS" ~ "Gao_Iet",
    subtype == "Lamassu-Hypothetical" ~ "Lamassu_Family",
    subtype == "shedu" ~ "Shedu",
    subtype == "hachiman_type_I" ~ "Hachiman_type_I",
    TRUE ~ subtype
    ),
    activity = ifelse(activity == "" & tool_final == "padloc", "Defense", activity),
    type = ifelse(name_of_profiles_in_sys == "acriia7", "Anti_CRISPR", type),
    type = ifelse(name_of_profiles_in_sys == "orf148", "Anti_defense_broad", type),
    type = ifelse(name_of_profiles_in_sys %in% c("gnarl1", "gnarl2"), "Anti_O-antigen-based_barrier", type)) %>%
  mutate(
    type = ifelse(subtype == "Uzume" & type == "", "Uzume", type),
    type = ifelse(subtype %in% c("RM_type_II", "RM_type_IV", "RM_type_IIG", "RM_type_HNH") & type == "", "RM", type),
    type = ifelse(subtype == "PD-Lambda-1" & type == "", "PD-Lambda-1", type),
    type = ifelse(subtype == "AbiE" & type == "", "AbiE", type),
    type = ifelse(subtype == "AbiL" & type == "", "AbiL", type),
    type = ifelse(subtype == "SEFIR" & type == "", "SEFIR", type),
    type = ifelse(subtype == "SoFic" & type == "", "SoFIC", type),
    type = ifelse(subtype == "HEC-06" & type == "", "HEC-06", type),
    type = ifelse(subtype == "Gao_Iet" & type == "", "Gao_Iet", type),
    type = ifelse(subtype == "AbiU" & type == "", "AbiU", type),
    type = ifelse(subtype == "Tiamat" & type == "", "Tiamat", type),
    type = ifelse(subtype == "retron_VII-A2" & type == "", "Retron", type),
    type = ifelse(subtype %in% c("DRT_other", "DRT_class_III") & type == "", "DRT", type),
    type = ifelse(subtype == "Menshen_other" & type == "", "Menshen", type),
    type = ifelse(subtype == "cas_type_IV-B" & type == "", "Cas", type),
    type = ifelse(subtype == "Paris" & type == "", "Paris", type),
    type = ifelse(subtype == "PD-T4-6" & type == "", "PD-T4-6", type),
    type = ifelse(subtype == "HEC-08" & type == "", "HEC-08", type),
    type = ifelse(subtype == "Shedu" & type == "", "Shedu", type),
    type = ifelse(subtype == "PifA" & type == "", "Pif", type),
    type = ifelse(subtype == "Mokosh_TypeII" & type == "", "Mokosh", type),
    type = ifelse(subtype == "PD-T7-2" & type == "Gao_Her", "PD-T7-2", type),
    type = ifelse(subtype == "CAS_Class1-Type-IV-B" & type == "", "Cas", type),
    type = ifelse(subtype == "HEC-04" & type == "", "HEC-04", type),
    type = ifelse(subtype == "septu_other" & type == "", "Septu", type),
    type = ifelse(subtype == "mza_other" & type == "", "Gao_Mza", type),
    type = ifelse(subtype == "gabija" & type == "", "Gabija", type),
    type = ifelse(subtype == "Hachiman_type_I" & type == "", "Hachiman", type),
    type = ifelse(subtype == "Lamassu_Family" & type == "", "Lamassu-Fam", type)
    )

# Merging & dereplicating the outputs
DF_PDC_DBA_to_derep <- rbind(DF_PDC, DBA_preped)

all_proteins_antidefense <- unlist(strsplit(DF_PDC_DBA_to_derep$protein_in_syst[DF_PDC_DBA_to_derep$activity == "Antidefense"], split = ","))
duplicated_proteins <- all_proteins_antidefense[duplicated(all_proteins_antidefense)]

DF_PDC_DBA_to_derep <- DF_PDC_DBA_to_derep %>%
  mutate(
    duplicates = map_lgl(protein_in_syst, ~ {
      prots <- str_split(.x, ",")[[1]]
      any(prots %in% duplicated_proteins)
    }),
    duplicates = if_else(activity == "Defense", FALSE, duplicates)
  )

DF_PDC_DBA <- DF_PDC_DBA_to_derep %>%
  filter(duplicates == FALSE)

DF_PDC_DBA_to_derep <- DF_PDC_DBA_to_derep %>%
  filter(duplicates == TRUE)

DF_PDC_DBA_to_derep <- DF_PDC_DBA_to_derep %>%
  group_by(protein_in_syst) %>%
  mutate(
    same_type = n() > 1 & n_distinct(type) == 1
  ) %>%
  filter(!same_type | tool_original == "defensefinder") %>%
  mutate(
    tool_final = if_else(same_type & tool_original == "defensefinder", "both", tool_final),
    duplicates = if_else(same_type & tool_original == "defensefinder", FALSE, duplicates)
  ) %>%
  ungroup() %>%
  select(-same_type)

DF_PDC_DBA <- rbind(DF_PDC_DBA, DF_PDC_DBA_to_derep[DF_PDC_DBA_to_derep$duplicates == FALSE, ])

DF_PDC_DBA_to_derep <- DF_PDC_DBA_to_derep %>%
  filter(duplicates == TRUE)

DF_PDC_DBA_to_derep <- DF_PDC_DBA_to_derep %>%
  group_by(seqid) %>%
  mutate(
    same_type = n() == 2 & n_distinct(type) == 1
  ) %>%
  filter(!same_type | genes_count == max(genes_count)) %>%
  mutate(duplicates = if_else(same_type, FALSE, duplicates)) %>%
  ungroup() %>%
  select(-same_type)
  
DF_PDC_DBA <- rbind(DF_PDC_DBA, DF_PDC_DBA_to_derep[DF_PDC_DBA_to_derep$duplicates == FALSE, ])

DF_PDC_DBA_to_derep <- DF_PDC_DBA_to_derep %>%
  filter(duplicates == TRUE)


# This one should be checked manually prior running!
# NADP system sometimes has gaps between the genes, which is not hard-coded for DBA when combining into systems
# Thus, where defensefinder finds it as a two-gene system, DBA will give 2 different outs
DF_PDC_DBA_to_derep <- DF_PDC_DBA_to_derep %>%
  group_by(seqid) %>%
  mutate(
    same_type = n() == 3 & type == "NADP"
  ) %>%
  filter(!same_type | tool_original == "defensefinder") %>%
  mutate(tool_final = if_else(same_type & tool_original == "defensefinder", "both", tool_final),
         duplicates = if_else(same_type & tool_original == "defensefinder", FALSE, duplicates)) %>%
  ungroup() %>%
  select(-same_type)

DF_PDC_DBA <- rbind(DF_PDC_DBA, DF_PDC_DBA_to_derep[DF_PDC_DBA_to_derep$duplicates == FALSE, ])

DF_PDC_DBA_to_derep <- DF_PDC_DBA_to_derep %>%
  filter(duplicates == TRUE)

# Manual checking discrepancies between the anti_crispr and anti_thoeris assignments
# Only leaving DBA assignments where discrepancies happen, and there's only one protein per system

proteins_df <- DF_PDC_DBA_to_derep$protein_in_syst[DF_PDC_DBA_to_derep$subtype == "tad2_acriia7"]
DF_PDC_DBA_to_derep <- DF_PDC_DBA_to_derep[!(DF_PDC_DBA_to_derep$subtype == "tad2_acriia7"),]
DF_PDC_DBA_to_derep$duplicates[DF_PDC_DBA_to_derep$protein_in_syst %in% proteins_df] <- FALSE

DF_PDC_DBA <- rbind(DF_PDC_DBA, DF_PDC_DBA_to_derep[DF_PDC_DBA_to_derep$duplicates == FALSE, ])

DF_PDC_DBA_to_derep <- DF_PDC_DBA_to_derep %>%
  filter(duplicates == TRUE)

# Manual decision-making for the duplicates that left
DF_PDC_DBA_to_derep <- DF_PDC_DBA_to_derep[DF_PDC_DBA_to_derep$sys_id %in% c("NEXT_MC03X035201F7_N131_L93942_K16.5_E1_P0_F0_dar_ddr_hdf_ulx_1945",
                                                                             "NEXT_V0107_N31_L43417_K51.9_E0_P0_F0_gam_2294",
                                                                             "NEXT_V0127_N160_L38840_K35_E0_P0_F0_acriia25_1811",
                                                                             "NEXT_V0127_N160_L38840_K35_E0_P0_F0_acriia3_1812",
                                                                             "NEXT_V0272_N29_L45614_K433.2_E0_P0_F0_acriia3_1819",
                                                                             "NEXT_V0272_N29_L45614_K433.2_E0_P0_F0_acriia25_1818",
                                                                             "NEXT_V0291_N27_L39952_K38.2_E0_P0_F0_acriia3_2209",
                                                                             "NEXT_V0291_N27_L39952_K38.2_E0_P0_F0_acriia25_2208",
                                                                             "NEXT_V0530_N1982_L5157_K0.9_E0_P0_F0_acric5_2989",
                                                                             "NEXT_V1124_N89_L37865_K41.8_E0_P0_F0_acriia3_1869",
                                                                             "NEXT_V1124_N89_L37865_K41.8_E0_P0_F0_acriia25_1868",
                                                                             "NEXT_V3010_N40_L43254_K1.2_E0_P0_F0_acriia3_1877",
                                                                             "NEXT_V3010_N40_L43254_K1.2_E0_P0_F0_acriia25_1876",
                                                                             "VREF_NC_049445.1_atd1_2192"),]



DF_PDC_DBA_to_derep$protein_in_syst[DF_PDC_DBA_to_derep$sys_id == "NEXT_MC03X035201F7_N131_L93942_K16.5_E1_P0_F0_dar_ddr_hdf_ulx_1945"] <- 
  "NEXT_MC03X035201F7_N131_L93942_K16.5_E1_P0_F0_8,NEXT_MC03X035201F7_N131_L93942_K16.5_E1_P0_F0_9,NEXT_MC03X035201F7_N131_L93942_K16.5_E1_P0_F0_10,NEXT_MC03X035201F7_N131_L93942_K16.5_E1_P0_F0_11,NEXT_MC03X035201F7_N131_L93942_K16.5_E1_P0_F0_13"

DF_PDC_DBA_to_derep$genes_count[DF_PDC_DBA_to_derep$sys_id == "NEXT_MC03X035201F7_N131_L93942_K16.5_E1_P0_F0_dar_ddr_hdf_ulx_1945"] <- 5
DF_PDC_DBA_to_derep$name_of_profiles_in_sys[DF_PDC_DBA_to_derep$sys_id == "NEXT_MC03X035201F7_N131_L93942_K16.5_E1_P0_F0_dar_ddr_hdf_ulx_1945"] <- "dara,ddra,ddrb,ddrb,hdf"


DF_PDC_DBA_to_derep$protein_in_syst[DF_PDC_DBA_to_derep$sys_id == "NEXT_V0107_N31_L43417_K51.9_E0_P0_F0_gam_2294"] <- 
  "NEXT_V0107_N31_L43417_K51.9_E0_P0_F0_1,NEXT_V0107_N31_L43417_K51.9_E0_P0_F0_2,NEXT_V0107_N31_L43417_K51.9_E0_P0_F0_70"

DF_PDC_DBA_to_derep$genes_count[DF_PDC_DBA_to_derep$sys_id == "NEXT_V0107_N31_L43417_K51.9_E0_P0_F0_gam_2294"] <- 3
DF_PDC_DBA_to_derep$name_of_profiles_in_sys[DF_PDC_DBA_to_derep$sys_id == "NEXT_V0107_N31_L43417_K51.9_E0_P0_F0_gam_2294"] <- "gam,gam,gam"

DF_PDC_DBA_to_derep$tool_final <- "both"
DF_PDC_DBA_to_derep$tool_final[DF_PDC_DBA_to_derep$sys_id == "VREF_NC_049445.1_atd1_2192"] <- "defensefinder"

DF_PDC_DBA <- rbind(DF_PDC_DBA, DF_PDC_DBA_to_derep)
DF_PDC_DBA$duplicates <- NULL

write.table(DF_PDC_DBA, "../results/Defense_antidefense_syss_combined_v2.tsv", sep='\t', row.names=F, col.names=T, quote=F)  # final table output