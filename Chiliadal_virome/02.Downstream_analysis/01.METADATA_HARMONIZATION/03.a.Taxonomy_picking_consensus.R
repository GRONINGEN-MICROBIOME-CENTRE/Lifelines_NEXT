setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

#############################################################
# Merging taxonomy assignment by geNomad, VITAP & RefSeq
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
# specific to geNomad fix, padding to N ;
pad_to_n <- function(x, N) {
  n <- str_count(x, fixed(";"))
  n[is.na(n)] <- 0L
  paste0(x, strrep(";", pmax(0L, N - n)))
}

# specific to geNomad fix, padding
fix_geNomad_taxonomy <- function(tax, source, pad_fun = pad_to_n) {
  stopifnot(length(tax) == length(source))
  x <- tax
  g <- source == "geNomad"
  
  # 1) Pre-fixing:  Family present but no Order → insert missing ';' after Class
  noorder <- g &
    grepl("viricetes", x) &
    !grepl("viricetes;;", x) &
    !grepl("ales;", x) &
    grepl("viridae", x)
  
  if (any(noorder)) {
    x[noorder] <- sub("viricetes", "viricetes;", x[noorder], fixed = TRUE)
  }
  
  # 2) Patterns that should be padded to 8 ranks
  pad_patterns <- c(
    "Viruses;*$",   # only Life
    "viria;*$",     # realm-level only (…viria)
    "virae;*$",     # nothing below Kingdom
    "viricota;*$",  # nothing below Phylum
    "viricetes;*$", # nothing below Class
    "ales;*$",      # nothing below Order
    "dae;*$"        # nothing below Family
  )
  
  for (p in pad_patterns) {
    idx <- g & grepl(p, x, perl = TRUE)
    if (any(idx)) x[idx] <- pad_fun(x[idx], 8)
  }
  
  # 3) Specific fixup
  fusello <- grepl("Viruses;Fuselloviridae;;", x)
  if (any(fusello)) {
    x[fusello] <- "Viruses;;;;;;Fuselloviridae;;"
  }
  
  x
}

# specific to RefSeq fix:
fix_RefSeq_taxonomy <- function(tax, source) {
  
  stopifnot(length(tax) == length(source))
  x <- tax
  r <- source == "RefSeq"
  
  # 1) Drop sub-ranks available only for well-studied human viruses then tidy separators
  drop_tokens_end <- function(s) {
    s <- gsub(
      "(^|;)[^;]*?(?:virinae|virineae|viricotina|Sarbecovirus)(?=;|$)",
      "", s, ignore.case = TRUE, perl = TRUE
    )
    s
  }
  x <- vapply(x, drop_tokens_end, character(1))
  
  gi   <- function(p) grepl(p, x, ignore.case = TRUE, perl = TRUE)
  subi <- function(p, r, i) sub(p, r, x[i], ignore.case = TRUE, perl = TRUE)
  gsubi <- function(p, r, i) gsub(p, r, x[i], ignore.case = TRUE, perl = TRUE)
  
  #2) Family present but no order → insert empty order after class
  noorder <- r & gi("viricetes") & !gi("viricetes;;") & !gi("ales;")
  if (any(noorder)) x[noorder] <- subi("viricetes", "viricetes;", noorder)
  
  # stricter fix when class is immediately followed by family (no order)
  class_then_family <- r & gi("viricetes;[^;]*viridae")
  if (any(class_then_family)) {
    x[class_then_family] <- gsubi("viricetes;", "viricetes;;", class_then_family)
  }
  
  # 3) No family but already has empty after class marker → add one more empty after class
  no_fam <- r & gi("viricetes;;") & !gi("viridae")
  if (any(no_fam)) x[no_fam] <- subi("viricetes", "viricetes;", no_fam)
  
  # 4) No genus:
  no_gen <- r & !gi("virus;")
  w_fam <- no_gen & gi("viridae")
  if (any(w_fam)) x[w_fam] <- gsubi("viridae", "viridae;", w_fam)
  wo_fam <- no_gen & !gi("viridae")
  if (any(wo_fam)) x[wo_fam] <- gsubi("viricetes", "viricetes;", wo_fam)
  
  # 5) Trailing-';' heuristics
  ends_with_patterns <- c(
    "Bocaparvovirus primate1",
    "Alphapapillomavirus 4",
    "nis",
    "Xuanwuvirus P884B11",
    "Peduovirus P24E6b",
    "Peduovirus P22H1"
  )
  ends_rx <- paste0("(", paste(ends_with_patterns, collapse = "|"), ")$")
  ends_idx <- grepl(ends_rx, x, perl = TRUE)
  if (any(ends_idx)) x[ends_idx] <- paste0(x[ends_idx], ";")
  
  # 6) Example-specific gap after class when 'phage' present
  phage_gap_idx <- grepl("phage", x, ignore.case = TRUE) & grepl("Caudoviricetes;;;;", x, perl = TRUE)
  if (any(phage_gap_idx)) x[phage_gap_idx] <- gsubi("viricetes", "viricetes;", phage_gap_idx)
  
  # 7) semicolon count → if exactly 8, double the last separator
  nsemi <- lengths(regmatches(x, gregexpr(";", x, perl = TRUE)))
  nsemi[is.na(nsemi)] <- 0L
  to_pad9 <- r & nsemi == 8L
  if (any(to_pad9)) x[to_pad9] <- sub(";(?!.*;)", ";;", x[to_pad9], perl = TRUE)
  
  # Clean-up
  tar_idx <- grepl("_Taranis;", x, perl = TRUE)
  if (any(tar_idx)) x[tar_idx] <- gsub("_Taranis;", "_Taranis", x[tar_idx], perl = TRUE)
  
  x
}
#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(stringr)
#############################################################
# 2. Load Input Data
#############################################################
ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Basic_ETOF_121062vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

# VITAP as the main source:
VITAP <- read.table('06.CLEAN_DATA/best_determined_lineages.tsv', sep='\t', header=T)
VITAP$VITAP_taxonomy <- paste0('Viruses;', sapply(VITAP$lineage, function(x) {
  paste(rev(strsplit(x, ";")[[1]]), collapse = ";")
}))

# RefSeq to rescue some unassigned:
refseq_taxonomy <- read.delim('06.CLEAN_DATA/viral_refseq_232_taxo.txt', sep='\t', header=T)

# here, vOTU clustering informaion is used to rescue vOTUs with unasigned taxonomy that were clustered w RefSeq sequences
vOTU_clustering <- read.table('06.CLEAN_DATA/NEXT_viral_clusters_MGS_VLP_long_format.txt', sep='\t', header=T)
vOTU_clustering <- vOTU_clustering[vOTU_clustering$Representative %in% ETOF_vOTUr$New_CID,]
#############################################################
# 3.1 Analysis: choosing the taxonomy assignemnt source
#############################################################
# geNomad:
taxa_master <- ETOF_vOTUr[,c("New_CID", "POST_CHV_length", "taxonomy")]
colnames(taxa_master)[grep('taxonomy', colnames(taxa_master))] <- "geNomad_taxonomy"
taxa_master[is.na(taxa_master$geNomad_taxonomy),]$geNomad_taxonomy <- "Unclassified"

# checking there are no conflicts among vOTU member assignments:
VREF <- vOTU_clustering[grepl('VREF', vOTU_clustering$Cluster_member),] # 251 VREF genomes clustered (with other genomes) to 195 vOTUs
VREF$RefSeq_OID <- gsub('VREF_', '', VREF$Cluster_member) # removing DB-specific tag (added during the DB processing)
VREF$RefSeq_taxonomy <- refseq_taxonomy$taxonomy[match(VREF$RefSeq_OID, refseq_taxonomy$genome_id)]

# checking that refseq viruses clustering to vOTUs indeed belong to same species:
mm_vref <- as.data.frame(table(VREF$Representative)) 
mm_vref <- mm_vref[mm_vref$Freq > 1,] # 32 vOTUs that have multiple members from Viral RefSeq
mm_vref <- VREF[VREF$Representative %in% mm_vref$Var1,]
# getting the species level
mm_vref$species <- gsub('^(([^ ]+ ){2}).*$', '\\1', sapply(strsplit(mm_vref$RefSeq_taxonomy, ';'), function(x) {
  if (length(x) >= 9) {
    x[9]           # 9th element when available
  } else {
    x[length(x)]   # last element otherwise
  }
}))

# calculating the number of different spp. per vOTU
refseq_check <- mm_vref %>%
  distinct(Representative, species) %>%
  count(Representative, name = "species_count") %>%
  filter(species_count > 1) # only a few instances, and all are "potato - potato"

# even though there are no conflicts, let's keep the original assignment to the V REF vOTUrs:
VREFrep <- VREF[grepl('VREF', VREF$Representative) & VREF$Representative == VREF$Cluster_member,]

# since no conflicts among v refseq members are spotted, a representative taxonomy can be chosen randomly:
VREF <- VREF %>%
  filter(!grepl("VREF", Representative)) %>%
  group_by(Representative) %>%
  slice_head(n=1)

VREF <- rbind(VREF, VREFrep)

# removing data from GE that I do not need:
rm(list = c("ETOF_vOTUr", "mm_vref", "refseq_check",
            "refseq_taxonomy", "vOTU_clustering", "VREFrep"))

# merging taxa annotations from different sources
taxa_master <- merge(taxa_master, 
                     VITAP[,c("Genome_ID", "VITAP_taxonomy")], 
                     by.x="New_CID",
                     by.y="Genome_ID",
                     all.x=T)
taxa_master[is.na(taxa_master$VITAP_taxonomy),]$VITAP_taxonomy <- "Unclassified"

taxa_master <- merge(taxa_master,
                     VREF[,c("Representative", "RefSeq_taxonomy")],
                     by.x = "New_CID",
                     by.y = "Representative",
                     all.x=T)

taxa_master[is.na(taxa_master$RefSeq_taxonomy),]$RefSeq_taxonomy <- "NOT_VREF"

# initial curation to match same taxonomy assignments
taxa_master$VITAP_mod <- gsub(';-;-;-;-','',taxa_master$VITAP_taxonomy)

taxa_master$geNomad_taxonomy_new <- gsub('Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;$', 
                                     'Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes', taxa_master$geNomad_taxonomy)

# parsing genera:
taxa_master$geNomad_family <- sapply(strsplit(taxa_master$geNomad_taxonomy_new, ";"), `[`, 7)
taxa_master[is.na(taxa_master$geNomad_family),]$geNomad_family <- "Unclassified"
taxa_master$geNomad_missed_family <- sapply(strsplit(taxa_master$geNomad_taxonomy_new, ";"), `[`, 6)
taxa_master[is.na(taxa_master$geNomad_missed_family),]$geNomad_missed_family <- "Unclassified"

taxa_master$VITAP_family <- sapply(strsplit(taxa_master$VITAP_taxonomy, ";"), `[`, 7)
taxa_master[is.na(taxa_master$VITAP_family),]$VITAP_family <- "-"

refseq_rename <- c(
  "Viruses;Human feces pecovirus",
  "Viruses;Thika virus",
  "Viruses;Chicken stool associated circular virus 2",
  "Viruses;Chicken stool associated circular virus 1"
)

# choosing the source:
taxa_master <- taxa_master %>%
  mutate(
    consensus_source = NA_character_,
    consensus_source = if_else(RefSeq_taxonomy != "NOT_VREF", "RefSeq", consensus_source),
    consensus_source = if_else(RefSeq_taxonomy %in% refseq_rename, "VITAP", consensus_source)
  ) %>%
  # handy aliases (only for matching; dropped at the end)
  mutate(g = geNomad_taxonomy_new, v = VITAP_taxonomy) %>%
  mutate(
    consensus_source = coalesce(
      consensus_source,
      case_when(
        # both unclassified / one unclassified cases
        g == "Unclassified" & v == "Unclassified" ~ "VITAP",
        g != "Unclassified" & v == "Unclassified" ~ "geNomad",
        g == "Unclassified" & v != "Unclassified" ~ "VITAP",
        
        # exact same classification
        g == VITAP_mod ~ "VITAP",
        
        # family/order-specific tweaks
        g == "Viruses;Anelloviridae" ~ "VITAP",
        str_detect(g, "Microviridae") ~ "VITAP",
        
        # Crassvirales / Caudoviricetes interplay
        str_detect(g, "Crassvirales") & str_detect(v, "Crassvirales") ~ "VITAP",
        str_detect(g, "Crassvirales") & str_detect(v, "Caudoviricetes;-") ~ "geNomad",
        str_detect(g, "Caudoviricetes$") & str_detect(v, "Crassvirales") ~ "VITAP",
        str_detect(g, "Caudoviricetes$") & !str_detect(v, "Caudoviricetes;-") ~ "VITAP",
        str_detect(g, "Caudoviricetes;") & str_detect(v, "Caudoviricetes;-") ~ "geNomad",
        
        # same genera / family match
        (geNomad_family == VITAP_family | geNomad_missed_family == VITAP_family) ~ "VITAP",
        
        # unassigned below Caudo in VITAP
        str_detect(v, ";-;-;-;-") ~ "geNomad",
        
        # some RNA families/orders
        str_detect(v, "Martellivirales|Negarnaviricota") ~ "VITAP",
        str_detect(g, "Picornavirales") ~ "geNomad",
        str_detect(v, "Tombusviridae") ~ "VITAP",
        str_detect(g, "Leviviricetes") ~ "geNomad",
        str_detect(v, "_Shirahamavirus") ~ "VITAP",
        
        # other families
        str_detect(g, "Aliceevansviridae|Demerecviridae|Solemoviridae") ~ "geNomad",
        str_detect(v, "Circoviridae|Betaflexiviridae") ~ "VITAP",
        
        # default fallthrough
        TRUE ~ "VITAP"
      )
    )
  ) %>%
  select(-g, -v)

table(taxa_master$consensus_source) # geNomad: 19789 vOTUs, RefSeq: 191, VITAP: 101082
# cleaning:
taxa_master$geNomad_family <- NULL
taxa_master$geNomad_missed_family <- NULL
taxa_master$VITAP_family <- NULL
rm(list = c("VITAP", "VREF", "refseq_rename"))

#############################################################
# 3.2 Analysis: unification of taxonomy entries from three sources
#############################################################

# preparing geNomad taxonomy:
taxa_master$geNomad_taxonomy_new <-
  fix_geNomad_taxonomy(taxa_master$geNomad_taxonomy_new, taxa_master$consensus_source)

# preparing RefSeq taxonomy:
taxa_master$RefSeq_taxonomy_new <- fix_RefSeq_taxonomy(taxa_master$RefSeq_taxonomy, taxa_master$consensus_source)

# final polishing:
taxa_master <- taxa_master %>%
  mutate(
    across(
      c(RefSeq_taxonomy_new, geNomad_taxonomy_new, VITAP_taxonomy),
      ~ str_count(., fixed(";")),
      .names = "{.col}_nsemi"
    )
  )

# padding (since RefSeq often has strains & species)
taxa_master$geNomad_taxonomy_new <- pad_to_n(taxa_master$geNomad_taxonomy_new, 9)
taxa_master$VITAP_taxonomy_new <- pad_to_n(taxa_master$VITAP_taxonomy, 9)
taxa_master$VITAP_taxonomy_new <- gsub('-;', ';', taxa_master$VITAP_taxonomy_new)

# masterfile
taxa_master <- taxa_master[,c("New_CID", "POST_CHV_length", "geNomad_taxonomy",
                        "VITAP_taxonomy", "RefSeq_taxonomy", "consensus_source",
                        "geNomad_taxonomy_new", "RefSeq_taxonomy_new", "VITAP_taxonomy_new")]

# to check that now all have 9 semicolons:
taxa_master <- taxa_master %>%
  mutate(
    across(
      c(RefSeq_taxonomy_new, geNomad_taxonomy_new, VITAP_taxonomy_new),
      ~ str_count(., fixed(";")),
      .names = "{.col}_nsemi"
    )
  )

# consensus unified lineage
taxa_master$consensus_lineage <- NA

taxa_master[taxa_master$consensus_source=="VITAP",]$consensus_lineage <- 
  taxa_master[taxa_master$consensus_source=="VITAP",]$VITAP_taxonomy_new

taxa_master[taxa_master$consensus_source=="geNomad",]$consensus_lineage <- 
  taxa_master[taxa_master$consensus_source=="geNomad",]$geNomad_taxonomy_new

taxa_master[taxa_master$consensus_source=="RefSeq",]$consensus_lineage <- 
  taxa_master[taxa_master$consensus_source=="RefSeq",]$RefSeq_taxonomy_new

# original source lineage
taxa_master$source_lineage <- NA

taxa_master[taxa_master$consensus_source=="VITAP",]$source_lineage <- 
  taxa_master[taxa_master$consensus_source=="VITAP",]$VITAP_taxonomy

taxa_master[taxa_master$consensus_source=="geNomad",]$source_lineage <- 
  taxa_master[taxa_master$consensus_source=="geNomad",]$geNomad_taxonomy

taxa_master[taxa_master$consensus_source=="RefSeq",]$source_lineage <- 
  taxa_master[taxa_master$consensus_source=="RefSeq",]$RefSeq_taxonomy

taxa_working <- taxa_master[,c("New_CID", "POST_CHV_length", "consensus_source", "source_lineage", "consensus_lineage")]
#############################################################
# 4. OUTPUT
#############################################################
# masterfile containing untouched taxonomy assignments etc:
write.table(taxa_master, '06.CLEAN_DATA/Master_taxonomy_121062vOTUs_not_curated.txt', sep='\t', quote=F, row.names = F)

# working file to run through ICTV & taxonomy expansion:
write.table(taxa_working, '06.CLEAN_DATA/Consensus_taxonomy_121062_not_curated.txt', sep='\t', quote=F, row.names = F)

