setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')
#############################################################
# Filtering RPKM table & supplementing metadata with 
# virome-based metrics (diversity, richness, temperate phage
# abundance)
#############################################################

#############################################################
# 0. Used files source
#############################################################
# Chiliadal_meta_ExtFiltered_v04.txt - long metadata for all successfully sequenced unique fecal samples from Chiliadal-included individuals (both available VLP and MGS samples)
# BC_PILOT_new_ids_chiliadal_format.txt - txt-format of Samples_BaseClear_pilot_sequencing.xlsx,  list of NEXT samples included in the Big Gut NEXT paper that were extracted and sequenced to test the new VLP protocol
# counts_filtered.txt - VLP MGS counts of reads aligned to DETECTED contigs (i.e., read counts of contigs not passing 75%-rule are excluded; the table is NOT decontaminated)
# VLP_MGS_decontaminated_RPKM_overlap.txt
# clean_RPKM_UPD.txt
# NEXT_viral_clusters_MGS_VLP_size.txt
#############################################################
# 1. Functions
#############################################################
# renaming VNP* columns in RPKM tables to CHV3*-format based on prechili:
rename_columns_prechili <- function(df, lookup, from_col = "V1", to_col = "V2") {
  matched <- colnames(df) %in% lookup[[from_col]]
  colnames(df)[matched] <- lookup[[to_col]][match(colnames(df)[matched], lookup[[from_col]])]
  return(df)
}

# counting the N of vOTUs in tables:
count_vOTUs <- function(df, samples, filter_rows = NULL) {
  df %>%
    select(all_of(intersect(colnames(df), samples))) %>%
    {
      if (!is.null(filter_rows)) filter(., rownames(.) %in% filter_rows) else .
    } %>%
    filter(rowSums(.) > 0) %>%
    nrow()
}

# Helper function to calculate % aligned reads
calc_perc_aligned <- function(count_mat, sample_ids, clean_reads_vec, filtered = NULL) {
  if (!is.null(filtered)) {
    count_mat <- count_mat[row.names(count_mat) %in% filtered, ]
  }
  aligned <- colSums(count_mat)
  sample_index <- match(sample_ids, names(aligned))
  aligned_reads <- aligned[sample_index]
  perc <- aligned_reads / clean_reads_vec * 100
  return(perc)
}

# Function to calculate per-sample richness
get_richness <- function(mat, sample_ids) {
  richness <- colSums(mat > 0, na.rm = TRUE)
  richness[match(sample_ids, names(richness))]
}

# Function to calculate shared vOTUs with controls
get_shared_with_controls <- function(mat, sample_ids, control_id = "controls") {
  sapply(sample_ids, function(i) {
    if (i %in% colnames(mat) && control_id %in% colnames(mat)) {
      sum(rowSums(mat[, c(i, control_id)] > 0) == 2)
    } else {
      NA
    }
  })
}

# get iqrs:
get_iqr <- function(vec) {
  
  SMR <- summary(vec)
  
  print(paste0(round(SMR[3], 2), ' (', round(SMR[2], 2), ' - ', round(SMR[5], 2), ')')) 
  
}

# seeks genome info from supplied ictv release table
# goes from Species to Realm
genome_info_seeker <- function(df, ictv) {
  
  # unsalvagable (all unclassified)
  uncl_idx <- rowSums(df[ ranks[c(-length(ranks), -1)] ] == "Unclassified") == length(ranks[c(-length(ranks), -1)])
  df$genome[uncl_idx] <- "Unclassified"
  
  for (rnk in rev(ranks) ) {
    
    if(!rnk %in% colnames(ictv)) next
    
    classified_rnk <- df %>%
      filter(! (!!sym(rnk)) %in% c("Unclassified", "Unassigned")) %>%
      filter(is.na(genome)) %>%
      pull(!!sym(rnk))
    
    ictv_lookup <- ictv %>%
      filter(!is.na(!!sym(rnk))) %>%
      filter(!!sym(rnk) %in% unique(classified_rnk)) %>%
      select(!!sym(rnk), genome_simple) %>%
      distinct()
    
    df <- left_join(df, ictv_lookup, by = rnk)
    df$genome[is.na(df$genome)] <- df$genome_simple[is.na(df$genome)]
    df$genome_simple <- NULL
  }
  
  return(df)
}

# seeks Host info from supplied ictv species release table
# goes from Species to Realm
assign_host_domain <- function(df, ictv_w_host) {
  
  uncl_idx <- rowSums(df[ ranks[c(-length(ranks), -1)] ] == "Unclassified") == length(ranks[c(-length(ranks), -1)])
  df$Host[uncl_idx] <- "Unknown"
  
  for (rnk in rev(ranks)) {
    
    if (! rnk %in% colnames(ictv_w_host)) next
    
    classified_rnk <- df %>%
      filter(is.na(Host)) %>% # since we go from species to realm, some ambigious can be filled up already, no need to reqrite them
      distinct(!!sym(rnk)) %>%
      pull()
    
    ictv_rnk <- ictv_w_host %>%
      filter(! is.na(!!sym(rnk))) %>%
      filter(!!sym(rnk) %in% classified_rnk) %>%
      select(!!sym(rnk), Host_simple) %>%
      count(across(everything())) %>% # choosing the host based on majority
      group_by(!!sym(rnk)) %>%
      arrange(desc(n), Host_simple) %>%
      slice_head(n=1)
    
    df <- df %>%
      left_join(ictv_rnk %>% select(-n), by=rnk)
    df$Host[is.na(df$Host)] <- df$Host_simple[is.na(df$Host)]
    df$Host_simple <- NULL
    
  }
  
  df$Host_simple <- df$Host
  df$Host_simple[grep("bacteria|archaea", df$Host_simple)] <- "Prokaryote"
  df$Host_simple[!df$Host_simple %in% c("Prokaryote", "Unknown")] <- "Eukaryote"
  return(df)
}

# identifies the DB origin of the vOTU cluster member
# writes it to column "DB_member"
get_member_origin <- function(df, member_col, dbs, class_col="DB_member"){
  
  df[class_col] <- NA
  
  for (db in dbs) {
    
    idx <- grepl(db, df[[member_col]], fixed = T)
    
    df[[class_col]][idx] <- gsub('-', '', db, fixed = T)
    
  }
  df[[class_col]][grepl('NEXT_V', df[[member_col]], fixed = T)] <- "NEXT_VLP"
  df[[class_col]][grepl('NEXT_M', df[[member_col]], fixed = T)] <- "NEXT_MGS"
  
  df
}

# gets relative abundance summed according to aggregator
get_RAb <- function(aggregator, mat, sample_ids) {
  
  RAb <- colSums(mat[row.names(mat) %in% aggregator,])/colSums(mat)*100
  RAb[match(sample_ids, names(RAb))]
}

# pulls RPKM table at the given rank
pull_rank <- function(taxa_table, rank, samples) {
  
  by_rank <- taxa_table %>%
    select(all_of(samples), rank) %>%
    group_by(!!sym(rank)) %>%
    summarise(across(everything(), sum), .groups = "drop") %>%
    mutate(across(all_of(samples), ~ .x / sum(.x)) * 100) %>%
    column_to_rownames(var = rank)
  
  return(by_rank)
}
#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(tidyr)
library(tidyverse)
library(vegan)
#############################################################
# 2. Load Input Data
#############################################################

# calculating the decontamination efficacy:
esmeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_ExtFiltered_v05.txt', sep='\t', header=T) 
esmeta$Timepoint_new <- esmeta$Timepoint_original
esmeta[esmeta$Type=="M",]$Timepoint_new <- 'Mother'
esmeta$Timepoint_new <- factor(esmeta$Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)

prechili <- read.table('~/Desktop/Projects_2021/Baseclear/01.TEST_22_samples/BC_PILOT_new_ids_chiliadal_format.txt', sep='\t')
prechili <- prechili[,c('V1', 'V2')]

# count table:
counts <- read.table('04.RAW_DATA/05.Abundance_tables/VLP_MGS_overlap/counts_filtered_UPD_w_mgc_nc.txt', sep='\t', header=T)
counts <- rename_columns_prechili(counts, prechili)

# filtering fecal samples not included in VLP MGS comparison:
esmeta <- esmeta[esmeta$Sequencing_ID %in% colnames(counts),]

controls_ids <- esmeta[!esmeta$Type %in% c('M', 'K'),]$Sequencing_ID

full_overlap_sids <- esmeta[esmeta$full_overlap==T,]$Sequencing_ID

# raw RPKM table
dirty <- read.table("04.RAW_DATA/05.Abundance_tables/VLP_MGS_overlap/RPKM_counts_VLP.txt", sep='\t', header=T)
dirty <- rename_columns_prechili(dirty, prechili)
dirty$controls <- rowSums(dirty[,colnames(dirty) %in% controls_ids] > 0)

# keeping only detected vOTUs in the ETOF:
ETOF <- data.table::fread('06.CLEAN_DATA/VLP_MGS_ETOF_full_rep.txt', sep='\t', header=T)
ETOFint <- ETOF[ETOF$New_CID %in% row.names(dirty),]

ab3kb <- ETOFint[ETOFint$POST_CHV_length >= 3000,]$New_CID

# decontaminated RPKM table:
clean <- read.table('04.RAW_DATA/05.Abundance_tables/VLP_MGS_overlap/clean_RPKM_UPD_w_mgc_nc.txt', sep='\t', header=T)
clean <- rename_columns_prechili(clean, prechili)
clean$controls <- rowSums(clean[,colnames(clean) %in% controls_ids] > 0)

# decontaminated count table:
cleaned_counts <- counts
cleaned_counts[clean[,-ncol(clean)]==0] <- 0

# CheckV clustering output
vOTU_cluster_size <- read.table('06.CLEAN_DATA/NEXT_viral_clusters_MGS_VLP_size.txt', sep='\t', header=T)
vOTU_clustering <- read.table('06.CLEAN_DATA/NEXT_viral_clusters_MGS_VLP_long_format.txt', sep='\t', header=T)

# all DBs used for the dereplication:
dbs <- c('NEXT_V', 'NEXT_M', 'Guerin', 
         'Yutin', 'NCBI_CrAss', 'NL_crAss',
         'Benler', 'COPSAC', 'GVD', 
         'IMGVR', 'MGV-', 'GPD', 'VREF')

vOTU_clustering <- get_member_origin(vOTU_clustering, "Cluster_member", dbs)

# BACPHLIP output
lifestyle_raw <- read.table('06.CLEAN_DATA/Lifestyle_prediction_bacphlip', sep='\t', header=T)

# curated & expanded taxonomy:
ranks <- c("Strain","Species","Genus","Family","Order","Class","Phylum","Kingdom","Realm", "Life")
ranks <- rev(ranks)

# ictv release:
ictv <- readxl::read_xlsx('06.CLEAN_DATA/ICTV_Master_Species_List_2024_MSL40.v2.xlsx', sheet=2, col_types = "text")
# simplification of genome entries:
ictv$genome_simple <- ictv$Genome
ictv$genome_simple[grep("ssDNA", ictv$genome_simple)] <- "ssDNA"
ictv$genome_simple[grep("dsDNA", ictv$genome_simple)] <- "dsDNA" # because we expect dsDNA-RT to be caught more frequently in their dsDNA form
ictv$genome_simple[grep("RNA", ictv$genome_simple)] <- "RNA"
ictv$genome_simple[!is.na(ictv$Genus) & ictv$Species == "Begomovirus chayotis"] <- "ssDNA" # I guess I googled it
# ictv release w host domain:
ictv_VMR <- readxl::read_xlsx('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/06.CLEAN_DATA/VMR_MSL40.v2.20251013.xlsx', sheet=2, col_types = "text")
ictv_VMR$Host_simple <- ictv_VMR$`Host source`
ictv_VMR$Host_simple[grep("(S)", ictv_VMR$`Host source`)] <- "Unknown"

# ICTV - curated taxonomy:
taxonomy_ictv <- read.table("06.CLEAN_DATA/Consensus_taxonomy_121062_ICTV_curated.txt", sep='\t', header=T)
taxonomy_ictv$tax_ictv <- apply(taxonomy_ictv[ranks], 1, function(rank)
  paste(rank, collapse = ";")
)

# ICTV curated and expanded taxonomy:
taxonomy_ictv_aai <- read.table("06.CLEAN_DATA/02.FINAL/Consensus_taxonomy_121062_ICTV_curated_AAI_expanded.txt", sep='\t', header=T)
taxonomy_ictv_aai$tax_ictv_aai <- apply(taxonomy_ictv_aai[ranks], 1, function(rank)
  paste(rank, collapse = ";")
)
# adding genome info:
taxonomy_ictv_aai$genome[taxonomy_ictv_aai$New_CID == "NEXT_V1007_N66_L4814_K1.7_E0_P0_F0"] <- "RNA" # because only Realm is available and most Riboviria are RNA
taxonomy_ictv_aai <- genome_info_seeker(taxonomy_ictv_aai, ictv)
# adding host info:
taxonomy_ictv_aai$Host <- NA
taxonomy_ictv_aai <- assign_host_domain(taxonomy_ictv_aai, ictv_VMR)

# creating merged table from two taxonomy options:
taxonomy <- taxonomy_ictv %>%
  select(New_CID, consensus_source, source_lineage, tax_ictv ) %>%
  left_join(taxonomy_ictv_aai %>%
              select(New_CID, tax_ictv_aai, Genus_OTUr, Family_OTUr, genome, Host, Host_simple),
            by = "New_CID")

# prokaryotic host prediction:
prok_host <- read.table('06.CLEAN_DATA/Host_prediction_to_genus_filtered_v2.tsv', sep='\t', header=T)

# PPV size and abundance:
VLP_PPV_size <- read.table('06.CLEAN_DATA/Intermediate/VLP_PPV_size_1110samples_120997vOTUrs.txt', sep='\t', header=T)
persistent_VLP_vOTUs <- read.table('06.CLEAN_DATA/Intermediate/List_VLP_PPV_vOTUrs_from_120997vOTUrs.txt', sep='\t', header=F)
persistent_MGS_vOTUs <- read.table('06.CLEAN_DATA/Intermediate/List_MGS_PPV_vOTUrs_from_120997vOTUrs.txt', sep='\t', header=F)

rm(list = c("ranks", "taxonomy_ictv", "taxonomy_ictv_aai", "ictv", "ictv_VMR"))
#############################################################
# 3.1 Analysis (collecting decontamination and filtering statistics)
#############################################################
# N vOTUs detected in 2,220 VLP and MGS samples depending on decontamination & filter length

N_vOTUs <- c(
  "dirty" = count_vOTUs(dirty, full_overlap_sids),
  "dirty filtered" = count_vOTUs(dirty, full_overlap_sids, ab3kb),
  "clean" = count_vOTUs(clean, full_overlap_sids),
  "clean filtered" = count_vOTUs(clean, full_overlap_sids, ab3kb)
)

# calculating the % aligned reads
esmeta$perc_aligned_d  <- calc_perc_aligned(counts, esmeta$Sequencing_ID, esmeta$clean_reads)

esmeta$perc_aligned_df <- calc_perc_aligned(counts, esmeta$Sequencing_ID, esmeta$clean_reads, ab3kb)

esmeta$perc_aligned_c  <- calc_perc_aligned(cleaned_counts, esmeta$Sequencing_ID, esmeta$clean_reads)

esmeta$perc_aligned_cf <- calc_perc_aligned(cleaned_counts, esmeta$Sequencing_ID, esmeta$clean_reads, ab3kb)

# calculating viral richness
esmeta$vir_richness_d <- get_richness(dirty, esmeta$Sequencing_ID)
esmeta$N_shared_cs_d <- get_shared_with_controls(dirty, esmeta$Sequencing_ID, control_id = "controls")

esmeta$vir_richness_df <- get_richness(dirty[row.names(dirty) %in% ab3kb,], esmeta$Sequencing_ID)
esmeta$N_shared_cs_df <- get_shared_with_controls(dirty[row.names(dirty) %in% ab3kb,], esmeta$Sequencing_ID, control_id = "controls")

esmeta$vir_richness_c <- get_richness(clean, esmeta$Sequencing_ID)
esmeta$N_shared_cs_c <- get_shared_with_controls(clean, esmeta$Sequencing_ID, control_id = "controls")

esmeta$vir_richness_cf <- get_richness(clean[row.names(clean) %in% ab3kb,], esmeta$Sequencing_ID)
esmeta$N_shared_cs_cf <- get_shared_with_controls(clean[row.names(clean) %in% ab3kb,], esmeta$Sequencing_ID, control_id = "controls")

# summary table for decontamination and filtering statistics:

# percentage of clean reads aligned to the vOTUs

perc_aligned_vlp <- c("dirty" = get_iqr(esmeta[esmeta$seq_type=="VLP" & esmeta$full_overlap==T,]$perc_aligned_d),
                    "dirty filtered" =  get_iqr(esmeta[esmeta$seq_type=="VLP" & esmeta$full_overlap==T,]$perc_aligned_df),
                    "clean" = get_iqr(esmeta[esmeta$seq_type=="VLP" & esmeta$full_overlap==T,]$perc_aligned_c),
                    "clean filtered" = get_iqr(esmeta[esmeta$seq_type=="VLP" & esmeta$full_overlap==T,]$perc_aligned_cf)
                  )

perc_aligned_mgs <- c("dirty" = get_iqr(esmeta[esmeta$seq_type=="MGS" & esmeta$full_overlap==T,]$perc_aligned_d),
                      "dirty filtered" =  get_iqr(esmeta[esmeta$seq_type=="MGS" & esmeta$full_overlap==T,]$perc_aligned_df),
                      "clean" = get_iqr(esmeta[esmeta$seq_type=="MGS" & esmeta$full_overlap==T,]$perc_aligned_c),
                      "clean filtered" = get_iqr(esmeta[esmeta$seq_type=="MGS" & esmeta$full_overlap==T,]$perc_aligned_cf)
)


# metrics to summarize:
cols_to_summarize <- c("vir_richness_d", "vir_richness_df", "vir_richness_c", "vir_richness_cf",
                       "N_shared_cs_d", "N_shared_cs_df", "N_shared_cs_c", "N_shared_cs_cf")

# IQR summarization
iqr_summary <- esmeta %>%
  filter(full_overlap == TRUE) %>%
  pivot_longer(cols = all_of(cols_to_summarize), names_to = "metric", values_to = "value") %>%
  mutate(
    category = case_when(
      grepl("vir_richness", metric) ~ "richness",
      grepl("N_shared", metric) ~ "shared"
    ),
    preprocessing = case_when(
      grepl("_d$", metric) ~ "dirty",
      grepl("_df$", metric) ~ "dirty filtered",
      grepl("_c$", metric) ~ "clean",
      grepl("_cf$", metric) ~ "clean filtered"
    )
  ) %>%
  group_by(Timepoint_new, seq_type, category, preprocessing) %>%
  summarise(
    IQR = get_iqr(value),
    .groups = "drop"
  )

iqr_summary_wide <- iqr_summary %>%
  unite("column_name", Timepoint_new, seq_type, category, sep = "_") %>%
  pivot_wider(
    names_from = column_name,
    values_from = IQR
  )

iqr_summary_wide <- as.data.frame(iqr_summary_wide)
row.names(iqr_summary_wide) <- iqr_summary_wide$preprocessing
iqr_summary_wide$preprocessing <- NULL

decontam_stat <- cbind(N_vOTUs,
  perc_aligned_vlp, perc_aligned_mgs)

decontam_stat <- merge(decontam_stat, iqr_summary_wide, by="row.names")
colnames(decontam_stat)[1] <- "RPKM_table_processing"

# removing bulky used dfs:
rm(list = c("cleaned_counts", "counts", 
            "dirty", "cols_to_summarize", "iqr_summary","N_vOTUs",
            "perc_aligned_mgs", "perc_aligned_vlp", "ETOFint"))
#############################################################
# 3.2 Analysis (generating clean RPKM tables)
#############################################################

# clean (and filtered) RPKM table
clean$controls <- NULL
clean <- clean[row.names(clean) %in% ab3kb,]
clean <- clean[rowSums(clean)>0,]

# cleanest RPKM table to save (only 2,220 VLP and MGS metaviromes)
cleanest <- clean[, colnames(clean) %in% esmeta[esmeta$full_overlap==T,]$Sequencing_ID]
cleanest <- cleanest[rowSums(cleanest)>0,]

# cleanest RPKM table for VLPs only:
cleanest_VLP <- cleanest[, colnames(cleanest) %in% esmeta[esmeta$seq_type=="VLP",]$Sequencing_ID]
cleanest_VLP <- cleanest_VLP[rowSums(cleanest_VLP)>0,]

# cleanest RPKM table for VLPs only:
cleanest_MGS <- cleanest[, colnames(cleanest) %in% esmeta[esmeta$seq_type=="MGS",]$Sequencing_ID]
cleanest_MGS <- cleanest_MGS[rowSums(cleanest_MGS)>0,]

#############################################################
# 3.2 Analysis (generating basic and extended ETOFs)
#############################################################
# basic ETOF:
baseETOF_vOTUr <- ETOF[ETOF$New_CID %in% row.names(cleanest),]

# extended ETOF:

# removing useless all-NA columns
keeper <- colnames(baseETOF_vOTUr)[colSums(is.na(baseETOF_vOTUr))!=nrow(baseETOF_vOTUr), drop = F]
ETOF_vOTUr <- baseETOF_vOTUr[, ..keeper] # ..keeper is because ETOF was read as data.table, not data.frame

#adding dereplication information to the ETOF:
ETOF_vOTUr$vOTU_size <- vOTU_cluster_size$Cluster_size[match(ETOF_vOTUr$New_CID, vOTU_cluster_size$Representative)]

# adding lifestyle assignment:
ETOF_vOTUr <- merge(ETOF_vOTUr, lifestyle_raw, by.x ="New_CID", by.y = "X", all.x=T)
ETOF_vOTUr$lifestyle <- NA
ETOF_vOTUr$lifestyle <- as.character(ETOF_vOTUr$lifestyle)
ETOF_vOTUr[ETOF_vOTUr$Temperate >=0.5,]$lifestyle <- "Temperate"
ETOF_vOTUr[ETOF_vOTUr$Virulent >0.5 & ETOF_vOTUr$completeness >= 90 & !is.na(ETOF_vOTUr$completeness),]$lifestyle <- "Virulent"
ETOF_vOTUr[is.na(ETOF_vOTUr$lifestyle),]$lifestyle <- "Unknown"

# adding taxonomy:
ETOF_vOTUr <- ETOF_vOTUr %>%
  select(-taxonomy) %>%
  left_join(taxonomy, by = "New_CID") %>%
  left_join(prok_host %>%
              rename(Host_taxonomy=Host.genus), by = c("New_CID" = "Virus")) %>%
  mutate(Host_taxonomy = ifelse(is.na(Host_taxonomy), "d__Unclassified;p__Unclassified;c__Unclassified;o__Unclassified;f__Unclassified;g__Unclassified", Host_taxonomy ))

# refining the lifestyle assignment:
ETOF_vOTUr[ETOF_vOTUr$Host_simple == "Eukaryote",]$lifestyle <- "Unknown"
ETOF_vOTUr[ETOF_vOTUr$genome == "RNA",]$lifestyle <- "Unknown"

# adding novelty column:
vOTU_clustering <- vOTU_clustering[vOTU_clustering$Representative %in% ETOF_vOTUr$New_CID,]

# calculating db contributions:
sources <- c("NEXT_VLP","NEXT_MGS","DB_sum")

vOTU_by_source <- vOTU_clustering %>%
  mutate(DB_member = factor(DB_member)) %>%
  count(Representative, DB_member) %>%
  pivot_wider(names_from = DB_member, values_from = n, values_fill = 0) %>%
  left_join(vOTU_cluster_size, by = "Representative") %>%
  rename(N_genomes = Cluster_size, vOTU_representative  = Representative) %>%
  mutate(DB_sum = N_genomes - NEXT_MGS - NEXT_VLP) %>%
  mutate(
    vOTU_cluster_type = apply(across(all_of(sources)) > 0, 1,
                              function(r) { x <- sources[r]; if (length(x)==0) "None" else paste(x, collapse = "+") })
  ) %>%
  mutate(vOTU_novelty = if_else(vOTU_cluster_type %in% c("NEXT_MGS", "NEXT_VLP+NEXT_MGS", "NEXT_VLP"), "novel", "described"))
  
ETOF_vOTUr <- ETOF_vOTUr %>%
  left_join(vOTU_by_source %>% select(vOTU_representative, vOTU_cluster_type, vOTU_novelty), by = c("New_CID" = "vOTU_representative"))

# adding prevalence in VLP column:
ETOF_vOTUr$prev_VLP <- rowSums(cleanest_VLP > 0)[match(ETOF_vOTUr$New_CID, row.names(cleanest_VLP))]
ETOF_vOTUr$prev_VLP[is.na(ETOF_vOTUr$prev_VLP)] <- 0
# adding prevalence in MGS column:
ETOF_vOTUr$prev_MGS <- rowSums(cleanest_MGS > 0)[match(ETOF_vOTUr$New_CID, row.names(cleanest_MGS))]
ETOF_vOTUr$prev_MGS[is.na(ETOF_vOTUr$prev_MGS)] <- 0

# # adding PPV column:
# it was supposed to be computed here as well, but due to the memory size limitation,
# and since I did not want to learn about extensive data.table solution, 
# I have computed it outside with:

#
# script name place holder (currently "./02.SCRIPTS/06.DYNAMICS/01.Draft_explore.R")
#

# persistent in VLP

ETOF_vOTUr$PPV_VLP <- if_else(ETOF_vOTUr$New_CID %in% persistent_VLP_vOTUs$V1, "yes", "no")
ETOF_vOTUr$PPV_MGS <- if_else(ETOF_vOTUr$New_CID %in% persistent_MGS_vOTUs$V1, "yes", "no")

# persistent_holo
ETOF_vOTUr$PPV_holo <- if_else(ETOF_vOTUr$New_CID %in% persistent_VLP_vOTUs$V1 | 
                                 ETOF_vOTUr$New_CID %in% persistent_MGS_vOTUs$V1, "yes", "no")

skimr::skim(ETOF_vOTUr)
# select columns and release to Nataliia

keep_work_ETOF <- c("New_CID", "provirus", "POST_CHV_length", 
                    "proviral_length", "gene_count", "viral_genes",
                    "host_genes", "checkv_quality", "miuvig_quality",
                    "completeness", "contamination", "kmer_freq",
                    "warnings", "Original_length", "PRU_status",
                    "vOTU_size", "lifestyle", "consensus_source",
                    "source_lineage", "tax_ictv", "tax_ictv_aai", "Genus_OTUr",
                    "Family_OTUr", "genome", "Host",
                    "Host_simple", "Host_taxonomy", "vOTU_cluster_type",
                    "vOTU_novelty", "PPV_VLP", "PPV_MGS",
                    "PPV_holo", "prev_VLP", "prev_MGS")

working_ETOF <- ETOF_vOTUr[, ..keep_work_ETOF]

skimr::skim(working_ETOF)

fs_etof <- c("provirus", "checkv_quality", "miuvig_quality",
             "PRU_status", "lifestyle", "consensus_source",
             "genome", "Host", "Host_simple", 
             "vOTU_cluster_type", "vOTU_novelty", 
             "PPV_VLP", "PPV_MGS", "PPV_holo")

wETOF_sumstat <- working_ETOF %>%
  mutate(across(all_of(fs_etof), as.factor)) %>%
  skimr::skim() %>%
  select(-all_of(c("character.min", "character.max", "character.empty", "character.whitespace")))
#############################################################
# 3.3 Analysis (supplementing esmeta w viral metrics)

# !! Here, in 3.3, I use both clean and cleanest tables,
# because I also want to know phage lifestyle metrics in NCs
#############################################################
# viral alpha-diversity:
esmeta$vir_diversity <- diversity(t(clean))[match(esmeta$Sequencing_ID, colnames(clean))]

# phage lifestyle dynamics:
temperate_phages <- working_ETOF$New_CID[working_ETOF$lifestyle=="Temperate"]
lytic_phages <- working_ETOF$New_CID[working_ETOF$lifestyle=="Virulent"]

# richness of temperate phages
esmeta$temperate_richness <- get_richness(clean[row.names(clean) %in% temperate_phages,], esmeta$Sequencing_ID)

# richness of lytic phages
esmeta$lytic_richness <- get_richness(clean[row.names(clean) %in% lytic_phages,], esmeta$Sequencing_ID)

# fraction of temperate phages from the total number of detected viruses
esmeta$temperate_perc <- esmeta$temperate_richness/esmeta$vir_richness_cf*100

# fraction of lytic phages from the total number of detected viruses
esmeta$lytic_perc <- esmeta$lytic_richness/esmeta$vir_richness_cf*100

# ratio between temperate and lytic phages per sample:
esmeta$temp_to_lytic_ratio <- (esmeta$temperate_richness + 1)/(esmeta$lytic_richness + 1 )

# relative abundance of temperate phages
esmeta$temperate_RAb <- get_RAb(temperate_phages, clean, esmeta$Sequencing_ID)

# relative abundance of lytic phages
esmeta$lytic_RAb <- get_RAb(lytic_phages, clean, esmeta$Sequencing_ID)

# virome enrichment (ViromeQC) based on single-copy genes only:

vqc <- esmeta %>%
  filter(Type %in% c('M', 'K')) %>%
  select(Universal_ID, bacterial_markers_alignment_rate, seq_type) %>%
  pivot_wider(id_cols = Universal_ID,
              names_from = seq_type,
              values_from = bacterial_markers_alignment_rate) %>%
  mutate(sc_enrichment_vlp = pmin(MGS / VLP, 100)) %>%
  mutate(sc_enrichment_mgs = pmin(median(MGS)/MGS, 100))

esmeta$sc_enrichment <- NA
esmeta$sc_enrichment[esmeta$seq_type == "VLP"] <- vqc$sc_enrichment_vlp[match(esmeta$Universal_ID[esmeta$seq_type == "VLP"], vqc$Universal_ID)]
esmeta$sc_enrichment[esmeta$seq_type == "MGS"] <- vqc$sc_enrichment_mgs[match(esmeta$Universal_ID[esmeta$seq_type == "MGS"], vqc$Universal_ID)]

# same, but for NCs (not sure if it makes for all of them, especially for no-template blanks)
vqc_nc <- esmeta %>%
  filter(!Type %in% c('M', 'K')) %>%
  select(Sequencing_ID, bacterial_markers_alignment_rate) %>%
  mutate(sc_enrichment = pmin(median(vqc$MGS)/bacterial_markers_alignment_rate, 100))

esmeta$sc_enrichment[esmeta$Sequencing_ID %in% controls_ids] <- vqc_nc$sc_enrichment[match(esmeta$Sequencing_ID[esmeta$Sequencing_ID %in% controls_ids], vqc_nc$Sequencing_ID)]

# PPV size
# it was supposed to be computed here as well, but due to the memory size limitation,
# and since I did not want to learn about extensive data.table solution, 
# I have computed it outside with:

#
# script name place holder (currently "./02.SCRIPTS/06.DYNAMICS/01.Draft_explore.R")
#

esmeta$PPV_fraction <- NA
esmeta$PPV_abundance <- NA

esmeta[esmeta$seq_type == "VLP",]$PPV_fraction <- VLP_PPV_size$PPV_fraction[match(esmeta[esmeta$seq_type == "VLP",]$Universal_ID, VLP_PPV_size$Universal_ID)]
esmeta[esmeta$seq_type == "VLP",]$PPV_abundance <- VLP_PPV_size$PPV_abundance[match(esmeta[esmeta$seq_type == "VLP",]$Universal_ID, VLP_PPV_size$Universal_ID)]

#############################################################
# 3.4 Analysis (supplementing esmeta w metrics based on virus 
# taxonomy / genome characteristics)
#############################################################

# calculating the percentage of viruses by genome type
vOTU_by_feature <- cleanest %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(working_ETOF) %>%
  column_to_rownames(var = "New_CID")

by_genome <- vOTU_by_feature %>%
  select(all_of(esmeta$Sequencing_ID[esmeta$full_overlap == T]), genome) %>%
  group_by(genome) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  mutate(across(all_of(esmeta$Sequencing_ID[esmeta$full_overlap == T]), ~ .x / sum(.x)) * 100) %>%
  pivot_longer(
    cols = all_of(esmeta$Sequencing_ID[esmeta$full_overlap == T]),
    names_to = "Sequencing_ID",
    values_to = "abundance"
  ) %>%
  pivot_wider(
    names_from = genome,
    values_from = abundance
  ) 

esmeta <- esmeta %>%
  left_join(by_genome) %>%
  rename("genome_unclassified" = "Unclassified")

# calculating the percentage of viruses by simplified host

by_simple_host <- vOTU_by_feature %>%
  select(all_of(esmeta$Sequencing_ID[esmeta$full_overlap == T]), Host_simple) %>%
  group_by(Host_simple) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  mutate(across(all_of(esmeta$Sequencing_ID[esmeta$full_overlap == T]), ~ .x / sum(.x)) * 100) %>%
  pivot_longer(
    cols = all_of(esmeta$Sequencing_ID[esmeta$full_overlap == T]),
    names_to = "Sequencing_ID",
    values_to = "abundance"
  ) %>%
  pivot_wider(
    names_from = Host_simple,
    values_from = abundance
  )

esmeta <- esmeta %>%
  left_join(by_simple_host) %>%
  rename("host_unknown" = "Unknown")


by_host <- vOTU_by_feature %>%
  select(all_of(esmeta$Sequencing_ID[esmeta$full_overlap == T]), Host) %>%
  group_by(Host) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  mutate(across(all_of(esmeta$Sequencing_ID[esmeta$full_overlap == T]), ~ .x / sum(.x)) * 100) %>%
  column_to_rownames(var = "Host") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sequencing_ID")

esmeta <- esmeta %>%
  left_join(by_host %>% select(-Unknown))

#############################################################
# 3.5 Analysis (creating RPKM tables aggregated by virus rank)
#############################################################

# creating species, genera, families, class RPKM tables?
taxa <- vOTU_by_feature %>%
  select( all_of(esmeta$Sequencing_ID[esmeta$full_overlap == T]), tax_ictv_aai) %>%
  separate(col = tax_ictv_aai,
           into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
           sep = ";",
           fill   = "right",               
           remove = FALSE)

# species
by_sp <-  pull_rank(taxa, "Species", esmeta$Sequencing_ID[esmeta$full_overlap == T])

# genera
by_genus <- pull_rank(taxa, "Genus", esmeta$Sequencing_ID[esmeta$full_overlap == T])

# family
by_family <- pull_rank(taxa, "Family", esmeta$Sequencing_ID[esmeta$full_overlap == T])

# class
by_class <- pull_rank(taxa, "Class", esmeta$Sequencing_ID[esmeta$full_overlap == T])

#############################################################
# 3.6 Analysis (creating working esmeta and its summary stat)
#############################################################
# working esmeta:
remove_work_esmeta <- c("replicate_type", 
                        "perc_aligned_df", "perc_aligned_c", 
                        "vir_richness_df", "N_shared_cs_df", "vir_richness_c", "N_shared_cs_c")

working_esmeta <- esmeta %>%
  filter(full_overlap == T) %>%
  select( -(all_of(remove_work_esmeta)))

# working VLP-only esmeta:
VLP_working_esmeta <- esmeta %>%
  filter(full_overlap == T & seq_type == "VLP") %>%
  select( -(all_of(remove_work_esmeta))) %>%
  select_if(~ !all(is.na(.)))

# preparing for summary stat:
fs <- c("Type", "Timepoint_original", "seq_type", 
        "isolation_method", "replicate_type", "bgnp_status", 
        "metaphlan4_unclassified_high_contaminants_factor_75", 
        "Family_structure", "non_dyads", "Timepoint_categorical")

esmeta_sumstat <- esmeta %>%
  mutate(across(all_of(fs), as.factor)) %>%
  skimr::skim() %>%
  select(-all_of(c("character.min", "character.max", "character.empty", "character.whitespace")))

work_VLP_esmeta_sumstat <- VLP_working_esmeta %>%
  mutate(across(all_of( fs[fs %in% colnames(VLP_working_esmeta)] ), as.factor)) %>%
  skimr::skim() %>%
  select(-all_of(c("character.min", "character.max", "character.empty", "character.whitespace")))

#############################################################
# 4. Output
#############################################################

# decontamination efficacy summary:
write.table(decontam_stat, "./07.RESULTS/01.Decontamination_efficacy_VLP_MGS_overlap.txt", sep='\t', quote=F, row.names=F)

# RPKM table for VLP and MGS samples:
write.table(cleanest, '06.CLEAN_DATA/02.FINAL/RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_2220_samples.txt', sep='\t', quote=F)

# IDs of 121,062 vOTUrs detected in 2,220 VLP and MGS samples:
#write.table(row.names(cleanest), '06.CLEAN_DATA/IDs_vOTUrs_VLP_MGS_dec99ANI_ab3kbp_2220_samples.txt', col.names=F, row.names=F, sep='\t', quote=F)

# RPKM table for VLP samples only:
write.table(cleanest_VLP, '06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', quote=F)

# RPKM table for MGS samples only:
write.table(cleanest_MGS, '06.CLEAN_DATA/02.FINAL/MGS_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', quote=F)

# basic ETOF:
write.table(baseETOF_vOTUr, '06.CLEAN_DATA/02.FINAL/Basic_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', quote=F, row.names = F)

# working ETOF (w 1 temporary column):
write.table(working_ETOF, "06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt", sep='\t', quote=F, row.names=F)

# saving full esmeta:
write.table(esmeta,"./06.CLEAN_DATA/02.FINAL/Chiliadal_meta_Ext_v05_suppl_w_virmetrics.txt", sep='\t', quote=F)

# saving succinct esmeta:
write.table(working_esmeta,"./06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt", sep='\t', quote=F)

# saving succinct VLP esmeta:
write.table(VLP_working_esmeta,"./06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_virmetrics.txt", sep='\t', quote=F)

# virus species table:
write.table(by_sp, "./06.CLEAN_DATA/02.FINAL/by_virus_rank/virus_species_table_2220samples_120997vOTUr.txt", sep='\t', quote=F)

# virus genus table:
write.table(by_genus, "./06.CLEAN_DATA/02.FINAL/by_virus_rank/virus_genus_table_2220samples_120997vOTUr.txt", sep='\t', quote=F)

# virus family table:
write.table(by_family, "./06.CLEAN_DATA/02.FINAL/by_virus_rank/virus_family_table_2220samples_120997vOTUr.txt", sep='\t', quote=F)

# virus class table:
write.table(by_class, "./06.CLEAN_DATA/02.FINAL/by_virus_rank/virus_class_table_2220samples_120997vOTUr.txt", sep='\t', quote=F)

# summary stat for the full esmeta:
write.table(esmeta_sumstat, "./06.CLEAN_DATA/02.FINAL/Summary_stat_Chiliadal_meta_Ext_v05_suppl.txt", sep='\t', quote=F, row.names = F)

# summary stat for the full esmeta:
write.table(work_VLP_esmeta_sumstat, "./06.CLEAN_DATA/02.FINAL/Summary_stat_Chiliadal_meta_VLP_Ext_v05_suppl.txt", sep='\t', quote=F, row.names = F)

# summary stat for the working esmeta:
write.table(wETOF_sumstat, "./06.CLEAN_DATA/02.FINAL/Summary_stat_Working_ETOF_120997vOTUr.txt", sep='\t', quote=F, row.names = F)

