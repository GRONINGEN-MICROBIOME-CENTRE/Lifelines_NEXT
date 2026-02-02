setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

#############################################################
# 1. Functions
#############################################################
# saturation_stat_fast is a function for creating "plot_data" 
# necessary for the viral entities saturation plot
# plot_data is a data frame with X number of rows, where X is 
# the number of samples and 3 columns: Samples (1 - X), 
# Mean_Detected_ent, standard deviation
# input: count_table (rows - entities, columns - samples, row.names - entities names)
# output: list with plot data, summarized permutation data, stats
saturation_stat_fast <- function(count_table, n_permutations) {
  
  bin_mat <- count_table > 0
  
  # 2) Helper to get a single richness curve
  richness_curve <- function(bm) {
    cum_pres <- t(apply(bm, 1, cumsum))
    colSums(cum_pres > 0)
  }
  
  # 3) Compute original
  original <- richness_curve(bin_mat)
  
  # 4) Compute all perms at once
  permuted_matrix <- replicate(
    n_permutations,
    richness_curve(bin_mat[, sample(ncol(bin_mat))])
  )
  
  # 5) Summarize
  mean_votus <- rowMeans(permuted_matrix)
  sd_votus   <- apply(permuted_matrix, 1, sd)
  
  plot_data <- data.frame(
    Samples             = seq_along(original),
    Mean_Detected_vOTUs = mean_votus,
    SD                  = sd_votus
  )
  
  # 6) Marginal‐gain stats on the flattened diffs
  marginal_gains <- apply(permuted_matrix, 2, diff)      # matrix: (n_samples-1) × n_permutations
  all_gains      <- as.vector(marginal_gains)           # flatten to one vector
  
  stat_tab <- data.frame(
    Metric = c(names(summary(all_gains)), "SD"),
    Value  = c(as.numeric(summary(all_gains)), sd(all_gains))
  )
  
  list(
    original = data.frame(Samples = seq_along(original),
                          Detected_vOTUs = original),
    permuted = plot_data,
    stat     = stat_tab
  )
}


library(tidyverse)
library(ggplot2)
library(dplyr)
library(see)
library(UpSetR)
#############################################################
# 2. Load Input Data
#############################################################

clean_RPKM <- data.table::fread('06.CLEAN_DATA/VLP_MGS_decontaminated_RPKM_overlap.txt', sep='\t', header=T)
row.names(clean_RPKM) <- clean_RPKM$vOTUs
clean_RPKM$vOTUs <- NULL

full_overlap <- read.table('06.CLEAN_DATA/Intermediate/MGS_VLP_samples_full_overlap.txt', sep='\t', header=T)

all_overlap <- c(full_overlap$VLP, full_overlap$MGS)

setdiff(colnames(clean_RPKM), all_overlap)

smeta <- read.delim('06.CLEAN_DATA/Intermediate/Chiliadal_metadata_ver_05_25052025.txt', sep='\t', header=T)
smeta$Timepoint_new <- factor(smeta$Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)

prechili <- read.table('~/Desktop/Projects_2021/Baseclear/01.TEST_22_samples/BC_PILOT_new_ids_chiliadal_format.txt', sep='\t')
prechili <- prechili[,c('V1', 'V2')]

smeta$VLP_ID_UPD <- smeta$Sequencing_ID_VLP
smeta$VLP_ID_UPD <- prechili$V2[match(smeta$Sequencing_ID_VLP, prechili$V1)]
smeta[is.na(smeta$VLP_ID_UPD),"VLP_ID_UPD"] <- smeta[is.na(smeta$VLP_ID_UPD),"Sequencing_ID_VLP"]
unique(smeta$Timepoint_new)

ETOF <- data.table::fread('06.CLEAN_DATA/VLP_MGS_ETOF_full_rep.txt', sep='\t', header=T)
ETOF_only <- ETOF[grep('NEXT_', ETOF$New_CID),]
ETOF_only$sample <- gsub('_.*', '', ETOF_only$Original_CID)

ETOF_only$method <- ifelse(grepl('NEXT_V', ETOF_only$New_CID), 'VLP', 'MGS')

# cleanest RPKM (to filter vOTUs etc):
clean_RPKM <- as.data.frame(clean_RPKM, row.names = row.names(clean_RPKM))
cleanest <- clean_RPKM[, colnames(clean_RPKM) %in% all_overlap]
ETOF_vOTUr <- ETOF[ETOF$New_CID %in% row.names(cleanest),]
cleanest <- cleanest[row.names(cleanest) %in% ETOF_vOTUr[ETOF_vOTUr$POST_CHV_length >=3000,]$New_CID,]

cleanest <- cleanest[rowSums(cleanest) >0,]
cleanest <- cleanest[,colSums(cleanest) >0]

# clean RPKM table only for VLP samples from full overlap:
cleanest_VLP <- cleanest[,colnames(cleanest) %in% full_overlap$VLP]
cleanest_VLP <- cleanest_VLP[rowSums(cleanest_VLP) > 0,]
colnames(cleanest_VLP)[colnames(cleanest_VLP) %in% prechili$V1] <- prechili$V2[match(colnames(cleanest_VLP)[colnames(cleanest_VLP) %in% prechili$V1],
                                                                                     prechili$V1)]

# saving point for cleanest_VLP

cleanest_MGS <- cleanest[,colnames(cleanest) %in% full_overlap$MGS]
cleanest_MGS <- cleanest_MGS[rowSums(cleanest_MGS) > 0,]

# saving point for cleanest_MGS


ETOF_vOTUr <- ETOF[ETOF$New_CID %in% row.names(cleanest),]
baseETOF_vOTUr <- ETOF_vOTUr

vOTU_clustering <- read.table('06.CLEAN_DATA/NEXT_viral_clusters_MGS_VLP_long_format.txt', sep='\t', header=T)
vOTU_clustering <- vOTU_clustering[vOTU_clustering$Representative %in% ETOF_vOTUr$New_CID,]

vOTU_cluster_size <- read.table('06.CLEAN_DATA/NEXT_viral_clusters_MGS_VLP_size.txt', sep='\t', header=T)
vOTU_cluster_size <- vOTU_cluster_size[vOTU_cluster_size$Representative %in% ETOF_vOTUr$New_CID,]

ETOF_vOTUr$vOTU_size <- vOTU_cluster_size$Cluster_size[match(ETOF_vOTUr$New_CID, vOTU_cluster_size$Representative)]

lifestyle_raw <- read.table('06.CLEAN_DATA/Lifestyle_prediction_bacphlip', sep='\t', header=T)
colnames(lifestyle_raw)[1] <- 'New_CID'
ETOF_vOTUr <- merge(ETOF_vOTUr, lifestyle_raw, by="New_CID")
ETOF_vOTUr$lifestyle <- NA
ETOF_vOTUr$lifestyle <- as.character(ETOF_vOTUr$lifestyle)
ETOF_vOTUr[ETOF_vOTUr$Temperate >=0.5,]$lifestyle <- "Temperate"
ETOF_vOTUr[ETOF_vOTUr$Virulent >0.5 & ETOF_vOTUr$completeness >= 90 & !is.na(ETOF_vOTUr$completeness),]$lifestyle <- "Virulent"
ETOF_vOTUr[is.na(ETOF_vOTUr$lifestyle),]$lifestyle <- "Unknown"

newtax <- read.table('06.CLEAN_DATA/MergedTaxonomy_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

ETOF_vOTUr <- merge(ETOF_vOTUr, newtax, 
                    by.x="New_CID", by.y="Genome_ID", all.x=T)

dim(cleanest)

VLP_metadata_to_update <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLPmatched_v02.txt', sep='\t', header=T)
vlpdiv <- data.frame(vegan::diversity(t(cleanest_VLP), "shannon")) %>%
  set_names("diversity")

temperates <- ETOF_vOTUr[ETOF_vOTUr$lifestyle=="Temperate",]$New_CID

VLP_metadata_to_update$virshannon <- vlpdiv$diversity[match(VLP_metadata_to_update$Sequencing_ID, row.names(vlpdiv))]
VLP_metadata_to_update$N_temperate <- colSums(cleanest_VLP[row.names(cleanest_VLP) %in% temperates,] > 0)[match(VLP_metadata_to_update$Sequencing_ID, colnames(cleanest_VLP))]
VLP_metadata_to_update$RAb_temperate <- 100*(colSums(cleanest_VLP[row.names(cleanest_VLP) %in% temperates,])/colSums(cleanest_VLP))[match(VLP_metadata_to_update$Sequencing_ID, colnames(cleanest_VLP))]

# for Cyrus

VLP_MGS_metadata_to_UPD <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v02.txt', sep='\t', header=T)
vlpdivall <- data.frame(vegan::diversity(t(cleanest), "shannon")) %>%
  set_names("diversity")

VLP_MGS_metadata_to_UPD$virshannon <- vlpdivall$diversity[match(VLP_MGS_metadata_to_UPD$Sequencing_ID, row.names(vlpdivall))]
colnames(cleanest)[colnames(cleanest) %in% prechili$V1] <- prechili$V2[match(colnames(cleanest)[colnames(cleanest) %in% prechili$V1],
                                                                                     prechili$V1)]

VLP_MGS_metadata_to_UPD$N_temperate <- colSums(cleanest[row.names(cleanest) %in% temperates,] > 0)[match(VLP_MGS_metadata_to_UPD$Sequencing_ID, colnames(cleanest))]
VLP_MGS_metadata_to_UPD$RAb_temperate <- 100*(colSums(cleanest[row.names(cleanest) %in% temperates,])/colSums(cleanest))[match(VLP_MGS_metadata_to_UPD$Sequencing_ID, colnames(cleanest))]


#############################################################
# 2. Analysis: redundant & complete genomes per method
#############################################################

N_disc_total <- as.vector.data.frame(table(ETOF_only$method))

N_disc_HQ <- as.vector.data.frame(table(ETOF_only[ETOF_only$checkv_quality %in% c('High-quality', 'Complete'),]$method))

N_disc_HQ10 <- table(ETOF_only[ETOF_only$checkv_quality %in% c('High-quality', 'Complete') & 
                  ETOF_only$POST_CHV_length >= 10000,]$method)

MedLen_dischq <- as.vector.data.frame(c())

for (i in c("MGS", "VLP")) {
  
  SM <- summary(ETOF_only[ETOF_only$checkv_quality %in% c('High-quality', 'Complete') &
                            ETOF_only$method == i,]$POST_CHV_length)
  
  char_method <- paste0(SM[3], ", IQR: (", SM[2], ", ", SM[5], ")")
  
  MedLen_dischq[i] <- char_method
  
}

# sample-wise:

Med_disc <- as.vector.data.frame(c())
Med_disc_hq <- as.vector.data.frame(c())
Med_disc_hq10 <- as.vector.data.frame(c())

for (i in c("MGS", "VLP")) {
  
  disc <- data.frame(table(ETOF_only[ETOF_only$method==i,]$sample))
  SM <- summary(disc$Freq)
  Med_disc[i] <- paste0(SM[3], ", IQR: (", SM[2], ", ", SM[5], ")")
  
  disc_hq <- data.frame(table(ETOF_only[ETOF_only$method==i & 
                                        ETOF_only$checkv_quality %in% c('High-quality', 'Complete'),]$sample))
  SM <- summary(disc_hq$Freq)
  Med_disc_hq[i] <- paste0(SM[3], ", IQR: (", SM[2], ", ", SM[5], ")")
  
  disc_hq10 <- data.frame(table(ETOF_only[ETOF_only$method==i & 
                                          ETOF_only$checkv_quality %in% c('High-quality', 'Complete') &
                                          ETOF_only$POST_CHV_length >= 10000,]$sample))
  SM <- summary(disc_hq10$Freq)
  Med_disc_hq10[i] <- paste0(SM[3], ", IQR: (", SM[2], ", ", SM[5], ")")
  
  if (i=="MGS") {
    smeta[, paste0("disc_", i)] <- disc$Freq[match(smeta$NG_ID, disc$Var1)]
    
    smeta[, paste0("disc_hq_", i)] <- disc_hq$Freq[match(smeta$NG_ID, disc_hq$Var1)]
    
    smeta[, paste0("disc_hq10_", i)] <- disc_hq10$Freq[match(smeta$NG_ID, disc_hq10$Var1)]
    
  } else {
    smeta[, paste0("disc_", i)] <- disc$Freq[match(smeta$VLP_ID_UPD, disc$Var1)]
    smeta[, paste0("disc_hq_", i)] <- disc_hq$Freq[match(smeta$VLP_ID_UPD, disc_hq$Var1)]
    smeta[, paste0("disc_hq10_", i)] <- disc_hq10$Freq[match(smeta$VLP_ID_UPD, disc_hq10$Var1)]
  }

  
}

per_m <- rbind(N_disc_total, 
               N_disc_HQ, 
               N_disc_HQ10, 
               MedLen_dischq,
               Med_disc,
               Med_disc_hq,
               Med_disc_hq10)

#############################################################
# 2. Analysis: Unique & novel genomes recovery
#############################################################

#ETOF <- ETOF[ETOF$New_CID %in% vOTU_clustering$Cluster_member,] #613,189
table(ETOF_vOTUr$miuvig_quality) # 16,208 complete or nearly complete (>90%)
min(ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality",]$POST_CHV_length)
max(ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality",]$POST_CHV_length)

vOTU_clustering$DB_member <- NA
vOTU_clustering[grep('NEXT_V', vOTU_clustering$Cluster_member),]$DB_member <- 'NEXT_VLP'
vOTU_clustering[grep('NEXT_M', vOTU_clustering$Cluster_member),]$DB_member <- 'NEXT_MGS'
vOTU_clustering[grep('Guerin', vOTU_clustering$Cluster_member),]$DB_member <- 'Guerin'
vOTU_clustering[grep('Yutin', vOTU_clustering$Cluster_member),]$DB_member <- 'Yutin'
vOTU_clustering[grep('NCBI_CrAss', vOTU_clustering$Cluster_member),]$DB_member <- 'NCBI_CrAss'
vOTU_clustering[grep('NL_crAss', vOTU_clustering$Cluster_member),]$DB_member <- 'NL_crAss'
vOTU_clustering[grep('Benler', vOTU_clustering$Cluster_member),]$DB_member <- 'Benler'
vOTU_clustering[grep('COPSAC', vOTU_clustering$Cluster_member),]$DB_member <- 'COPSAC'
vOTU_clustering[grep('GVD', vOTU_clustering$Cluster_member),]$DB_member <- 'GVD'
vOTU_clustering[grep('IMGVR', vOTU_clustering$Cluster_member),]$DB_member <- 'IMGVR'
vOTU_clustering[grep('MGV-', vOTU_clustering$Cluster_member),]$DB_member <- 'MGV'
vOTU_clustering[grep('GPD', vOTU_clustering$Cluster_member),]$DB_member <- 'GPD'
vOTU_clustering[grep('VREF', vOTU_clustering$Cluster_member),]$DB_member <- 'VREF'

vOTU_clustering$DB_member <- factor(vOTU_clustering$DB_member)

vOTU_by_source <- vOTU_clustering %>%
  mutate(DB_member = factor(DB_member)) %>%
  count(Representative, DB_member) %>%
  pivot_wider(names_from = DB_member, values_from = n, values_fill = 0)

merged <- vOTU_by_source %>%
  
  left_join(vOTU_cluster_size, by = "Representative") %>%
  
  rename(N_genomes = Cluster_size, vOTU_representative  = Representative) %>%
  
  #left_join(genera, by = c("vOTU_representative" = "Cluster_member")) %>%
  
  #rename(Genus = Representative) %>%
  
  left_join(ETOF %>% select(New_CID, POST_CHV_length, miuvig_quality),
            by = c("vOTU_representative" = "New_CID")) #%>%
  
  #left_join(families,  by = c("vOTU_representative" = "Cluster_member")) %>%
  
  #rename(Family = Representative) %>%
  
  # mutate(
  #   vOTU_cluster_type = if_else(N_genomes == NEXT, "NEXT", "Mixed"),
  #   miuvig_quality = factor(miuvig_quality,
  #                           levels = c("High-quality", "Genome-fragment"),
  #                           ordered = TRUE)
  
merged$DB_sum <- merged$N_genomes - merged$NEXT_MGS - merged$NEXT_VLP

sources <- c("NEXT_VLP","NEXT_MGS","DB_sum")

merged <- merged %>%
  mutate(
    vOTU_cluster_type = apply(across(all_of(sources)) > 0, 1,
                         function(r) { x <- sources[r]; if (length(x)==0) "None" else paste(x, collapse = "+") })
  )

merged_HQ <- merged[merged$miuvig_quality=="High-quality",] %>%
  mutate(
    vOTU_cluster_type = apply(across(all_of(sources)) > 0, 1,
                              function(r) { x <- sources[r]; if (length(x)==0) "None" else paste(x, collapse = "+") })
  )

all <- table(merged$vOTU_cluster_type)
hq_only <- table(merged_HQ$vOTU_cluster_type)

vOTUstat <- cbind(all, hq_only)

merged$vOTUr_source <- NA
merged[grep('NEXT_V', merged$vOTU_representative),]$vOTUr_source <- 'NEXT_VLP'
merged[grep('NEXT_M', merged$vOTU_representative),]$vOTUr_source <- 'NEXT_MGS'
merged[is.na(merged$vOTUr_source),"vOTUr_source"] <- "DB"

##### adding info to ETOF_vOTUr:
ETOF_vOTUr <- merge(ETOF_vOTUr, merged[,c("vOTU_representative", "vOTU_cluster_type", "vOTUr_source")], by.x="New_CID", by.y="vOTU_representative", all = T)

# UpSet plot (all vOTUs)
listInput <- list(DB=merged[merged$DB_sum>0,]$vOTU_representative,
                  MGS=merged[merged$NEXT_MGS>0,]$vOTU_representative,
                  VLP=merged[merged$NEXT_VLP>0,]$vOTU_representative)


upset_all <- upset(fromList(listInput), order.by = "freq", sets.bar.color = "#C00000", 
      number.angles = 20,
      sets.x.label = "N vOTUs", scale.sets = "identity",
      text.scale = c(1, 1, 1, 0.7, 1, 1))
 
png('05.PLOTS/05.VLP_MGS/UpSet_plot_all_vOTUs.png', width=10, height=7, units="cm", res = 300)
upset_all
dev.off()

# UpSet plot (HQ vOTUs)
listInputHQ <- list(DB=merged_HQ[merged_HQ$DB_sum>0,]$vOTU_representative,
                  MGS=merged_HQ[merged_HQ$NEXT_MGS>0,]$vOTU_representative,
                  VLP=merged_HQ[merged_HQ$NEXT_VLP>0,]$vOTU_representative)


upset_allhq <- upset(fromList(listInputHQ), order.by = "freq", sets.bar.color = "#C00000", 
      number.angles = 20,
      sets.x.label = "N vOTUs", scale.sets = "identity",
      text.scale = c(1, 1, 1, 0.7, 1, 1))

png('05.PLOTS/05.VLP_MGS/UpSet_plot_HQ_vOTUs.png', width=10, height=7, units="cm", res = 300)
upset_allhq
dev.off()

### CALCUATING TRUE NOVELS taking into account the QoL:
# derive the Quality of the Longest (QoL) contig per contributing source per vOTUr,
# to truly understand which source is indespensible
vOTU_clustering$length <- ETOF$POST_CHV_length[match(vOTU_clustering$Cluster_member, ETOF$New_CID)]
vOTU_clustering$quality <- ETOF$miuvig_quality[match(vOTU_clustering$Cluster_member, ETOF$New_CID)]

summary_df <- vOTU_clustering %>%
  filter(Representative %in% hq_votus) %>%
  mutate(DB_member = if_else(DB_member %in% c("NEXT_MGS", "NEXT_VLP"), DB_member, "DB")) %>%
  arrange(desc(length), desc(quality)) %>%
  group_by(Representative, DB_member) %>%
  slice_head(n = 1) %>%    ungroup() %>%
  select(Representative, DB_member, quality) %>%
  pivot_wider(
    names_from = DB_member,
    values_from = quality,
    names_prefix = ""
  ) %>%
  rename(DB_QoL = DB, MGS_QoL=NEXT_MGS, VLP_QoL=NEXT_VLP) %>%
  replace_na(list(
    DB_QoL = "Absent",
    MGS_QoL = "Absent",
    VLP_QoL = "Absent"
  ))

summary_df <- merge(summary_df, ETOF_vOTUr[,c("New_CID", "vOTU_cluster_type", "vOTUr_source")],
                    by.x="Representative", by.y="New_CID")

summary_table <- summary_df %>%
  group_by(vOTU_cluster_type) %>%
  summarise(
    N_HQ = n(),  # Total number of vOTUs in this cluster type
    HQ_MGS = sum(MGS_QoL == "High-quality"),
    HQ_VLP = sum(VLP_QoL == "High-quality"),
    HQ_DB  = sum(DB_QoL  == "High-quality"),
    HQ_MGS_VLP = sum(MGS_QoL == "High-quality" & VLP_QoL == "High-quality"),
    HQ_DB_VLP = sum(DB_QoL == "High-quality" & VLP_QoL == "High-quality"),
    HQ_DB_MGS = sum(DB_QoL == "High-quality" & MGS_QoL == "High-quality"),
    HQ_DB_MGS_VLP = sum(DB_QoL == "High-quality" & MGS_QoL == "High-quality" & VLP_QoL == "High-quality")
  ) %>%
  arrange(desc(N_HQ))


# total VLP HQ:
3374 + 4045 + 3399 # hq but derep w MGS & DB + unique + hq but derep w MGS

# prop novel then:
(4045 + 3399)/(3374 + 4045 + 3399) #68.8


# 86.8% of all novels would have been recovered by VLPs alone 

# total MGS HQ:
2766 + 439 + 1846 # hq but derep w VLP & DB + unique + hq but derep w VLP
(439 + 1846)/(2766 + 439 + 1846) # 45.2

# 26.4% of novels would have been recovered by MGS alone

# does every VLP sample just carries more HQs to start with? (judge based on ETOF)

# also solves problem with saturation curve -> use this method to recalculate the saturation curve p method
# but not sure what to do about the stacked bar plot

#############################################################
# 2. Analysis: testing how frequent it is for a genome-fragment
# quality vOTU to have high quality members
#############################################################
# vOTU_clustering$quality <- ETOF$miuvig_quality[match(vOTU_clustering$Cluster_member, 
#                                                      ETOF$New_CID)]
# vOTU_clustering$length <- ETOF$POST_CHV_length[match(vOTU_clustering$Cluster_member, 
#                                                      ETOF$New_CID)]
# 
# cluster_summary <- vOTU_clustering %>%
#   group_by(Representative) %>%  
#   summarise(
#     N_high = sum(quality == "High-quality"),
#     N_low  = sum(quality == "Genome-fragment")
#   )
# 
# cluster_summary$rmiuvig_quality <- ETOF_vOTUr$miuvig_quality[match(df_summary$Representative, 
#                                                                    ETOF_vOTUr$New_CID)]
# 
# dim(df_summary[df_summary$N_high > 0 & df_summary$rmiuvig_quality=="Genome-fragment",])
# there are only 248 cases like that (1.5%), some of them are caused by CheckV erroneously 
# considering some short genomes HQ (as w vOTU==NEXT_V0156_N33_L31446_K2.8_E0_P0_F0)
# DECISION: IGNORE IT
#############################################################
# 2. Analysis: every sample contributing to novel discovery
#############################################################
vOTUs_uniq_to_vlp <- merged_HQ[merged_HQ$vOTU_cluster_type %in% c("NEXT_VLP", "NEXT_VLP+NEXT_MGS"),]$vOTU_representative
vOTUs_uniq_to_mgs <- merged_HQ[merged_HQ$vOTU_cluster_type %in% c("NEXT_MGS", "NEXT_VLP+NEXT_MGS"),]$vOTU_representative

VLP_hq_uniq <- cleanest[row.names(cleanest) %in% vOTUs_uniq_to_vlp, colnames(cleanest) %in% full_overlap$VLP]
MGS_hq_uniq <- cleanest[row.names(cleanest) %in% vOTUs_uniq_to_mgs, colnames(cleanest) %in% full_overlap$MGS]

VLP_vOTU_HQ_new_cumulative <- saturation_stat_fast(VLP_hq_uniq, 100)
MGS_vOTU_HQ_new_cumulative <- saturation_stat_fast(MGS_hq_uniq, 100)

tmp <- VLP_vOTU_HQ_new_cumulative[["permuted"]]
tmp$Method <- "VLP"

tmp2 <- MGS_vOTU_HQ_new_cumulative[["permuted"]]
tmp2$Method <- "MGS"

TMP <- rbind(tmp, tmp2)

novel_satu <- ggplot(TMP, aes(x = Samples, y = Mean_Detected_vOTUs, color=Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = Mean_Detected_vOTUs - SD, ymax = Mean_Detected_vOTUs + SD), width = 0.1) +
  labs(title = "N novel vOTUs discovered with increasing number of samples",
       x = "Number of fecal samples", y = "N vOTUs") +
  theme_minimal()

ggsave('05.PLOTS/05.VLP_MGS/Novel_detected_saturation.png',
       novel_satu,  "png", width=16, height=12, units="cm", dpi = 300)
#############################################################
# 2. Analysis: sample-based method assessment
#############################################################
smeta$clean_richness_VLP <- colSums(cleanest > 0)[match(smeta$Sequencing_ID_VLP, colnames(cleanest))]
smeta$clean_richness_MGS <- colSums(cleanest > 0)[match(smeta$NG_ID, colnames(cleanest))]

wilcox.test(smeta$clean_richness_VLP, smeta$clean_richness_MGS, paired=T)


p_VLP_MGS_richness <- smeta %>%
  filter(Universal_ID %in% full_overlap$Universal_ID) %>%
  select(Universal_ID, clean_richness_MGS, clean_richness_VLP) %>%
  pivot_longer(!Universal_ID, 
               names_to = "Source",
               values_to = "Richness") %>%
  ggplot(aes(Source, Richness)) +
  geom_violinhalf(aes(Source, Richness, fill = Source, color=Source), flip=1) +
  scale_color_manual(values=c("darkred", "#234C6A")) +
  ggnewscale::new_scale_color() +
  geom_boxplot(aes(Source, Richness, fill = Source, color=Source),width = 0.1, outlier.shape = NA) +
  scale_color_manual(values=c("darkred", "#234C6A")) +
  geom_rect(xmin=1, xmax=1.1, ymin=10, ymax=2000, fill="white") +
  geom_rect(xmin=1.9, xmax=2, ymin=10, ymax=2000, fill="white") +
  geom_line(aes(group = Universal_ID),alpha=0.1, color="darkgrey") +
  geom_point(aes(Source, Richness, color=Source), size=1.5, alpha=0.3) + 
  scale_x_discrete(labels = c("MGS", "VLP")) +
  theme_bw() +
  theme(legend.position = "none")

ggsave('05.PLOTS/05.VLP_MGS/Overall_richness_comparison.png',
       p_VLP_MGS_richness,  "png", width=10, height=10, units="cm", dpi = 300)

# For HQs:

hq_votus <- ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality",]$New_CID
smeta$clean_richness_VLPHQ <- colSums(cleanest[row.names(cleanest) %in% hq_votus,] > 0)[match(smeta$Sequencing_ID_VLP, colnames(cleanest))]
smeta$clean_richness_MGSHQ <- colSums(cleanest[row.names(cleanest) %in% hq_votus,] > 0)[match(smeta$NG_ID, colnames(cleanest))]

p_VLP_MGS_richnessHQ <- smeta %>%
  filter(Universal_ID %in% full_overlap$Universal_ID) %>%
  select(Universal_ID, clean_richness_MGSHQ, clean_richness_VLPHQ) %>%
  pivot_longer(!Universal_ID, 
               names_to = "Source",
               values_to = "Richness") %>%
  ggplot(aes(Source, Richness)) +
  geom_violinhalf(aes(Source, Richness, fill = Source, color=Source), flip=1) +
  scale_color_manual(values=c("darkred", "#234C6A")) +
  ggnewscale::new_scale_color() +
  geom_boxplot(aes(Source, Richness, fill = Source, color=Source),width = 0.1, outlier.shape = NA) +
  scale_color_manual(values=c("darkred", "#234C6A")) +
  geom_rect(xmin=1, xmax=1.1, ymin=10, ymax=50, fill="white") +
  geom_rect(xmin=1.9, xmax=2, ymin=10, ymax=100, fill="white") +
  geom_line(aes(group = Universal_ID),alpha=0.1, color="darkgrey") +
  geom_point(aes(Source, Richness, color=Source), size=1.5, alpha=0.3) + 
  labs(y="Richness of HQ vOTUs") +
  scale_x_discrete(labels = c("MGS", "VLP")) +
  theme_bw() +
  theme(legend.position = "none")

ggsave('05.PLOTS/05.VLP_MGS/HQ_richness_comparison.png',
       p_VLP_MGS_richnessHQ,  "png", width=10, height=10, units="cm", dpi = 300)

Richness_MGS_VLP <- cor.test(smeta$clean_richness_MGS, smeta$clean_richness_VLP, method = "spearman")
Richness_MGS_VLP$p.value

# explanation why HQs are higher in VLP vs MGS:
hq_votus_ssDNA <- ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality" &
                               ETOF_vOTUr$genome_simple=="ssDNA",]$New_CID
smeta$clean_richness_VLPHQss <- colSums(cleanest[row.names(cleanest) %in% hq_votus_ssDNA,] > 0)[match(smeta$Sequencing_ID_VLP, colnames(cleanest))]
smeta$clean_richness_MGSHQss <- colSums(cleanest[row.names(cleanest) %in% hq_votus_ssDNA,] > 0)[match(smeta$NG_ID, colnames(cleanest))]

HQ_VLP_VLPssDNA <- cor.test(smeta$clean_richness_VLPHQ, smeta$clean_richness_VLPHQss, method = "spearman")
HQ_VLP_VLPssDNA$p.value

table(ETOF_vOTUr[ETOF_vOTUr$vOTU_cluster_type=="NEXT_VLP" &
                   ETOF_vOTUr$miuvig_quality=="High-quality" &
                   ETOF_vOTUr$genome_simple %in% c('ssDNA', 'RNA'), genome_simple]) #123 RNA, 2825 ssDNA


#############################################################
# 3. Intermediate output: contig IDs for further analysis
#############################################################
write.table(row.names(cleanest), "./06.CLEAN_DATA/VLP_MGS_vOTUr_dec99ANI_ab3kbp_2220samples.txt", sep='\t', row.names=F, col.names = F, quote=F)
write.table(cleanest, "./06.CLEAN_DATA/RPKM_table_VLP_MGS_dec99ANI_ab3kbp_2220_samples.txt", sep='\t', quote=F)
write.table(baseETOF_vOTUr, "./06.CLEAN_DATA/ETOF_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt", sep='\t', quote=F, row.names=F)
write.table(ETOF_vOTUr, "./06.CLEAN_DATA/Extended_ETOF_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt", sep='\t', quote=F)


# clean VLP RPKM table based on all 127553 vOTUs discovered in MGS-VLP full overlap samples:
write.table(cleanest_VLP, './06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ab3kbp_1110_samples.txt', sep='\t', quote=F)

# Updated VLP only metadata now containing temperate phage relative abundance and richness of temperate phages
write.table(VLP_metadata_to_update, '06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLPmatched_v02.1.txt', sep='\t', quote=F, row.names=F)

# Updated VLP-MGS metadata now containing temperate phage relative abundance and richness of temperate phages & virshannin
write.table(VLP_MGS_metadata_to_UPD, '06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v02.1.txt', sep='\t', quote=F, row.names=F)

# clean MGS RPKM table based on all 127553 vOTUs discovered in MGS-VLP full overlap samples:
write.table(cleanest_MGS, './06.CLEAN_DATA/02.FINAL/MGS_only_RPKM_table_VLP_MGS_dec99ANI_ab3kbp_1110_samples.txt', sep='\t', quote=F)

# table with per-contig HQ vOTU recovery info by source
write.table(summary_df, './06.CLEAN_DATA/03.RESULTS_to_plot/HQ_vOTU_recoverable_by_source.txt', sep='\t', quote=F, row.names=F)

# summary stats of HQ vOTU recovery across sources
write.table(summary_table, './06.CLEAN_DATA/03.RESULTS_to_plot/HQ_vOTU_recovery_stat.txt', sep='\t', quote=F, row.names=F)

