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

# get iqrs:
get_iqr <- function(vec) {
  
  SMR <- summary(vec)
  
  print(paste0(round(SMR[3], 2), ' (', round(SMR[2], 2), ' - ', round(SMR[5], 2), ')')) 
  
}

# gets the summaries of different curve fittings
podgonian <- function(DF, x, y){
  
  logger <- list()
  
  F1 <- as.formula(paste0("log(", y, ") ~ log(", x, ")"))
  
  fit_power <- lm(F1, data = DF)
  logger[["lm_fit_power"]] <- summary(fit_power)
  
  b <- unname(fit_power$coefficients[2])
  median_index <- round(median(DF[,x]))
  a <- sqrt(DF[median_index,y])
  
  # power function
  F2 <- as.formula(paste(y, " ~ a * ", x, "^b"))
  
  fit_root <- nls(F2, 
                  data = DF, 
                  start = list(a=a,b=b))
  logger[["fit_root"]] <- summary(fit_root)
  
  root_params <- coef(fit_root)
  MDR <- root_params["a"] * root_params["b"] * (nrow(DF)^(root_params["b"]-1))
  logger[["MDR"]] <- MDR
  
  # plateau:
  F3 <- as.formula(paste(y, " ~ ", x))
  fit_plat <- drc::drm(F3, # drc is called here in the isolated way, otherwise it masks "select" from dplyr
                       data = DF, 
                       fct = drc::L.4())
  
  logger[["fit_plat"]] <- summary(fit_plat)
  
  logger[["plat_vs_root"]] <- AIC(fit_plat, fit_root)
  
  # asymptotic model
  F4 <- as.formula( paste(y, " ~ SSasymp(", x, ", Asym, R0, lrc)"))
  
  fit_asym <- nls(F4, 
                  data = DF)
  
  logger[["fit_asym"]] <- summary(fit_asym)
  
  logger[["asym_vs_root"]] <- AIC(fit_asym, fit_root)
  logger[["asym_vs_plat"]] <- AIC(fit_asym, fit_plat)
  
  
  params <- coef(fit_asym)
  Asym <- params["Asym"]
  R0   <- params["R0"]
  lrc  <- params["lrc"]
  
  rate <- exp(lrc)
  # first derivative
  slope <- (Asym - R0) * rate * exp(-rate * nrow(DF))
  
  logger[["slope"]] <- slope
  
  return(logger)
  
  
}

smart_round <- function(x){
  
  ifelse(x < 0.01,
        as.numeric(formatC(x, format = "e", digits = 2)),
                round(x, 2))

}

library(tidyverse)
library(ggplot2)
library(dplyr)
library(see)
library(UpSetR)
library(lme4)
library(lmerTest)
library(purrr)
library(broom)
#library(drc) # for fitting curves
#############################################################
# 2. Load Input Data
#############################################################

clean_RPKM <- read.table('06.CLEAN_DATA/02.FINAL/RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_2220_samples.txt', sep='\t', header=T)

smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)
smeta$Timepoint_new <- factor(smeta$Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)

ETOF <- read.table('06.CLEAN_DATA/VLP_MGS_ETOF_full_rep.txt', sep='\t', header = T)

ETOF_only <- ETOF %>%
  filter(grepl('NEXT_', ETOF$New_CID)) %>%
  mutate(sample = gsub('_.*', '', Original_CID)) %>%
  mutate( method = ifelse(grepl('NEXT_V', New_CID), 'VLP', 'MGS') )

ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

vOTU_clustering <- read.table('06.CLEAN_DATA/NEXT_viral_clusters_MGS_VLP_long_format.txt', sep='\t', header=T)
vOTU_clustering <- vOTU_clustering[vOTU_clustering$Representative %in% ETOF_vOTUr$New_CID,]

vOTU_cluster_size <- read.table('06.CLEAN_DATA/NEXT_viral_clusters_MGS_VLP_size.txt', sep='\t', header=T)
vOTU_cluster_size <- vOTU_cluster_size[vOTU_cluster_size$Representative %in% ETOF_vOTUr$New_CID,]

VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)

MGS <- read.table('06.CLEAN_DATA/02.FINAL/MGS_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)
#############################################################
# 3.1 Analysis: redundant & hq genomes per method
#############################################################

N_disc_total <- as.vector.data.frame(table(ETOF_only$method))

N_disc_HQ <- as.vector.data.frame(table(ETOF_only[ETOF_only$checkv_quality %in% c('High-quality', 'Complete'),]$method))

N_disc_HQ10 <- table(ETOF_only[ETOF_only$checkv_quality %in% c('High-quality', 'Complete') & 
                  ETOF_only$POST_CHV_length >= 10000,]$method)

disc <- data.frame(table(ETOF_only$sample))

disc_hq <- data.frame(table(ETOF_only$sample[ETOF_only$checkv_quality %in% c('High-quality', 'Complete')]))

disc_hq10 <- data.frame(table(ETOF_only[ETOF_only$checkv_quality %in% c('High-quality', 'Complete') &
                                          ETOF_only$POST_CHV_length >= 10000,]$sample))

smeta[, "N_discovered"] <- disc$Freq[match(smeta$Sequencing_ID, disc$Var1)]

smeta[, "N_discovered_HQ"] <- disc_hq$Freq[match(smeta$Sequencing_ID, disc_hq$Var1)]

smeta[, "N_discovered_HQ10"] <- disc_hq10$Freq[match(smeta$Sequencing_ID, disc_hq10$Var1)]

# sample-wise:
MedLen_dischq <- as.vector.data.frame(c())
Med_disc <- as.vector.data.frame(c())
Med_disc_hq <- as.vector.data.frame(c())
Med_disc_hq10 <- as.vector.data.frame(c())

for (i in c("MGS", "VLP")) {
  
  MedLen_dischq[i] <- get_iqr(ETOF_only[ETOF_only$checkv_quality %in% c('High-quality', 'Complete') &
                                          ETOF_only$method == i,]$POST_CHV_length)
  
  Med_disc[i] <- get_iqr(smeta$N_discovered[smeta$seq_type == i])
  
  Med_disc_hq[i] <- get_iqr(smeta$N_discovered_HQ[smeta$seq_type == i])
  
  Med_disc_hq10[i] <- get_iqr(smeta$N_discovered_HQ10[smeta$seq_type == i])
  
}

per_m <- rbind(N_disc_total, 
               N_disc_HQ, 
               N_disc_HQ10, 
               MedLen_dischq,
               Med_disc,
               Med_disc_hq,
               Med_disc_hq10)

# tidying up
rm(list=setdiff(ls(), c("clean_RPKM", "ETOF", "ETOF_only", "ETOF_vOTUr", "per_m", "smeta", "get_iqr", "saturation_stat_fast", "VLP", "MGS")))
#############################################################
# 3.2 Analysis: N detected vOTUs saturation curves
#############################################################
# since there is a lot of interdependency (paired samples, 
# longitudinal samples) no sense to make saturation over 2,220
# -> merging to holovirome table

##### creating aholovirome binary table
colnames(VLP) <- smeta$Universal_ID[match(colnames(VLP), smeta$Sequencing_ID)]

colnames(MGS) <- smeta$Universal_ID[match(colnames(MGS), smeta$Sequencing_ID)]

## ordering
MGS <- MGS[, colnames(VLP)]  # if needed

all_votus <- ETOF_vOTUr$New_CID

# initiate matrix
m1 <- matrix(0,
             nrow = length(all_votus),
             ncol = ncol(VLP),
             dimnames = list(all_votus, colnames(VLP)))
m2 <- m1

# populating matrices
m1[rownames(VLP), ] <- (as.matrix(VLP) > 0) * 1
m2[rownames(MGS), ] <- (as.matrix(MGS) > 0) * 1

## unified presence/absence
unified <- ((m1 + m2) > 0) * 1

holo <- as.data.frame(unified)
rm(list = c("m1", "m2", "all_votus", "unified"))

gen_holo <- holo %>%
  rownames_to_column("New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, Genus_OTUr)) %>% # library(drc) masks "select", wtf!  
  group_by(Genus_OTUr) %>%
  summarise(across(where(is.numeric), sum)) %>% # 27,017 genera detected
  column_to_rownames("Genus_OTUr")

fam_holo <- holo %>%
  rownames_to_column("New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, Family_OTUr)) %>% # library(drc) masks "select", wtf!  
  group_by(Family_OTUr) %>%
  summarise(across(where(is.numeric), sum)) %>% # 7,596 families detected
  column_to_rownames("Family_OTUr")

all_vOTU_cumulative <- saturation_stat_fast(holo, 1000)
all_gOTU_cumulative <- saturation_stat_fast(gen_holo, 1000)
all_fOTU_cumulative <- saturation_stat_fast(fam_holo, 1000)

all_cumulative <- all_vOTU_cumulative[["permuted"]] %>%
  mutate(Rank = "Species") %>%
  bind_rows(all_gOTU_cumulative[["permuted"]] %>% mutate(Rank = "Genus")) %>%
  bind_rows(all_fOTU_cumulative[["permuted"]] %>% mutate(Rank = "Family"))

all_satu <- ggplot(all_cumulative, aes(x = Samples, y = Mean_Detected_vOTUs, color=Rank)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = Mean_Detected_vOTUs - SD, ymax = Mean_Detected_vOTUs + SD), width = 0.1) +
  labs(x = "Number of fecal samples", y = "Number of OTUs") +
  scale_color_manual(values = c("#7C3E66", "#A5BECC", "#243A73")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave('05.PLOTS/05.VLP_MGS/all_OTUs_detected_holovirome.png',
       all_satu,  "png", width=10, height=8, units="cm", dpi = 300)

# looking for approximation:
# Species looks like something * sqrt()

all_cumulative %>%
  filter(Rank == "Species") %>%
  ggplot(aes(log(Samples), log(Mean_Detected_vOTUs))) +
  geom_point() +
  theme_minimal() +
  labs(y="log(mean detected species)") # ehm kind of straight line?

sp_stat <- podgonian(all_cumulative[all_cumulative$Rank == "Species",],
                  "Samples",
                  "Mean_Detected_vOTUs") # best is root square

gn_stat <- podgonian(all_cumulative[all_cumulative$Rank == "Genus",],
                     "Samples",
                     "Mean_Detected_vOTUs") # best is root square, but gain is only 2 genera ps -> plateau

fm_stat <- podgonian(all_cumulative[all_cumulative$Rank == "Family",],
                     "Samples",
                     "Mean_Detected_vOTUs") # best is root square, but gain is only 1 (0.5) ps -> def plateau

#############################################################
# 3.3 Analysis: descriptive stat NEXT virome catalog
#############################################################
# % genomes, % temperate etc was derived from summary stat of ETOF
virclass <- ETOF_vOTUr  %>%
  separate(
    col    = tax_ictv_aai,
    into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  group_by(Class) %>%
  summarise(sum=n()) %>%
  mutate(perc = round(sum / nrow(ETOF_vOTUr) * 100, 1))

#############################################################
# 3.4 Analysis: staton on MGS vs VLP genomics
#############################################################
genomics_compare <- c('raw_reads', 'human_reads', 'clean_reads',
                      'contigs_0_bp', 'contigs_1000_bp', 'sc_enrichment')

results_genomics <- map_dfr(genomics_compare, function(genomics) {
    
    print(paste("Data:", genomics))
  
  f1 <- paste("log10(", genomics, ")", "~ seq_type")
    
    formula <- as.formula(paste0(f1, "+ (1 | NEXT_ID)"))
    
    model <- lmer(
      formula,
      REML = FALSE,
      data = smeta
    )
   
    eff_res <- effsize::cohen.d(as.formula(f1), data = smeta) 
    
    model_summary <- summary(model)$coefficients %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      filter(rowname != "(Intercept)") %>%
      mutate(Cohens_D = abs(eff_res$estimate)) %>%
      mutate(rowname = gsub('seq_type', '', rowname)) %>%
      rename(`Metavirome type` = rowname) %>%
      mutate(`Genomic metric` = genomics) %>%
      mutate(across(where(is.numeric), ~ smart_round(.)) ) %>%
      relocate(`Genomic metric`)
    
    model_summary
    
  }) 

write.table(results_genomics, '07.RESULTS/Compare_genomic_metrics_MGS_VLP.txt', sep='\t', quote=F, row.names=F)
#############################################################
# 3.3 Analysis: N novel discovered vOTUs saturation curves
#############################################################
# vOTUs_uniq_to_vlp <- ETOF_vOTUr %>%
#   filter(vOTU_cluster_type %in% c("NEXT_VLP", "NEXT_VLP+NEXT_MGS", "")) %>%
#   pull(New_CID)
# 
# vOTUs_uniq_to_mgs <- ETOF_vOTUr %>%
#   filter(vOTU_cluster_type %in% c("NEXT_MGS", "NEXT_VLP+NEXT_MGS")) %>%
#   pull(New_CID)
# 
# VLP_uniq <- clean_RPKM[, colnames(clean_RPKM) %in% smeta$Sequencing_ID[smeta$seq_type == "VLP"]]
# VLP_uniq <- VLP_uniq [rowSums(VLP_uniq) >0,]
# MGS_uniq <- clean_RPKM[row.names(clean_RPKM) %in% vOTUs_uniq_to_mgs, colnames(clean_RPKM) %in% smeta$Sequencing_ID[smeta$seq_type == "MGS"]]
# MGS_uniq <- MGS_uniq [rowSums(MGS_uniq) >0,]
# 
# VLP_vOTU_HQ_new_cumulative <- saturation_stat_fast(VLP_uniq, 100)
# MGS_vOTU_HQ_new_cumulative <- saturation_stat_fast(MGS_uniq, 100)
# 
# novel_satu <- ggplot(TMP, aes(x = Samples, y = Mean_Detected_vOTUs, color=Method)) +
#   geom_line() +
#   geom_point() +
#   geom_errorbar(aes(ymin = Mean_Detected_vOTUs - SD, ymax = Mean_Detected_vOTUs + SD), width = 0.1) +
#   labs(title = "N vOTUs detected with increasing number of samples",
#        x = "Number of fecal samples", y = "N vOTUs") +
#   theme_minimal()
# 
# ggsave('05.PLOTS/05.VLP_MGS/vOTUs_detected_MGS_VLP.png',
#        novel_satu,  "png", width=16, height=12, units="cm", dpi = 300)

#############################################################
# 2. Analysis: Unique & novel genomes recovery
#############################################################

table(ETOF_vOTUr$miuvig_quality) # 15,969 complete or nearly complete (>90%)
min(ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality",]$POST_CHV_length) #3,000
max(ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality",]$POST_CHV_length) #396,781 

# vOTU_clustering$DB_member <- NA
# vOTU_clustering[grep('NEXT_V', vOTU_clustering$Cluster_member),]$DB_member <- 'NEXT_VLP'
# vOTU_clustering[grep('NEXT_M', vOTU_clustering$Cluster_member),]$DB_member <- 'NEXT_MGS'
# vOTU_clustering[grep('Guerin', vOTU_clustering$Cluster_member),]$DB_member <- 'Guerin'
# vOTU_clustering[grep('Yutin', vOTU_clustering$Cluster_member),]$DB_member <- 'Yutin'
# vOTU_clustering[grep('NCBI_CrAss', vOTU_clustering$Cluster_member),]$DB_member <- 'NCBI_CrAss'
# vOTU_clustering[grep('NL_crAss', vOTU_clustering$Cluster_member),]$DB_member <- 'NL_crAss'
# vOTU_clustering[grep('Benler', vOTU_clustering$Cluster_member),]$DB_member <- 'Benler'
# vOTU_clustering[grep('COPSAC', vOTU_clustering$Cluster_member),]$DB_member <- 'COPSAC'
# vOTU_clustering[grep('GVD', vOTU_clustering$Cluster_member),]$DB_member <- 'GVD'
# vOTU_clustering[grep('IMGVR', vOTU_clustering$Cluster_member),]$DB_member <- 'IMGVR'
# vOTU_clustering[grep('MGV-', vOTU_clustering$Cluster_member),]$DB_member <- 'MGV'
# vOTU_clustering[grep('GPD', vOTU_clustering$Cluster_member),]$DB_member <- 'GPD'
# vOTU_clustering[grep('VREF', vOTU_clustering$Cluster_member),]$DB_member <- 'VREF'

#vOTU_clustering$DB_member <- factor(vOTU_clustering$DB_member)

#vOTU_by_source <- vOTU_clustering %>%
#  mutate(DB_member = factor(DB_member)) %>%
#  count(Representative, DB_member) %>%
#  pivot_wider(names_from = DB_member, values_from = n, values_fill = 0)

#merged <- vOTU_by_source %>%
  
#  left_join(vOTU_cluster_size, by = "Representative") %>%
  
#  rename(N_genomes = Cluster_size, vOTU_representative  = Representative) %>%
  
 # #left_join(genera, by = c("vOTU_representative" = "Cluster_member")) %>%
  
#  #rename(Genus = Representative) %>%
  
 # left_join(ETOF %>% select(New_CID, POST_CHV_length, miuvig_quality),
 #           by = c("vOTU_representative" = "New_CID")) #%>%
  
#  #left_join(families,  by = c("vOTU_representative" = "Cluster_member")) %>%
  
#  #rename(Family = Representative) %>%
  
#  # mutate(
#  #   vOTU_cluster_type = if_else(N_genomes == NEXT, "NEXT", "Mixed"),
#  #   miuvig_quality = factor(miuvig_quality,
#  #                           levels = c("High-quality", "Genome-fragment"),
#  #                           ordered = TRUE)
  
#merged$DB_sum <- merged$N_genomes - merged$NEXT_MGS - merged$NEXT_VLP

#sources <- c("NEXT_VLP","NEXT_MGS","DB_sum")

# merged <- merged %>%
#   mutate(
#     vOTU_cluster_type = apply(across(all_of(sources)) > 0, 1,
#                          function(r) { x <- sources[r]; if (length(x)==0) "None" else paste(x, collapse = "+") })
#   )
# 
# merged_HQ <- merged[merged$miuvig_quality=="High-quality",] %>%
#   mutate(
#     vOTU_cluster_type = apply(across(all_of(sources)) > 0, 1,
#                               function(r) { x <- sources[r]; if (length(x)==0) "None" else paste(x, collapse = "+") })
#   )

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

