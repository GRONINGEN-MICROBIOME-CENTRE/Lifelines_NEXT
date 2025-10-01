setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

#############################################################
# Functions
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
ETOF_vOTUr[ETOF_vOTUr$Temperate >=0.5,]$lifestyle <- "Temperate"
ETOF_vOTUr[ETOF_vOTUr$Virulent >0.5 & ETOF_vOTUr$completeness >= 90 & !is.na(ETOF_vOTUr$completeness),]$lifestyle <- "Virulent"
ETOF_vOTUr[is.na(ETOF_vOTUr$lifestyle),]$lifestyle <- "Unknown"

newtax <- read.table('06.CLEAN_DATA/MergedTaxonomy_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

ETOF_vOTUr <- merge(ETOF_vOTUr, newtax, 
                    by.x="New_CID", by.y="Genome_ID", all.x=T)

dim(cleanest)


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

ggplot(TMP, aes(x = Samples, y = Mean_Detected_vOTUs, color=Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = Mean_Detected_vOTUs - SD, ymax = Mean_Detected_vOTUs + SD), width = 0.1) +
  labs(title = "N vOTUs detected with increasing number of samples",
       x = "Number of fecal VLP samples", y = "Number of vOTUs detected in the NEXT virome") +
  theme_minimal()


smeta$clean_richness_VLP <- colSums(clean_RPKM > 0)[match(smeta$Sequencing_ID_VLP, colnames(clean_RPKM))]
smeta$clean_richness_MGS <- colSums(clean_RPKM > 0)[match(smeta$NG_ID, colnames(clean_RPKM))]

wilcox.test(smeta$clean_richness_VLP, smeta$clean_richness_MGS, paired=T)




#############################################################
# 3. Intermediate output: contig IDs for further analysis
#############################################################
write.table(row.names(cleanest), "./06.CLEAN_DATA/VLP_MGS_vOTUr_dec99ANI_ab3kbp_2220samples.txt", sep='\t', row.names=F, col.names = F, quote=F)
write.table(cleanest, "./06.CLEAN_DATA/RPKM_table_VLP_MGS_dec99ANI_ab3kbp_2220_samples.txt", sep='\t', quote=F)
write.table(baseETOF_vOTUr, "./06.CLEAN_DATA/ETOF_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt", sep='\t', quote=F)
write.table(ETOF_vOTUr, "./06.CLEAN_DATA/Extended_ETOF_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt", sep='\t', quote=F)
