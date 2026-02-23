setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# holovirome dynamics
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
NMDS_maker <- function(data, ID, metadata, label_dist){
  
  scores_df <- as.data.frame(scores(data
                                    , display = "sites")) %>%
    rownames_to_column(var = paste0(ID)) %>%
    left_join(metadata) %>%
    filter( !(!!sym(ID) %in% c('CHV107612E10', 'CHV004401C06')))
  
  centroids <- scores_df %>%
    group_by(Timepoint_new) %>%
    summarise(V1 = median(NMDS1),
              V2 = median(NMDS2), .groups = "drop") %>%
    as.data.frame()
  
  plotik <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2, color = Timepoint_new)) +
    ggrastr::rasterise(geom_point(size = 1, alpha = 0.6), dpi = 600) +
    geom_point(data=centroids, aes(V1, V2, fill = Timepoint_new),  shape = 23, size = 3, color = "black") +
    theme_bw() +
    labs(
      x = "NMDS1",
      y = "NMDS2",
      color = "Timepoint",
      fill = "Timepoint") + 
    ggtitle(label = label_dist, subtitle = paste("stress value is", round(data$stress, 4))) +
    theme(axis.title = element_text(size=8),
          axis.text = element_text(size=7),
          title = element_text(size=8),
          plot.subtitle = element_text(size = 7),
          legend.text = element_text(size=7),
          legend.title = element_text(size = 8))
  
  return(plotik)
}


PCOA_maker <- function(data, ID, metadata, label_dist){
  
  meh <- data$points
  meh <- meh %>%
    as.data.frame() %>%
    rownames_to_column(var = paste0(ID)) %>%
    left_join(metadata)
  
  centroids <- meh %>%
    group_by(Timepoint_new) %>%
    summarise(V1 = median(V1),
              V2 = median(V2), .groups = "drop") %>%
    as.data.frame()
  
  pvar <- round(data$eig / sum(data$eig) * 100, 1)
  
  plotik <- ggplot(meh, aes(x = V1, y = V2, color = Timepoint_new)) +
    ggrastr::rasterise(geom_point(size = 1, alpha = 0.6), dpi = 600) +
    geom_point(data=centroids, aes(V1, V2, fill = Timepoint_new),  shape = 23, size = 3, color = "black") +
    theme_bw() +
    labs(
      x = paste0("PC1: ", pvar[1], "%"),
      y = paste0("PC1: ", pvar[2], "%"),
      color = "Timepoint",
      fill = "Timepoint") + 
    ggtitle(label = label_dist) +
    theme(axis.title = element_text(size=8),
          axis.text = element_text(size=7),
          title = element_text(size=8),
          plot.subtitle = element_text(size = 7),
          legend.text = element_text(size=7),
          legend.title = element_text(size = 8))
  
  return(plotik)
}

# to report the N comparisons & median of jaccard similarity
# used in 3.3 Analysis
stat_box_data <- function(y) {
  return( 
    data.frame(
      y = -0.05,  # change to rise up
      label = paste('n =', length(y), '\n',
                    'med =', round(median(y), 2), '\n')
    )
  )
}
#############################################################
# 1. Loading libraries
#############################################################
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(MetBrewer)
library(vegan)

#############################################################
# 2. Load Input Data
#############################################################
smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)
smeta$Timepoint_new <- factor(smeta$Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)

long_smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)
long_smeta$Timepoint_new <- factor(smeta$Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)

holovirome <- read.table('06.CLEAN_DATA/Intermediate/Holovirome_RPKM_1110samples_120997vOTUs.txt', sep='\t', header=T)

VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)

ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)
temperate_list <- ETOF_vOTUr$New_CID[ETOF_vOTUr$lifestyle == "Temperate"]
#############################################################
# 3.1 Analysis: holo vs VLP diversity/richness
#############################################################
smeta$vir_richness_holo <- colSums(holovirome)[match(smeta$Universal_ID, colnames(holovirome))]
smeta$temperate_richness_holo <- colSums(holovirome[row.names(holovirome) %in% temperate_list,])[match(smeta$Universal_ID, colnames(holovirome))]

# maybe holovirome grows faster / more substantially than VLP because of the temperate phage/prophage contribution
# maybe association w bacterial diversity/richness
# so more bacteria come but their phages do not really get out there or what?

holo_vs_vlp <- smeta %>%
  select(Timepoint_new, vir_richness_cf, vir_richness_holo) %>%
  pivot_longer(cols = !Timepoint_new) %>%
ggplot(aes(Timepoint_new, value, fill = name)) +
  geom_jitter(aes(color = name), position = position_jitterdodge(), size=0.1, alpha=1)+
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  theme_bw() +
  labs(x = "Timepoint", y = "N viruses", fill = "Metavirome\nfraction", color = "Metavirome\nfraction") +
  scale_fill_manual(labels = c("Active", "Holovirome"), values = c(met.brewer("Kandinsky")[2], "#547792")) +
  scale_color_manual(labels = c("Active", "Holovirome"), values = c(met.brewer("Kandinsky")[2], "#547792"))

ggsave('05.PLOTS/06.DYNAMICS/Richness_over_time_holo_VLP.pdf',
       holo_vs_vlp,  "pdf", width=10, height=8, units="cm", dpi = 300)

rich_growth <- map_dfr(c("vir_richness_cf", "vir_richness_holo"), function(richness) {
  
  formula <- as.formula(paste0(richness, "~ exact_age_months_at_collection + (1|NEXT_ID)"))
  
  model <- lmer(
    formula,
    REML = FALSE,
    data = smeta[smeta$Type == "K",]
  )
  
  model_summary <- summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(rowname != "(Intercept)") %>%
    mutate(richness_type = richness)
})

rich_growth$p_adj <- p.adjust(rich_growth$`Pr(>|t|)`, "BH")

# does holo grows faster? (yes)
llka <- smeta %>%
  select(Timepoint_new, exact_age_months_at_collection, NEXT_ID, vir_richness_cf, vir_richness_holo) %>%
  pivot_longer(cols = !c(Timepoint_new, exact_age_months_at_collection, NEXT_ID))

model <- lmer(value ~ exact_age_months_at_collection*name + (1|NEXT_ID), data = llka[llka$Timepoint_new!="Mother",])
summary(model)

#############################################################
# 3.2 Analysis: beta-diversity
#############################################################
set.seed(444) 

###### VLP
VLP_dist <- list()
pcoa_res <- list()
nmds_res <- list()

VLP_RAB <- as.data.frame(t(VLP)/colSums(VLP))
VLP_dist[["BC"]] <- vegdist(VLP_RAB, method = "bray", parallel = 6)

VLP_HELL <-as.data.frame( decostand(t(VLP), method = "hellinger"))
VLP_dist[["EUC"]] <- vegdist(VLP_HELL, method = "euclidean", parallel = 6)

VLP_rCLR <- as.data.frame( decostand(t(VLP), method = "rclr"))
VLP_dist[["AIT"]] <- vegdist(VLP_rCLR, method = "euclidean", parallel = 6)

VLP_dist[["JAC"]] <- vegdist(t(VLP)>0, method = "jaccard", parallel = 6)

dist_ord_chooser <- map_dfr(c("BC", "EUC", "AIT", "JAC"), function(distance){
  
  #PCoA
  pcoa_res[[distance]] <<- cmdscale(VLP_dist[[distance]], k = 2, eig = TRUE)
  
  pvar <- round(pcoa_res[[distance]]$eig / sum(pcoa_res[[distance]]$eig) * 100, 1)
  
  #NMDS
  nmds_res[[distance]] <<- metaMDS(
    VLP_dist[[distance]], autotransform = F, 
    k = 2,
    parallel = 6,
    trymax = 200)
  
  data.frame(ord = c('PCoA', 'NMDS'),
             var_exp = c(paste(pvar[1:3], collapse = ", "), NA),
             stress_value = c(NA, nmds_res[[distance]]$stress),
             converged = c(NA, nmds_res[[distance]]$converged),
             iterations = c(NA, nmds_res[[distance]]$iters),
             dist_method = rep(distance, 2))
  
})

###### HOLO

holo_dist <- list()
holo_pcoa_res <- list()
holo_nmds_res <- list()

holo_dist[["JAC"]] <- vegdist(t(holovirome), method = "jaccard", parallel = 6)

dist_ord_chooser_holo <- map_dfr(c("JAC"), function(distance){
  
  #PCoA
  holo_pcoa_res[[distance]] <<- cmdscale(holo_dist[[distance]], k = 2, eig = TRUE)
  
  pvar <- round(holo_pcoa_res[[distance]]$eig / sum(holo_pcoa_res[[distance]]$eig) * 100, 1)
  
  #NMDS
  holo_nmds_res[[distance]] <<- metaMDS(
    holo_dist[[distance]], autotransform = F, 
    k = 2,
    parallel = 6,
    trymax = 200)
  
  data.frame(ord = c('PCoA', 'NMDS'),
             var_exp = c(paste(pvar[1:3], collapse = ", "), NA),
             stress_value = c(NA, holo_nmds_res[[distance]]$stress),
             converged = c(NA, holo_nmds_res[[distance]]$converged),
             iterations = c(NA, holo_nmds_res[[distance]]$iters),
             dist_method = rep(distance, 2))
  
})

# VLP, BC:

# NMDS
ggsave("05.PLOTS/06.DYNAMICS/ord_plots/VLP_BC_NMDS.pdf",
       NMDS_maker(nmds_res[["BC"]], "Sequencing_ID", smeta, "VLP, Bray-Curtis"), "pdf", width=10, height=8, units="cm", dpi = 300)

# PCOA
ggsave("05.PLOTS/06.DYNAMICS/ord_plots/VLP_BC_PCoA.pdf",
       PCOA_maker(pcoa_res[["BC"]], "Sequencing_ID", smeta, "VLP, Bray-Curtis"), "pdf", width=10, height=8, units="cm", dpi = 300)

# VLP, AIT:
# NMDS
ggsave("05.PLOTS/06.DYNAMICS/ord_plots/VLP_AIT_NMDS.pdf",
       NMDS_maker(nmds_res[["AIT"]], "Sequencing_ID", smeta, "VLP, rCLR & Aitchison"), "pdf", width=10, height=8, units="cm", dpi = 300)

# PCOA
ggsave("05.PLOTS/06.DYNAMICS/ord_plots/VLP_AIT_PCoA.pdf",
       PCOA_maker(pcoa_res[["AIT"]], "Sequencing_ID", smeta, "VLP, rCLR & Aitchison"), "pdf", width=10, height=8, units="cm", dpi = 300)

# HOLO, JAC:

# NMDS
ggsave("05.PLOTS/06.DYNAMICS/ord_plots/HOLO_JAC_NMDS.pdf",
       NMDS_maker(holo_nmds_res[["JAC"]], "Universal_ID", smeta, "Holovirome, Jaccard"), "pdf", width=10, height=8, units="cm", dpi = 300)

# PCOA
ggsave("05.PLOTS/06.DYNAMICS/ord_plots/HOLO_JAC_PCoA.pdf",
       PCOA_maker(holo_pcoa_res[["JAC"]], "Universal_ID", smeta, "Holovirome, Jaccard"), "pdf", width=10, height=8, units="cm", dpi = 300)


#############################################################
# 3.3 Analysis: convergence w maternal virome over time
#############################################################
jaccard_similarity <- 1 - as.matrix(VLP_dist$JAC)

upper_tri_idx <- which(upper.tri(jaccard_similarity, diag = FALSE), arr.ind = TRUE)

esmeta <- smeta %>%
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) # making the 2nd pregnancy as separate fams


dist_list <- data.frame(
  sample_1 = rownames(jaccard_similarity)[upper_tri_idx[, 1]],
  sample_2 = colnames(jaccard_similarity)[upper_tri_idx[, 2]],
  similarity = jaccard_similarity[upper_tri_idx]
) %>%
  left_join(esmeta %>% select(Sequencing_ID, 
                              Universal_ID,
                              Modified_NEXT_ID_without_preg_number,
                              FAMILY,
                              FAMILYupd, 
                              Type, 
                              Timepoint_new), by = c("sample_1" = "Sequencing_ID" )) %>%
  left_join(esmeta %>% select(Sequencing_ID, 
                              Universal_ID,
                              Modified_NEXT_ID_without_preg_number,
                              FAMILY,
                              FAMILYupd,
                              Type,
                              Timepoint_new), by = c("sample_2" = "Sequencing_ID"), suffix = c("_1", "_2")) %>%
  mutate(same_dyad = ifelse(FAMILYupd_1 == FAMILYupd_2 &
                              Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                              Type_1 != Type_2, T, F)) %>%
  mutate(unrelated = ifelse(FAMILY_1 != FAMILY_2, T, F)) %>%
  mutate(same_infant = ifelse(Universal_ID_1 != Universal_ID_2 &
                                      Type_1 == "K" &
                                      Type_2 == "K" &
                                      Modified_NEXT_ID_without_preg_number_1 == Modified_NEXT_ID_without_preg_number_2, T, F)) %>%
  mutate(unrelated_dyad = ifelse(FAMILYupd_1 != FAMILYupd_2 &
                                   Type_1 != Type_2, T, F))

plot_data <- dist_list %>%
  select(similarity, same_dyad:same_infant) %>%
  pivot_longer(cols = !similarity, 
               names_to = "Category", 
               values_to = "Is_True") %>%
  filter(Is_True == TRUE) %>%
  mutate(method = gsub(".*_", "", Category)) %>%
  mutate(method = if_else(method %in% c('inf', 'mom'), "inter", method) )

ggplot(plot_data, aes(x = Category, y = similarity, fill = method)) +
  geom_violin(alpha = 0.7, linewidth = 0.4, scale = "width") + 
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, linewidth = 0.4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        title = element_text(size=10),
        axis.title = element_text(size=7)) +
  labs(title = "Similarity by method and kinship",
       x = "",
       y = "Jaccard similarity") +
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9, size = 2
  )

# virome of dyads is more similar than virome of non-dyads and infant virome resembles more of its mother virome?


baba_data <- dist_list %>%
  filter(same_dyad == TRUE) %>%
  # Use a temporary column to swap
  mutate(
    temp_s1 = sample_1,
    sample_1 = ifelse(Timepoint_new_1 != "Mother", sample_2, sample_1),
    sample_2 = ifelse(Timepoint_new_1 != "Mother", temp_s1, sample_2)
  ) %>%
  select(-temp_s1)

baba_data <- baba_data %>%
  select(sample_1, sample_2, similarity) %>%
  left_join(esmeta %>% select(Sequencing_ID, 
                              Universal_ID,
                              Modified_NEXT_ID_without_preg_number,
                              FAMILY,
                              FAMILYupd, 
                              Type, 
                              Timepoint_new), by = c("sample_1" = "Sequencing_ID" )) %>%
  left_join(esmeta %>% select(Sequencing_ID, 
                              Universal_ID,
                              Modified_NEXT_ID_without_preg_number,
                              FAMILY,
                              FAMILYupd,
                              Type,
                              Timepoint_new, 
                              exact_age_months_at_collection), by = c("sample_2" = "Sequencing_ID"), suffix = c("_1", "_2"))


ggplot(baba_data, aes(x = Timepoint_new_2, y = similarity, fill = Timepoint_new_2)) +
  geom_violin(alpha = 0.7, linewidth = 0.4, scale = "width") + 
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, linewidth = 0.4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        title = element_text(size=10),
        axis.title = element_text(size=7)) +
  labs(title = "Similarity by method and kinship",
       x = "",
       y = "Jaccard similarity") +
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9, size = 2
  )

model <- lmer(similarity ~ exact_age_months_at_collection + (1|Modified_NEXT_ID_without_preg_number_2), data = baba_data)
summary(model)$coefficients


vlp_converge <- ggplot(baba_data, aes(x = exact_age_months_at_collection, y = similarity)) +
  geom_line(aes(group = Modified_NEXT_ID_without_preg_number_2), alpha = 1, linewidth = 0.2, color = "gray60") +  
  geom_smooth(method = "lm", color = "firebrick", linewidth = 1, se = TRUE, fill = "firebrick", alpha = 0.3) +
  stat_summary(aes(x = as.numeric(gsub("M", "", Timepoint_new_2)), 
                   group = Timepoint_new_2), 
               fun.data = mean_cl_boot, 
               geom = "errorbar", 
               width = 0.5, 
               color = "darkred") +
  stat_summary(aes(x = as.numeric(gsub("M", "", Timepoint_new_2)),
                   group = Timepoint_new_2), 
               fun = mean, 
               geom = "point", 
               size = 1.2, 
               color = "darkred") +
  scale_x_continuous(breaks = c(0, 1, 3, 6, 12)) +
  theme_minimal() + 
  labs(x = "Infant age, months", y = "Jaccard similarity to maternal virome, VLP") + 
  annotate(geom="text", x = 3, y=0.3, label = "beta = 6.2e-3\np-value=5.9e-45", size = 2) +
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=7))
  

ggsave("05.PLOTS/06.DYNAMICS/VLP_infant_converge_mom_trajectory.pdf", 
       vlp_converge, "pdf", width=8, height=8, units="cm", dpi = 300)


nebaba_data <- baba_data %>%
  left_join(smeta %>% select(Sequencing_ID, vir_richness_cf, vir_diversity, temperate_richness, lytic_richness, bacShannon, contigs_1000_bp, perc_aligned_cf, reads_lost_QC, sc_enrichment, sequencing_batch), by = c("sample_2" = "Sequencing_ID"))

modelka <- lmer(similarity ~ exact_age_months_at_collection*bacShannon + (1|Modified_NEXT_ID_without_preg_number_2), data = nebaba_data)
summary(modelka)$coefficients


# deda data:
deda_data <- dist_list %>%
  filter(unrelated_dyad == TRUE) %>%
  # Use a temporary column to swap
  mutate(
    temp_s1 = sample_1,
    sample_1 = ifelse(Timepoint_new_1 != "Mother", sample_2, sample_1),
    sample_2 = ifelse(Timepoint_new_1 != "Mother", temp_s1, sample_2)
  ) %>%
  select(-temp_s1)

deda_data <- deda_data %>%
  select(sample_1, sample_2, similarity) %>%
  group_by(sample_2) %>%
  summarise(mean_dist = mean(similarity)) %>%
  ungroup() %>%
  left_join(esmeta %>% select(Sequencing_ID, 
                              Universal_ID,
                              exact_age_months_at_collection,
                              Modified_NEXT_ID_without_preg_number,
                              FAMILY,
                              FAMILYupd, 
                              Type, 
                              Timepoint_new), by = c("sample_2" = "Sequencing_ID" )) %>%
  select(sample_2, mean_dist, Timepoint_new, exact_age_months_at_collection, Modified_NEXT_ID_without_preg_number) %>%
  rename(similarity = mean_dist, Timepoint_new_2 = Timepoint_new, Modified_NEXT_ID_without_preg_number_2 = Modified_NEXT_ID_without_preg_number)

better_together <- baba_data %>%
  select(sample_2, similarity, Timepoint_new_2, exact_age_months_at_collection, Modified_NEXT_ID_without_preg_number_2) %>%
  mutate(Mother = "Own") %>%
  bind_rows(deda_data %>% mutate(Mother = "Unrelated")) %>%
  mutate(Line_ID = paste0(Modified_NEXT_ID_without_preg_number_2, "_", Mother)) %>%
  mutate(Timepoint_group = paste0(Timepoint_new_2, "_", Mother))

vlp_converge_moms <- ggplot(better_together, aes(x = exact_age_months_at_collection, y = similarity)) +
  geom_line(aes(group = Line_ID, color = Mother), alpha = 0.7, linewidth = 0.2) +  
  scale_color_manual(values = c("grey60", "grey90")) +
  ggnewscale::new_scale_color() +
  geom_smooth(aes(color = Mother, fill = Mother), method = "lm",  linewidth = 1, se = TRUE, alpha = 0.3) +
  scale_color_manual(values = c("firebrick", "#67B2D8")) +
  scale_fill_manual(values = c("firebrick", "#67B2D8")) +
  ggnewscale::new_scale_color() +
  stat_summary(aes(x = as.numeric(gsub("M", "", Timepoint_new_2)), 
                   group = Timepoint_group, color = Mother), 
               fun.data = mean_cl_boot, 
               geom = "errorbar", 
               width = 0.5) +
  stat_summary(aes(x = as.numeric(gsub("M", "", Timepoint_new_2)),
                   group = Timepoint_group, color = Mother), 
               fun = mean, 
               geom = "point", 
               size = 1.2) +
  scale_color_manual(values = c("darkred", "darkblue")) +
  scale_x_continuous(breaks = c(0, 1, 3, 6, 12)) +
  theme_minimal() + 
  labs(x = "Infant age, months", y = "Jaccard similarity to maternal virome, VLP") + 
  annotate("rect", xmin = 1, xmax = 5, ymin = 0.28, ymax = 0.32, 
           alpha = 0.2, fill = "firebrick") +
  annotate(geom="text", x = 3, y=0.3, label = "beta = 6.2e-3\np-value=5.9e-45", size = 2) +
  annotate("rect", xmin = 1, xmax = 5, ymin = 0.22, ymax = 0.26, 
           alpha = 0.2, fill = "#67B2D8") +
  annotate(geom="text", x = 3, y=0.24, label = "beta = 3.7e-3\np-value=2.5e-87", size = 2) +
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=7),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7),
        legend.position = "right")

ggsave("05.PLOTS/06.DYNAMICS/VLP_infant_converge_mom_unrelated_trajectory.pdf", 
       vlp_converge_moms, "pdf", width=10, height=8, units="cm", dpi = 300)

model <- lmer(similarity ~ exact_age_months_at_collection + (1|Modified_NEXT_ID_without_preg_number_2), data = baba_data)
summary(model)$coefficients

model2 <- lmer(similarity ~ exact_age_months_at_collection + (1|Modified_NEXT_ID_without_preg_number_2), data = deda_data)
summary(model2)$coefficients

# now I forgot how it looked in holovirome.. actually useful to have it anyways so make it for holo as well anyways
#############################################################
# 3.4 Analysis: Which viruses drive convergence?
#############################################################
infant_VLP <- VLP %>%
  select(all_of(smeta %>% filter(Type == "K") %>% pull(Sequencing_ID)))

infant_VLP <- infant_VLP[rowSums(infant_VLP > 0) > 0.05*ncol(infant_VLP),] # 10% prevalence
infant_VLP <- infant_VLP[,colSums(infant_VLP) > 0]

infant_VLP <- as.data.frame(t(infant_VLP))

adamdriver <- baba_data %>% 
  left_join(infant_VLP %>% rownames_to_column(var = "sample_2"))

adtest <- map_dfr(colnames(infant_VLP), function(vOTU){
  
  formula <- as.formula(paste0("similarity ~ ", "log(`",vOTU, "` + 1 )", " + exact_age_months_at_collection + (1|Modified_NEXT_ID_without_preg_number_2)"))
  
  model <- lmer(
    formula,
    REML = FALSE,
    data = adamdriver
  )
  
 summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(!rowname %in% c("(Intercept)", "exact_age_months_at_collection")) %>%
    mutate(vOTU = vOTU)
  
})

adtest <- adtest %>%
  mutate(p_adjusted = p.adjust(`Pr(>|t|)`, "BH")) %>%
  filter(p_adjusted < 0.05) %>%
  left_join(ETOF_vOTUr, by = c("vOTU" = "New_CID")) 

summary_peaks <- adamdriver %>%
  select(Timepoint_new_2, all_of(adtest$vOTU)) %>%
  pivot_longer(!Timepoint_new_2) %>%
  group_by(name, Timepoint_new_2) %>%
  summarise(mean_abund = mean(value)) %>%
  group_by(name) %>%
  summarise(peak_time = Timepoint_new_2[which.max(mean_abund)])

babydrive <- infant_VLP %>%
  select(all_of(adtest$vOTU)) %>%
  rownames_to_column(var = "Sequencing_ID") %>%
  left_join(smeta %>% select(Sequencing_ID, Timepoint_new)) %>%
  select(-Sequencing_ID) %>%
  pivot_longer(!Timepoint_new) %>%
  group_by(name,Timepoint_new) %>%
  summarise(mean =  mean(value))

# 1. Calculate Z-scores within each vOTU
babydrive_scaled <- babydrive %>%
  group_by(name) %>%
  # Z-score = (value - mean) / standard deviation
  mutate(z_score = scale(mean)) %>% 
  ungroup()

# 2. Identify the peak timepoint for each vOTU (to use for sorting)
votu_order <- babydrive_scaled %>%
  group_by(name) %>%
  slice_max(z_score, n = 1, with_ties = FALSE) %>%
  arrange(Timepoint_new, desc(z_score)) %>%
  pull(name)

# 3. Pivot to Matrix format for the heatmap
heatmap_matrix <- babydrive_scaled %>%
  select(name, Timepoint_new, z_score) %>%
  pivot_wider(names_from = Timepoint_new, values_from = z_score) %>%
  # Turn into a matrix and order by the peak time we found
  column_to_rownames("name") %>%
  .[votu_order, ] %>%
  as.matrix()

library(ComplexHeatmap)
library(circlize)

# 1. 'mat_scaled' is your matrix of Z-scored mean abundances (vOTUs x Timepoints)
# 2. 'anno_df' contains your LifeCycle, MaternalMatch, and HostGenus

badtest <- adtest  %>%
  separate(
    col    = Host_taxonomy,
    into = c('Domain', 'Phylum', 'Class', 'Order', 'Family',  'Genus'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  group_by(Genus) %>%
  summarise(sum = n()) %>%
  filter(sum < 12) %>%
  pull(Genus)

anno_df <- adtest  %>%
  separate(
    col    = Host_taxonomy,
    into = c('Domain', 'Phylum', 'Class', 'Order', 'Family',  'Genus'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  mutate(Host_genus = case_when(Genus == "g__Phocaeicola" ~ "g__Bacteroides",
                                Genus %in% c(badtest, "g__Unclassified") ~ "Other",
                                .default = Genus)) %>%
  arrange(match(vOTU, votu_order)) %>%
  column_to_rownames("vOTU") %>%
  select(Host_genus)


unique_hosts <- gsub("g__", "", unique(anno_df$Host_genus))
host_cols <- setNames(as.character(MetBrewer::met.brewer("Renoir", n = 6)), 
                      unique_hosts)


# Define the annotation
row_anno <- rowAnnotation(
  Host = gsub("g__", "", anno_df$Host_genus),
  col = list(Host = host_cols),
  show_annotation_name = TRUE
)


plotka <- Heatmap(heatmap_matrix, 
        name = "Z-score", 
        cluster_columns = FALSE, 
        cluster_rows = FALSE,      # Important: keep your manual "tornado" sort
        right_annotation = row_anno,
        show_row_names = FALSE, row_title = "Covergence-driving vOTUs", column_title = "Timepoint", column_names_rot = 0 
        )     #

pdf("05.PLOTS/06.DYNAMICS/busya.pdf", width=10/2.5, height=8/2.5)
Heatmap(heatmap_matrix, 
        name = "Z-score", 
        cluster_columns = FALSE, 
        cluster_rows = FALSE,      # Important: keep your manual "tornado" sort
        right_annotation = row_anno,
        show_row_names = FALSE, row_title = "Covergence-driving vOTUs", column_title = "Timepoint", column_names_rot = 0 
) 
dev.off()

# Create color mappings
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

Heatmap(heatmap_matrix, 
        name = "Z-score Abundance",
        col = col_fun,
        cluster_columns = FALSE,     # Keep time chronological
        cluster_rows = FALSE,        # We manually sorted them by peak time
        row_split = badtest$Genus, # Optional: split the heatmap by life cycle
        show_row_names = FALSE)

#############################################################
# 3.5 Analysis: genome-type sample composition over time
#############################################################
holovirome %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, genome)) %>%
  group_by(genome) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(!genome) %>%
  mutate(Timepoint_new = metadata$Timepoint_new[match(name, metadata$Universal_ID)]) %>%
  mutate(genome = factor(genome, levels = c("dsDNA", "ssDNA", "RNA", "Unclassified"), ordered = T)) %>%
  ggplot(aes(Timepoint_new, value, fill=genome)) +
  geom_bar(position="fill", stat="identity")

buz <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, genome)) %>%
  group_by(genome) %>%
  summarise(across(where(is.numeric), ~sum(. > 0, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(!genome) %>%
  left_join(metadata %>% select(Sequencing_ID, Timepoint_new, NEXT_ID, vir_richness_cf, vir_diversity, bacShannon), by = c("name" = "Sequencing_ID"))


buzy <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, genome)) %>%
  group_by(genome) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x, na.rm = TRUE))) %>%
  pivot_longer(!genome) %>%
  left_join(metadata %>% select(Sequencing_ID, Timepoint_new, NEXT_ID, vir_richness_cf, vir_diversity, bacShannon), by = c("name" = "Sequencing_ID"))


relab <- buzy%>%
  filter(genome!="Unclassified") %>%
  mutate(genome = factor(genome, levels = c("dsDNA", "ssDNA", "RNA", "Unclassified"), ordered = T)) %>%
  ggplot(aes(Timepoint_new, value, fill=genome)) +
  ggrastr::rasterise(geom_jitter(aes(color = genome), position = position_jitterdodge(), size=0.1, alpha=1), dpi = 300)+
  geom_boxplot(aes(fill = genome), alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~genome, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_rect(NA),
        legend.position = "none",
        axis.title = element_text(size=8),
        axis.text = element_text(size=7),
        strip.text = element_text(size = 8)) + 
  labs(x = "Timepoint", y = "Abundance", fill = "Viral genome", color = "Viral genome") + 
  scale_fill_manual(values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                               "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                               "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2])) +
  scale_color_manual(values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                                "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                                "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2]))

ggsave("05.PLOTS/06.DYNAMICS/VLP_abundance_by_genome.pdf",
       relab, "pdf", width = 14, height = 7, units = "cm")

rich_prop <- buz%>%
  filter(genome!="Unclassified") %>%
  mutate(genome = factor(genome, levels = c("dsDNA", "ssDNA", "RNA", "Unclassified"), ordered = T)) %>%
  ggplot(aes(Timepoint_new, value/vir_richness_cf, fill=genome)) +
  ggrastr::rasterise(geom_jitter(aes(color = genome), position = position_jitterdodge(), size=0.1, alpha=1), dpi = 300)+
  geom_boxplot(aes(fill = genome), alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~genome, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_rect(NA),
        legend.position = "none",
        axis.title = element_text(size=8),
        axis.text = element_text(size=7),
        strip.text = element_text(size = 8)) + 
  labs(x = "Timepoint", y = "Richness proportion", fill = "Viral genome", color = "Viral genome") + 
  scale_fill_manual(values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                               "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                               "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2])) +
  scale_color_manual(values = c("dsDNA"=MetBrewer::met.brewer("VanGogh2")[4],
                                "RNA"=MetBrewer::met.brewer("VanGogh2")[1], 
                                "ssDNA"=MetBrewer::met.brewer("VanGogh2")[2]))

ggsave("05.PLOTS/06.DYNAMICS/VLP_richness_by_genome.pdf",
       rich_prop, "pdf", width = 14, height = 7, units = "cm")


by_genome_test <- buz %>%
  select(name, genome, value) %>%
  pivot_wider( names_from = genome, id_cols = name) %>%
  select(-Unclassified) %>%
  mutate(info_type = "Richness") %>%
  bind_rows(buzy %>%  select(name, genome, value) %>%
              pivot_wider( names_from = genome, id_cols = name) %>%
              select(-Unclassified) %>%
              mutate(info_type = "Abundance") ) %>%
  left_join(metadata %>% select(Sequencing_ID, NEXT_ID, Timepoint_new, exact_age_months_at_collection), by = c("name" = "Sequencing_ID")) %>%
  mutate(ssDNAtodsDNA = log((ssDNA + 1)/(dsDNA + 1)),
         RNAtodsDNA = log((RNA + 1)/(dsDNA + 1)))

coda_genome <- map_dfr(c('Richness', 'Abundance'), function(var_type){
  
  map_dfr(c('ssDNAtodsDNA', 'RNAtodsDNA'), function(prop_type) {
    
    formula <- as.formula(paste0(prop_type, " ~ ", " exact_age_months_at_collection + (1|NEXT_ID)"))
    
    model <- lmer(
      formula,
      REML = FALSE,
      data = by_genome_test[by_genome_test$info_type == var_type & by_genome_test$Timepoint_new != "Mother",]
    )
    
    summary(model)$coefficients %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      filter(!rowname %in% c("(Intercept)")) %>%
      mutate(variable = var_type,
             genome = prop_type)
  })
  
  
  
})

coda_genome$p_adjust <- p.adjust(coda_genome$`Pr(>|t|)`, "BH")

#############################################################
# 3.6 Analysis: genome-type sample composition over time
#############################################################


