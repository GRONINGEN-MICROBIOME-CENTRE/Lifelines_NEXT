setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# Composition dynamics of VLP data
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
    ggrastr::rasterise(geom_point(size = 1, alpha = 0.8), dpi = 600) +
    geom_point(data=centroids, aes(V1, V2, fill = Timepoint_new), alpha = 0.8, shape = 23, size = 3, color = "black") +
    geom_path(
      data = centroids %>% filter(Timepoint_new != "Mother"),
      aes(x = V1, y = V2),
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      colour = "firebrick",
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    theme_bw() +
    labs(
      x = paste0("PC1: ", pvar[1], "%"),
      y = paste0("PC1: ", pvar[2], "%"),
      color = "Timepoint",
      fill = "Timepoint") + 
    scale_fill_manual(values = met.brewer("Cassatt1")[c(4:1, 8)]) +
    scale_color_manual(values = met.brewer("Cassatt1")[c(4:1, 8)]) +
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
library(vegan)
library(MetBrewer)
#############################################################
# 2. Load Input Data
#############################################################
# abundance tables etc:
VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)

# metadatas:
smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)

smeta <- smeta %>%
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) %>% # treating 2nd pregnancy as a separate family:
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T))

#############################################################
# 3.0 Analysis: calculating distances
#############################################################
set.seed(444)

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

################ saving ordination plots ################ 
# VLP, BC:

#NMDS
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

# VLP, JAC:

# NMDS
ggsave("05.PLOTS/06.DYNAMICS/ord_plots/VLP_JAC_NMDS.pdf",
      NMDS_maker(nmds_res[["JAC"]], "Sequencing_ID", smeta, "VLP, Jaccard"), "pdf", width=10, height=8, units="cm", dpi = 300)

# PCOA
ggsave("05.PLOTS/06.DYNAMICS/ord_plots/VLP_JAC_PCoA.pdf",
      PCOA_maker(pcoa_res[["JAC"]], "Sequencing_ID", smeta, "VLP, Jaccard"), "pdf", width=10, height=8, units="cm", dpi = 300)

######## best solution #########
PCOA_maker(pcoa_res[["BC"]], "Sequencing_ID", smeta, "VLP, Bray-Curtis")
#############################################################
# 3.1 Analysis: individual specificity
#############################################################
# recalculate stat

similarity <- 1 - as.matrix(VLP_dist$BC)

upper_tri_idx <- which(upper.tri(similarity, diag = FALSE), arr.ind = TRUE)

# exact copy of a part from 03.VLP_MGS_sample_comparison.R
dist_list <- data.frame(
  sample_1 = rownames(similarity)[upper_tri_idx[, 1]],
  sample_2 = colnames(similarity)[upper_tri_idx[, 2]],
  similarity = similarity[upper_tri_idx]
) %>%
  left_join(smeta %>% select(Sequencing_ID, 
                             Universal_ID,
                             Modified_NEXT_ID_without_preg_number,
                             FAMILY,
                             FAMILYupd, 
                             Type, 
                             Timepoint_new), by = c("sample_1" = "Sequencing_ID" )) %>%
  left_join(smeta %>% select(Sequencing_ID, 
                             Universal_ID,
                             Modified_NEXT_ID_without_preg_number,
                             FAMILY,
                             FAMILYupd,
                             Type,
                             Timepoint_new), by = c("sample_2" = "Sequencing_ID"), suffix = c("_1", "_2")) %>%
  # mother and infant from the same family (pregnancy)
  mutate(same_dyad = ifelse(FAMILYupd_1 == FAMILYupd_2 &
                              Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                              Type_1 != Type_2, T, F)) %>%
  # twins:
  mutate(twins = ifelse(FAMILYupd_1 == FAMILYupd_2 &
                          Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                          Type_1 == "K" &
                          Type_2 == "K",
                        T, F)) %>%
  # unrelated
  mutate(unrelated = ifelse(FAMILY_1 != FAMILY_2, T, F)) %>%
  # unrelated infants:
  mutate(unrelated_inf = ifelse(FAMILY_1 != FAMILY_2 & 
                                  Type_1 == "K" &
                                  Type_2 == "K", 
                                T, F)) %>%
  # unrelated mothers:
  mutate(unrelated_mom = ifelse(FAMILY_1 != FAMILY_2 & 
                                  Type_1 == "M" &
                                  Type_2 == "M", 
                                T, F)) %>%
  # within-infant longitudinal similarities
  mutate(same_infant = ifelse(Universal_ID_1 != Universal_ID_2 &
                                Type_1 == "K" &
                                Type_2 == "K" &
                                Modified_NEXT_ID_without_preg_number_1 == Modified_NEXT_ID_without_preg_number_2, T, F)) %>%
  # similarity between unrelated mother-infant pair
  mutate(unrelated_dyad = ifelse(FAMILY_1 != FAMILY_2 &
                                   Type_1 != Type_2, T, F))

plot_data <- dist_list %>%
  select(similarity, same_dyad:same_infant) %>%
  pivot_longer(cols = !similarity, 
               names_to = "Category", 
               values_to = "Is_True") %>%
  filter(Is_True == TRUE) %>%
  mutate(method = gsub(".*_", "", Category)) %>%
  mutate(method = if_else(method %in% c('inf', 'mom'), "inter", method) )

my_comparisons <- list( c("same_infant", "twins"), 
                        c("same_infant", "unrelated_inf"), 
                        c("twins", "unrelated_inf"),
                        c("unrelated_inf", "unrelated_mom"),
                        c("same_infant", "unrelated_mom"))

long_vs_unrelated <- plot_data %>%
  filter(Category %in% c('same_infant', 'twins', 'unrelated_inf', 'unrelated_mom')) %>%
  mutate(Kinship = if_else(grepl('unrelated', Category), 'Unrelated', 'Related')) %>%
  ggplot(aes(x = Category, y = similarity, fill = Kinship)) +
  geom_violin(alpha = 0.7, linewidth = 0.4, scale = "width") + 
  scale_fill_manual(labels = c("Related", "Unrelated"), values = c("#FFCAD4", 
                                                                   "#B7C9F2")) +
  ggnewscale::new_scale_fill() +
  geom_boxplot(aes(fill = Kinship), alpha = 0.7, width = 0.1, color = "black", outlier.shape = NA, linewidth = 0.4) +
  scale_fill_manual(labels = c("Related", "Unrelated"), values = c("firebrick", "#67B2D8")) +
  theme_bw() +
  scale_x_discrete(labels = c("Within-\ninfant", "Between\ntwins", "Between\ninfants", "Between\nmothers")) +
  theme(axis.text.x = element_text(angle = 0),
        legend.position = "bottom",
        axis.title = element_text(size=10),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  labs(x = "",
       y = "Similarity (1 - BC dissimilarity)") +
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.9, size = 2
  ) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size=3, p.adjust.method = "BH")

ggsave('05.PLOTS/06.DYNAMICS/Within_vs_between.pdf', 
       long_vs_unrelated, "pdf", width=8, height=11, units="cm", dpi = 300) # be aware that ggsave conflicts a lot w ggpubr here, so first delete the previous version of the figure and then try to save

#############################################################
# 3.2 Permutation analysis for BC similarity
#############################################################
# partially copied and pasted from 03.VLP_MGS_sample_comparison.R

# n diff possible (and sensible) combos:

busy <- c('same_infant',
          'twins', 
          'unrelated_inf', 
          'unrelated_mom')

combos <- as.data.frame(t(combn(busy, 2))) %>%
  mutate(p_value_real = NA) %>%
  mutate(pval_perm_based = NA)

start.time <- Sys.time()

for (combo in 1:nrow(combos)) {
  
  feature1 <- combos[combo, 1]
  feature2 <- combos[combo, 2]
  
  print(paste("Testing", feature1, "versus", feature2))
  
  # n permutations
  n_perm <- 1:1000
  
  # subset rows only where features are true, otherwise comparison is not fair
  df_coi <- dist_list %>%
    filter(!!sym(feature1) == T | !!sym(feature2) == T) %>%
    mutate(new_factor = case_when(
      !!sym(feature1) ~ feature1,
      !!sym(feature2) ~ feature2
    )) # creates new factor to make it easier to calculate wilcoxon
  
  if (any(df_coi[,feature1] & df_coi[, feature2], na.rm = TRUE)) {
    warning("Overlapping between ", feature1, "and", feature2)
  } # should not be the case if no same_feces_inf / same_feces_mom are chosen together w same_feces_inter
  
  original <- wilcox.test(similarity ~ new_factor, data = df_coi)
  
  # store perm pvalues
  p_perm <- numeric()
  
  for (i in n_perm) {
    
    print(i)
    
    # if comparing unrelated vs unrelated, no permutation restricting by ID:
    if ( (feature1 == "unrelated_mom" & feature2 == "unrelated_inf") | (feature1 == "unrelated_inf" & feature2 == "unrelated_mom") ) {
      perm_table <- df_coi %>%
        mutate(new_factor = sample(new_factor)) # permuting assignment
      
    } else {
      # else restricting because virome is individual specific
      perm_table <- df_coi %>%
        group_by(Modified_NEXT_ID_without_preg_number_1) %>% 
        mutate(new_factor = sample(new_factor)) %>% #  permuting assignment
        ungroup()
    }
    
    res_perm <- wilcox.test(similarity ~ new_factor, data = perm_table)
    
    p_perm[i] <- res_perm$p.value
    
  }
  
  combos[combo, "p_value_real"] <- original$p.value
  combos[combo, "pval_perm_based"] <- sum(p_perm <= original$p.value)/length(n_perm)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

combos$pval_adj <- p.adjust(combos$pval_perm_based, method = "BH")
combos$n_perm <- i

writexl::write_xlsx(combos, '07.RESULTS/Compare_composition_similarity_BC_1000perm.xlsx')
#############################################################
# 3.3 Analysis: infant virome composition dynamics
#############################################################
metadata <- smeta[match(labels(VLP_dist[["BC"]]), smeta$Sequencing_ID), ]

metadata <- metadata %>%
  filter(Type == "K") %>%
  droplevels()

IBC <- as.matrix(VLP_dist[["BC"]])
IBC <- IBC[row.names(IBC) %in% metadata$Sequencing_ID,
                 colnames(IBC) %in% metadata$Sequencing_ID]

###### TIME ######
result <- adonis2(
  IBC ~ Timepoint_new,
  data = metadata,
  permutations = 999,
  strata = metadata$NEXT_ID
) # Model      3     9.10 0.0232 6.5232  0.001 ***
result <- as.data.frame(result)
writexl::write_xlsx(combos, '07.RESULTS/Adonis_infant_VLP_virome_time_1000perm.xlsx')

###### COMPOSITION SUBJECT-SPECIFCITY ######
result_inf_spec <- adonis2(IBC ~ NEXT_ID, data = metadata, permutations = 999)
result_inf_spec <- as.data.frame(result_inf_spec)
writexl::write_xlsx(combos, '07.RESULTS/Adonis_infant_VLP_virome_individual_specificity_1000perm.xlsx')

###### CONVERGENCE TO INFANT CENTROID ######
bd <- betadisper(as.dist(IBC), group = metadata$Timepoint_new)

# differs between timepoints?
inf_centr <- anova(bd) %>%
  as.data.frame()
writexl::write_xlsx(inf_centr, '07.RESULTS/Betadisper_infant_VLP_virome.xlsx')

move_to_centroid <-as.data.frame( bd$distances) %>%
  rename("Distance_to_centroid" = 1) %>%
  mutate(similarity_to_centroid = 1 - Distance_to_centroid) %>%
  rownames_to_column(var = "Sequencing_ID") %>%
  left_join(metadata) %>%
  ggplot(aes(Timepoint_new, similarity_to_centroid, color = Timepoint_new, fill = Timepoint_new)) +
  ggrastr::rasterise(geom_jitter(position = position_jitterdodge(), size=0.1), dpi = 600) +
  geom_boxplot(alpha = 0.3, color = "black", outlier.shape = NA, width= 0.3) +
  scale_fill_manual(values = met.brewer("Cassatt1")[c(4:1)]) +
  scale_color_manual(values = met.brewer("Cassatt1")[c(4:1)]) +
  labs(x = "Timepoint", y = "Similarity to centroid (1 - BC dissimilarity)", fill = "Timepoint", color = "Timepoint") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=10),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('05.PLOTS/06.DYNAMICS/Betadisper_infant_VLP_similarity.pdf', 
       move_to_centroid, "pdf", width=8, height=11, units="cm", dpi = 300)
#############################################################
# 3.3 Analysis: convergence w maternal virome over time
#############################################################
plot_data <- dist_list %>%
  select(similarity, same_dyad:same_infant) %>%
  pivot_longer(cols = !similarity, 
               names_to = "Category", 
               values_to = "Is_True") %>%
  filter(Is_True == TRUE) %>%
  mutate(method = gsub(".*_", "", Category)) %>%
  mutate(method = if_else(method %in% c('inf', 'mom'), "inter", method) )

# convergence to the virome of own mother/unrelated adult
inftomom <- dist_list %>%
  filter(same_dyad == TRUE | unrelated_dyad == TRUE) %>%
  mutate(
    temp_s1 = sample_1,
    sample_1 = ifelse(Timepoint_new_1 != "Mother", sample_2, sample_1),
    sample_2 = ifelse(Timepoint_new_1 != "Mother", temp_s1, sample_2)
  ) %>%
  select(-temp_s1) %>%
  rename(mother=sample_1, infant=sample_2) %>%
  select(mother, infant, similarity, same_dyad, unrelated_dyad) %>%
  left_join(smeta, by = c("infant" = "Sequencing_ID"))

# selecting a mean of similarities to unrelated adult, otherwise difficult trajectories
inftomom_unrelated <- inftomom %>%
  filter(unrelated_dyad == T) %>%
  select(mother, infant, similarity) %>%
  group_by(infant) %>%
  summarise(mean_dist = mean(similarity)) %>%
  ungroup() %>%
  left_join(smeta, by = c("infant" = "Sequencing_ID" )) %>%
  rename(similarity = mean_dist)

# combining data for infant vs own and infant vs unrelated as a long data:
better_together <- inftomom %>%
  filter(same_dyad == T) %>%
  mutate(Mother = "Own") %>%
  select(-c(same_dyad, mother, unrelated_dyad)) %>%
  bind_rows(inftomom_unrelated %>% mutate(Mother = "Unrelated")) %>%
  mutate(Line_ID = paste0(Modified_NEXT_ID_without_preg_number, "_", Mother)) %>%
  mutate(Timepoint_group = paste0(Timepoint_new, "_", Mother))

# does similarity to own mom grows over time?

# save for the supplementary table & test if any virmetrics support this convergence such as viral diversity

convergence_stat <- map_dfr(c('Own', 'Unrelated'), function(mom){
  
  F1 <- as.formula(paste0("similarity ~ exact_age_months_at_collection + (1|Modified_NEXT_ID_without_preg_number)"))
  
  model <- lmer(F1, data = better_together[better_together$Mother == mom,])
  
  summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column(var = "rowname") %>%
    filter(rowname != "(Intercept)") %>%
    mutate(Mother_relation = mom)
  
})

writexl::write_xlsx(convergence_stat, '07.RESULTS/Convergence_to_mom.xlsx')
# is the slope for related vs unrelated higher?
model3 <- lmer(similarity ~ exact_age_months_at_collection*Mother + (1|Modified_NEXT_ID_without_preg_number), data = better_together)
summary(model3)$coefficients

vlp_converge_moms <- ggplot(better_together, aes(x = exact_age_months_at_collection, y = similarity)) +
  geom_line(aes(group = Line_ID, color = Mother), alpha = 0.7, linewidth = 0.2) +  
  scale_color_manual(values = c("#FFCAD4", "#B7C9F2")) +
  ggnewscale::new_scale_color() +
  geom_smooth(aes(color = Mother, fill = Mother), method = "lm",  linewidth = 1, se = TRUE, alpha = 0.3) +
  scale_color_manual(values = c("firebrick", "#67B2D8")) +
  scale_fill_manual(values = c("firebrick", "#67B2D8")) +
  ggnewscale::new_scale_color() +
  stat_summary(aes(x = as.numeric(gsub("M", "", Timepoint_new)), 
                   group = Timepoint_group, color = Mother), 
               fun.data = mean_cl_boot, 
               geom = "errorbar", 
               width = 0.5) +
  stat_summary(aes(x = as.numeric(gsub("M", "", Timepoint_new)),
                   group = Timepoint_group, color = Mother), 
               fun = mean, 
               geom = "point", 
               size = 1.2) +
  scale_color_manual(values = c("darkred", "darkblue")) +
  scale_x_continuous(breaks = c(0, 1, 3, 6, 12)) +
  theme_minimal() + 
  labs(x = "Infant age, months", y = "Similarity to maternal virome (1 - BC dissimilarity)") + 
  annotate("rect", xmin = 6, xmax = 12, ymin = 0.7, ymax = 0.8, 
           alpha = 0.2, fill = "firebrick") +
  annotate(geom="text", x = 9, y=0.75, label = "beta = 6.5e-3\np-value=1.1e-16", size = 2.5) +
  annotate("rect", xmin = 6, xmax = 12, ymin = 0.6, ymax = 0.7, 
           alpha = 0.2, fill = "#67B2D8") +
  annotate(geom="text", x = 9, y=0.65, label = "beta = 3.9e-3\np-value=2.0e-75", size = 2.5) +
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.position = "right")

ggsave("05.PLOTS/06.DYNAMICS/VLP_infant_converge_mom_unrelated_trajectory.pdf",
      vlp_converge_moms, "pdf", width=10, height=9, units="cm", dpi = 300)


