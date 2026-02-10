setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')
#############################################################
# Checking within- and between-distances for MGS and VLP samples
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
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
library(dplyr)
library(purrr)
library(tidyverse)
library(broom)
library(skimr)
library(ggplot2)
library(UpSetR)
library(ggforce)
library(Matrix)
library(lme4)
library(lmerTest)
#############################################################
# 2. Load Input Data
#############################################################
# cleaner:
rm(list=setdiff(ls(), "jaccard_mat"))

esmeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T, check.names = F)

# creating new Timepoint factor:
esmeta$Timepoint_new <- factor(esmeta$Timepoint_new, levels = c('M1', 'M3', 'M6', 'M12', 'Mother'), ordered = T)

rpkm <- read.table('06.CLEAN_DATA/02.FINAL/RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_2220_samples.txt', sep='\t', header=T)
#############################################################
# 3.1 Analysis comparing infant vs mother metrics
#############################################################

# editing family structure to make parsing swaps easier:
dyad_esmeta <- esmeta %>%
  rename("inv_ver" = "invertebrates, vertebrates") %>%
  bind_rows(esmeta %>% filter(non_dyads == "Twin family" & Type == "M") %>% mutate(Family_structure = "M2P1")) %>% # setting mothers of twins as Mother2
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILY_2 = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) %>% # making the 2nd pregnancy as separate fams
  mutate(FAMILY_2 = if_else( grepl( "M2|K2", Family_structure ), paste0(FAMILY_2, "_T2"), FAMILY_2)) # making twin families as separate fams

# omitting some columns so that it is easy to check:
comfy <-  dyad_esmeta %>%
  select(NEXT_ID, Type, Timepoint_original, FAMILY, FAMILY_2, Universal_ID,  # Invidivual / fecal sample data
         Sequencing_ID, seq_type, dna_conc, isolation_method, # fecal sample processing details
         raw_reads, human_reads, clean_reads, reads_lost_QC, bacterial_markers_alignment_rate, sc_enrichment, # sequencing results
         contigs_0_bp, contigs_1000_bp, # assembly characteristics
         check_VLP_for_mixup, # inferred from NEXT paper
         perc_aligned_cf, vir_richness_d, N_shared_cs_d, vir_richness_cf, N_shared_cs_cf, # filtration and decontamination stats
         vir_diversity, temperate_richness, lytic_richness, temperate_perc, lytic_perc,
         temp_to_lytic_ratio, temperate_RAb, lytic_RAb, # virus sample characteristics
         PPV_fraction, PPV_abundance, dsDNA, ssDNA, RNA, genome_unclassified # taxonomy / composition
  ) %>%
  mutate(Timepoint_original = factor(Timepoint_original, levels = c("P12", "P28", "B", "M1", "M3", "M6", "M12"), ordered = T)) 

# 1. infants that have more viruses than their mothers:
rich_by_fam <- dyad_esmeta %>%
  filter(seq_type == "VLP") %>% # rare in MGS: 4 cases when an infant has more viruses than their mother
  arrange(Timepoint_new) %>%
  pivot_wider(id_cols = FAMILY_2,
              names_from = Timepoint_new,
              values_from = vir_richness_cf) %>%
  filter(!if_all(c("M1", "M3", "M6"), is.na)) %>%
  mutate(richness_discordance = if_any(c("M1", "M3", "M6"), ~!is.na(.x) & .x >= Mother)) %>%
  mutate(richness_discordance = ifelse(!is.na(richness_discordance), richness_discordance, FALSE )) %>%
  mutate(too_high = if_any(c("M1", "M3"), ~!is.na(.x) & .x >= median(esmeta$vir_richness_cf[esmeta$seq_type == "VLP" & esmeta$Type == "M"]))) # higher than maternal median

# 2. infants that have lower temperate phage abundance compared to their mothers:
temp_by_fam <- dyad_esmeta %>%
  filter(seq_type == "VLP") %>%
  arrange(Timepoint_new) %>%
  pivot_wider(id_cols = FAMILY_2,
              names_from = Timepoint_new,
              values_from = temperate_RAb) %>%
  filter(!if_all(c("M1", "M3"), is.na)) %>%
  mutate(temperate_discordance = if_any(c("M1", "M3"), ~!is.na(.x) & .x < Mother)) %>%
  mutate(temperate_discordance = ifelse(!is.na(temperate_discordance), temperate_discordance, FALSE ))

# 3. sus fams that wher einfant have too high richness & too low temperate phage abundance
check_both <- rich_by_fam %>%
  select(c(FAMILY_2, richness_discordance, too_high)) %>%
  left_join(temp_by_fam %>% select(c(FAMILY_2, temperate_discordance))) %>%
  filter(richness_discordance == T & temperate_discordance == T) # 5 fams

########################### sample swap!
# in the first iteration (w previous version of esmeta, Chiliadal_meta_VLP_MGS_matched_v04_suppl_w_virmetrics.txt):
# after manual exploration of candidates, FAM0941 needs to be swapped (also confirmed by NEXT paper)

########################### observation:
# so far, an excessive number of contigs is a great predictor for the higher sharedness of viruses to NCs


# 4. infants have more viruses at M1 than at M3 or other

# MGS
rich_by_infant_MGS <- dyad_esmeta %>%
  filter(seq_type == "MGS") %>%
  arrange(Timepoint_new) %>%
  pivot_wider(id_cols = FAMILY_2,
              names_from = Timepoint_new,
              values_from = vir_richness_cf) %>%
  filter(!if_all(c("M1", "M3"), is.na)) %>%
  mutate(M1_too_high = if_any(c("M3", "M6", "M12"), ~!is.na(.x) & .x <= M1)) %>% # 72 cases
  mutate(M3_too_high = if_any(c("M6", "M12"), ~!is.na(.x) & .x <= M3)) %>% # 47 cases
  mutate(M1_too_high = ifelse(is.na(M1_too_high), FALSE, M1_too_high)) %>%
  mutate(M3_too_high = ifelse(is.na(M3_too_high), FALSE, M3_too_high)) 

# VLP
rich_by_infant_VLP <- dyad_esmeta %>%
  filter(seq_type == "VLP") %>%
  arrange(Timepoint_new) %>%
  pivot_wider(id_cols = FAMILY_2,
              names_from = Timepoint_new,
              values_from = vir_richness_cf) %>%
  filter(!if_all(c("M1", "M3"), is.na)) %>%
  mutate(M1_too_high = if_any(c("M3", "M6", "M12"), ~!is.na(.x) & .x <= M1)) %>%
  mutate(M3_too_high = if_any(c("M6", "M12"), ~!is.na(.x) & .x <= M3)) %>%
  mutate(M1_too_high = ifelse(is.na(M1_too_high), FALSE, M1_too_high)) %>% # 90 cases (excluding cases when M6 or M12 too low)
  mutate(M3_too_high = ifelse(is.na(M3_too_high), FALSE, M3_too_high)) %>% # 55 cases (excluding cases when M6 or M12 too low)
  mutate(M6_or_M12_to_low = ifelse(M6 < 100 | M12 < 100, T, F)) %>% # 43 cases
  mutate(M6_or_M12_to_low = ifelse(!is.na(M6_or_M12_to_low), M6_or_M12_to_low, F))

# let's check if there is some biology:
rich_by_infant_VLP <- rich_by_infant_VLP %>%
  mutate(M1_too_high_MGS = ifelse(FAMILY_2 %in% rich_by_infant_MGS$FAMILY_2[rich_by_infant_MGS$M1_too_high == T], T, F)) %>% # 36 out of 90 M1 VLP cases
  mutate(M3_too_high_MGS = ifelse(FAMILY_2 %in% rich_by_infant_MGS$FAMILY_2[rich_by_infant_MGS$M3_too_high == T], T, F)) # 13 out of 55 cases

# Number of M1 & Number of M3:
dim(rich_by_infant_VLP[!is.na(rich_by_infant_VLP$M1), ]) # 179
dim(rich_by_infant_VLP[!is.na(rich_by_infant_VLP$M3), ]) # 300

# UpSet plot showing mainly cases when no explanation of too high M1 or M3:
upset_data <- rich_by_infant_VLP %>%
  select(all_of(c("M1_too_high", "M3_too_high", "M6_or_M12_to_low", "M1_too_high_MGS",  "M3_too_high_MGS"))) %>%
  mutate(across(everything(), as.integer)) %>%
  as.data.frame()


upset_all <- upset(upset_data, order.by = "freq", sets.bar.color = "#C00000", 
                   number.angles = 20,
                   sets.x.label = "N vOTUs", scale.sets = "identity",
                   text.scale = c(1, 1, 1, 0.7, 1, 1))

png('05.PLOTS/08.SAMPLE_SWAP_EXPLORE/VLP_M1_and_M3_explore.png', width=16, height=12, units="cm", res = 300)
upset_all
dev.off()


#############################################################
# 3.2 Analysis of samples excluded from big gut next paper
#############################################################

# 1. manual check for mix up

potmixups <- unique(comfy$FAMILY[comfy$check_VLP_for_mixup == T])

check_bgnp <- comfy %>%
  filter(FAMILY %in% potmixups)

# none too suspicious
# clean up:
rm(list = c("upset_data", "upset_all", 
            "rich_by_infant_VLP", "rich_by_infant_MGS", 
            "potmixups", "check_bgnp"))
#############################################################
# 3.3 Analysis of composition similarity
#############################################################

# dissimilarity between samples (binary jaccard is probably best here)
jaccard_mat <- vegan::vegdist(t(rpkm), method = "jaccard", binary = T)

######## saving point for jaccard DISSIMILARITY
jaccard_saver_point <- as.matrix(jaccard_mat)
######## saving point for jaccard DISSIMILARITY

# similarity
jaccard_similarity <- 1 - as.matrix(jaccard_mat)

upper_tri_idx <- which(upper.tri(jaccard_similarity, diag = FALSE), arr.ind = TRUE)

esmeta <- esmeta %>%
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) %>% # making the 2nd pregnancy as separate fams
  mutate(easy_ID = paste0(gsub('GS|LP', '', seq_type), '_', 
                          gsub('AM', '', FAMILYupd), '_', 
                          Family_structure, '_', 
                          Timepoint_original)) # creating an easy ID since Chiliadal samples were named sequentially

row.names(jaccard_similarity) <- esmeta$easy_ID[match(row.names(jaccard_similarity) , esmeta$Sequencing_ID)]
colnames(jaccard_similarity) <- esmeta$easy_ID[match(colnames(jaccard_similarity) , esmeta$Sequencing_ID)]

# 1. paired distances between different samples:

dist_list <- data.frame(
  sample_1 = rownames(jaccard_similarity)[upper_tri_idx[, 1]],
  sample_2 = colnames(jaccard_similarity)[upper_tri_idx[, 2]],
  similarity = jaccard_similarity[upper_tri_idx]
) %>%
  # has richness 0:
  filter(sample_1 != "M_F0257_K1P1_M6") %>%
  filter(sample_2 != "M_F0257_K1P1_M6") %>%
  left_join(esmeta %>% select(easy_ID, 
                               Universal_ID,
                               Modified_NEXT_ID_without_preg_number,
                               FAMILY, 
                               FAMILYupd,
                               Type, 
                              seq_type,
                              Timepoint_new), by = c("sample_1" = "easy_ID" )) %>%
  left_join(esmeta %>% select(easy_ID, 
                               Universal_ID,
                               Modified_NEXT_ID_without_preg_number,
                               FAMILY,
                               FAMILYupd,
                               Type, seq_type,
                              Timepoint_new), by = c("sample_2" = "easy_ID"), suffix = c("_1", "_2")) %>%
  select(sample_1, sample_2, similarity, sort(colnames(.))) %>%
  # 1. Same infant fecal aliquote, different sequencing types
  mutate(same_feces_inf = ifelse(Universal_ID_1 == Universal_ID_2 &
                                   Type_1 == "K" &
                                   Type_2 == "K", T, F)) %>% 
  # 2. Same maternal fecal aliquote, different sequencing types
  mutate(same_feces_mom = ifelse(Universal_ID_1 == Universal_ID_2 &
                                   Type_1 == "M" &
                                   Type_2 == "M", T, F)) %>%
  # 3. Within-infant similarities in VLP
  mutate(same_infant_VLP = ifelse(seq_type_1 == "VLP" & 
                                    seq_type_2 == "VLP" & 
                                    Type_1 == "K" &
                                    Type_2 == "K" &
                                    Modified_NEXT_ID_without_preg_number_1 == Modified_NEXT_ID_without_preg_number_2, T, F)) %>%
  # 4. Within-mom first vs second pregnancy comparison in VLP
  mutate(same_mom_VLP = ifelse(seq_type_1 == "VLP" & 
                                 seq_type_2 == "VLP" & 
                                 Type_1 == "M" &
                                 Type_2  == "M" &
                                 Modified_NEXT_ID_without_preg_number_1 == Modified_NEXT_ID_without_preg_number_2, T, F)) %>%
  # 5. Within-infant similarities in MGS
  mutate(same_infant_MGS = ifelse(seq_type_1 == "MGS" & 
                                    seq_type_2 == "MGS" & 
                                    Type_1 == "K" &
                                    Type_2 == "K" &
                                    Modified_NEXT_ID_without_preg_number_1 == Modified_NEXT_ID_without_preg_number_2, T, F)) %>%
  # 6. Within-mom first vs second pregnancy comparison in MGS
  mutate(same_mom_MGS = ifelse(seq_type_1 == "MGS" & 
                                 seq_type_2 == "MGS" & 
                                 Type_1 == "M" &
                                 Type_2 == "M" &
                                 Modified_NEXT_ID_without_preg_number_1 == Modified_NEXT_ID_without_preg_number_2, T, F)) %>%
  # 7. Between-twins in VLP
  mutate(twins_VLP = ifelse(seq_type_1 == "VLP" & 
                              seq_type_2 == "VLP" & 
                              FAMILYupd_1 == FAMILYupd_2 &
                              Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                              Type_1 == "K" &
                              Type_2 == "K", T, F)) %>%
  # 8. Between-twins in MGS
  mutate(twins_MGS = ifelse(seq_type_1 == "MGS" & 
                              seq_type_2 == "MGS" & 
                              FAMILYupd_1 == FAMILYupd_2 &
                              Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                              Type_1 == "K" &
                              Type_2 == "K", T, F)) %>%
  # 9. Between-siblings in VLP
  mutate(siblings_VLP = ifelse(seq_type_1 == "VLP" & 
                                 seq_type_2 == "VLP" & 
                                 FAMILY_1 == FAMILY_2 &
                                 FAMILYupd_1 != FAMILYupd_2 &
                                 Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                                 Type_1 == "K" &
                                 Type_2 == "K", T, F)) %>%
  # 10. Between-siblings in MGS
  mutate(siblings_MGS = ifelse(seq_type_1 == "MGS" & 
                                 seq_type_2 == "MGS" & 
                                 FAMILY_1 == FAMILY_2 &
                                 FAMILYupd_1 != FAMILYupd_2 &
                                 Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                                 Type_1 == "K" &
                                 Type_2 == "K", T, F)) %>%
  # 11. Mother-infant in VLP
  mutate(same_dyad_VLP = ifelse(seq_type_1 == "VLP" & 
                                    seq_type_2 == "VLP" & 
                                    FAMILYupd_1 == FAMILYupd_2 &
                                    Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                                    twins_VLP != T, T, F)) %>%
  # 12. Mother-infant in MGS
  mutate(same_dyad_MGS = ifelse(seq_type_1 == "MGS" & 
                                    seq_type_2 == "MGS" & 
                                    FAMILYupd_1 == FAMILYupd_2 &
                                    Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                                    twins_MGS != T, T, F)) %>%
  # 13. Unrelated in VLP
  mutate(unrelated_VLP = ifelse(seq_type_1 == "VLP" & 
                                  seq_type_2 == "VLP" & 
                                  FAMILY_1 != FAMILY_2, T, F)) %>%
  # 14. Unrelated in MGS
  mutate(unrelated_MGS = ifelse(seq_type_1 == "MGS" & 
                                  seq_type_2 == "MGS" & 
                                  FAMILY_1 != FAMILY_2, T, F)) %>%
  # 15. Within-infant similarities between MGS and VLP
  mutate(same_infant_inter = ifelse(seq_type_1 != seq_type_2 & 
                                      Universal_ID_1 != Universal_ID_2 &
                                      Type_1 == "K" &
                                      Type_2 == "K" &
                                      Modified_NEXT_ID_without_preg_number_1 == Modified_NEXT_ID_without_preg_number_2, T, F)) %>%
  # 16. Between-twins similarities between MGS and VLP
  mutate(twins_inter = ifelse(seq_type_1 != seq_type_2 &
                                FAMILYupd_1 == FAMILYupd_2 &
                                Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                                Type_1 == "K" &
                                Type_2 == "K", T, F)) %>%
  # 17. Mother-infant between MGS and VLP
  mutate(same_dyad_inter = ifelse(seq_type_1 != seq_type_2 & 
                                      FAMILYupd_1 == FAMILYupd_2 &
                                      Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                                      twins_inter != T
                                    , T, F)) %>%
  # 18. Unrelated between VLP and MGS
  mutate(unrelated_inter = ifelse(seq_type_1 != seq_type_2 & 
                                    FAMILY_1 != FAMILY_2, T, F)) %>%
  # 19. Same fecal samples, diff methods (same_feces_inf + same_feces_mom)
  mutate(same_feces_inter = ifelse(Universal_ID_1 == Universal_ID_2, T, F)) 

# preparing data for plotting (elongating & categorizing)
plot_data <- dist_list %>%
  select(similarity, same_feces_inf:same_feces_inter) %>%
  pivot_longer(cols = !similarity, 
               names_to = "Category", 
               values_to = "Is_True") %>%
  filter(Is_True == TRUE) %>%
  mutate(method = gsub(".*_", "", Category)) %>%
  mutate(method = if_else(method %in% c('inf', 'mom'), "inter", method) )

# Visualize all:
p_all <- ggplot(plot_data, aes(x = Category, y = similarity, fill = method)) +
  geom_violin(alpha = 0.7, linewidth = 0.4) + 
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, linewidth = 0.4) +
  theme_minimal() +
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

png('05.PLOTS/08.SAMPLE_SWAP_EXPLORE/VLP_MGS_all_jaccard.png', width=27, height=15, units="cm", res = 300)
p_all
dev.off()

# visualize sensible:

my_comparisons <- list( c("unrelated_VLP", "unrelated_MGS"), 
                        c("unrelated_VLP", "unrelated_inter"), 
                        c("unrelated_MGS", "unrelated_inter"),
                        c("same_feces_inter", "unrelated_VLP"),
                        c("same_feces_inter", "unrelated_MGS"),
                        c("same_feces_inter", "unrelated_inter")
                        )


p_select <- plot_data %>%
  filter(Category %in% c('same_feces_inter',
                         'unrelated_MGS', 
                         'unrelated_VLP',
                         'unrelated_inter'
                         )) %>%
  ggplot(aes(x = Category, y = similarity, fill = method)) +
  geom_violin(alpha = 0.7, linewidth = 0.4, scale = "width") + 
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, linewidth = 0.4) +
  theme_bw() +
  scale_x_discrete(labels = c("Paired\nMGS and VLP", "Unrelated\nMGS and VLP", "Unrelated\nMGS", "Unrelated\nVLP")) +
  scale_fill_manual(labels = c("INTER", "MGS", "VLP"), values = c(MetBrewer::met.brewer("VanGogh2")[1], 
                                                                  MetBrewer::met.brewer("Kandinsky")[1], 
                                                                  MetBrewer::met.brewer("Kandinsky")[2])) +
  theme(legend.position = "right") +
  labs(x = "",
       y = "Jaccard binary similarity",
       fill = "Method") +
  stat_summary(
    fun.data = stat_box_data,
    geom = "text",
    hjust = 0.5,
    vjust = 0.8, size = 2.5
  ) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size=3, p.adjust.method = "BH")

ggsave('05.PLOTS/05.VLP_MGS/VLP_MGS_select_jaccard.png', 
       p_select,  "png", width=17, height=12, units="cm", dpi = 300)


#############################################################
# 3.4 Permutation analysis for jaccard similarity
#############################################################
# reprod
set.seed(444)

# n diff possible (and sensible) combos:

busy <- c('same_feces_inter',
          'same_infant_MGS', 
          'same_infant_VLP', 
          'unrelated_MGS', 
          'unrelated_VLP',
          'same_dyad_MGS',
          'same_dyad_VLP',
          'unrelated_inter')

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
    if ( (feature1 == "unrelated_MGS" & feature2 == "unrelated_VLP") | (feature2 == "unrelated_MGS" & feature1 == "unrelated_VLP") ) {
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

#############################################################
# 3.5 Exploring why same_feces_inter is sometimes very low
#############################################################
# max unrelated distance
noise_lookup <- dist_list %>%
  filter(unrelated_inter == TRUE) %>%
  bind_rows(
    dist_list %>%
      filter(unrelated_inter == TRUE) %>%
      rename(sample_1 = sample_2, sample_2 = sample_1)
  ) %>%
  group_by(sample_1) %>%
  slice_max(similarity, n = 1, with_ties = FALSE) %>%
  select(sample_1, closest_unrelated = sample_2, max_unrelated_sim = similarity)

# distance between paired samples  
flagged_cases <-   dist_list %>%
    filter(same_feces_inter == TRUE) %>%
    left_join(noise_lookup) %>%
    select(sample_1, sample_2, similarity, closest_unrelated, max_unrelated_sim, Universal_ID_1)

# separate VLP and MGS are built here as well because I do not bother
stranger_test <- ggplot(flagged_cases, aes(x = reorder(sample_1, similarity), y = similarity)) +
  geom_point(aes(color = "Paired")) + 
  geom_errorbar(aes(ymin = similarity, ymax = max_unrelated_sim), linetype = "dotted") +
  geom_point(aes(y = max_unrelated_sim, color = "Max unrelated"), shape = 4) +
  scale_color_manual(breaks = c("Paired", "Max unrelated"), 
                     values = c("red", "blue"),
                     labels = c("Paired", "Max unrelated")) +
  coord_flip() +
  labs(title = "Jaccard similarity between paired fecal samples vs closest unrelated (VLP) sample",
       x = "Sample ID",
       y = "Jaccard Similarity") +
  theme_bw() + 
  theme(axis.text.y = element_text(size=1),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("05.PLOTS/08.SAMPLE_SWAP_EXPLORE/Stranger_test_VLP_MGS.png", stranger_test, width = 8 , height = 13, dpi = 300)

# sorting by most:
flagged_cases <- flagged_cases %>%
  arrange(similarity, -max_unrelated_sim) #%>%
  #filter(similarity < max_unrelated_sim) #considering only these cases

# manual checking, how to make it better?

# manual inspection: 
fam1 <- esmeta$easy_ID[esmeta$FAMILYupd %in% c('FAM0566', 'FAM0386')]
becker <- jaccard_similarity[#row.names(jaccard_similarity) %in% fam1,
  , colnames(jaccard_similarity) %in% fam1]

#### conclusion these two MGS samples were swapped during extracted DNA aliquoting (see box position)
#esmeta$easy_ID[esmeta$Sequencing_ID == "C06F056645I3"] <- "M_F0386_K1P1_M6"
#esmeta$easy_ID[esmeta$Sequencing_ID == "C06F038645I4"] <- "M_F0566_K1P1_M6"

### no time to follow it up further

#############################################################
# 3.6 Does within-sample similarity depend on tech specs?
#############################################################
potimpact <- c("temperate_RAb", 
                     "lytic_RAb", 
                     "ssDNA", 
                     "sc_enrichment",
                     "vir_diversity",
                     "contigs_1000_bp",
                     "perc_aligned_cf",
                     "N_shared_cs_cf")

pon <- flagged_cases %>%
  left_join(esmeta %>% 
              filter(seq_type == "VLP") %>% 
              select(all_of(potimpact), 
                     NEXT_ID, 
                     Universal_ID, 
                     exact_age_days_at_collection,
                     Timepoint_new), by = c("Universal_ID_1" = "Universal_ID"))



results <- map_dfr(c('similarity', 'max_unrelated_sim'), function(sim) {
  map_dfr(potimpact, function(virvar) {
    
    print(paste("METRIC:", sim, "VIRVAR:", virvar))
    
    formula <- as.formula(paste(sim, "~",  virvar, "+ exact_age_days_at_collection + (1 | NEXT_ID)"))
    
    model <- lmer(
      formula,
      REML = FALSE,
      data = pon
    )
    
    anova_res <- car::Anova(model, type = "II", test.statistic = "Chisq") %>%
      tidy() %>%
      filter(term == virvar)
    
    eff_res <- effectsize::eta_squared(model, partial = F) %>%
      as.data.frame() %>%
      filter(Parameter == virvar)
    
    data.frame(
      virvar_VLP = virvar,
      sim_metric = sim,
      p_value = anova_res$p.value,
      eta_sq = eff_res$Eta2
    )
    
  }) 
}) 

results <- results %>%
  mutate(
    neglog10_p = -log10(p_value + min(results[results$p_value!=0,]$p_value))
  ) %>%
  mutate(neglog10_p = ifelse(p_value >= 0.05, NA, neglog10_p))

png("05.PLOTS/08.SAMPLE_SWAP_EXPLORE/Paired_vs_unpaired_Jaccard_and_viral_features.png", width = 12, height = 10, units = "cm", res = 300)

results %>%
ggplot(aes(x = virvar_VLP, y = sim_metric)) +
  geom_tile(color = "grey", aes(fill = eta_sq)) +
  geom_point(aes(size = neglog10_p), shape = 21, colour = "black", alpha = 1, fill = NA) +
  scale_fill_gradient2(name = "Eta sq",
                       low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  scale_size_continuous(name = "-log10(p)", range = c(0, 6)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) + 
   labs(x = "VLP sample feature", y = "Jaccard similarity to") +
   scale_y_discrete(labels = c("closest\nunrelated", "paired"))
dev.off() 


png('05.PLOTS/08.SAMPLE_SWAP_EXPLORE/VLP_MGS_jaccard.png', width=16, height=12, units="cm", res = 300)
ggplot(pon, aes(Timepoint_new, similarity)) +
  geom_violin(fill = "#DBE2EF", color = "#112D4E") +
  geom_sina(color = "#3F72AF", alpha = 0.7, size = 1) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = NA, color = "#112D4E") +
  theme_bw() +
  labs(x = "Timepoint", y = "Jaccard index", title = "Pairwise comparison of VLP and MGS samples")
dev.off()

# to save the flagged cases: 
flagged_cases_saver <- flagged_cases %>%
  left_join(esmeta %>% filter(seq_type == "VLP") %>% select(Universal_ID, Sequencing_ID), by = c("Universal_ID_1" = "Universal_ID")) %>%
  rename("Sequencing_ID_VLP" = "Sequencing_ID") %>%
  left_join(esmeta %>% filter(seq_type == "MGS") %>% select(Universal_ID, Sequencing_ID), by = c("Universal_ID_1" = "Universal_ID")) %>%
  rename("Sequencing_ID_MGS" = "Sequencing_ID")
#############################################################
# X. OUTPUT
#############################################################
# Jaccard dist matrix (otherwise too long to calculate)
write.table(jaccard_saver_point, "06.CLEAN_DATA/Intermediate/Jaccard_dissimilarity_RPKM_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_2220_samples.txt", sep='\t', quote=F)

# save flagged cases just in case:
write.table(flagged_cases_saver, "06.CLEAN_DATA/Intermediate/Jaccard_simmilarity_MGS_vs_VLP.txt", sep='\t', quote = F, row.names = F)

# save results of permutations:
write.table(combos, "07.RESULTS/Jaccard_similarity_perm_test_combinations.txt", sep='\t', quote = F, row.names = F)

writexl::write_xlsx(combos, '07.RESULTS/Jaccard_similarity_perm_test_combinations.xlsx')

