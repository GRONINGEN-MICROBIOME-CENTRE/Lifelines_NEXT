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
                    'med =', round(median(y), 3), '\n')
    )
  )
}
#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(tidyverse)
library(ggplot2)
#############################################################
# 2. Load Input Data
#############################################################
# cleaner:
rm(list=setdiff(ls(), "jaccard_mat"))

esmeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T, check.names = F)
esmeta <- esmeta %>%
  mutate(Timepoint_new = factor(Timepoint_new, levels = c('M1', 'M3', 'M6', 'M12', 'Mother'), ordered = T),
    secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) %>% # making the 2nd pregnancy as separate fams
  mutate(easy_ID = paste0(gsub('GS|LP', '', seq_type), '_', 
                          gsub('AM', '', FAMILYupd), '_', 
                          Family_structure, '_', 
                          Timepoint_original)) # creating an easy ID since Chiliadal samples were named sequentially

rpkm <- read.table('06.CLEAN_DATA/02.FINAL/RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_2220_samples.txt', sep='\t', header=T)
#############################################################
# 3.1 Analysis of composition similarity
#############################################################
# dissimilarity between samples (binary jaccard is probably best here)
jaccard_mat <- vegan::vegdist(t(rpkm), method = "jaccard", binary = T)

######## saving point for jaccard DISSIMILARITY
jaccard_saver_point <- as.matrix(jaccard_mat)
######## saving point for jaccard DISSIMILARITY

# similarity
jaccard_similarity <- 1 - as.matrix(jaccard_mat)

upper_tri_idx <- which(upper.tri(jaccard_similarity, diag = FALSE), arr.ind = TRUE)

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
  # 1. Unrelated in VLP
  mutate(unrelated_VLP = ifelse(seq_type_1 == "VLP" & 
                                  seq_type_2 == "VLP" & 
                                  FAMILY_1 != FAMILY_2, T, F)) %>%
  # 2. Unrelated in MGS
  mutate(unrelated_MGS = ifelse(seq_type_1 == "MGS" & 
                                  seq_type_2 == "MGS" & 
                                  FAMILY_1 != FAMILY_2, T, F)) %>%
  # 3. Unrelated between VLP and MGS
  mutate(unrelated_inter = ifelse(seq_type_1 != seq_type_2 & 
                                    FAMILY_1 != FAMILY_2, T, F)) %>%
  # 4. Same fecal samples, diff methods (same_feces_inf + same_feces_mom)
  mutate(same_feces_inter = ifelse(Universal_ID_1 == Universal_ID_2, T, F)) 

#############################################################
# 3.2 Jaccard similarity plot
#############################################################

# preparing data for plotting (elongating & categorizing)
plot_data <- dist_list %>%
  select(similarity, unrelated_VLP:same_feces_inter) %>%
  pivot_longer(cols = !similarity, 
               names_to = "Category", 
               values_to = "Is_True") %>%
  filter(Is_True == TRUE) %>%
  mutate(method = gsub(".*_", "", Category))

# Jaccard index inter:
plot_data %>%
  filter(Category == "same_feces_inter") %>%
  pull(similarity) %>%
  summary() # rounded values: 19 (IQR: 9-33)


my_comparisons <- list( c("unrelated_VLP", "unrelated_MGS"), 
                        c("same_feces_inter", "unrelated_inter"),
                        c("same_feces_inter", "unrelated_VLP"),
                        c("same_feces_inter", "unrelated_MGS")
)


p_select <- ggplot(plot_data, aes(x = Category, y = similarity, fill = method)) +
  geom_violin(alpha = 0.7, linewidth = 0.4, scale = "width") + 
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, linewidth = 0.4) +
  theme_bw() +
  scale_x_discrete(labels = c("Paired\nMGS and VLP", "Unrelated\nMGS and VLP", "Unrelated\nMGS", "Unrelated\nVLP")) +
  scale_fill_manual(labels = c("INTER", "MGS", "VLP"), values = c(MetBrewer::met.brewer("VanGogh2")[1], 
                                                                  MetBrewer::met.brewer("Kandinsky")[1], 
                                                                  MetBrewer::met.brewer("Kandinsky")[2])) +
  theme(legend.position = "right",
        legend.text = element_text(size=8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        legend.key.size = unit(0.7, 'line'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-10)) +
  labs(x = "",
       y = "Jaccard binary similarity",
       fill = "Method") +
  stat_summary(fun.data = stat_box_data,
               geom = "text",
               hjust = 0.5,
               vjust = 0.8, size = 2.5) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, 
                             label = "p.signif", 
                             method = "wilcox.test", # just for visualization
                             size=3, 
                             tip.length = 0.00,
                             step.increase = c(0.00, 0.00, 0.02, 0.04),
                             vjust = 0.001,
                             p.adjust.method = "BH")

ggsave('05.PLOTS/05.VLP_MGS/VLP_MGS_select_jaccard.pdf', 
       p_select, "pdf", width=11, height=7, units="cm", dpi = 300)

#############################################################
# 3.2 Permutation analysis for jaccard similarity
#############################################################
# reprod
set.seed(444)

# combos:

busy <- c('same_feces_inter',
          'unrelated_MGS', 
          'unrelated_VLP',
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

colnames(combos)[1:2] <- c("Comparison group 1", "Comparison group 2")

combos <- combos %>%
  mutate(across(
    where(is.character), ~ recode(.x,
             "same_feces_inter" = "Paired_MGS_and_VLP",
             "unrelated_inter" = "Unrelated_MGS_and_VLP",
             "unrelated_MGS" = "Unrelated_MGS",
             "unrelated_VLP" = "Unrelated_VLP")
  ))

# save results of permutations:
writexl::write_xlsx(combos, '07.RESULTS/Jaccard_similarity_perm_test_combinations.xlsx')



