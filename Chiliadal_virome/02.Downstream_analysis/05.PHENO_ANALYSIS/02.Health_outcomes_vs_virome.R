setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")
#############################################################
# Here we explore whether virome is associated with health outcomes
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
virome_shaper <- function(pheno, df, dist){
  
  pheno_df <- df %>%
    filter(!is.na(!!sym(pheno)))
  
  dister <- dist[colnames(dist) %in% pheno_df$Sequencing_ID,
                 row.names(dist) %in% pheno_df$Sequencing_ID]
  
  pcoa_res <- cmdscale(as.dist(dister), k = 2, eig = TRUE)
  
  meh <- pcoa_res$points
  pcoa_plot_data <- meh %>%
    as.data.frame() %>%
    rownames_to_column(var = "Sequencing_ID") %>%
    left_join(pheno_df)
  
  centroids <- pcoa_plot_data %>%
    group_by(!!sym(pheno), Timepoint_new) %>%
    summarise(V1 = median(V1),
              V2 = median(V2), .groups = "drop") %>%
    as.data.frame()
  
  pvar <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)
  
  return(list("pcoa_plot_data" = pcoa_plot_data, 
              "centroids" = centroids, 
              "pvar" = pvar))
  
}

# 0. Updated from 01.Phenotypes_vs_viral_features.R
# 1. function is rigid to the name of the column containing sample ID:
# it has to be Sequencing_ID
# 2. check how it pulls the beta and p-value IF you decide to correct for
# anything else besides time
running_LMM_upd <- function(phenos_list, # vector of phenotypes to test
                            metadata_w_phenos, # metadata containing these phenos as columns
                            RPKM_table, # filtered by prevalence table
                            pseudocount, # calculated pseudocount, no default!
                            correct_for = NULL # vector of factors to correct for, NULL by default
) { 
  
  map_dfr(phenos_list, function(pheno){
    
    if (is.null(correct_for)) {
      
      pheno_df <- metadata_w_phenos %>%
        filter(!is.na(!!sym(pheno)))
      
    } else {
      
      pheno_df <- metadata_w_phenos %>%
        filter(!is.na(!!sym(pheno))) %>%
        filter(if_all(all_of(correct_for), ~ !is.na(.x)))
      
    }
    
    smaller <- RPKM_table[,colnames(RPKM_table) %in% pheno_df$Sequencing_ID]
    smaller <- smaller[rowSums(smaller) > 0,]
    
    feature_set <- row.names(smaller)
    
    smaller <- smaller %>%
      rownames_to_column(var = "New_CID") %>%
      pivot_longer(!New_CID) %>%
      pivot_wider(names_from = New_CID) %>%
      left_join(pheno_df, by = c("name" = "Sequencing_ID"))
    
    right_side <- c("Timepoint_new", correct_for, pheno, "(1|NEXT_ID)")
    
    map_dfr(feature_set, function(feature) {
      
      #F1 <- as.formula(paste0("log(", feature, "+", pseudocount/2, ")", " ~ ", pheno, " + Timepoint_new + (1|NEXT_ID)"))
      
      F1 <- as.formula(paste0("log(", feature, "+", pseudocount/2, ")", " ~ ", paste(right_side, collapse = " + ")))
      
      model <- lmer(F1, 
                    data = smaller,
                    REML = F)  
      
      summary(model)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "rowname") %>%
        filter(grepl(pheno, rowname)
        ) %>%
        mutate(feature = feature,
               n_nzero_feature = sum(smaller %>% 
                                       filter(!!sym(feature) != 0) %>%
                                       pull(!!sym(feature)) > 0),
               n_non_NA_pheno = nrow(smaller)  )
      
    })
  })
}

# very rigid function for everything; I made it to clean up this Rscript
running_adonis <- function(phenotypes, # character vector of phenotypes
                           metadata_w_phenos, # df, metadata containing NEXT_ID, phenotype columns, etc
                           distance_matrix, # df, distance matrix (if Sequencing_ID and matrix names are disorderd, returns NULL)
                           run_type, # character, arguments: "cross" or anything else; cross-sectional or longitudinal type of phenotypes
                           run_by, # arguments: "terms", "margin"
                           perms = 999, # number of permutations, 999 by default
                           correct_for = NULL # vector of factors to correct for, NULL by default
) {
  map_dfr(phenotypes, 
          function(pheno){
            
            if (is.null(correct_for)) {
              
              pheno_df <- metadata_w_phenos %>%
                filter(!is.na(!!sym(pheno)))
              
            } else {
              
              pheno_df <- metadata_w_phenos %>%
                filter(!is.na(!!sym(pheno))) %>%
                filter(if_all(all_of(correct_for), ~ !is.na(.x)))
              
            }
            
            dister <- distance_matrix[colnames(distance_matrix) %in% pheno_df$Sequencing_ID,
                                      row.names(distance_matrix) %in% pheno_df$Sequencing_ID]
            
            right_side <- c("Timepoint_new", correct_for, pheno)
            
            F1 <- as.formula(paste0("dister ~ ", paste(right_side, collapse = " + ")))
            
            if (!identical(pheno_df$Sequencing_ID, colnames(dister))) return(NULL)
            
            res <-  if (run_type == "cross") {
              adonis2(F1,
                      data=pheno_df,
                      by = run_by,
                      permutations = how(blocks = pheno_df$NEXT_ID, nperm = perms),
                      parallel = 4) 
            } else {
              adonis2(F1,
                      data=pheno_df,
                      by = run_by,
                      strata = pheno_df$NEXT_ID,
                      permutations = perms,
                      parallel = 4) 
            }
            
            res %>%
              as.data.frame() %>%
              rownames_to_column(var = "rowname") %>%
              filter(rowname == pheno) %>%
              mutate(N_permutations = perms) %>%
              mutate(N_samples = nrow(pheno_df)) %>%
              mutate(corrected_for = ifelse(is.null(correct_for), "none", paste(correct_for, collapse = ", "))  )
            
          })
  
}

# another very rigid function for everything; I also made it to clean up this Rscript
# always corrected for perc_aligned_cf, reads_lost_QC, sequencing_batch, sc_enrichment
running_virvar_corrected <- function(phenos_to_run, # character vector of phenotypes
                                     metadata_w_phenos, # df, metadata containing NEXT_ID, phenotype columns, virvars
                                     vivars_to_run, # virvars to run
                                     correct_for = NULL # vector of factors to correct for, NULL by default
){
  
  map_dfr(phenos_to_run,
          function(pheno){
            
            if (is.null(correct_for)) {
              
              pheno_df <- metadata_w_phenos %>%
                filter(!is.na(!!sym(pheno)))
              
            } else {
              
              pheno_df <- metadata_w_phenos %>%
                filter(!is.na(!!sym(pheno))) %>%
                filter(if_all(all_of(correct_for), ~ !is.na(.x)))
              
            }
            
            right_side <- c("Timepoint_new", 
                            correct_for, 
                            c("perc_aligned_cf", "reads_lost_QC", "sequencing_batch", "sc_enrichment"),
                            pheno,
                            "(1|NEXT_ID)")
            
            map_dfr(vivars_to_run, function(virvar) {
              
              pseudocount <- virvar_pseudocounts$value[virvar_pseudocounts$name == virvar]
              
              if (virvar == "vir_diversity") {
                F1 <- as.formula(paste0(virvar, " ~ ", paste(right_side, collapse = " + ")))
              } else {
                F1 <- as.formula(paste0("log(", virvar, " + ", pseudocount/2, ") ~ ", paste(right_side, collapse = " + ")))
              }
              
              model <- lmer(F1, 
                            data = pheno_df,
                            REML = F)  
              
              summary(model)$coefficients %>%
                as.data.frame() %>%
                rownames_to_column(var = "rowname") %>%
                filter(grepl(pheno, rowname) 
                ) %>%
                mutate(virvar = virvar,
                       n_nzero_virvar = sum(pheno_df %>% 
                                              filter(!!sym(virvar) != 0) %>%
                                              pull(!!sym(virvar)) > 0),
                       n_non_NA_pheno = nrow(pheno_df) ) %>%
                mutate(corrected_for = ifelse(is.null(correct_for), "none", paste(correct_for, collapse = ", ")) )
            })
          })
  
}

# another very rigid function for everything; I also made it to clean up this Rscript
# always corrected for perc_aligned_cf, reads_lost_QC, sequencing_batch, sc_enrichment
running_virvar_corrected_timepoint <- function(phenos_to_run, # character vector of phenotypes
                                               metadata_w_phenos, # df, metadata containing NEXT_ID, phenotype columns, virvars
                                               vivars_to_run, # virvars to run
                                               correct_for = NULL # vector of factors to correct for, NULL by default
){
  
  map_dfr(phenos_to_run,
          function(pheno){
            
            if (is.null(correct_for)) {
              
              pheno_df <- metadata_w_phenos %>%
                filter(!is.na(!!sym(pheno)))
              
            } else {
              
              pheno_df <- metadata_w_phenos %>%
                filter(!is.na(!!sym(pheno))) %>%
                filter(if_all(all_of(correct_for), ~ !is.na(.x)))
              
            }
            
            right_side <- c(correct_for, 
                            c("perc_aligned_cf", "reads_lost_QC", "sequencing_batch", "sc_enrichment"),
                            pheno)
            
            map_dfr(unique(pheno_df$Timepoint_new), function(timepoint) {
              
              pheno_df <- pheno_df %>%
                filter(Timepoint_new == timepoint)
              
              map_dfr(vivars_to_run, function(virvar) {
                
                pseudocount <- virvar_pseudocounts$value[virvar_pseudocounts$name == virvar]
                
                if (length(pseudocount) == 0) {
                  pseudocount <- min_rpkm
                }
                
                if (virvar == "vir_diversity") {
                  F1 <- as.formula(paste0(virvar, " ~ ", paste(right_side, collapse = " + ")))
                } else {
                  F1 <- as.formula(paste0("log(", virvar, " + ", pseudocount/2, ") ~ ", paste(right_side, collapse = " + ")))
                }
                
                model <- lm(F1, 
                            data = pheno_df)  
                
                summary(model)$coefficients %>%
                  as.data.frame() %>%
                  rownames_to_column(var = "rowname") %>%
                  filter(grepl(pheno, rowname) 
                  ) %>%
                  mutate(virvar = virvar,
                         timepoint = timepoint,
                         n_nzero_virvar = sum(pheno_df %>% 
                                                filter(!!sym(virvar) != 0) %>%
                                                pull(!!sym(virvar)) > 0),
                         n_non_NA_pheno = nrow(pheno_df) ) %>%
                  mutate(corrected_for = ifelse(is.null(correct_for), "none", paste(correct_for, collapse = ", ")) )
              })
            })
          })
  
}
#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(MetBrewer)
library(vegan)
library(patchwork)
#############################################################
# 2. Load Input Data
#############################################################
# phenotype selection:
phenos_selected <- readxl::read_xlsx("07.RESULTS/Adonis2_betadisper_VLP_virome.xlsx")

# metadatas:
smeta_w_phenos <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_phenotypes.txt', sep='\t', header=T)
smeta_w_phenos <- smeta_w_phenos %>%
  mutate(Timepoint_new = factor(Timepoint_new, levels = c("M1", "M3", "M6", "M12"), ordered = T))

# abundance tables etc:
# VLP
VLP <- read.table("06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt", sep='\t', header=T)
VLP <- VLP[,colnames(VLP) %in% smeta_w_phenos$Sequencing_ID]
VLP <- VLP[rowSums(VLP) > 0,]

# ordering smeta same way as VLP:
smeta_w_phenos <- smeta_w_phenos[match(colnames(VLP), smeta_w_phenos$Sequencing_ID), ]

# transform:
VLP_RAB <- as.data.frame(t(as.data.frame(t(VLP)/colSums(VLP))))
VLP_dist <- vegdist(t(VLP_RAB), method = "bray", parallel = 6)
dist_matrix <- as.matrix(VLP_dist)

# # virus metadata
ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)
ETOF_vOTUr <- ETOF_vOTUr %>%
  separate(
    col    = tax_ictv_aai,
    into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
    sep    = ";",
    fill   = "right",
    remove = FALSE
  ) %>%
  separate(
    col    = Host_taxonomy,
    into = paste0('Host_', c('Domain', 'Phylum',  'Class', 'Order', 'Family',  'Genus')),
    sep    = ";",
    fill   = "right",
    remove = FALSE
  )

# results of running virvars against shaping phenos:
virvars_to_correct_for <- readxl::read_xlsx("07.RESULTS/LMM_virVars_shaping_phenos.xlsx")
virvars_to_correct_for <- virvars_to_correct_for %>%
  filter(FDR < 0.05)
#############################################################
# 3.0 Analysis: composition vs health outcomes
#############################################################
adonis_cross <- running_adonis(phenos_selected %>%
                            filter(type == "cross-sectional" & analysis == "health_outcomes") %>%
                            pull(new_name),
                          smeta_w_phenos,
                          dist_matrix,
                          "cross",
                          "margin")

adonis_long <- running_adonis(phenos_selected %>%
                            filter(type == "longitudinal" & analysis == "health_outcomes") %>%
                            pull(new_name),
                          smeta_w_phenos,
                          dist_matrix,
                          "long",
                          "margin")

adonis_results <- adonis_cross %>%
  bind_rows(adonis_long) %>%
  mutate(FDR = p.adjust(`Pr(>F)`, "BH")) %>%
  ungroup()

## corrected for feeding and delivery modes
correct_factors <- c("infant_ffq_feeding_mode_simple", "birth_deliverybirthcard_mode_binary")

adonis_cross_corrected <- running_adonis(phenos_selected %>%
                                                   filter(type == "cross-sectional" & analysis == "health_outcomes") %>%
                                                   pull(new_name),
                                                 smeta_w_phenos,
                                                 dist_matrix,
                                                 "cross",
                                                 "margin",
                                                 correct_for = correct_factors)

adonis_long_corrected <- running_adonis(phenos_selected %>%
                                                   filter(type == "longitudinal" & analysis == "health_outcomes") %>%
                                                   pull(new_name),
                                                 smeta_w_phenos,
                                                 dist_matrix,
                                                 "long",
                                                 "margin",
                                                 correct_for = correct_factors)

adonis_results_corrected <- adonis_cross_corrected %>%
  bind_rows(adonis_long_corrected) %>%
  mutate(FDR = p.adjust(`Pr(>F)`, "BH")) %>%
  ungroup()

adonis_together <- adonis_results %>%
  mutate(corrected = "notcorr") %>%
  bind_rows(adonis_results_corrected %>%
              mutate(corrected = "corr")) %>%
  pivot_wider(names_from = corrected,
              values_from = !c(rowname, corrected),
              names_sep = "_"
  )
#############################################################
# 3.1 Analysis: composition vs health outcomes
#############################################################

R2viz <- adonis_together %>%
  filter(FDR_notcorr < 0.05) %>%
  select(rowname, R2_corr, FDR_corr) %>%
  mutate(sign_after_corr = ifelse(FDR_corr < 0.05, "yes", "no")) %>%
  arrange(-desc(R2_corr)) %>%
  mutate(plot_name = case_when(rowname == "infant_health_eczema_diagnosis_strict" ~ "Eczema, strict",
                               rowname == "infant_health_wheeze_cross" ~ "Wheeze",
                               rowname == "infant_health_eczema_diagnosis_relaxed" ~ "Eczema, relaxed",
                               rowname == "infant_health_fever_cross" ~ "Fever",
                               rowname == "infant_growth_standardized_weight_slope_kg" ~ "Weight-gain slope",
                               .default = rowname)) %>%
  mutate(plot_name = factor(plot_name, levels = plot_name)) %>%
  ggplot(aes(x = plot_name, y = R2_corr, fill = sign_after_corr)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  geom_text(aes(label=round(R2_corr, 3)), hjust = -0.1, size = 2.5) +
  scale_fill_manual(values = c("grey", "firebrick")) +
  theme_bw() +
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = "bottom") + 
  labs(y = "Adonis R-sqaured", x="Health outcomes", fill = "Significant after correction") +
  ylim(0,0.008)

ggsave("05.PLOTS/07.Health_outcomes/VLP_health_outcomes_adonis.pdf",
       R2viz, "pdf", width=8, height=8, units="cm", dpi = 300)

#############################################################
# 3.1.2 Analysis: visualization of the health phenotypes
#############################################################

smeta_w_phenos %>%
  filter(!is.na(birth_deliverybirthcard_mode_binary) & 
           !is.na(infant_ffq_feeding_mode_simple) &
           !is.na(infant_health_eczema_diagnosis_relaxed)) %>%
  count(infant_health_eczema_diagnosis_relaxed, Timepoint_new) # choosing relax, because there centroides are stable

plot_smeta <- smeta_w_phenos %>%
  filter(
    !is.na(infant_health_eczema_diagnosis_relaxed) & 
           !is.na(infant_ffq_feeding_mode_simple) & 
           !is.na(birth_deliverybirthcard_mode_binary) ) %>%
 rename("Eczema, relaxed" = "infant_health_eczema_diagnosis_relaxed",
          "Wheeze" = "infant_health_wheeze_cross")
  

# eczema
ES_list <- virome_shaper("Eczema, relaxed",
                         plot_smeta,
                         dist_matrix)

pheno <- "Eczema, relaxed"
ES <- ggplot(ES_list[["pcoa_plot_data"]], aes(x = V1, y = V2, color = !!sym(pheno))) +
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.8), dpi = 600) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = !!sym(pheno), fill=!!sym(pheno)), linetype = 2) +
  geom_point(data=ES_list[["centroids"]], aes(V1, V2, fill = !!sym(pheno)), alpha = 0.7, shape = 23, size = 2, color = "black") +
  scale_fill_manual(values = c("grey", "#BF092F")) +
  scale_color_manual(values = c("grey", "#BF092F")) +
  ggnewscale::new_scale_color() +
  geom_path(
    data = ES_list[["centroids"]],
    aes(x = V1, y = V2, group = !!sym(pheno),  colour = !!sym(pheno)),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("black", "#9B0F06")) +
  theme_bw() +
  labs(
    x = paste0("PC1: ", ES_list[["pvar"]][1], "%"),
    y = paste0("PC1: ", ES_list[["pvar"]][2], "%")) + 
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9),
        legend.position = "right")

ggsave("05.PLOTS/07.Health_outcomes/VLP_eczema_PCoA_corrected.pdf",
       ES, "pdf", width=12, height=8, units="cm", dpi = 300)

#############################################################
# 3.3 Analysis: differentially abundant vOTUs
#############################################################

min_rpkm <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  pivot_longer(!New_CID) %>%
  filter(value != 0) %>%
  pull(value) %>%
  min()

VLP_filtered_10 <- VLP[rowSums(VLP > 0) >= round(0.1 * ncol(VLP), 0), ]
VLP_filtered_10 <- VLP_filtered_10[,colSums(VLP_filtered_10) > 0]

correct_factors <- c("infant_ffq_feeding_mode_simple", "birth_deliverybirthcard_mode_binary")

vOTU_vs_pehnos <- running_LMM_upd(phenos_list = "infant_health_eczema_diagnosis_relaxed",
                              metadata_w_phenos = smeta_w_phenos,
                              RPKM_table = VLP_filtered_10,
                              pseudocount = min_rpkm,
                              correct_for = correct_factors)

vOTU_vs_pehnos$FDR <- p.adjust(vOTU_vs_pehnos$`Pr(>|t|)`, 'BH') # no significant results
#############################################################
# 3.3 Analysis: trying family again
#############################################################

by_family <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr) %>%
  group_by(Family) %>%
  summarise(across(all_of(colnames(VLP)), ~ sum(.x))) %>%
  filter( rowSums( across(all_of(colnames(VLP))) > 0) >= 0.1*ncol(VLP)) %>%
  filter(!Family %in% c("Unclassified", "Unassigned")) %>%
  column_to_rownames(var = "Family")

family_vs_phenos <- running_LMM_upd(phenos_list = "infant_health_eczema_diagnosis_relaxed",
                                metadata_w_phenos = smeta_w_phenos,
                                RPKM_table = by_family,
                                pseudocount = min_rpkm,
                                correct_for = correct_factors)

family_vs_phenos$FDR <- p.adjust(family_vs_phenos$`Pr(>|t|)`, 'BH') # no significant results
#############################################################
# 3.4 Analysis: genome composition
#############################################################
by_genome_RPKM <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(genome, New_CID)) %>%
  group_by(genome) %>%
  summarise(across(all_of(colnames(VLP)), ~ sum(.x))) %>%
  filter( rowSums( across(all_of(colnames(VLP))) > 0) >= 0.1*ncol(VLP)) %>%
  filter(genome != "Unclassified") %>%
  mutate(genome = paste0(genome, "_RPKM")) %>%
  column_to_rownames(var = "genome")

genome_vs_phenos <- running_LMM_upd(phenos_list = "infant_health_eczema_diagnosis_relaxed",
                                metadata_w_phenos = smeta_w_phenos,
                                RPKM_table = by_genome_RPKM,
                                pseudocount = min_rpkm,
                                correct_for = correct_factors)

genome_vs_phenos$FDR <- p.adjust(genome_vs_phenos$`Pr(>|t|)`, 'BH')
#############################################################
# 3.5 Analysis: host composition
#############################################################
by_host_RPKM <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, Host) %>% mutate(Host = ifelse(Host == "invertebrates, vertebrates", "inver_ver", Host))) %>%
  group_by(Host) %>%
  summarise(across(all_of(colnames(VLP)), ~ sum(.x))) %>%
  filter( rowSums(across(all_of(colnames(VLP))) > 0) >= 0.1*ncol(VLP) ) %>%
  filter(!Host %in% c("Unknown")) %>%
  mutate(Host = paste0(Host, "_RPKM")) %>%
  column_to_rownames(var = "Host")

host_vs_phenos <- running_LMM_upd(phenos_list = "infant_health_eczema_diagnosis_relaxed",
                              metadata_w_phenos = smeta_w_phenos,
                              RPKM_table = by_host_RPKM,
                              pseudocount = min_rpkm,
                              correct_for = correct_factors)

host_vs_phenos$FDR <- p.adjust(host_vs_phenos$`Pr(>|t|)`, 'BH') 


by_genhost_RPKM <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, Host_simple)) %>%
  group_by(Host_simple) %>%
  summarise(across(all_of(colnames(VLP)), ~ sum(.x))) %>%
  filter( rowSums(across(all_of(colnames(VLP))) > 0) >= 0.1*ncol(VLP) ) %>%
  filter(!Host_simple %in% c("Unknown")) %>%
  mutate(Host_simple = paste0(Host_simple, "_RPKM")) %>%
  column_to_rownames(var = "Host_simple")

genhost_vs_phenos <- running_LMM_upd(phenos_list = "infant_health_eczema_diagnosis_relaxed",
                              metadata_w_phenos = smeta_w_phenos,
                              RPKM_table = by_genhost_RPKM,
                              pseudocount = min_rpkm,
                              correct_for = correct_factors)

genhost_vs_phenos$FDR <- p.adjust(genhost_vs_phenos$`Pr(>|t|)`, 'BH') # abundance of phages is down?

#############################################################
# 3.5.1 Analysis: visualization: abundance by host
#############################################################

# Prokaryote_RPKM is lower in M1 and M3, but switches in M6 and M12

VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, Host_simple)) %>%
  group_by(Host_simple) %>%
  summarise(across(all_of(colnames(VLP)), ~ sum(.x))) %>%
  mutate(Host_simple = paste0(Host_simple, "_RPKM")) %>%
  pivot_longer(!Host_simple, names_to = "Sequencing_ID") %>%
  pivot_wider(names_from = Host_simple) %>%
  left_join(smeta_w_phenos) %>%
  filter(!is.na(infant_health_eczema_diagnosis_relaxed) #&
           #!is.na(infant_ffq_feeding_mode_simple) &
           #!is.na(birth_deliverybirthcard_mode_binary)
         ) %>%
  ggplot(aes(Timepoint_new, log(Prokaryote_RPKM + min_rpkm/2), fill = infant_health_eczema_diagnosis_relaxed)) +
  geom_boxplot()


#############################################################
# 3.6 Analysis: viral sample metrics
#############################################################

virvars <- c("vir_richness_cf", "vir_diversity", "temperate_RAb", "lytic_RAb",
             "temperate_richness", "lytic_richness")

pehnos_significant_before_correction <- c("infant_health_eczema_diagnosis_relaxed", "infant_health_fever_cross", "infant_health_wheeze_cross", "infant_growth_standardized_weight_slope_kg")

virvar_pseudocounts <- smeta_w_phenos %>%
  summarise(across(all_of(virvars), 
                   ~ min(.x[.x > 0], na.rm = TRUE))) %>%
  pivot_longer(cols = all_of(virvars))

virvar_corrected_1 <- running_virvar_corrected(phenos_to_run = pehnos_significant_before_correction,
                                               metadata_w_phenos = smeta_w_phenos,
                                               vivars_to_run = c("temperate_RAb"),
                                               correct_for = "birth_deliverybirthcard_mode_binary")

virvar_corrected_2 <- running_virvar_corrected(phenos_to_run = pehnos_significant_before_correction,
                                               metadata_w_phenos = smeta_w_phenos,
                                               vivars_to_run = c("vir_richness_cf", "lytic_RAb", "lytic_richness"),
                                               correct_for = "infant_ffq_feeding_mode_simple")

virvar_notcorrected <- running_virvar_corrected(phenos_to_run = pehnos_significant_before_correction,
                                                metadata_w_phenos = smeta_w_phenos,
                                                vivars_to_run = c("vir_diversity", "temperate_richness"),
                                                correct_for = NULL)


virvars_vs_phenos <- virvar_corrected_1 %>%
  bind_rows(virvar_corrected_2) %>%
  bind_rows(virvar_notcorrected) %>%
  mutate(cat_virvar = ifelse(grepl("vir_", virvar), "Richness", "Lifestyle")) %>%
  group_by(cat_virvar) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`, 'BH')) %>%
  ungroup()

# does eczema stay significant even after correction for family history of disease?

trab_corrected <- running_virvar_corrected(phenos_to_run = "infant_health_eczema_diagnosis_relaxed",
                                               metadata_w_phenos = smeta_w_phenos,
                                               vivars_to_run = c("temperate_RAb"),
                                               correct_for = c("birth_deliverybirthcard_mode_binary", 
                                                               "infant_scorad_family_history_allergic_disease"))
#############################################################
# 3.6.1 Analysis: visualization
#############################################################

eczema_trab <- smeta_w_phenos %>%
  filter(!is.na(infant_health_eczema_diagnosis_relaxed)
  ) %>%
  ggplot(aes(Timepoint_new, log(temperate_RAb + virvar_pseudocounts$value[virvar_pseudocounts$name == "temperate_RAb"]), 
            fill = infant_health_eczema_diagnosis_relaxed)) +
  ggrastr::rasterise(geom_jitter(aes(color = infant_health_eczema_diagnosis_relaxed), position = position_jitterdodge(), size=0.1, alpha=1), dpi = 300)+
  geom_boxplot(alpha=0.6, outlier.shape = NA) +
  theme_minimal() +
  theme(strip.background = element_rect(NA),
        legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm")) + 
  labs(y = "Temperate phage abundance, log-scaled", x="Timepoint", fill = "Eczema", color = "Eczema") +
  scale_fill_manual(values = c("grey", "#BF092F"), labels = c("no", "yes")) +
  scale_color_manual(values = c("grey", "#BF092F"), labels = c("no", "yes")) 

ggsave("05.PLOTS/07.Health_outcomes/VLP_eczema_temperate_RAb.pdf",
       eczema_trab, "pdf", width=6, height=8, units="cm", dpi = 300)
#############################################################
# 3.7 Analysis: viral sample metrics per timepoint
#############################################################
virvar_timepoint_1 <- running_virvar_corrected_timepoint(phenos_to_run = pehnos_significant_before_correction,
                                           metadata_w_phenos = smeta_w_phenos,
                                           vivars_to_run = c("temperate_RAb"),
                                           correct_for = "birth_deliverybirthcard_mode_binary")

virvar_timepoint_2 <- running_virvar_corrected_timepoint(phenos_to_run = pehnos_significant_before_correction,
                                              metadata_w_phenos = smeta_w_phenos,
                                              vivars_to_run = c("vir_richness_cf", "lytic_RAb", "lytic_richness"),
                                              correct_for = "infant_ffq_feeding_mode_simple")

virvar_timepoint_3 <- running_virvar_corrected_timepoint(phenos_to_run = pehnos_significant_before_correction,
                                              metadata_w_phenos = smeta_w_phenos,
                                              vivars_to_run = c("vir_diversity", "temperate_richness"),
                                              correct_for = NULL)

virvars_vs_phenos_timepoint <- virvar_timepoint_1 %>%
  bind_rows(virvar_timepoint_2) %>%
  bind_rows(virvar_timepoint_3) %>%
  mutate(cat_virvar = ifelse(grepl("vir_", virvar), "Richness", "Lifestyle")) %>%
  group_by(cat_virvar) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`, 'BH')) %>%
  ungroup()

#############################################################
# 3.8 Analysis: by genome per timepoint
#############################################################

smeta_w_genome <- by_genome_RPKM %>%
  rownames_to_column(var = "genome") %>%
  pivot_longer(!genome, names_to = "Sequencing_ID") %>%
  pivot_wider(names_from = genome) %>%
  left_join(smeta_w_phenos)


genomes_vs_phenos_timepoint <- running_virvar_corrected_timepoint(phenos_to_run = pehnos_significant_before_correction,
                                                                  metadata_w_phenos = smeta_w_genome,
                                                                  vivars_to_run = c("RNA_RPKM", "dsDNA_RPKM", "ssDNA_RPKM"),
                                                                  correct_for = NULL)

genomes_vs_phenos_timepoint <- genomes_vs_phenos_timepoint %>%
  group_by(timepoint) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`, 'BH')) %>%
  ungroup()

#############################################################
# 3.8.1 Analysis: by genome per timepoint visualization
#############################################################
# detection panel

p1 <- smeta_w_genome %>%
  filter(!is.na(infant_health_wheeze_cross)) %>%
  group_by(Timepoint_new, infant_health_wheeze_cross) %>%
  summarise(n_total = n(),
            n_present = sum(RNA_RPKM > 0),
            n_frac = sum(RNA_RPKM > 0)/n() * 100,
            .groups = "drop") %>%
  ggplot(aes(x = Timepoint_new, y = n_frac, fill = infant_health_wheeze_cross)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.6, color = "black") +
  labs(y = "Presence, %", x = NULL, fill = "Wheeze") +
  scale_fill_manual(values = met.brewer("Ingres")[c(4,5)]) +
  ylim(0, 85) +
  theme_minimal() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.5, units = "cm"))

# abundance
p2 <- smeta_w_genome %>%
  filter(!is.na(infant_health_wheeze_cross) & 
           RNA_RPKM > 0
         ) %>%
  ggplot(aes(x = Timepoint_new, y = log(RNA_RPKM), 
             fill = infant_health_wheeze_cross,
             color = infant_health_wheeze_cross)) +
  ggrastr::rasterise(geom_jitter(aes(color = infant_health_wheeze_cross), 
                                 position = position_jitterdodge(), 
                                 size=0.1, alpha=1), dpi = 600) +
  geom_boxplot(alpha=0.4, color = "black", outlier.shape = NA) +
  labs(y = "RNA viruses, log(RPKM > 0)", x = "Timepoint", fill = "Wheeze", color = "Wheeze") +
  scale_fill_manual(values = met.brewer("Ingres")[c(4,5)]) +
  scale_color_manual(values = met.brewer("Ingres")[c(4,5)]) +
  theme_minimal() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.5, units = "cm"))

plot_RNA_wheeze <- (p1 / p2) + plot_layout(guides = "collect", heights = c(3, 4)) & theme(legend.position = 'bottom')

ggsave("05.PLOTS/07.Health_outcomes/VLP_wheeze_RNA.pdf",
       plot_RNA_wheeze, "pdf", width=7, height=9, units="cm", dpi = 300)

#############################################################
# 3.9 Analysis: by host per timepoint
#############################################################

smeta_w_host <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, Host, genome) %>% mutate(Host = ifelse(Host == "invertebrates, vertebrates", "inver_ver", Host))) %>%
  filter(genome == "RNA") %>%
  group_by(Host) %>%
  summarise(across(all_of(colnames(VLP)), ~ sum(.x))) %>%
  filter( rowSums(across(all_of(colnames(VLP))) > 0) >= 0.1*ncol(VLP) ) %>%
  filter(!Host %in% c("Unknown")) %>%
  mutate(Host = paste0(Host, "_RPKM")) %>%
  pivot_longer(!Host, names_to = "Sequencing_ID") %>%
  pivot_wider(names_from = Host) %>%
  left_join(smeta_w_phenos) %>%
  filter(Timepoint_new == "M1")

smeta_w_host %>%
  filter(!is.na(infant_health_wheeze_cross)) %>%
  summarise(across(all_of(c("invertebrates_RPKM", "plants_RPKM", "vertebrates_RPKM")), ~ sum(.x > 0)))

bos <- running_virvar_corrected_timepoint(phenos_to_run = "infant_health_wheeze_cross",
                                   metadata_w_phenos = smeta_w_host,
                                   vivars_to_run = c("plants_RPKM", "vertebrates_RPKM"),
                                   correct_for = NULL)

bos <- bos %>%
  group_by(timepoint) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`, 'BH')) %>%
  ungroup()

smeta_w_host %>%
  filter(!is.na(infant_health_wheeze_cross)) %>%
  group_by(Timepoint_new, infant_health_wheeze_cross) %>%
  summarise(n_total = n(),
            n_present = sum(plants_RPKM > 0),
            n_frac = sum(plants_RPKM > 0)/n() * 100,
            .groups = "drop") %>%
  ggplot(aes(x = Timepoint_new, y = n_frac, fill = infant_health_wheeze_cross)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.6, color = "black") +
  labs(y = "Presence, %", x = NULL, fill = "Wheeze") +
  scale_fill_manual(values = met.brewer("Ingres")[c(4,5)]) +
  ylim(0, 100) +
  theme_minimal() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.5, units = "cm"))

smeta_w_host %>%
  filter(!is.na(infant_health_wheeze_cross) & 
           vertebrates_RPKM > 0
  ) %>%
  ggplot(aes(x = Timepoint_new, y = log(vertebrates_RPKM), 
             fill = infant_health_wheeze_cross,
             color = infant_health_wheeze_cross)) +
  ggrastr::rasterise(geom_jitter(aes(color = infant_health_wheeze_cross), 
                                 position = position_jitterdodge(), 
                                 size=0.1, alpha=1), dpi = 600) +
  geom_boxplot(alpha=0.4, color = "black", outlier.shape = NA) +
  labs(y = "RNA viruses, log(RPKM > 0)", x = "Timepoint", fill = "Wheeze", color = "Wheeze") +
  scale_fill_manual(values = met.brewer("Ingres")[c(4,5)]) +
  scale_color_manual(values = met.brewer("Ingres")[c(4,5)]) +
  theme_minimal() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.5, units = "cm"))

# follow-up (inn a separate file?)
