setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# Here we explore which phenotypes shape the virome
#############################################################

#############################################################
# 0. Used files source
#############################################################
selected_cross <- c("mother_health_smoked_one_whole_year_p18", 
                    "birth_deliverybirthcard_mode_binary",
                    "birth_deliverybirthcard_place_delivery_simple",
                    "family_pets_any", 
                    "infant_health_eczema_diagnosis_relaxed",
                    "infant_health_eczema_diagnosis_strict",
                    "infant_health_fever",
                    "infant_health_wheeze",
                    "mother_birthcardhealth_gravidity", # parity has more missing data
                    "mother_birthcardself_gestational_age_weeks",
                    "infant_misc_sex",
                    "infant_scorad_family_history_allergic_disease",
                    "infant_growth_standardized_weight_slope_kg",
                    "infant_health_cough" # check distribution
) 

selected_long <- c("infant_health_cough", # check distribution
                      "infant_health_eczema_questionnaire", # check distribution
                      "infant_health_fever", # check distribution
                      "infant_health_wheeze", # check distribution
                      "infant_ffq_feeding_mode_simple",
                      "infant_BITSS",
                      "infant_med_antibiotics_J01", # induction
                      "infant_scorad_measurement", # for eczema
                      "infant_sleep_day_hours"
)

phenos_select <- data.frame(phenotype = selected_cross, 
           type = rep("cross-sectional", length(selected_cross))) %>%
  bind_rows(data.frame(phenotype = selected_long, 
                       type = rep("longitudinal", length(selected_long)))) %>%
  mutate(analysis = ifelse(phenotype %in% c("mother_health_smoked_one_whole_year_p18",
                                            "birth_deliverybirthcard_mode_binary",
                                            "birth_deliverybirthcard_place_delivery_simple",
                                            "infant_ffq_feeding_mode_simple",
                                            "family_pets_any",
                                            "mother_birthcardhealth_gravidity", 
                                            "mother_birthcardself_gestational_age_weeks",
                                            "infant_misc_sex",
                                            "infant_scorad_family_history_allergic_disease",
                                            "infant_med_antibiotics_J01"), "shaping", "health_outcomes")) %>%
  mutate(new_name = case_when(phenotype %in% c("infant_health_fever", "infant_health_wheeze", "infant_health_cough") & type == "cross-sectional" ~ paste0(phenotype, "_cross"),
                              phenotype %in% c("infant_health_fever", "infant_health_wheeze", "infant_health_cough") & type == "longitudinal" ~ paste0(phenotype, "_long"),
                              .default = phenotype)) %>%
  mutate(cat_vs_num = ifelse(phenotype %in% c("mother_birthcardself_gestational_age_weeks", "infant_growth_standardized_weight_slope_kg",
                                              "infant_scorad_measurement", "infant_sleep_day_hours"), "numeric", "categorical"))

write.table(phenos_select, "06.CLEAN_DATA/Intermediate/Phenotype_selection.txt", sep='\t', row.names = F, quote=F)
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
#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(tidyverse)
library(MetBrewer)
library(vegan)
library(skimr)
#############################################################
# 2. Load Input Data
#############################################################
# metadatas:
smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)

smeta <- smeta %>%
  filter(Type == "K") %>%
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) %>% # treating 2nd pregnancy as a separate family:
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T))

# Longitudinal phenotypes:
L_phenos <- read.table('06.CLEAN_DATA/Phenotypes/masterfile_longitudinal_2023_09_29.txt', sep='\t', header=T)
L_phenos <- L_phenos[L_phenos$SAMPLE_ID %in% smeta$Universal_ID,]

# Cross-sectional phenotypes:
C_phenos <- read.table('06.CLEAN_DATA/Phenotypes/masterfile_cross_sectional_2023_11_15.txt', sep = '\t', header=T)
C_phenos <- C_phenos[C_phenos$FAMILY %in% smeta$FAMILY,]

# merging smeta w phenos:
smeta_w_phenos <- smeta %>%
  left_join(C_phenos %>% 
              select(all_of(selected_cross), next_id_infant), 
            by = c("NEXT_ID" = "next_id_infant")) %>%
  left_join(L_phenos %>% 
              select(all_of(selected_long), next_id_infant, timepoint), 
            by = c("NEXT_ID" = "next_id_infant", "Timepoint_new" = "timepoint"),
            suffix = c("_cross", "_long")) %>%
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12"), ordered = T))

# abundance tables etc:
VLP <- read.table("06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt", sep='\t', header=T)
VLP <- VLP[,colnames(VLP) %in% smeta$Sequencing_ID]
VLP <- VLP[rowSums(VLP) > 0,]

# ordering smeta same way as VLP:
smeta_w_phenos <- smeta_w_phenos[match(colnames(VLP), smeta_w_phenos$Sequencing_ID), ]
# saving it for further analyses:
write.table(smeta_w_phenos, "06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_phenotypes.txt", sep='\t', quote = F, row.names = F)

# BC dist:
VLP_RAB <- as.data.frame(t(VLP)/colSums(VLP))
VLP_dist <- vegdist(VLP_RAB, method = "bray", parallel = 6)
dist_matrix <- as.matrix(VLP_dist)

# virus metadata
ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)
ETOF_vOTUr <- ETOF_vOTUr %>%
  separate(
    col    = tax_ictv_aai,
    into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  )
#############################################################
# 3.1 Analysis: searching for the virome shaping phenotypes
#############################################################
set.seed(444)
cross_summary <- C_phenos %>%
  select(all_of(selected_cross)) %>% 
  mutate(across(where(is.character), as.factor)) %>%
  skim()

long_summary <- L_phenos %>%
  select(all_of(selected_long)) %>% 
  mutate(across(where(is.character), as.factor)) %>%
  skim() 

# what about tech features? they explain it by less than 1% -> skip correcting for (I guess it is also because vlp data is super sparse)

# virome composition vs phenotypes
adonis_results_cross <- map_dfr(phenos_select %>%
                                  filter(type == "cross-sectional") %>%
                                  pull(new_name), 
                                function(pheno){
  
  pheno_df <- smeta_w_phenos %>%
    filter(!is.na(!!sym(pheno)))
  
  dister <- dist_matrix[colnames(dist_matrix) %in% pheno_df$Sequencing_ID,
                     row.names(dist_matrix) %in% pheno_df$Sequencing_ID]
  
  F1 <- as.formula(paste0("dister ~ Timepoint_new + ", pheno))
  perm <- 999
  
  if ( identical(pheno_df$Sequencing_ID, colnames(dister)) ) {
    res <- adonis2(F1, 
            data=pheno_df,
            by = "terms",
            #strata = pheno_df$NEXT_ID,
            #permutations = perm,
            permutations = how(blocks = pheno_df$NEXT_ID, nperm = 999),
            parallel = 4) 
    
  } else {
    return(NULL)
  }
  
  res %>%
    as.data.frame() %>%
    rownames_to_column(var = "rowname") %>%
    filter(rowname == pheno) %>%
    mutate(N_permutations = perm) %>%
    mutate(N_samples = nrow(pheno_df))
  
})

adonis_results_long <- map_dfr(phenos_select %>%
                                 filter(type == "longitudinal") %>%
                                 pull(new_name), 
                               function(pheno){
  
  pheno_df <- smeta_w_phenos %>%
    filter(!is.na(!!sym(pheno)))
  
  dister <- dist_matrix[colnames(dist_matrix) %in% pheno_df$Sequencing_ID,
                        row.names(dist_matrix) %in% pheno_df$Sequencing_ID]
  
  F1 <- as.formula(paste0("dister ~ Timepoint_new + ", pheno))
  perm <- 999
  
  if ( identical(pheno_df$Sequencing_ID, colnames(dister)) ) {
    res <- adonis2(F1, 
                   data=pheno_df,
                   by = "margin",
                   strata = pheno_df$NEXT_ID,
                   permutations = perm,
                   parallel = 4) 
    
  } else {
    return(NULL)
  }
  
  res %>%
    as.data.frame() %>%
    rownames_to_column(var = "rowname") %>%
    filter(rowname == pheno) %>%
    mutate(N_permutations = perm) %>%
    mutate(N_samples = nrow(pheno_df))
  
})

adonis_results <- adonis_results_cross %>%
  bind_rows(adonis_results_long) %>%
  mutate(analysis = phenos_select$analysis[match(rowname, phenos_select$new_name)]) %>%
  group_by(analysis) %>%
  mutate(FDR = p.adjust(`Pr(>F)`, "BH")) %>%
  ungroup()

# testing how homogeneous is the variance between pheno categories
betadisper <- map_dfr(phenos_select %>%
                        filter(cat_vs_num == "categorical") %>%
                        pull(new_name), 
                      function(pheno){
  
  pheno_df <- smeta_w_phenos %>%
    filter(!is.na(!!sym(pheno)))
  
  dister <- dist_matrix[colnames(dist_matrix) %in% pheno_df$Sequencing_ID,
                        row.names(dist_matrix) %in% pheno_df$Sequencing_ID]
  
  if ( identical(pheno_df$Sequencing_ID, colnames(dister)) ) {
    bd <- betadisper(as.dist(dister), pheno_df[,pheno])
  } else {
    return(NULL)
  }
  
  permutest(bd)$tab %>%
    rownames_to_column(var = "rowname") %>%
    filter(rowname == "Groups") %>%
    mutate(pheno = pheno)
  
})


# results summary table:
results_summary <- phenos_select %>%
  left_join(cross_summary %>%
              mutate(type = "cross-sectional") %>%
              bind_rows(long_summary %>%
                          mutate(type = "longitudinal") ) %>%
              select(-skim_type), by = c("phenotype" = "skim_variable", "type" = "type")) %>%
  left_join(adonis_results %>%
              select(-SumOfSqs, -analysis) %>%
              rename(p_value_adonis=`Pr(>F)`), by = c("new_name" = "rowname")) %>%
  left_join(betadisper %>%
              select(`F`, `Pr(>F)`, pheno) %>%
              rename(p_value_disper = `Pr(>F)`),
              by = c("new_name" = "pheno"), suffix = c("", "_disper"))

writexl::write_xlsx(results_summary, '07.RESULTS/Adonis2_betadisper_VLP_virome.xlsx')

#############################################################
# 3.2.1 Analysis: virome shaping phenotype visualization
#############################################################

R2viz <- results_summary %>%
  filter(analysis == "shaping") %>%
  filter(FDR < 0.05 & (p_value_disper > 0.05 | is.na(p_value_disper))) %>%
  select(new_name, R2) %>%
  bind_rows(data.frame(new_name = c("NEXT_ID", "Timepoint"),
                       R2 = c(0.513, 0.023))) %>% # from 02.Virome_composition_dynamics.R; section 3.3
  arrange(-desc(R2)) %>%
  mutate(plot_name = case_when(new_name == "NEXT_ID" ~ "Infant ID",
                               new_name == "infant_ffq_feeding_mode_simple" ~ "Feeding mode",
                               new_name == "birth_deliverybirthcard_mode_binary" ~ "Delivery mode",
                               new_name == "infant_scorad_family_history_allergic_disease" ~ "Family allergies",
                               new_name == "birth_deliverybirthcard_place_delivery_simple" ~ "Place of delivery",
                               new_name == "mother_birthcardself_gestational_age_weeks" ~ "Gestational age",
                               new_name == "family_pets_any" ~ "Pets",
                               new_name == "infant_misc_sex" ~ "Sex",
                               .default = new_name)) %>%
  mutate(plot_name = factor(plot_name, levels = plot_name)) %>%
  ggplot(aes(x = plot_name, y = R2)) +
  geom_bar(stat = "identity", fill="#D25353") + 
  coord_flip() +
  geom_text(aes(label=round(R2, 3)), hjust = -0.1, size = 2.5) +
  theme_bw() +
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size=9)) + 
  labs(y = "Technical feature", x="Adonis R-sqaured") +
  ylim(0,0.54)

ggsave("05.PLOTS/07.Health_outcomes/VLP_shaping_factors_adonis.pdf",
       R2viz, "pdf", width=8, height=8, units="cm", dpi = 300)

#############################################################
# 3.2.2 Analysis: virome shaping phenotype visualization
#############################################################
plot_smeta <- smeta_w_phenos %>%
  mutate(infant_ffq_feeding_mode_simple = case_when(!is.na(infant_ffq_feeding_mode_simple) & infant_ffq_feeding_mode_simple == "excl_BF" ~ "Breastfeeding",
                                                    !is.na(infant_ffq_feeding_mode_simple) & infant_ffq_feeding_mode_simple == "excl_FF" ~ "Formula feeding",
                                                    .default = infant_ffq_feeding_mode_simple)) %>%
  mutate(birth_deliverybirthcard_mode_binary = case_when(!is.na(birth_deliverybirthcard_mode_binary) & birth_deliverybirthcard_mode_binary == "CS" ~ "C-section",
                                                    !is.na(birth_deliverybirthcard_mode_binary) & birth_deliverybirthcard_mode_binary == "VG" ~ "Vaginal delivery",
                                                    .default = birth_deliverybirthcard_mode_binary)) %>%
  rename("Feeding mode" = "infant_ffq_feeding_mode_simple",
         "Delivery mode" = "birth_deliverybirthcard_mode_binary")
  
# feeding mode:
FM_list <- virome_shaper("Feeding mode",
                            plot_smeta,
                            dist_matrix)

pheno <- "Feeding mode"
FM <- ggplot(FM_list[["pcoa_plot_data"]], aes(x = V1, y = V2, color = !!sym(pheno))) +
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.8), dpi = 600) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = !!sym(pheno), fill=!!sym(pheno)), linetype = 2) +
  geom_point(data=FM_list[["centroids"]], aes(V1, V2, fill = !!sym(pheno)), alpha = 0.8, shape = 23, size = 2, color = "black") +
  scale_fill_manual(values = met.brewer("Hiroshige")[c(4,6)]) +
  scale_color_manual(values = met.brewer("Hiroshige")[c(4,6)]) +
  ggnewscale::new_scale_color() +
  geom_path(
    data = FM_list[["centroids"]],
    aes(x = V1, y = V2, group = !!sym(pheno),  colour = !!sym(pheno)),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = met.brewer("Hiroshige")[c(1,9)]) +
  theme_bw() +
  labs(
    x = paste0("PC1: ", FM_list[["pvar"]][1], "%"),
    y = paste0("PC1: ", FM_list[["pvar"]][2], "%")) + 
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9),
        legend.position = "bottom")

ggsave("05.PLOTS/07.Health_outcomes/VLP_feeding_mode_PCoA.pdf",
       FM, "pdf", width=8, height=9, units="cm", dpi = 300)

# delivery mode:
DM_list <- virome_shaper("Delivery mode",
                         plot_smeta,
                         dist_matrix)

pheno <- "Delivery mode"
DM <- ggplot(DM_list[["pcoa_plot_data"]], aes(x = V1, y = V2, color = !!sym(pheno))) +
  ggrastr::rasterise(geom_point(size = 1, alpha = 0.8), dpi = 600) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = !!sym(pheno), fill=!!sym(pheno)), linetype = 2) +
  geom_point(data=DM_list[["centroids"]], aes(V1, V2, fill = !!sym(pheno)), alpha = 0.8, shape = 23, size = 2, color = "black") +
  scale_fill_manual(values = met.brewer("Morgenstern")[c(2,5)]) +
  scale_color_manual(values = met.brewer("Morgenstern")[c(2,5)]) +
  ggnewscale::new_scale_color() +
  geom_path(
    data = DM_list[["centroids"]],
    aes(x = V1, y = V2, group = !!sym(pheno),  colour = !!sym(pheno)),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = met.brewer("Morgenstern")[c(1,8)]) +
  theme_bw() +
  labs(
    x = paste0("PC1: ", DM_list[["pvar"]][1], "%"),
    y = paste0("PC1: ", DM_list[["pvar"]][2], "%")) + 
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9),
        legend.position = "bottom")

ggsave("05.PLOTS/07.Health_outcomes/VLP_delivery_mode_PCoA.pdf",
       DM, "pdf", width=8, height=9, units="cm", dpi = 300)

