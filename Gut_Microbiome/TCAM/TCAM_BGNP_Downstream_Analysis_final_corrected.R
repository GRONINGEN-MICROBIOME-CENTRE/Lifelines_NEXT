################################################################################
##### Big gut NEXT paper: TCAM analysis
### Author(s): Nataliia Kuzub, Asier Fernández-Pato, Trishla Sinha
### Last updated: 25th March, 2024
################################################################################

################################################################################
# Loading libraries
################################################################################
library(vegan)
library(dplyr)
library(ggplot2)
library(readr)
library(ggExtra)

# Set working directory 
setwd("C:\\Users\\Natal\\Documents\\UMCG\\TCAM_BGNP")

################################################################################
# Loading functions
################################################################################

do_preprocessing_infants <- function(TCAM, cross_phenotypes, cross_selection, timepoints){
  # TCAM table preprocesing: output of the TCAM.ipynb
  rownames(TCAM) <- TCAM$X
  TCAM <- subset(TCAM, select = -c(X))
  
  # Fitering out the samples excluded for this analysis
  cross_phenotypes <- cross_phenotypes[cross_phenotypes$next_id_infant %in% row.names(TCAM), ]
  cross_phenotypes$mother_birthcard_parity <- ifelse(is.na(cross_phenotypes$mother_birthcard_parity), NA, 
                                                     ifelse(cross_phenotypes$mother_birthcard_parity >= 2, 2, cross_phenotypes$mother_birthcard_parity))
  row.names(cross_phenotypes) <- cross_phenotypes$next_id_infant
  
  # Selecting data relevant for infant 
  infant_cross_selection <- cross_selection[cross_selection$TCAM == 1 & cross_selection$only_run_vaginal_delivery == 0 & !is.na(cross_selection$TCAM),]
  column_names <- infant_cross_selection[[1]]
  infant_cross_phenotypes <- cross_phenotypes %>% select(all_of(column_names))
  
  # infant_cross_phenotypes_filtered <-drop_phenotypes(infant_cross_phenotypes, non_NA_values, minimum_variance)
  
  # Preparing phenotypes for associations
  to_exclude <- c("FAMILY", "next_id_mother", "next_id_infant", "next_id_partner", "infant_relations", "sibling_number", "twin_pair")
  infant_cross_phenotypes_filtered <- infant_cross_phenotypes[,!(colnames(infant_cross_phenotypes) %in% to_exclude)]
  
  
  return(list(TCAM, infant_cross_phenotypes_filtered))
}



do_preprocessing_infants_vagonly <- function(TCAM, cross_phenotypes, cross_selection, timepoints){
  # TCAM table preprocesing: output of the TCAM.ipynb
  rownames(TCAM) <- TCAM$X
  TCAM <- subset(TCAM, select = -c(X))
  
  # Fitering out the samples excluded for this analysis
  cross_phenotypes <- cross_phenotypes[cross_phenotypes$next_id_infant %in% row.names(TCAM) & 
                                         cross_phenotypes$birth_deliverybirthcard_mode_binary == "VG" &
                                         !is.na(cross_phenotypes$birth_deliverybirthcard_mode_binary), ]
  row.names(cross_phenotypes) <- cross_phenotypes$next_id_infant
  
  # Remove CS infants from the TCAM table
  TCAM <- TCAM[row.names(TCAM) %in% row.names(cross_phenotypes), ]
  
  # Selecting data relevant for infant 
  infant_cross_selection <- cross_selection[cross_selection$TCAM == 1 & cross_selection$only_run_vaginal_delivery == 1 & !is.na(cross_selection$TCAM),]
  column_names <- infant_cross_selection[[1]]
  infant_cross_phenotypes <- cross_phenotypes %>% select(all_of(column_names))
  
  # infant_cross_phenotypes_filtered <-drop_phenotypes(infant_cross_phenotypes, non_NA_values, minimum_variance)
  
  # Preparing phenotypes for associations
  to_exclude <- c("FAMILY", "next_id_mother", "next_id_infant", "next_id_partner", "infant_relations", "sibling_number", "twin_pair")
  infant_cross_phenotypes_filtered <- infant_cross_phenotypes[,!(colnames(infant_cross_phenotypes) %in% to_exclude)]
  
  
  return(list(TCAM, infant_cross_phenotypes_filtered))
}


do_preprocessing_mothers <- function(TCAM, cross_phenotypes, cross_selection, timepoints){
  # TCAM table preprocesing: output of the TCAM.ipynb
  rownames(TCAM) <- TCAM$X
  TCAM <- subset(TCAM, select = -c(X))
  
  # Fitering out the samples excluded for this analysis
  cross_phenotypes <- cross_phenotypes[cross_phenotypes$next_id_mother %in% row.names(TCAM), ]
  cross_phenotypes$mother_birthcard_parity <- ifelse(is.na(cross_phenotypes$mother_birthcard_parity), NA, 
                                                     ifelse(cross_phenotypes$mother_birthcard_parity >= 2, 2, cross_phenotypes$mother_birthcard_parity))
  row.names(cross_phenotypes) <- cross_phenotypes$next_id_mother
  
  # Selecting data relevant for infant 
  mother_cross_selection <- cross_selection[cross_selection$TCAM_mothers == 1 & !is.na(cross_selection$TCAM_mothers),]
  column_names <- mother_cross_selection[[1]]
  mother_cross_phenotypes <- cross_phenotypes %>% select(all_of(column_names))
  
  
  # Preparing phenotypes for associations
  to_exclude <- c("FAMILY", "next_id_mother", "next_id_infant", "next_id_partner", "infant_relations", "sibling_number", "twin_pair")
  mother_cross_phenotypes_filtered <- mother_cross_phenotypes[,!(colnames(mother_cross_phenotypes) %in% to_exclude)]
  
  
  return(list(TCAM, mother_cross_phenotypes_filtered))
}


do_permanova_simple <- function(processed_data, permutation_num){
  # Generate the Aitchison distance matrices
  Distance <- vegdist(processed_data[[1]], method = "euclidean" )
  Distance <- data.frame(as.matrix((Distance)))
  
  # Run PERMANOVA to estimate the effect of phenotypes on overall bacterial composition
  adonis_results_raw <- data.frame(Df = numeric(ncol(processed_data[[2]])),
                                   F = numeric(ncol(processed_data[[2]])),
                                   R2 = numeric(ncol(processed_data[[2]])),
                                   p_value = numeric(ncol(processed_data[[2]])))
  
  rownames(adonis_results_raw) <- colnames(processed_data[[2]])
  
  for (i in 1:ncol(processed_data[[2]])) {
    # Subset distance matrix and phenotypes, removing rows with NAs
    na_rows <- which(is.na(processed_data[[2]][, i]))
    distmat_cleaned <- Distance[!(rownames(Distance) %in% row.names(processed_data[[2]])[na_rows]),
                                !(colnames(Distance) %in% row.names(processed_data[[2]])[na_rows])]
    phenos2 <- processed_data[[2]][row.names(processed_data[[2]]) %in% row.names(distmat_cleaned), ]
    phenos2 <- phenos2[ match(rownames(distmat_cleaned), rownames(phenos2) ), ]
    # Run adonis2
    ad1 <- adonis2(distmat_cleaned ~ phenos2[[i]],
                   permutations = permutation_num,
                   parallel = 8,
                   na.action = na.fail,
                   by = "margin")
    # Save results
    adonis_results_raw[i, ] <- c(ad1$Df[1], ad1$F[1], ad1$R2[1], ad1$"Pr(>F)"[1])
  }
  
  # Calculate FDR-adjusted p-values
  adonis_results_raw$FDR <- p.adjust(adonis_results_raw$p_value, method = "BH")
  
  return(adonis_results_raw)
}


do_cor_phenotype_preselection <- function(cross_phenotypes, correction_phenotype, next_id_col){
  cross_phenotypes_selection <- cross_phenotypes[, c(next_id_col, correction_phenotype)]
  cross_phenotypes_selection <- cross_phenotypes_selection[(!is.na(cross_phenotypes_selection[[next_id_col]])) &
                                                             (!is.na(cross_phenotypes_selection[[correction_phenotype]])), ]
  cross_phenotypes_selection <- cross_phenotypes_selection[!duplicated(cross_phenotypes_selection[[next_id_col]]),]
  row.names(cross_phenotypes_selection) <- cross_phenotypes_selection[[next_id_col]]
  return(cross_phenotypes_selection)
}


do_permanova_corrected <- function(processed_data, cross_phenotypes_selection, permutation_num, correction_phenotype, next_id_col){
  # Generate the Aitchison distance matrices
  Distance <- vegdist(processed_data[[1]], method = "euclidean" )
  Distance <- data.frame(as.matrix((Distance)))
  # Run PERMANOVA to estimate the effect of phenotypes on overall bacterial composition after correcting for mode of delivery
  
  phenos1 <- processed_data[[2]][,!(colnames(processed_data[[2]]) %in% correction_phenotype)]
  
  adonis_results_corrected <- data.frame(Df = numeric(ncol(phenos1)),
                                         F = numeric(ncol(phenos1)),
                                         R2 = numeric(ncol(phenos1)),
                                         p_value = numeric(ncol(phenos1)))
  
  rownames(adonis_results_corrected) <- colnames(phenos1)
  
  for (i in 1:ncol(phenos1)) {
    print(i)  # delete after the run
    # Subset distance matrix and phenotypes, removing rows with NAs
    na_rows <- which(is.na(phenos1[, i]))
    
    distmat_cleaned <- Distance[!(rownames(Distance) %in% row.names(phenos1)[na_rows]),
                                !(colnames(Distance) %in% row.names(phenos1)[na_rows])]
    
    distmat_cleaned <- distmat_cleaned[(rownames(distmat_cleaned) %in% cross_phenotypes_selection[[next_id_col]]),
                                       (colnames(distmat_cleaned) %in% cross_phenotypes_selection[[next_id_col]])]
    
    phenos2 <- phenos1[row.names(phenos1) %in% row.names(distmat_cleaned), ]
    phenos2 <- phenos2[ match(rownames(distmat_cleaned), rownames(phenos2) ), ]
    
    cor_phenotypes <- cross_phenotypes_selection[row.names(cross_phenotypes_selection) %in% row.names(distmat_cleaned), ]
    cor_phenotypes <- cor_phenotypes[ match(rownames(distmat_cleaned), rownames(cor_phenotypes) ), ]
    cor_phenotypes <- cor_phenotypes[, (colnames(cor_phenotypes) %in% correction_phenotype)]
    print("Running adonis")  # delete after the run
    # Run adonis2
    ad1 <- adonis2(distmat_cleaned ~ cor_phenotypes + phenos2[[i]],
                   permutations = permutation_num,
                   parallel = 8,
                   na.action = na.fail,
                   by = "margin")
    # Save results
    adonis_results_corrected[i, ] <- c(ad1$Df[2], ad1$F[2], ad1$R2[2], ad1$"Pr(>F)"[2])
  }
  # Calculate FDR-adjusted p-values
  adonis_results_corrected$FDR <- p.adjust(adonis_results_corrected$p_value, method = "BH")
  
  return(adonis_results_corrected)
}


overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

do_lm_factor_phenos <- function(processed_data){
  results <- data.frame(factor = character(), pheno = character(), p = numeric(), R = numeric())
  
  for (i in c(1:20)){
    for (k in (1:ncol(processed_data[[2]]))){
      lm_model <- lm(processed_data[[1]][order(rownames(processed_data[[1]])), ][[i]] ~ 
                       processed_data[[2]][order(rownames(processed_data[[2]])), ][[k]])
      R <- summary(lm_model)$r.squared
      p <- overall_p(lm_model)
      vector <- c(colnames(processed_data[[1]][i]), colnames(processed_data[[2]][k]), R, p)
      results <- rbind(results, vector)
    }
  }
  colnames(results) <- c("factor", "pheno", "R", "p")
  return(results)
}

plot_tcam_results <- function(processed_data, phenotype, phenotype_values, plot_colors, label_x, label_y, factor1, factor2){
  all <- merge(processed_data[[2]], processed_data[[1]], by = "row.names")
  all <- all[!is.na(all[[phenotype]]), ]
  all[[phenotype]] <- factor(all[[phenotype]], phenotype_values)
  
  
  TCAM_plot <- ggplot(all, aes(.data[[factor1]], .data[[factor2]], color = .data[[phenotype]])) +
    geom_point(size = 2.5, alpha = 0.7) +
    stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = .data[[phenotype]], color = .data[[phenotype]]), linetype = 2, linewidth = 0.8) +
    xlab(label_x) +
    ylab(label_y) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    theme_bw() +
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, linewidth = 1.2, colour = "grey30"),
          axis.title.y = element_text(margin = margin(r = 10)),
          axis.title.x = element_text(margin = margin(t = 10)),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, colour = "grey30"),
          legend.position = "bottom") +
    guides(fill = "none")
  
  return(ggMarginal(TCAM_plot, type = "densigram", groupFill = T))
}

plot_r2_results <- function(data_permanova){
  FDR_sign_data <- data_permanova[data_permanova$FDR < 0.05, ]
  
  plot_r2 <- ggplot(FDR_sign_data, aes(y=R2, x= reorder(row.names(FDR_sign_data), FDR_sign_data$R2 ) )) + 
    geom_bar(stat = "identity", fill = "#086788") +
    coord_flip() +
    labs(x = "Phenotype") +
    theme_bw() +
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, linewidth = 1.2, colour = "grey30"),
          axis.title.y = element_text(margin = margin(r = 10)),
          axis.title.x = element_text(margin = margin(t = 10)),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, colour = "grey30"),
          legend.position = "bottom")
  return(plot_r2)
}

plot_loadings <- function(TCAM_loadings, factor){
  first_rows <- head(TCAM_loadings[order(TCAM_loadings[[factor]]),], 10)
  last_rows <- tail(TCAM_loadings[order(TCAM_loadings[[factor]]),], 10)
  selected_loadings <- rbind(first_rows, last_rows)
  plot_fin <- ggplot(selected_loadings, aes(y=.data[[factor]], x=reorder(X, .data[[factor]])) ) + 
    geom_bar(stat = "identity", fill = "#086788") +
    coord_flip() +
    labs(x = "Bug") +
    theme_bw() +
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, linewidth = 1.2, colour = "grey30"),
          axis.title.y = element_text(margin = margin(r = 10)),
          axis.title.x = element_text(margin = margin(t = 10)),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, colour = "grey30"),
          legend.position = "bottom")
  return(plot_fin)
}


################################################################################
# Loading raw data
################################################################################

cross_selection_corrected <- read.delim("masterfile_cross_sectional_selection_AFTER_correlation_&_correction_03_04_2024.txt")
cross_selection_corrected[cross_selection_corrected$variable_name == "birth_delivery_mode_complex", "TCAM"] <- 0
cross_selection_corrected[cross_selection_corrected$variable_name == "birth_delivery_mode_simple", "TCAM"] <- 0

cross_phenotypes <- read.delim("masterfile_cross_sectional_2023_11_15.txt")
cross_selection_mothers <- read.delim("masterfile_cross_sectional_selection_AFTER_correlation_&_correction_09_05_2024.txt")

TCAM_infant_W2_M1_M2_M3 <- read.delim("TCAM_infant_W2_M1_M2_M3.txt", sep = ",")
TCAM_infant_M6_M9_M12 <- read.delim("TCAM_infant_M6_M9_M12.txt", sep = ",")
TCAM_infant_M1_M3_M6_M12 <- read.delim("TCAM_infant_M1_M3_M6_M12.txt", sep = ",")
TCAM_mother_P12_B_M3 <- read.delim("TCAM_mother_P12_B_M3.txt", sep = ",")

TCAM_var_expl_infant_W2_M1_M2_M3 <- read.delim("TCAM_variance_expl_infant_W2_M1_M2_M3.txt", sep = ",")
TCAM_var_expl_infant_M6_M9_M12 <- read.delim("TCAM_variance_expl_infant_M6_M9_M12.txt", sep = ",")
TCAM_var_expl_infant_M1_M3_M6_M12 <- read.delim("TCAM_variance_expl_infant_M1_M3_M6_M12.txt", sep = ",")
TCAM_var_expl_mother_P12_B_M3 <- read.delim("TCAM_variance_expl_mother_P12_B_M3.txt", sep = ",")

TCAM_loadings_infant_W2_M1_M2_M3 <- read.delim("TCAM_loadings_infant_W2_M1_M2_M3.txt", sep = ",")
TCAM_loadings_infant_M6_M9_M12 <- read.delim("TCAM_loadings_infant_M6_M9_M12.txt", sep = ",")
TCAM_loadings_infant_M1_M3_M6_M12 <- read.delim("TCAM_loadings_infant_M1_M3_M6_M12.txt", sep = ",")
TCAM_loadings_mother_P12_B_M3 <- read.delim("TCAM_loadings_mother_P12_B_M3.txt", sep = ",")


################################################################################
# Preprocessing raw data
################################################################################from the selection after the correlation

processed_infant_W2_M1_M2_M3 <- do_preprocessing_infants(TCAM_infant_W2_M1_M2_M3, cross_phenotypes, cross_selection_corrected, c('W2', 'M1', 'M2', 'M3'))
processed_infant_W2_M1_M2_M3_vag_only <- do_preprocessing_infants_vagonly(TCAM_infant_W2_M1_M2_M3, cross_phenotypes, cross_selection_corrected, c('W2', 'M1', 'M2', 'M3'))

processed_infant_M6_M9_M12 <- do_preprocessing_infants(TCAM_infant_M6_M9_M12, cross_phenotypes, cross_selection_corrected, c('M6', 'M9', 'M12'))
processed_infant_M6_M9_M12_vag_only <- do_preprocessing_infants_vagonly(TCAM_infant_M6_M9_M12, cross_phenotypes, cross_selection_corrected, c('M6', 'M9', 'M12'))

processed_infant_M1_M3_M6_M12 <- do_preprocessing_infants(TCAM_infant_M1_M3_M6_M12, cross_phenotypes, cross_selection_corrected, c('M1', 'M3', 'M6', 'M12'))
processed_infant_M1_M3_M6_M12_vag_only <- do_preprocessing_infants_vagonly(TCAM_infant_M1_M3_M6_M12, cross_phenotypes, cross_selection_corrected, c('M1', 'M3', 'M6', 'M12'))

processed_mother_P12_B_M3 <- do_preprocessing_mothers(TCAM_mother_P12_B_M3, cross_phenotypes, cross_selection_mothers, c('P12', 'B', 'M3'))

################################################################################
# PERMANOVA analysis
################################################################################

permanova_simple_infant_W2_M1_M2_M3 <- do_permanova_simple(processed_infant_W2_M1_M2_M3, 5000)
permanova_simple_infant_M6_M9_M12 <- do_permanova_simple(processed_infant_M6_M9_M12, 5000)
permanova_simple_infant_M1_M3_M6_M12 <- do_permanova_simple(processed_infant_M1_M3_M6_M12, 5000)

permanova_simple_mother_P12_B_M3 <- do_permanova_simple(processed_mother_P12_B_M3, 5000)

permanova_simple_infant_W2_M1_M2_M3_vag_only <- do_permanova_simple(processed_infant_W2_M1_M2_M3_vag_only, 5000)
permanova_simple_infant_M6_M9_M12_vag_only <- do_permanova_simple(processed_infant_M6_M9_M12_vag_only, 5000)
permanova_simple_infant_M1_M3_M6_M12_vag_only <- do_permanova_simple(processed_infant_M1_M3_M6_M12_vag_only, 5000)

cross_phenotypes_selection_infant <- do_cor_phenotype_preselection(cross_phenotypes, "birth_deliverybirthcard_mode_binary", "next_id_infant")
cross_phenotypes_selection_mother <- do_cor_phenotype_preselection(cross_phenotypes, "birth_deliverybirthcard_mode_binary", "next_id_mother")

permanova_corrected_infant_W2_M1_M2_M3 <- do_permanova_corrected(processed_infant_W2_M1_M2_M3, cross_phenotypes_selection_infant, 5000, "birth_deliverybirthcard_mode_binary", "next_id_infant")
permanova_corrected_infant_M6_M9_M12 <- do_permanova_corrected(processed_infant_M6_M9_M12, cross_phenotypes_selection_infant, 5000, "birth_deliverybirthcard_mode_binary", "next_id_infant")
permanova_corrected_infant_M1_M3_M6_M12 <- do_permanova_corrected(processed_infant_M1_M3_M6_M12, cross_phenotypes_selection_infant, 5000, "birth_deliverybirthcard_mode_binary", "next_id_infant")

permanova_corrected_mother_P12_B_M3 <- do_permanova_corrected(processed_mother_P12_B_M3, cross_phenotypes_selection_mother, 5000, "birth_deliverybirthcard_mode_binary", "next_id_mother")

################################################################################
# Getting the best factor combinations to plot results
################################################################################

factor_vs_pheno_infant_W2_M1_M2_M3 <- do_lm_factor_phenos(processed_infant_W2_M1_M2_M3)
factor_vs_pheno_infant_M6_M9_M12 <- do_lm_factor_phenos(processed_infant_M6_M9_M12)
factor_vs_pheno_infant_M1_M3_M6_M12 <- do_lm_factor_phenos(processed_infant_M1_M3_M6_M12)

factor_vs_pheno_infant_M1_M3_M6_M12 <- do_lm_factor_phenos(processed_mother_P12_B_M3)

################################################################################
# Plotting results: coordinate plots factor1 vs factor2
################################################################################

plot_infant_W2_M1_M2_M3_DM <- plot_tcam_results(processed_infant_W2_M1_M2_M3, "birth_deliverybirthcard_mode_binary", c("VG", "CS"), c("#0055AA","#C40003"), "Factor 1 (7.97%)", "Factor 2 (3.73%)", "X0", "X1")
plot_infant_W2_M1_M2_M3_FM <- plot_tcam_results(processed_infant_W2_M1_M2_M3, "infant_birthcard_feeding_mode_after_birth", c("BF", "MF", "FF"), c("#f61067","#5e239d", "#00f0b5"), "Factor 1 (7.97%)", "Factor 2 (3.73%)", "X0", "X1")
plot_infant_W2_M1_M2_M3_EN <- plot_tcam_results(processed_infant_W2_M1_M2_M3, "infant_ffq_ever_never_breastfed", c("ever_BF", "never_BF"), c("#f61067", "#00f0b5"), "Factor 1 (7.97%)", "Factor 2 (3.73%)", "X0", "X1")
plot_infant_W2_M1_M2_M3_PA <- plot_tcam_results(processed_infant_W2_M1_M2_M3, "mother_birthcard_parity", c(0, 1, 2), c("#fe8181", "#f01819", "#400810"), "Factor 1 (7.97%)", "Factor 2 (3.73%)", "X0", "X1")


plot_infant_M6_M9_M12_PA <- plot_tcam_results(processed_infant_M6_M9_M12, "mother_birthcard_parity", c(0, 1, 2), c("#fe8181", "#f01819", "#400810"), "Factor 1 (5.21%)", "Factor 2 (2.98%)", "X0", "X1")

plot_infant_M1_M3_M6_M12_FM <- plot_tcam_results(processed_infant_M1_M3_M6_M12, "infant_birthcard_feeding_mode_after_birth", c("BF", "MF", "FF"), c("#f61067","#5e239d", "#00f0b5"), "Factor 1 (4.45%)", "Factor 2 (1.93%)", "X0", "X1")
plot_infant_M1_M3_M6_M12_PA <- plot_tcam_results(processed_infant_M1_M3_M6_M12, "mother_birthcard_parity", c(0, 1, 2), c("#fe8181", "#f01819", "#400810"), "Factor 1 (4.45%)", "Factor 2 (1.93%)", "X0", "X1")
plot_infant_M1_M3_M6_M12_DM <- plot_tcam_results(processed_infant_M1_M3_M6_M12, "birth_deliverybirthcard_mode_binary", c("VG", "CS"), c("#0055AA","#C40003"), "Factor 1 (4.45%)", "Factor 2 (1.93%)", "X0", "X1")
plot_infant_M1_M3_M6_M12_EN <- plot_tcam_results(processed_infant_M1_M3_M6_M12, "infant_ffq_ever_never_breastfed", c("ever_BF", "never_BF"), c("#f61067","#00f0b5"), "Factor 1 (4.45%)", "Factor 2 (1.93%)", "X0", "X1")

plot_mother_P12_B_M3_S <- plot_tcam_results(processed_mother_P12_B_M3, "mother_health_smoked_one_whole_year_p18", c(0, 1), c("#f61067","#00f0b5"), "Factor 1 (9.1%)", "Factor 2 (3.4%)", "X0", "X1")

plot_infant_W2_M1_M2_M3_DM
plot_infant_W2_M1_M2_M3_FM
plot_infant_W2_M1_M2_M3_EN
plot_infant_W2_M1_M2_M3_PA

plot_infant_M6_M9_M12_PA

plot_infant_M1_M3_M6_M12_FM
plot_infant_M1_M3_M6_M12_PA
plot_infant_M1_M3_M6_M12_DM
plot_infant_M1_M3_M6_M12_EN

plot_mother_P12_B_M3_S


################################################################################
# Plotting results: coordinate plots with customized factor choice (_cfc)
################################################################################

plot_infant_W2_M1_M2_M3_DM_cfc <- plot_tcam_results(processed_infant_W2_M1_M2_M3, "birth_deliverybirthcard_mode_binary", c("VG", "CS"), c("#0055AA","#C40003"), "Factor 1 (7.97%)", "Factor 6 (2.15%)", "X0", "X5")
plot_infant_W2_M1_M2_M3_FM_cfc <- plot_tcam_results(processed_infant_W2_M1_M2_M3, "infant_birthcard_feeding_mode_after_birth", c("BF", "MF", "FF"), c("#f61067","#5e239d", "#00f0b5"), "Factor 3 (3.16%)", "Factor 13 (1.3%)", "X2", "X12")
plot_infant_W2_M1_M2_M3_EN_cfc <- plot_tcam_results(processed_infant_W2_M1_M2_M3, "infant_ffq_ever_never_breastfed", c("ever_BF", "never_BF"), c("#f61067", "#00f0b5"), "Factor 3 (3.16%)", "Factor 10 (1.57%)", "X2", "X9")
plot_infant_W2_M1_M2_M3_PA_cfc <- plot_tcam_results(processed_infant_W2_M1_M2_M3, "mother_birthcard_parity", c(0, 1, 2), c("#fe8181", "#f01819", "#400810"), "Factor 1 (7.97%)", "Factor 2 (3.73%)", "X0", "X1")  # it's not a mistake - the highest correlations I have is between them


plot_infant_M6_M9_M12_PA_cfc <- plot_tcam_results(processed_infant_M6_M9_M12, "mother_birthcard_parity", c(0, 1, 2), c("#fe8181", "#f01819", "#400810"), "Factor 1 (5.21%)", "Factor 18 (0.81%)", "X0", "X17")

plot_infant_M1_M3_M6_M12_FM_cfc <- plot_tcam_results(processed_infant_M1_M3_M6_M12, "infant_birthcard_feeding_mode_after_birth", c("BF", "MF", "FF"), c("#f61067","#5e239d", "#00f0b5"), "Factor 14 (0.91%)", "Factor 19 (0.79%)", "X13", "X18")
plot_infant_M1_M3_M6_M12_PA_cfc <- plot_tcam_results(processed_infant_M1_M3_M6_M12, "mother_birthcard_parity", c(0, 1, 2), c("#fe8181", "#f01819", "#400810"), "Factor 1 (4.45%)", "Factor 4 (1.54%)", "X0", "X3")
plot_infant_M1_M3_M6_M12_DM_cfc <- plot_tcam_results(processed_infant_M1_M3_M6_M12, "birth_deliverybirthcard_mode_binary", c("VG", "CS"), c("#0055AA","#C40003"), "Factor 1 (4.45%)", "Factor 8 (1.3%)", "X0", "X7")
plot_infant_M1_M3_M6_M12_EN_cfc <- plot_tcam_results(processed_infant_M1_M3_M6_M12, "infant_ffq_ever_never_breastfed", c("ever_BF", "never_BF"), c("#f61067","#00f0b5"), "Factor 8 (1.3%)", "Factor 19 (0.79%)", "X7", "X18")


plot_mother_P12_B_M3_S_cfc <- plot_tcam_results(processed_mother_P12_B_M3, "mother_health_smoked_one_whole_year_p18", c(0, 1), c("#f61067","#00f0b5"), "Factor 2 (3.4%)", "Factor 8 (1.4%)", "X1", "X7")


plot_infant_W2_M1_M2_M3_DM_cfc
plot_infant_W2_M1_M2_M3_FM_cfc
plot_infant_W2_M1_M2_M3_EN_cfc
plot_infant_W2_M1_M2_M3_PA_cfc

plot_infant_M6_M9_M12_PA_cfc

plot_infant_M1_M3_M6_M12_FM_cfc
plot_infant_M1_M3_M6_M12_PA_cfc
plot_infant_M1_M3_M6_M12_DM_cfc
plot_infant_M1_M3_M6_M12_EN_cfc

plot_mother_P12_B_M3_S_cfc

################################################################################
# Plotting results: R2 per analysis
################################################################################

plot_r2_PS_infant_W2_M1_M2_M3 <- plot_r2_results(permanova_simple_infant_W2_M1_M2_M3)
plot_r2_PC_infant_W2_M1_M2_M3 <- plot_r2_results(permanova_corrected_infant_W2_M1_M2_M3)

plot_r2_PS_infant_M6_M9_M12 <- plot_r2_results(permanova_simple_infant_M6_M9_M12)
plot_r2_PC_infant_M6_M9_M12 <- plot_r2_results(permanova_corrected_infant_M6_M9_M12)

plot_r2_PS_infant_M1_M3_M6_M12 <- plot_r2_results(permanova_simple_infant_M1_M3_M6_M12)
plot_r2_PC_infant_M1_M3_M6_M12 <- plot_r2_results(permanova_corrected_infant_M1_M3_M6_M12)

plot_r2_PS_infant_W2_M1_M2_M3
plot_r2_PC_infant_W2_M1_M2_M3

plot_r2_PS_infant_M6_M9_M12
plot_r2_PC_infant_M6_M9_M12

plot_r2_PS_infant_M1_M3_M6_M12
plot_r2_PC_infant_M1_M3_M6_M12


################################################################################
# Plotting results: R2 overall
################################################################################

top_loadings_infant_W2_M1_M2_M3_X0 <- plot_loadings(TCAM_loadings_infant_W2_M1_M2_M3, "X0")
top_loadings_infant_W2_M1_M2_M3_X1 <- plot_loadings(TCAM_loadings_infant_W2_M1_M2_M3, "X1")
top_loadings_infant_W2_M1_M2_M3_X2 <- plot_loadings(TCAM_loadings_infant_W2_M1_M2_M3, "X2")
top_loadings_infant_W2_M1_M2_M3_X5 <- plot_loadings(TCAM_loadings_infant_W2_M1_M2_M3, "X5")
top_loadings_infant_W2_M1_M2_M3_X9 <- plot_loadings(TCAM_loadings_infant_W2_M1_M2_M3, "X9")
top_loadings_infant_W2_M1_M2_M3_X12 <- plot_loadings(TCAM_loadings_infant_W2_M1_M2_M3, "X12")


top_loadings_infant_M6_M9_M12_X0 <- plot_loadings(TCAM_loadings_infant_M6_M9_M12, "X0")
top_loadings_infant_M6_M9_M12_X17 <- plot_loadings(TCAM_loadings_infant_M6_M9_M12, "X17")

top_loadings_infant_M1_M3_M6_M12_X13 <- plot_loadings(TCAM_loadings_infant_M1_M3_M6_M12, "X13")
top_loadings_infant_M1_M3_M6_M12_X18 <- plot_loadings(TCAM_loadings_infant_M1_M3_M6_M12, "X18")

top_loadings_infant_M1_M3_M6_M12_X7 <- plot_loadings(TCAM_loadings_infant_M1_M3_M6_M12, "X7")
top_loadings_infant_M1_M3_M6_M12_X3 <- plot_loadings(TCAM_loadings_infant_M1_M3_M6_M12, "X3")


################################################################################
# Saving results
################################################################################

write.table(permanova_simple_infant_W2_M1_M2_M3,"permanova_simple_infant_W2_M1_M2_M3_new_pheno_selection.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(permanova_simple_infant_M6_M9_M12,"permanova_simple_infant_M6_M9_M12_new_pheno_selection.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(permanova_simple_infant_M1_M3_M6_M12,"permanova_simple_infant_M1_M3_M6_M12_new_pheno_selection.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)

write.table(permanova_simple_infant_W2_M1_M2_M3_vag_only,"permanova_simple_infant_W2_M1_M2_M3_new_pheno_selection_vagonly.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(permanova_simple_infant_M6_M9_M12_vag_only,"permanova_simple_infant_M6_M9_M12_new_pheno_selection_vagonly.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(permanova_simple_infant_M1_M3_M6_M12_vag_only,"permanova_simple_infant_M1_M3_M6_M12_new_pheno_selection_vagonly.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)


write.table(permanova_corrected_infant_W2_M1_M2_M3,"permanova_corrected_infant_W2_M1_M2_M3_new_pheno_selection.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(permanova_corrected_infant_M6_M9_M12,"permanova_corrected_infant_M6_M9_M12_new_pheno_selection.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(permanova_corrected_infant_M1_M3_M6_M12,"permanova_corrected_infant_M1_M3_M6_M12_new_pheno_selection.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)



write.table(permanova_simple_mother_P12_B_M3,"permanova_simple_mother_P12_B_M3.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(permanova_corrected_mother_P12_B_M3,"permanova_corrected_mother_P12_B_M3.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)



























































