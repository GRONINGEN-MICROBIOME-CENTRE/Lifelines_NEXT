################################################################################
### TCAM analysis
### Author(s): Nataliia Kuzub, Asier Fernández-Pato, Trishla Sinha
### Last updated: 13th December, 2024
################################################################################

################################################################################
# Loading libraries
################################################################################
library(vegan)
library(dplyr)
library(ggplot2)
library(readr)
library(ggExtra)

################################################################################
# Loading functions
################################################################################

do_preprocessing_infants <- function(TCAM, cross_phenotypes, timepoints){
  # TCAM table preprocesing: output of the TCAM.ipynb
  rownames(TCAM) <- TCAM$X
  TCAM <- subset(TCAM, select = -c(X))
  
  # Fitering out the samples excluded for this analysis
  cross_phenotypes$next_id_infant <- gsub("^FAM", "NEXT", cross_phenotypes$Family)
  ## Dropping all the dynamic phenotypes from the df
  cross_phenotypes <- cross_phenotypes[cross_phenotypes$next_id_infant %in% row.names(TCAM), 
                                       !(colnames(cross_phenotypes) %in% c("Family", "Type", "Timepoint", "Binary_phenotype_dynamic", "Factor_phenotype_dynamic", "Quant_phenotype_dynamic", "Covariate1", "Covariate2"))]
  row.names(cross_phenotypes) <- NULL
  infant_cross_phenotypes <- cross_phenotypes[!duplicated(cross_phenotypes), ]
  row.names(infant_cross_phenotypes) <- infant_cross_phenotypes$next_id_infant
  
  # Preparing phenotypes for associations
  to_exclude <- c("FAMILY", "next_id_mother", "next_id_infant", "next_id_partner", "infant_relations", "sibling_number", "twin_pair")
  infant_cross_phenotypes_filtered <- infant_cross_phenotypes[,!(colnames(infant_cross_phenotypes) %in% to_exclude)]
  
  infant_cross_phenotypes_filtered <- infant_cross_phenotypes_filtered %>%
    mutate(across(where(is.character), as.factor))
  
  return(list(TCAM, infant_cross_phenotypes_filtered))
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
  
  plot_r2 <- ggplot(FDR_sign_data, aes(y=R2, x= reorder(row.names(FDR_sign_data), R2) )) + 
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

cross_phenotypes <- read.delim("simulated_metadata.txt")
TCAM_infant_W2_M1_M3 <- read.delim("TCAM_infant_W2_M1_M3.txt", sep = ",")
TCAM_var_expl_infant_W2_M1_M3 <- read.delim("TCAM_variance_expl_infant_W2_M1_M3.txt", sep = ",")
TCAM_loadings_infant_W2_M1_M3 <- read.delim("TCAM_loadings_infant_W2_M1_M3.txt", sep = ",")

################################################################################
# Preprocessing raw data
################################################################################from the selection after the correlation

processed_infant_W2_M1_M3 <- do_preprocessing_infants(TCAM_infant_W2_M1_M3, cross_phenotypes, c('W2', 'M1', 'M3'))

################################################################################
# PERMANOVA analysis
################################################################################

permanova_simple_infant_W2_M1_M3 <- do_permanova_simple(processed_infant_W2_M1_M3, 5000)

################################################################################
# Getting the best factor combinations to plot results
################################################################################

factor_vs_pheno_infant_W2_M1_M3 <- do_lm_factor_phenos(processed_infant_W2_M1_M3)

################################################################################
# Plotting results: coordinate plots factor1 vs factor2
################################################################################

plot_infant_W2_M1_M2_M3_BPS <- plot_tcam_results(processed_infant_W2_M1_M3, "Binary_phenotype_static", c("BinaryStatic1", "BinaryStatic2"), c("#0055AA","#C40003"), "Factor 1 (0.46%)", "Factor 2 (0.45%)", "X0", "X1")
plot_infant_W2_M1_M2_M3_BPS

################################################################################
# Plotting results: R2 per analysis
################################################################################

plot_r2_PS_infant_W2_M1_M3 <- plot_r2_results(permanova_simple_infant_W2_M1_M3)
plot_r2_PS_infant_W2_M1_M3  # Empty plot, since there was no FDR < 0.05 in the input data

################################################################################
# Plotting results: R2 overall
################################################################################

top_loadings_infant_W2_M1_M3_X0 <- plot_loadings(TCAM_loadings_infant_W2_M1_M3, "X0")
top_loadings_infant_W2_M1_M3_X0

################################################################################
# Saving results
################################################################################

write.table(permanova_simple_infant_W2_M1_M3,"permanova_simple_infant_W2_M1_M3_static_pheno_selection.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
