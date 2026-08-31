# ────────────────────────────────────────────────
# Libraries
# ────────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(pROC)
library(lme4)
library(foreach)
library(patchwork)
library(scales)
library(ggbeeswarm)
library(grid)

# ────────────────────────────────────────────────
# Load data
# ────────────────────────────────────────────────

foodseq_plant_infant <- read.csv("foodseq_plant_regrouped_binary.csv", row.names = 1) # only infants
foodseq_animal_all <- read.csv("foodseq_animal_asv_binary.csv", row.names = 1)
LLNEXT.metadata <- read.delim("LLNEXT_metadata_family_structure_28_08_2023.txt")
ffq.data <- read.delim("masterfile_longitudinal_2023_09_29.txt",header=T,sep="\t")

# ────────────────────────────────────────────────
# Keep only infant data and sample IDs
# ────────────────────────────────────────────────

# Keep only infant animal FoodSeq
foodseq_animal_all$NG_ID <- rownames(foodseq_animal_all)
foodseq_animal_all_meta <- left_join(foodseq_animal_all, LLNEXT.metadata[,c("NG_ID", "SAMPLE_ID", "Type")], by="NG_ID")
foodseq_animal_infant <- foodseq_animal_all_meta[foodseq_animal_all_meta$Type=="infant",]

# Keep only infant metadata
metadata_selected <- LLNEXT.metadata[LLNEXT.metadata$NG_ID %in% foodseq_animal_infant$NG_ID, c("NG_ID", "Type", "FAMILY", "Timepoint_categorical", "SAMPLE_ID")]

# Set row names and delete metadata columns
rownames(foodseq_animal_infant) <- foodseq_animal_infant$SAMPLE_ID
foodseq_animal_infant <- foodseq_animal_infant[, !colnames(foodseq_animal_infant) %in% c("SAMPLE_ID", "NG_ID", "Type")]

# Clean environment
rm(foodseq_animal_all, foodseq_animal_all_meta, LLNEXT.metadata)

# ────────────────────────────────────────────────
# Set prevalence filter
# ────────────────────────────────────────────────

# Plant
low_prev_columns <- names(foodseq_plant_infant)[colSums(foodseq_plant_infant) < 7]
foodseq_plant_infant_prev <- foodseq_plant_infant[ , !(colnames(foodseq_plant_infant) %in% low_prev_columns)]

# Animal
low_prev_columns <- names(foodseq_animal_infant)[colSums(foodseq_animal_infant) < 7]
foodseq_animal_prev <- foodseq_animal_infant[ , !(colnames(foodseq_animal_infant) %in% low_prev_columns)]

# Create "meat" variable
foodseq_animal_prev$Meat <- ifelse(rowSums(foodseq_animal_prev[, c("Cattle", "Chicken", "Pig", "Turkey")] == 1) >= 1, 1, 0)

# Clean environment
rm(foodseq_plant_infant, foodseq_animal_infant)

# ────────────────────────────────────────────────
# Combine plant and animal foodseq
# ────────────────────────────────────────────────

foodseq_plant_infant_prev$SAMPLE_ID <- rownames(foodseq_plant_infant_prev)
foodseq_animal_prev$SAMPLE_ID <- rownames(foodseq_animal_prev)

foodseq_infant <- left_join(foodseq_plant_infant_prev, foodseq_animal_prev, by="SAMPLE_ID")

rm(foodseq_plant_infant_prev, foodseq_animal_prev)

# ────────────────────────────────────────────────
# FFQ data processing
# ────────────────────────────────────────────────

# Keep only infant samples in the same order as FoodSeq
ffq_infant <- ffq.data[match(foodseq_infant$SAMPLE_ID, ffq.data$SAMPLE_ID),]

# Keep SAMPLE_IDs as rownames in foodseq and ffq
rownames(foodseq_infant) <- foodseq_infant$SAMPLE_ID
foodseq_infant$SAMPLE_ID <- NULL
rownames(ffq_infant) <- ffq_infant$SAMPLE_ID
ffq_infant$SAMPLE_ID <- NULL

# Recode to numeric
ffq_binary <- ffq_infant %>%
  mutate(across(everything(), ~ case_when(
    .x == "yes" ~ 1,
    .x == "no" ~ 0,
    .x == "excl_BF" ~ 0,
    .x == "excl_FF" ~ 1,
    .x == "development_delay" ~ 1,
    .x == "development_typical" ~ 0,
    .x == "1_medium_high_GI_distress" ~ 1,
    .x == "0_no_little_GI_distress" ~ 0,
    .x == "never_BF" ~ 2,
    TRUE ~ as.numeric(.x)
  )))

# Prevalence filtering
low_prev_columns <- names(ffq_binary)[colSums(ffq_binary > 0, na.rm = T) < 7]
ffq_filtered <- ffq_binary[, !(names(ffq_binary) %in% low_prev_columns)]

# Clean environment
rm(ffq.data, ffq_infant, ffq_binary)


# ────────────────────────────────────────────────
# Select FFQ items for correlations
# ────────────────────────────────────────────────

plant_patterns <- c(
  "apple", "banana",
  "legume", "citrus", "nut",
  "potato", "rice", "_breads")

meat_patterns <- c("fish", "dairy", "cream", "buttermilk", "quark", "meat", "egg", "cheese", "feeding", "formula")


select_by_patterns <- function(x, patterns, ignore.case = T) {
  pattern <- paste(patterns, collapse = "|")
  x[, grepl(pattern, names(x), ignore.case = ignore.case)]
}

ffq_plants <- select_by_patterns(ffq_filtered, plant_patterns)
ffq_meat <- select_by_patterns(ffq_filtered, meat_patterns)

animal_exclude_ffq <- c("infant_ffq_meat_substitute",
                        "infant_ffq_meat_substitute_freq_weighted",
                        "infant_ffq_group_meat_substitute_gday")

ffq_meat <- ffq_meat[, !(names(ffq_meat) %in% animal_exclude_ffq)]


# ────────────────────────────────────────────────
# Select FoodSeq items for correlations
# ────────────────────────────────────────────────

plant_patterns_asv <- c("apple", "banana", "beans", "citrus", "nut", "potato", "rice", "wheat")

meat_patterns_asv <- c("atlantic", "bison", "cattle", "chicken", "pig", "turkey", "meat")

foodseq_plants <- select_by_patterns(foodseq_infant, plant_patterns_asv)
foodseq_meat <- select_by_patterns(foodseq_infant, meat_patterns_asv)


# ────────────────────────────────────────────────
# New column names for foodseq
# ────────────────────────────────────────────────
#
new_names <- c(
  Apples.and.pears = "apple",
  banana..NA..fehi.banana..hairy.banana..plantain = "banana",
  Beans..peas..and.olives = "legume",
  Citrus = "citrus",
  Nuts = "nuts",
  Potatoes.and.sweet.potatoes = "potato",
  Rice.and.grains = "grains",
  Wheat = "Wheat"
)

names(foodseq_plants)[match(names(new_names), names(foodseq_plants))] <- new_names


# ────────────────────────────────────────────────
# Combine plant and animal FoodSeq
# ────────────────────────────────────────────────

foodseq_plants$SAMPLE_ID <- rownames(foodseq_plants)
foodseq_meat$SAMPLE_ID <- rownames(foodseq_meat)

foodseq_binary <- merge(foodseq_plants, foodseq_meat, by="SAMPLE_ID")


# ────────────────────────────────────────────────
# Combine plant and animal FFQ
# ────────────────────────────────────────────────

ffq_plants$SAMPLE_ID <- rownames(ffq_plants)
ffq_meat$SAMPLE_ID <- rownames(ffq_meat)

ffq_numeric <- merge(ffq_plants, ffq_meat, by="SAMPLE_ID")


# ────────────────────────────────────────────────
# Metadata
# ────────────────────────────────────────────────

metadata = data.frame(do.call(rbind,strsplit(split = "_",foodseq_binary[,"SAMPLE_ID"])))
colnames(metadata) = c("ID","Time")


# ────────────────────────────────────────────────
# Only retain relevant FoodSeq and FFQ
# ────────────────────────────────────────────────

foodseq_binary <- foodseq_binary[, !colnames(foodseq_binary) %in% c("SAMPLE_ID", "Bison")]

ffq_numeric <- ffq_numeric[, !colnames(ffq_numeric) %in% c("SAMPLE_ID", 
                                                           "infant_ffq_breastfeeding_freq_weighted",
                                                           "infant_ffq_feeding_mode_cups_numeric",
                                                           "infant_ffq_stopped_breastfeeding",
                                                           "infant_ffq_grains_to_formula_follow_on_milk",
                                                           "infant_ffq_grains_to_formula_follow_on_milk_freq_weighted")]


# ────────────────────────────────────────────────
# Association Foodseq-FFQ ------------------------
# ────────────────────────────────────────────────

invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))} 




glmer_output = foreach(i = 1:ncol(ffq_numeric),.combine = rbind) %:%
  foreach( j = 1:ncol(foodseq_binary),.combine = rbind) %do% {
    
    #glmer_output = foreach(i = 1,.combine = rbind) %:%
    #  foreach( j = 13:14,.combine = rbind) %do% {
    print (paste(i,j))
    g1 = try(glmer(foodseq_binary[,j] ~ invrank(ffq_numeric[,i])  + (1|ID),family = binomial, data = metadata,nAGQ = 5,
                   control = glmerControl(
                     optimizer = "bobyqa",
                     optCtrl = list(maxfun = 1e5)
                   )))
    if(class(g1)[1]!= "try-error") {s1 = summary(g1)
    if(length(g1@optinfo$conv$lme4)==0){
      message = NA
      code = NA
    } else  {
      message = g1@optinfo$conv$lme4[["messages"]]
      message = ifelse(is.null(message),NA,unlist(message))
      code = g1@optinfo$conv$lme4[["code"]]
      code = ifelse(is.null(code),NA,code)
    }
    OR = exp(s1$coef[2,1])
    confint1 = confint(g1,method = "Wald")["invrank(ffq_numeric[, i])",]
    data.frame(Foodseq = colnames(foodseq_binary)[j],
               FFQ = colnames(ffq_numeric)[i],
               ffq_type = ifelse(length(table(ffq_numeric[,1]))<3,"bool","numeric"),
               success = "glmer",
               message = message,
               code = code,
               OR = OR,
               OR2.5CI = exp(confint1[1]),
               OR97.5CI = exp(confint1[2]),
               Beta = s1$coef[2,1],
               SE = s1$coef[2,2],
               Z = s1$coef[2,3],
               P = s1$coef[2,4])
    } else {
      data.frame(Foodseq = colnames(foodseq_binary)[j],
                 FFQ = colnames(ffq_numeric)[i],
                 ffq_type = ifelse(length(table(ffq_numeric[,1]))<3,"bool","numeric"),
                 success = F,
                 message = NA,
                 OR =NA,
                 OR2.5CI = NA,
                 OR97.5CI = NA,
                 code = NA,
                 Beta = NA,
                 SE = NA,
                 Z = NA,
                 P = NA
      )
    }
  }


glm_output = foreach(i = 1:ncol(ffq_numeric),.combine = rbind) %:%
  foreach( j = 1:ncol(foodseq_binary),.combine = rbind) %do% {
    
    #glmer_output = foreach(i = 1,.combine = rbind) %:%
    #  foreach( j = 13:14,.combine = rbind) %do% {
    print (paste(i,j))
    w<-NULL
    g1 = withCallingHandlers(
      glm(foodseq_binary[,j] ~ invrank(ffq_numeric[,i]),family = binomial),
      warning = function(wrn) {
        w <<- c(w, conditionMessage(wrn))
        invokeRestart("muffleWarning")
      }
    )
    s1 = summary(g1)
    confint1 = confint(g1)["invrank(ffq_numeric[, i])",]
    data.frame(Foodseq = colnames(foodseq_binary)[j],
               FFQ = colnames(ffq_numeric)[i],
               ffq_type = ifelse(length(table(ffq_numeric[,1]))<3,"bool","numeric"),
               success = "glm",
               message = ifelse(is.null(w),NA,w[1]),
               code = code,
               OR =exp(s1$coef[2,1]),
               OR2.5CI = exp(confint1[1]),
               OR97.5CI = exp(confint1[2]),
               Beta = s1$coef[2,1],
               SE = s1$coef[2,2],
               Z = s1$coef[2,3],
               P = s1$coef[2,4])
  }

final_output = glmer_output
final_output[!is.na(glmer_output$message),] = glm_output[!is.na(glmer_output$message),]


# ────────────────────────────────────────────────
# Add AUC calculation
# ────────────────────────────────────────────────

auc_result <- matrix(
  NA,
  nrow = ncol(foodseq_binary),
  ncol = ncol(ffq_numeric),
  dimnames = list(colnames(foodseq_binary), colnames(ffq_numeric))
)

for (i in 1:ncol(foodseq_binary)) {
  for (j in 1:ncol(ffq_numeric)) {
    
    non_na_rows <- !is.na(foodseq_binary[[i]]) & !is.na(ffq_numeric[[j]])
    
    if (sum(non_na_rows) > 0 &&
        length(unique(foodseq_binary[[i]][non_na_rows])) == 2) {
      auc_result[i, j] <- as.numeric(
        auc(
          response = foodseq_binary[[i]][non_na_rows],
          predictor = ffq_numeric[[j]][non_na_rows],
          quiet = T)
        )
    }
  }
}

# add auc values to final_output
final_output$AUC <- mapply(
  function(foodseq, ffq) auc_result[foodseq, ffq],
  final_output$Foodseq,
  final_output$FFQ
)


# ────────────────────────────────────────────────
# Prior FDR select relevant associations
# ────────────────────────────────────────────────

associations_to_keep <- data.frame(
  Foodseq = c("apple", "apple",
              "banana",
              "legume", "legume", "legume",
              "citrus", "citrus",
              "nuts", "nuts", "nuts",
              "potato", "potato", "potato", "potato", "potato",
              "grains", "grains", "grains",
              "Wheat", "Wheat", "Wheat", "Wheat", "Wheat", "Wheat", "Wheat", "Wheat", "Wheat",
              
              "Atlantic.salmon", "Atlantic.salmon", "Atlantic.salmon",
              
              "Cattle", "Cattle", "Cattle", "Cattle", "Cattle", "Cattle",
              "Cattle", "Cattle", "Cattle", "Cattle", "Cattle", "Cattle",
              
              "Cattle", "Cattle",
              
              "Chicken", "Chicken", "Chicken", 
              
              "Chicken", "Chicken", "Chicken", "Chicken", "Chicken",
              
              "Pig", "Pig", "Pig", "Pig", "Pig",
              
              "Turkey", "Turkey", "Turkey", "Turkey", "Turkey",
              
              "Meat", "Meat", "Meat", "Meat", "Meat",
              
              "Cattle", "Cattle", "Cattle", "Cattle", "Cattle"),
  
  FFQ     = c("infant_ffq_fresh_fruits_apple",
              "infant_ffq_fresh_fruits_apple_weighted",
              
              "infant_ffq_fresh_fruits_banana_weighted",
              
              "infant_ffq_legumes",
              "infant_ffq_legumes_freq_weighted",
              "infant_ffq_group_legumes_gday",
              
              "infant_ffq_fresh_fruits_citrus",
              "infant_ffq_fresh_fruits_citrus_weighted",
              
              "infant_ffq_peanut_butter_and_other_nut_spreads",
              "infant_ffq_peanut_butter_and_other_nut_spreads_freq_weighted",
              "infant_ffq_group_nuts_gday",
              
              "infant_ffq_group_potatoes_gday",
              "infant_ffq_potato_boiled_or_mashed",
              "infant_ffq_potato_boiled_or_mashed_freq_weighted",
              "infant_ffq_potato_fried_or_baked",
              "infant_ffq_potato_fried_or_baked_freq_weighted",
              
              "infant_ffq_rice_and_grains",
              "infant_ffq_rice_and_grains_freq_weighted",
              "infant_ffq_group_grain_products_and_rice_gday",
              
              "infant_ffq_group_breads_gday",
              "infant_ffq_breads",
              "infant_ffq_breads_freq_weighted",
              "infant_ffq_breads_brown",
              "infant_ffq_breads_brown_weighted",
              "infant_ffq_breads_white",
              "infant_ffq_breads_white_weighted",
              "infant_ffq_breads_whole_grain",
              "infant_ffq_breads_whole_grain_weighted",
              
              "infant_ffq_group_fish_gday",
              "infant_ffq_fish_freq_weighted",
              "infant_ffq_fish",
              
              "infant_ffq_group_dairy_gday",
              "infant_ffq_cheese",
              "infant_ffq_group_cheese_gday",
              "infant_ffq_cheese_freq_weighted",
              "infant_ffq_milk_buttermilk",
              "infant_ffq_milk_buttermilk_freq_weighted",
              "infant_ffq_whipped_cream",
              "infant_ffq_whipped_cream_freq_weighted",
              "infant_ffq_yoghurt_or_quark",
              "infant_ffq_yoghurt_or_quark_freq_weighted",
              "infant_ffq_ice_cream",
              "infant_ffq_ice_cream_freq_weighted",
              
              "infant_ffq_feeding_mode_simple",
              "infant_ffq_formula_follow_on_freq_weighted",
              
              "infant_ffq_eggs",
              "infant_ffq_group_egg_gday",
              "infant_ffq_eggs_freq_weighted",
              
              "infant_ffq_group_meat_gday",
              "infant_ffq_meat",
              "infant_ffq_meat_freq_weighted",
              "infant_ffq_meat_cold_cuts",
              "infant_ffq_meat_cold_cuts_freq_weighted",
              
              "infant_ffq_group_meat_gday",
              "infant_ffq_meat",
              "infant_ffq_meat_freq_weighted",
              "infant_ffq_meat_cold_cuts",
              "infant_ffq_meat_cold_cuts_freq_weighted",
              
              "infant_ffq_group_meat_gday",
              "infant_ffq_meat",
              "infant_ffq_meat_freq_weighted",
              "infant_ffq_meat_cold_cuts",
              "infant_ffq_meat_cold_cuts_freq_weighted",
              
              "infant_ffq_group_meat_gday",
              "infant_ffq_meat",
              "infant_ffq_meat_freq_weighted",
              "infant_ffq_meat_cold_cuts",
              "infant_ffq_meat_cold_cuts_freq_weighted",
              
              "infant_ffq_group_meat_gday",
              "infant_ffq_meat",
              "infant_ffq_meat_freq_weighted",
              "infant_ffq_meat_cold_cuts",
              "infant_ffq_meat_cold_cuts_freq_weighted")
)

final_output_subset <- final_output %>%
  semi_join(associations_to_keep, by = c("Foodseq", "FFQ"))


# ────────────────────────────────────────────────
# Permutation FDR
# ────────────────────────────────────────────────

set.seed(143)

perm_res2 = lapply(1:1000,\(x) {
  foodseq_perm = foodseq_binary[sample(1:nrow(foodseq_binary)),]
  foreach (i = 1:nrow(final_output_subset),.combine = c) %do%{
    g1 = glm(foodseq_perm[,final_output_subset[i,1]] ~ invrank(ffq_numeric[,final_output_subset[i,2]]),family = binomial)
    summary(g1)$coef[2,4]
    
  }
})
s1 = do.call(cbind,perm_res2)
final_output_subset$Perm = NA
for(i in 1:nrow(final_output_subset)){
  P = sum(final_output_subset$P<= final_output_subset$P[i])
  FP = sum(s1 <=final_output_subset$P[i]) / 1000
  fdr = FP/P
  if(fdr>1) fdr = 1
  final_output_subset$Perm[i] = fdr
}

# If FFQ-FoodSeq food pairs are biologically plausible, assign "matched" 
final_output$Paired = ifelse (paste(final_output$Foodseq,final_output$FFQ) %in% paste(final_output_subset$Foodseq,final_output_subset$FFQ),"matched","random")



# ────────────────────────────────────────────────
# general plot
# ────────────────────────────────────────────────

cols <- c(
  "random" = "#4E79A7",
  "matched" = "#E15759"
)

w.test.stat = wilcox.test(final_output$Z ~ final_output$Paired)

p_density <- ggplot(final_output, aes(x = Z, fill = Paired, color = Paired)) +
  geom_density(alpha = 0.35, linewidth = 1) +
  theme_classic(base_size = 14) +
  labs(x = NULL, y = "Density") + 
  scale_color_manual(name = "FFQ-Foodseq",values = cols) + 
  scale_fill_manual(name = "FFQ-Foodseq",values = cols) +
  theme(legend.position = "top")
p_density <- p_density + 
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0("Wilcoxon p = ", signif(w.test.stat$p.value, 3)),
    hjust = 1.1,
    vjust = 1.5,
    size = 4.5
  )

p_box <- ggplot(final_output, aes(x = Z, y = Paired, fill = Paired, color = Paired)) +
  geom_boxplot(width = 0.45, alpha = 0.7, outlier.alpha = 0.4,outlier.shape = NA) +
  geom_jitter(
    height = 0.12,
    width = 0,
    size = 0.7,
    alpha = 0.3
  ) + 
  theme_classic(base_size = 14) +
  labs(x = "T statistics", y = NULL) +
  scale_color_manual(name = "FFQ-Foodseq",values = cols) + 
  scale_fill_manual(name = "FFQ-Foodseq",values = cols) + 
  theme(legend.position = "none")

p_out <- p_density / p_box + plot_layout(heights = c(3, 1))
p_out

ggsave("foodseq_ffq_pairs.tiff", plot = p_out, dpi = 600, width = 5.8, height = 4.5, units = "in", compression = "lzw", bg = "white")



# ────────────────────────────────────────────────
# Associations shown individually
# ────────────────────────────────────────────────

fill_pal <- c("Not detected" = "#C4C4C4", "Detected" = "#91D0FA")
line_pal <- c("Not detected" = "#6D6D6D", "Detected" = "#1C6EA7")

cat_pal <- c(
  "No" = "#C4C4C4",
  "Yes" = "#91D0FA",
  "Exclusively human milk-fed" = "#C4C4C4",
  "Exclusively formula-fed" = "#91D0FA"
)

# Theme
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#EFEAE3", linewidth = 0.45),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size),
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1.2, lineheight = 0.95),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1.3),
      legend.key.height = unit(0.42, "cm"),
      legend.key.width = unit(0.42, "cm"),
      plot.margin = margin(14, 14, 14, 10)
    )
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

assoc_label <- function(foodseq_var, ffq_var) {
  x <- final_output_subset %>%
    filter(Foodseq == foodseq_var, FFQ == ffq_var) %>%
    slice(1)
  
  paste0(
    "FDR = ", format_p(x$Perm)
  )
}

num_plot <- function(foodseq_var, ffq_var, title, ylab) {
  df <- data.frame(
    FoodSeq = factor(
      foodseq_binary[[foodseq_var]],
      levels = c(0, 1),
      labels = c("Not detected", "Detected")
    ),
    FFQ = ffq_numeric[[ffq_var]]
  ) %>%
    filter(!is.na(FFQ))
  
  ggplot(df, aes(FoodSeq, FFQ)) +
    geom_violin(aes(fill = FoodSeq), trim = FALSE, alpha = 0.45,
                color = NA, width = 0.88) +
    geom_boxplot(aes(fill = FoodSeq, color = FoodSeq), width = 0.55,
                 outlier.shape = NA, size = 0.7,
                 median.linewidth = 1.2, alpha = 0.92) +
    geom_quasirandom(aes(fill = FoodSeq, color = FoodSeq), groupOnX = TRUE,
                     width = 0.11, varwidth = FALSE, shape = 21,
                     size = 1.4, stroke = 0.55, alpha = 0.85) +
    scale_fill_manual(values = fill_pal) +
    scale_color_manual(values = line_pal) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
    labs(
      x = "DNA detection",
      y = ylab,
      title = title,
      subtitle = assoc_label(foodseq_var, ffq_var)
    ) +
    theme_pub() +
    theme(legend.position = "none")
}

cat_plot <- function(foodseq_var, ffq_var, title, fill_lab, labels) {
  df <- data.frame(
    FoodSeq = factor(
      foodseq_binary[[foodseq_var]],
      levels = c(0, 1),
      labels = c("Not detected", "Detected")
    ),
    FFQ = factor(
      ffq_numeric[[ffq_var]],
      levels = c(0, 1),
      labels = labels
    )
  ) %>%
    filter(!is.na(FFQ))
  
  df_sum <- df %>%
    count(FoodSeq, FFQ, name = "n") %>%
    group_by(FoodSeq) %>%
    mutate(
      prop = n / sum(n),
      label = ifelse(prop >= 0.08, paste0(round(prop * 100), "%"), "")
    ) %>%
    ungroup()
  
  ggplot(df_sum, aes(FoodSeq, prop, fill = FFQ)) +
    geom_col(width = 0.62, color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5),
              size = 3, color = "black") +
    scale_fill_manual(values = cat_pal, name = fill_lab) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.03))) +
    labs(
      x = "DNA detection",
      y = "Infants, %",
      title = title,
      subtitle = assoc_label(foodseq_var, ffq_var)
    ) +
    theme_pub() +
    theme(legend.position = "right")
}


p1 <- num_plot(
  "Cattle",
  "infant_ffq_group_dairy_gday",
  "Cattle DNA vs Dairy Intake (g/day)",
  "Dairy intake, g/day"
)

p2 <- num_plot(
  "Cattle",
  "infant_ffq_formula_follow_on_freq_weighted",
  "Cattle DNA vs Follow-on formula (freq)",
  "Follow-on formula intake,\nfrequency"
)

p3 <- num_plot(
  "potato",
  "infant_ffq_potato_boiled_or_mashed_freq_weighted",
  "Potato DNA vs Boiled/mashed Potatoes (freq)",
  "Boiled/mashed potato intake,\nfrequency"
)

p4 <- cat_plot(
  "Cattle",
  "infant_ffq_feeding_mode_simple",
  "Cattle DNA vs Feeding Mode",
  "Feeding mode",
  c("Exclusively human milk-fed", "Exclusively formula-fed")
)

p5 <- cat_plot(
  "Pig",
  "infant_ffq_meat_cold_cuts",
  "Pig DNA vs Meat Cold Cuts",
  "Cold cuts",
  c("No", "Yes")
)

p6 <- num_plot(
  "nuts",
  "infant_ffq_peanut_butter_and_other_nut_spreads_freq_weighted",
  "Nut DNA vs Nut Spreads (freq)",
  "Nut spread intake,\nfrequency"
)

# Additional associations
p7 <- num_plot(
  "potato",
  "infant_ffq_group_potatoes_gday",
  "Potato DNA vs Potato Intake (g/day)",
  "Potato intake, g/day"
)

p8 <- num_plot(
  "Pig",
  "infant_ffq_meat_cold_cuts_freq_weighted",
  "Pig DNA vs Meat Cold Cuts (freq)",
  "Cold cuts intake,\nfrequency"
)

p9 <- num_plot(
  "nuts",
  "infant_ffq_group_nuts_gday",
  "Nut DNA vs Nut Intake (g/day)",
  "Nut intake, g/day"
)

fig_associations <- (
  (p1 | p2 | p3) /
    (p4 | p5 | p6) /
    (p7 | p8 | p9)
) +
  plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 14, color = "black"),
    plot.tag.position = c(0.02, 0.98)
  )

fig_associations


ggsave(filename = "foodseq_ffq_associations.tiff",
       plot = fig_associations,
       width = 300, height = 300, units = "mm",
       dpi = 600, compression = "lzw")



