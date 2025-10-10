# Load scripts 
source("~/Desktop/LLNEXT/Analysis/scripts/new_scripts_2025/Basic_functions_XGBoost.R")

# Load packages 
library(coda4microbiome)
library(caret)
library(pROC)
library(xgboost)
library(tidyverse)


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis")
metadata<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
phenos <-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_2023_11_15.txt")
eczema<-phenos 
names (eczema)[2]<-"NEXT_ID"
eczema<-eczema[!duplicated(eczema$NEXT_ID), ]
eczema$FAMILY=NULL
eczema$infant_relations=NULL
metadata<-left_join(metadata, eczema)
row.names(metadata)<-metadata$NG_ID

# Only birth 
mother_metadata_B<-metadata[metadata$Timepoint_categorical=="B",]
mother_metadata_B<- mother_metadata_B%>% drop_na(infant_health_eczema_diagnosis_strict)
table (mother_metadata_B$infant_health_eczema_diagnosis_strict)
table (mother_metadata_B$infant_scorad_family_history_allergic_disease)


# Assign metadata to final_data
final_data <- mother_metadata_B
final_data$infant_health_eczema_diagnosis_strict <- as.factor(final_data$infant_health_eczema_diagnosis_strict)
final_data$infant_scorad_family_history_allergic_disease <- as.factor(final_data$infant_scorad_family_history_allergic_disease)


# Taxa at SGB level 
taxa <- read.delim("taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
SGB <- taxa[, grep("t__", colnames(taxa))]
colnames(SGB) <- gsub(".*s__", "", colnames(SGB))
SGB_clean <- SGB[rownames(SGB) %in% rownames(final_data), ]
SGB_clean <- SGB_clean[rownames(final_data), ]
Prev_analysis <- apply(SGB_clean, 2, function(x) sum(x != 0) / length(x))
Prev_analysis <- names(Prev_analysis[Prev_analysis > 0.3])
SGB_filtered_prev <- SGB_clean[, Prev_analysis]
SGB_filtered_prev_mat <- as.matrix(SGB_filtered_prev)
eczema <-final_data$infant_health_eczema_diagnosis_strict

# Pathways 
pathways <- read.delim("pathways/metacyc_mothers_ALR_zeroTreated_20250722.tsv")
path_clean <- pathways [rownames(pathways ) %in% rownames(final_data), ]
path_clean <- path_clean[rownames(final_data), ]
Prev_analysis <- apply(path_clean, 2, function(x) sum(x != 0) / length(x))
Prev_analysis <- names(Prev_analysis[Prev_analysis > 0.5]) # Chosen this filter to make features in species and pathways comparable 
path_filtered_prev <- path_clean[, Prev_analysis]
path_filtered_prev_mat <- as.matrix(path_filtered_prev)
eczema <-final_data$infant_health_eczema_diagnosis_strict

#Pred_probabilities <- Fig_xgboost_LOO(SGB_filtered_prev_mat, as.numeric(final_data$infant_health_eczema_diagnosis_strict) - 1)
#No ratios and species 
Pred_probabilities <- Fig_xgboost_LOO(SGB_filtered_prev_mat, as.numeric(final_data$infant_health_eczema_diagnosis_strict) - 1, ratios = F)

# Pathways
Pred_probabilities_path <- Fig_xgboost_LOO(path_filtered_prev_mat, as.numeric(final_data$infant_health_eczema_diagnosis_strict) - 1, ratios = F)

#Alpha diversity
Pred_probabilities_alpha <- Fig_xgboost_LOO(as.matrix( select(final_data, shannon)) , as.numeric(final_data$infant_health_eczema_diagnosis_strict) - 1, ratios = F, Do_logistic=T)
#Smoking
final_data %>% select(mother_health_smoked_one_whole_year_p18, infant_health_eczema_diagnosis_strict) %>% drop_na() -> input_smoke
Pred_probabilities_smoke <- Fig_xgboost_LOO( select(input_smoke, -infant_health_eczema_diagnosis_strict) %>% as.matrix()  , as.numeric(input_smoke$infant_health_eczema_diagnosis_strict) - 1, ratios = F, Do_logistic=T)
#Family history
final_data %>% select(infant_scorad_family_history_allergic_disease, infant_health_eczema_diagnosis_strict) %>% drop_na() -> input_familyhistory
Pred_probabilities_familyhistory <- Fig_xgboost_LOO( select(input_familyhistory, -infant_health_eczema_diagnosis_strict) %>% as.matrix()  , as.numeric(input_familyhistory$infant_health_eczema_diagnosis_strict) - 1, ratios = F, Do_logistic=T)
#Both
Pred_probabilities_all <- Fig_xgboost_LOO(merge(SGB_filtered_prev_mat, select(final_data, shannon), by='row.names' ) %>% select(-`Row.names`) , as.numeric(final_data$infant_health_eczema_diagnosis_strict) - 1, ratios = F, Do_logistic=F)


Y <- final_data$infant_health_eczema_diagnosis_strict == "yes"
Scores = Get_scores_model( Pred_probabilities, Y)
Scores_path=Get_scores_model( Pred_probabilities_path, Y)
Scores_alpha = Get_scores_model( Pred_probabilities_alpha, Y)
Scores_smoke = Get_scores_model( Pred_probabilities_smoke, input_smoke$infant_health_eczema_diagnosis_strict=='yes' )
Scores_familyhistory = Get_scores_model( Pred_probabilities_familyhistory, input_familyhistory$infant_health_eczema_diagnosis_strict=='yes' )
Scores_all = Get_scores_model( Pred_probabilities_all, Y)

tibble(Y_hat = Pred_probabilities_familyhistory, Y = input_familyhistory$infant_health_eczema_diagnosis_strict ) %>%
ggplot(aes(x=Y, y=Y_hat)) + geom_boxplot() + theme_bw() + geom_jitter() 


###Association
final_data %>% glm(infant_health_eczema_diagnosis_strict ~ shannon + infant_scorad_family_history_allergic_disease + mother_health_smoked_one_whole_year_p18, . , family='binomial') %>% summary()
###Include top features model species
library("shapviz")
model <- xgboost(
  data = SGB_filtered_prev_mat,
  label = as.numeric(final_data$infant_health_eczema_diagnosis_strict) - 1,
  nrounds = 100,
  objective = "binary:logistic",
  verbose = 0
)
shap_values = shapviz(model, X_pred = SGB_filtered_prev_mat )
sv_importance(shap_values) -> Plot
sv_importance(
  shap_values,
  kind = "bee",
  max_display = 10,
  number_size = 10,
)
#include
output_model<-merge(SGB_filtered_prev_mat, select(final_data, shannon, infant_health_eczema_diagnosis_strict,infant_scorad_family_history_allergic_disease, mother_health_smoked_one_whole_year_p18  ), by='row.names' ) %>% select(-`Row.names`) %>% 
  glm(infant_health_eczema_diagnosis_strict ~   infant_scorad_family_history_allergic_disease +mother_health_smoked_one_whole_year_p18 +shannon , . , family='binomial') %>% summary()
output_model

output_model<-merge(SGB_filtered_prev_mat, select(final_data, shannon, infant_health_eczema_diagnosis_strict ), by='row.names' ) %>% select(-`Row.names`) %>% 
  glm(infant_health_eczema_diagnosis_strict ~ Dialister_invisus.t__SGB5825_group  + GGB9730_SGB15291.t__SGB15291 + Coprococcus_comes.t__SGB4577_group   + shannon , . , family='binomial') %>% summary()


df <- merge(
  SGB_filtered_prev_mat,
  select(final_data, shannon, infant_health_eczema_diagnosis_strict),
  by = "row.names"
) %>%
  select(-Row.names)

ggplot(df, aes(x = infant_health_eczema_diagnosis_strict, y = Dialister_invisus.t__SGB5825_group, fill = infant_health_eczema_diagnosis_strict)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # No boxplot outliers
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) + # Show individual points
  theme_minimal() +
  labs(
    title = "Dialister invisus (SGB5825) vs Eczema Diagnosis",
    x = "Eczema Diagnosis (Strict)",
    y = "Dialister invisus (SGB5825) Group"
  ) +
  theme(legend.position = "none")


library(pROC)
library(ggplot2)
library(dplyr)

# Compute ROC objects
roc1 <- roc(controls = Pred_probabilities[Y == FALSE],
            cases    = Pred_probabilities[Y == TRUE])
roc2 <- roc(controls = Pred_probabilities_path[Y == FALSE],
            cases    = Pred_probabilities_path[Y == TRUE])
roc3 <- roc(controls = Pred_probabilities_alpha[Y == FALSE],
            cases    = Pred_probabilities_alpha[Y == TRUE])

# Create data frames from ROC objects
#df1 <- data.frame(
  #FPR = 1 - roc1$specificities,
 # TPR = roc1$sensitivities,
 # Model = sprintf("Taxanomies (AUC = %.2f)", auc(roc1))
#)

#df2 <- data.frame(
  #FPR = 1 - roc2$specificities,
  #TPR = roc2$sensitivities,
 # Model = sprintf("Pathways (AUC = %.2f)", auc(roc2))
#)

#df3 <- data.frame(
  #FPR = 1 - roc3$specificities,
  #TPR = roc3$sensitivities,
 # Model = sprintf("Alpha Diversity (AUC = %.2f)", auc(roc3))
#)

# Combine all into one dataframe
#roc_df <- bind_rows(df1, df2, df3)

#ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
 # geom_smooth(se = FALSE, method = "loess", span = 0.5) +
  #geom_abline(linetype = "dashed", color = "gray") +
  #labs(title = "Smoothed ROC Curves", x = "Sensitivity", y = "Specificity") +
 # theme_minimal() +
  #theme(legend.title = element_blank())


# ROC's 
roc1 <- roc(controls = Pred_probabilities[Y == FALSE],
            cases    = Pred_probabilities[Y == TRUE],
            smooth = F)

roc2 <- roc(controls = Pred_probabilities_path[Y == FALSE],
            cases    = Pred_probabilities_path[Y == TRUE],
            smooth = F)

roc3 <- roc(controls = Pred_probabilities_alpha[Y == FALSE],
            cases    = Pred_probabilities_alpha[Y == TRUE],
            smooth = F)

# Compute AUCs
auc1 <- auc(roc1)
auc2 <- auc(roc2)
auc3 <- auc(roc3)


col1 <- "#E63946"     # Bright coral red
col2 <- "#F4A261"     # Soft orange
col3 <- "#2A9D8F"     # Bright teal


plot(roc1, col = col1, lwd = 2, main = "", xlim = c(1, 0), ylim = c(0, 1))


lines(roc2, col = col2, lwd = 2)
lines(roc3, col = col3, lwd = 2)


legend("bottomright",
       legend = c(
         sprintf("SGBs: %.2f", auc1),
         sprintf("Pathways: %.2f", auc2),
         sprintf("Alpha diversity: %.2f", auc3)
       ),
       col = c(col1, col2, col3),
       lwd = 2,
       cex = 0.7) 
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_4")


