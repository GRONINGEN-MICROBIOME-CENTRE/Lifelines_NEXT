##################### Figure 4 ############################
# Author: Trishla Sinha,
# Last update: 5th October, 2025


library(wesanderson)

### Figure 2A 
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/alpha_diversity/old/mothers/")
alpha_dynamic<-read.delim("alpha_diversity_GAM_mothers_dynamic_phenotypes.txt")
alpha_cross<-read.delim("alpha_diversity_mmrm_trait_cross_sectional_mothers.txt")

sig_alpha_dynamic<-alpha_dynamic[alpha_dynamic$FDR<0.05,]
sig_alpha_dynamic$converged=NULL
sig_alpha_dynamic$t.value <- sig_alpha_dynamic$beta/ sig_alpha_dynamic$SE
names (sig_alpha_dynamic)[1]<-"outcome"
names (sig_alpha_dynamic)[2]<-"phenotype"
names (sig_alpha_dynamic)[3]<-"levels"
names (sig_alpha_dynamic)[4]<-"Estimate"
sig_alpha_cross<-alpha_cross[alpha_cross$FDR<0.05,]
sig_alpha_cross$type=NULL
sig_alpha_cross$Covar.type=NULL
names (sig_alpha_cross)[1]<-"outcome"
names (sig_alpha_cross)[2]<-"phenotype"

merged_alpha_sig <- bind_rows(sig_alpha_dynamic, sig_alpha_cross)
merged_alpha_sig$variable <- paste0(merged_alpha_sig$phenotype,"_", merged_alpha_sig$levels)
merged_alpha_sig$variable <- sub("_$", "", merged_alpha_sig$variable)

merged_alpha_sig <- merged_alpha_sig %>% filter(
  !grepl("complex|mode_simple|eczema",phenotype))

ggplot(merged_alpha_sig, aes(x = reorder(variable, t.value), y = t.value, fill = t.value > 0)) +
  geom_bar(stat = "identity", color = "black") +
  xlab("") +
  ylab("T-Value") +
  ggtitle("") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.position = "top",  
    legend.justification = c(1.2, 0.1)   
  ) +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("Decreased", "Increased")) +
  coord_flip() +
  guides(fill = guide_legend(title = ""))

ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_4/figure_4_a_maternal_alpha_diversity.pdf", 
       width = 6, height = 5)

# Figure 2B
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint/old/")
ResultsAdonis<-read.delim("result_adonis_per_timepoint_mothers_10000_perm.txt")
ReSultsTCAM<-read.delim("result_adonis_tcam_mothers_10000_perm.txt")
ReSultsTCAM$Timepoint <-"Overall"
names (ReSultsTCAM) [4]<-"p.value"
merged_adonis <- bind_rows(ResultsAdonis, ReSultsTCAM)
significant<-merged_adonis [merged_adonis $p.value<0.005,]
significant$Timepoint <- factor(significant$Timepoint, 
                                                 levels = c("P12", "P28", "B", "M3", "Overall"))

significant <- significant %>%
  mutate(R2 = as.numeric(as.character(R2)),  
         R2 = round(R2, 3))  

significant2 <- significant %>% filter(
  !grepl("complex|vitaminK|mode_simple",Phenotype))
significant2<-significant2[significant2$Df<3,]
significant2$Timepoint <- factor(as.character(significant2$Timepoint),
                                 levels = levels(significant2$Timepoint),
                                 ordered = T)
segments <- significant2 %>% group_by(Phenotype) %>% 
  summarise(Timepoint2=min(Timepoint),
            xend=max(Timepoint)) %>% rename(Timepoint=Timepoint2)


adonis_mother_per_timepoint<-ggplot()+
  geom_segment(data=segments, aes(x=Timepoint, y=Phenotype, 
                                  xend=xend,yend=Phenotype),color="black",linewidth=0.2)+
  geom_point(data = significant2,
             aes(x = Timepoint, y=Phenotype,
                 size = R2, fill = Timepoint),shape=21,color="black") +
  # facet_wrap(~ Timepoint) + 
  # coord_flip() +  
  xlab("Timepoint") +  
  ylab("") +  
  theme_bw() +
  theme(axis.text.y=element_text(size=8),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())+
  scale_size_continuous(range = c(1,4))+
  scale_fill_manual(values = c( "#f90404", "#f78310", "#fbd123", "#ECA4C9", "black"))

adonis_mother_per_timepoint

ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/figure_4/mother_adonis.pdf", 
       width = 8, height = 8)

# Figure 4C 

setwd("~/Desktop/LLNEXT/Analysis/results/alpha_diversity/old/mothers/")
results_mother_alpha_diversity<-read.delim("alpha_diversity_mmrm_trait_cross_sectional_mothers.txt")

# Load metadata and phenotypes
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
cross_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_2023_11_15.txt")
cross_selection <-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_selection_AFTER_correlation_&_correction_10_07_2024.txt")

# Merging files and selecting mother relevant data
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_mothers<-metadata[metadata$Type=="mother", ]
metadata_mothers <- subset(metadata_mothers, Timepoint_categorical != "M1" & Timepoint_categorical != "M2")
metadata_mothers$Timepoint_categorical <- factor(metadata_mothers$Timepoint_categorical, 
                                levels = c("P12", "P28", "B", "M3"))
names (metadata_mothers)[2]<-"next_id_mother"

# Selecting data relevant for mother 
mother_cross_selection<-cross_selection[cross_selection$mother_microbiome_selection==1,]
column_names <- mother_cross_selection[[1]]
mother_cross_phenotypes <- cross_phenotypes %>% select(all_of(column_names))
mother_cross_phenotypes <- subset(mother_cross_phenotypes, twin_pair != "pair_2")
mother_cross_phenotypes$FAMILY=NULL
mother_cross_phenotypes$infant_relations=NULL

mother_metadata_cross_phenotypes<-left_join(metadata_mothers, mother_cross_phenotypes, by="next_id_mother")
row.names(mother_metadata_cross_phenotypes)<-mother_metadata_cross_phenotypes$NG_ID
mother_metadata_cross_phenotypes_1<-mother_metadata_cross_phenotypes
# Remove phenotypes with too many NA's 
mother_metadata_cross_phenotypes_2 <-mother_metadata_cross_phenotypes_1
names (mother_metadata_cross_phenotypes_2)[2]<-"NEXT_ID"
mother_metadata_cross_phenotypes_2$NEXT_ID=as.factor(mother_metadata_cross_phenotypes_2$NEXT_ID)

# Select outcome 
shannon<-mother_metadata_cross_phenotypes_2 %>% select(shannon)
# Select and normalize phenotypes 
metadata_columns <- c("NG_ID", "NEXT_ID","Modified_NEXT_ID_without_preg_number", "days_from_first_collection","FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                      "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                      "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                      "exact_age_days_at_collection", "exact_age_months_at_collection",
                      "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                      "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                      "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                      "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint")


phenos <- mother_metadata_cross_phenotypes_2[, !(colnames(mother_metadata_cross_phenotypes_2) %in% metadata_columns)] 
phenos$next_id_infant=NULL

eczema<-mother_metadata_cross_phenotypes_2 %>%
  filter(!is.na(infant_health_eczema_diagnosis_strict))

ggplot(eczema, aes(Timepoint_categorical, y=shannon, fill=infant_health_eczema_diagnosis_strict, color=infant_health_eczema_diagnosis_strict)) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 2))+
  scale_color_manual(values = wes_palette("GrandBudapest1", n = 2))+
  #scale_fill_manual(values=c("#e60049", "#0bb4ff", "#50e991", "#f46a9b", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#b30000", "#7c1158", "#4421af"))+
  #scale_color_manual(values=c("#e60049", "#0bb4ff", "#50e991", "#f46a9b", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#b30000", "#7c1158", "#4421af"))+
  geom_boxplot(alpha=0.4, outlier.colour = NA)+
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0))+
  # geom_jitter() +
  ggtitle("")+
  theme_bw()+labs(x="", y = "Maternal Shannon Diversity", fill="Infant eczema", color="Infant eczema")+
  theme(
    plot.title = element_text(color="black", size=16, face="bold"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=16, angle = 60, hjust = 1))
#strip.text.x = element_text(size = 10))

# FDR 0.02 
ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_4/eczema_diversity.pdf", 
       width = 6, height = 5)

delivery_place_simple<-mother_metadata_cross_phenotypes_2 %>%
  filter(!is.na(birth_deliverybirthcard_place_delivery_simple))


ggplot(delivery_place_simple, aes(Timepoint_categorical, y=shannon, fill=birth_deliverybirthcard_place_delivery_simple, color=birth_deliverybirthcard_place_delivery_simple)) +
  #scale_fill_manual(values = wes_palette("GrandBudapest1", n = 2))+
  #scale_color_manual(values = wes_palette("GrandBudapest1", n = 2))+
  scale_fill_manual(values=c("#50e991", "#66D4FF"))+
  scale_color_manual(values=c("#50e991", "#66D4FF"))+
  geom_boxplot(alpha=0.2, outlier.colour = NA)+
  geom_point(alpha=0.2,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0))+
  ggtitle("")+
  theme_bw()+labs(x="", y = "Maternal Shannon Diversity", fill="Place of delivery", color="Place of delivery")+
  theme(
    plot.title = element_text(color="black", size=16, face="bold"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    axis.text.y = element_text(face="bold", size=16),
    axis.text.x = element_text(size=16, angle = 60, hjust = 1))

# FDR 0.006
ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_4/mother_shannon_place_delivery.pdf", 
       width = 6, height = 5)



################# Figure 4 D and E 

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
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) + 
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



pdf('Figure_auc_prediction_eczema.pdf', width = 6, height = 6)
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_4")

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

dev.off()
