####################### Figure 2 ########################

# Author: Trishla Sinha, Cyrus Mallon
# Last update: 1st October, 2025


library (ggplot2)
library(dplyr)
library(stringr)
library(MOFA2)
library(BiocManager)
library(tidyverse)
library(patchwork)
library(magrittr)


# Figure 2 A, B

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/NEXT_MEFISTO")
load("feature_matrix_26_01_2024.rda")
load("data_for_MOFA_26_01_2024.rda")
MOFAobject<-load("MOFAobject_untrained_26_01_2024.rda")

MOFAobject <- load_model("LLNEXT_microbiome_model_26_01_2024.hdf5", load_interpol_Z = FALSE)

res_variance <- calculate_variance_explained(MOFAobject)
start_factor <- 1
end_factor <- 10

# create a vector of factor names
factor_names <- paste0("factor_", seq(start_factor, end_factor))

# plot as bar charts
p_variance_bar_charts <-
  as_tibble(res_variance$r2_per_factor$single_group) %>% 
  mutate(factor = all_of(factor_names)) %>%
  mutate(factor = factor(factor, levels = factor_names )) %>%
  ggplot() +
  geom_bar(aes(x = factor, y = microbiome), stat = "identity") +
  ylab("variance explained (%)") +
  xlab("factor") +
  ggtitle("Variance Explained by Factor") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 


# examine variance total
p <- plot_variance_explained(MOFAobject, plot_total = T)

# variance decomposition by factor
p_variance_by_factor <-
  p[[1]] + theme(axis.text.x = element_text(angle = 90))

# total variance explained per view 
p_variance_total <-
  p[[2]] + theme(axis.text.x = element_text(angle = 90))

# plot projections
p_dimensions_reduced <-
  plot_factors(MOFAobject, color_by = "time")



mat <- MOFAobject@covariates$single_group

# convert to matrix
tib <- as_tibble(t(mat))

# set names of matrix
tib$sample <- attr(mat, "dimnames")[[2]]

# get factor scores
z_scores <- as_tibble(MOFAobject@expectations$Z$single_group)
z_scores$sample <- rownames(MOFAobject@expectations$Z$single_group)

# Phenotypes that were significant after post-hoc phenotype analysis
phenotype_data <- MOFAobject@samples_metadata[,c("NG_ID_short", "birth_deliverybirthcard_mode_binary",
                                                 "birth_deliverybirthcard_place_delivery_simple",
                                                 "infant_birthcard_feeding_mode_after_birth",
                                                 "infant_birthcard_feeding_mode_after_birth",
                                                 "mother_health_smoked_one_whole_year_p18",
                                                 "mother_birthcard_parity",
                                                 "infant_misc_sex")]

data <-
  tib %>%
  left_join(., z_scores, by = "sample") %>%
  mutate(NG_ID_short = str_sub(sample, 1, 8)) %>%
  left_join(., phenotype_data, by = c("NG_ID_short")) %>%
  mutate(across(where(is.character), ~na_if(., "<NA>"))) %>%
  mutate(birth_deliverybirthcard_mode_binary = as.factor(birth_deliverybirthcard_mode_binary))


delivery_colors <- c( "#ac2120", "#65a03f")

fac1_trajectory_sp <-
  data %>% 
  filter(birth_deliverybirthcard_mode_binary != "nan") %>%
  mutate(
    time = factor(time, levels = c(0.5, 1, 2, 3, 6, 9, 12)) 
  ) %>%
  ggplot(aes(x = time, y = Factor1, group = interaction(time, birth_deliverybirthcard_mode_binary), color = birth_deliverybirthcard_mode_binary, fill = birth_deliverybirthcard_mode_binary)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 0.3) +
  geom_boxplot(alpha = 0.4, position = position_dodge(width = 0.7), width = 0.6,outlier.shape = NA) +
  scale_color_manual(values = delivery_colors) +
  scale_fill_manual(values = delivery_colors) +
  stat_summary(
    fun = mean,   
    geom = "line",  
    aes(group = birth_deliverybirthcard_mode_binary),  # draw lines for each group
    size = 1.5,
    alpha = 0.6
  ) +
  scale_x_discrete(
    labels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12")  
  ) +
  theme_classic()+labs(x="Timepoint Categorical", y = "Factor 1 score", fill="Mode of Delivery", color="Mode of Delivery")


fac1_trajectory_sp

# fac 3 ----
feeding_colors <- c("#1B9E77", "#D95F02","#7570B3")

fac3_trajectory_sp <-
  data %>%
  filter(infant_birthcard_feeding_mode_after_birth != "nan") %>%
  mutate(
    infant_birthcard_feeding_mode_after_birth = factor(infant_birthcard_feeding_mode_after_birth, levels = c("BF", "MF", "FF")),
    time = factor(time, levels = c(0.5, 1, 2, 3, 6, 9, 12)) 
  ) %>%
  mutate(infant_birthcard_feeding_mode_after_birth = factor(infant_birthcard_feeding_mode_after_birth, levels = c("BF", "MF", "FF") )) %>%
  ggplot(aes(x = time, y = Factor3, group = interaction(time, infant_birthcard_feeding_mode_after_birth), color = infant_birthcard_feeding_mode_after_birth, fill  = infant_birthcard_feeding_mode_after_birth)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 0.3) +
  geom_boxplot(alpha = 0.4, position = position_dodge(width = 0.7), width = 0.6, outlier.shape = NA) +
  scale_color_manual(values = feeding_colors ) +
  scale_fill_manual(values = feeding_colors ) +
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = infant_birthcard_feeding_mode_after_birth),
    size = 1.5,
    alpha = 0.6
  ) +
  scale_x_discrete(
    labels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12")  
  ) +
  labs(fill = "feeding mode",
       y = "Factor 3 Scores",
       x = "time (month)") +
  theme_classic()+labs(x="Timepoint Categorical", y = "Factor 3 score", fill="Feeding Mode", color="Feeding Mode")


fac3_trajectory_sp

feature_matrix <-
  feature_matrix %>% rename(feature_description = taxon, feature = taxonID)

features_metadata(MOFAobject) <- features_metadata(MOFAobject) %>%
  left_join(.,feature_matrix, by = c("feature")) %>%
  separate(col = feature_description, sep = "\\.",
           into = c("kingdom", "phylum", "class",
                    "order", "family", "genus", "species","taxon"),
           remove = FALSE)

# Fac 1 weights----
fac1_colors_weights  <- c(CS = "#ac2120", VG = "#65a03f")


df_weights1 <- get_weights(MOFAobject, as.data.frame = T, factors = 1) %>%
  left_join(., features_metadata(MOFAobject), by = c("feature", "view"))
df_weights1$species_new<-paste0(df_weights1$species, " ", df_weights1$taxon)
df_weights1 %<>%
  filter(!is.na(taxon)) %>%
  mutate(species_new = gsub("s__","", species_new)) %>%
  mutate(species_new = gsub("t__","", species_new)) %>%
mutate(species_new = gsub("_"," ", species_new)) 

# filter to top 10 positive and negative ones
df_top <- df_weights1 %>% group_by(species_new) %>% select(-feature_description) %>% 
  summarize(mean_weight = mean(value), n_spec = n()) %>% 
  ungroup() %>%
  # cesarean is linked to higher latent values
  mutate(type = ifelse(mean_weight < 0, "Negative Weights", "Positive Weights")) %>%
  group_by(type) %>%
  slice_max(abs(mean_weight), n= 5) %>% ungroup() %>% arrange(mean_weight) %>%
  mutate(species_new = factor(species_new, levels = .$species_new))


fac1_weights <-
  df_top %>% 
  mutate(birthmode = case_when(
    type == "Negative Weights" ~ "CS",
    type == "Positive Weights" ~ "VG",
    TRUE ~ "undetermined"  
  )) %>%
  ggplot(aes(x= species_new, y = mean_weight, fill = birthmode)) + 
  geom_bar(stat="identity") +
  coord_flip() + 
  xlab("") + 
  ggtitle("Influence of Species Across Factor") +
  geom_point(data = filter(df_weights1, species_new %in% df_top$species_new),
             aes(x = species_new, y = value), inherit.aes = FALSE, alpha = 0.3)  +
  scale_fill_manual(values = fac1_colors_weights) +
  ylab("Weight (Factor 1)") +
  theme_classic() +
  theme(
    text = element_text(size = 10),
    plot.title = element_blank()
  ) +
  guides(fill = guide_legend(title = "Mode of Delivery"))

fac1_weights



# fac3 weights ----
df_weights1 <- get_weights(MOFAobject, as.data.frame = T, factors = 3) %>%
  left_join(., features_metadata(MOFAobject), by = c("feature", "view"))

df_weights1$species_new<-paste0(df_weights1$species, " ", df_weights1$taxon)
df_weights1 %<>%
  filter(!is.na(taxon)) %>%
  mutate(species_new = gsub("s__","", species_new)) %>%
  mutate(species_new = gsub("t__","", species_new)) %>%
  mutate(species_new = gsub("_"," ", species_new)) 

# summarize by mean weights across all species in the genus and 
# filter to top 10 positive and negative ones
df_top <- df_weights1 %>% group_by(species_new) %>% select(-feature_description) %>% 
  summarize(mean_weight = mean(value), n_spec = n()) %>% 
  ungroup() %>%
  mutate(type = ifelse(mean_weight < 0, "Negative Weights", "Positive Weights")) %>%
  group_by(type) %>%
  slice_max(abs(mean_weight), n= 5) %>% ungroup() %>% arrange(mean_weight) %>%
  mutate(species_new = factor(species_new, levels = .$species_new))


fac3_colors_weights  <- c(BF = "#1B9E77", FF = "#7570B3")

fac3_weights <-
  df_top %>% 
  mutate(feedingmode = case_when(
    type == "Negative Weights" ~ "FF",
    type == "Positive Weights" ~ "BF",
    TRUE ~ "undetermined"  
  )) %>%
  ggplot(aes(x= species_new, y = mean_weight, fill = feedingmode)) + 
  geom_bar(stat="identity") +
  coord_flip() + 
  #theme_bw()  + 
  #theme(legend.position = "top") + 
  xlab("") + 
  #guides(fill = guide_legend(title="")) +
  ggtitle("Influence of Species Across Factor") +
  geom_point(data = filter(df_weights1, species_new %in% df_top$species_new),
             aes(x = species_new, y = value), inherit.aes = FALSE, alpha = 0.3)  +
  scale_fill_manual(values = fac3_colors_weights) +
  ylab("Weight (Factor 3)") +
  theme_classic() +
  theme(#plot.title = element_text(hjust = 0.5),
    text = element_text(size = 10),
    plot.title = element_blank())+
  guides(fill = guide_legend(title = "Feeding Mode"))

fac3_weights

# Figure 2a, b, c, d 

a<- (fac1_trajectory_sp) 

a
ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/figure_2/Figure_2a_mfisto_birth_mode.pdf", 
       width = 6, height = 5)

b <- (fac1_weights) 

b
ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/figure_2/Figure_2b_mfisto_birth_mode_weights.pdf", 
       width = 6, height = 5)

c <- fac3_trajectory_sp
c
ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/figure_2/Figure_2c_mfisto_feeding_mode.pdf", 
       width = 6, height = 5)
 
d <- fac3_weights
d
ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/figure_2/Figure_2d_mfisto_feeding_mode_weights.pdf", 
       width = 6, height = 5)
  

# Figure 2 E
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/alpha_diversity")
alpha_dynamic<-read.delim("Alpha_diversity_dynamic_phenotypes_results_all_23_09_2025.txt")
alpha_dynamic <- alpha_dynamic %>% filter(!is.na(FDR))
sig_alpha_dynamic<-alpha_dynamic[alpha_dynamic$FDR<0.05,]
sig_alpha_dynamic$converged=NULL
sig_alpha_dynamic$t.value <- sig_alpha_dynamic$Estimate_cor_base/ sig_alpha_dynamic$SE_cor_base
names (sig_alpha_dynamic)[1]<-"outcome"
names (sig_alpha_dynamic)[2]<-"phenotype"
names (sig_alpha_dynamic)[3]<-"levels"
names (sig_alpha_dynamic)[4]<-"Estimate"

alpha_cross<-read.delim("Alpha_diversity_cross_phenotypes_results_23_09.txt")
alpha_cross <- alpha_cross %>% filter(!is.na(FDR))
sig_alpha_cross<-alpha_cross[alpha_cross$FDR<0.05,]
sig_alpha_cross$model_sucess_cor_base=NULL
sig_alpha_cross$Covar.type_cor_base=NULL
sig_alpha_cross$t.value<-sig_alpha_cross$t.value_cor_base
names (sig_alpha_cross)[1]<-"outcome"
names (sig_alpha_cross)[2]<-"phenotype"

merged_alpha_sig <- bind_rows(sig_alpha_dynamic, sig_alpha_cross)
merged_alpha_sig$variable <- paste0(merged_alpha_sig$phenotype,"_", merged_alpha_sig$levels)
merged_alpha_sig$variable <- sub("_$", "", merged_alpha_sig$variable)

# Removing those phenotypes with an insignificant p value after correction for mode of delivery & feeding mode
# No factors here become insignificant appart from mode of delivery & feeding mode (which we want to retain in the plot)
merged_alpha_sig$p_cor_delivery_breastfed

merged_alpha <- ggplot(merged_alpha_sig, aes(x = reorder(variable, t.value), y = t.value, fill = t.value > 0)) + 
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
merged_alpha

ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/figure_2/shannon_diversity_infants_all.pdf", 
       width = 6, height = 5)

# Figure 2e
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/adonis_per_timepoint")
ResultsAdonis<-read.delim("combined_results_timepoint_adonis.txt")
ReSultsTCAM<-read.delim("combined_results_TCAM_adonis.txt")
ReSultsTCAM$Timepoint <-"Overall"
merged_adonis <- bind_rows(ResultsAdonis, ReSultsTCAM)
significant<-merged_adonis[merged_adonis$p.value_cor_delivery_breastfed<0.001,]
significant$Timepoint <- factor(significant$Timepoint, levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12", "Overall"))
significant <- significant %>%
  mutate(R2 = as.numeric(as.character(R2_cor_delivery_breastfed)),  
         R2 = round(R2, 3))  

significant2 <- significant %>% filter(
  !grepl("vitaminK|after_birth|mode_simple",Phenotype))

significant2$Timepoint <- factor(as.character(significant2$Timepoint),
                                 levels = levels(significant2$Timepoint),
                                 ordered = T)
segments <- significant2 %>% group_by(Phenotype) %>% 
  summarise(Timepoint2=min(Timepoint),
            xend=max(Timepoint)) %>% rename(Timepoint=Timepoint2)


segments <- rbind(segments,data.frame(Phenotype="birth_deliverybirthcard_mode_binary",
                                      Timepoint="M6",xend="M6"))
segments <- rbind(segments,data.frame(Phenotype="birth_deliverybirthcard_mode_binary",
                                      Timepoint="M9",xend="M9"))


adonis_infant_per_timepoint<-ggplot()+
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
  scale_fill_manual(values = c( "#b5dd88","#41c0b4", "#4397bb", "#eca4c9","#a42097", "#390962", "black"))
adonis_infant_per_timepoint

ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/main_figures/Figure_2_new/infant_adonis.pdf", 
       width = 8, height = 4)
