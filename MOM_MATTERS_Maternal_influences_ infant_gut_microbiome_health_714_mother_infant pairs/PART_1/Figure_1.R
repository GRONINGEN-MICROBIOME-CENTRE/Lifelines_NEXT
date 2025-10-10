################# Figure 1 ##################

# Load functions 
mixed_models_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
    # print (c(Bug, Prevalence))
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      Model0 = as.formula(paste( c(Bug2,  " ~ dna_conc + clean_reads_FQ_1 + BATCH_NUMBER+(1|NEXT_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ dna_conc  + clean_reads_FQ_1 + BATCH_NUMBER+",pheno2, "+ (1|NEXT_ID)"), collapse="" ))
      lmer(Model2, df_pheno, REML = F) -> resultmodel2
      M = "Mixed"
      as.data.frame(anova(resultmodel2, base_model))['resultmodel2','Pr(>Chisq)']->p_simp
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
  
}

# Load libraries 
library(tidyverse)
library(lmerTest)

# Figure 1A Cohort overview

# Mother alpha diversity + infant alpha diversity
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata")
metadata<-read.delim("LLNEXT_metadata_15_04_2024.txt")

metadata<- metadata %>%
  filter(!(Type == "mother" & Timepoint_categorical %in% c("M1", "M2")))
metadata$Timepoint_categorical=factor(metadata$Timepoint_categorical, levels = c("P12","P28","B", "W2", "M1", "M2", "M3", "M6", 'M9', "M12"))
metadata$Type=factor(metadata$Type, levels = c("mother", "infant"))
timepoint_colors <- c("#f90404", "#f78310", "#fbd123","#b5dd88","#41c0b4", "#4397bb", "#eca4c9", "#cb4563","#a42097", "#390962" )
summary (metadata$Timepoint_categorical)

# Testing significance in mothers 
mothers <-metadata[metadata$Type=="mother",]
results_mothers<-lmer(shannon ~ clean_reads_FQ_1 + dna_conc + Timepoint_categorical + (1|NEXT_ID), mothers)
summary (results_mothers)
#Fixed effects:
#Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)               3.701e+00  6.551e-02  1.401e+03  56.505  < 2e-16 ***
#clean_reads_FQ_1          1.195e-08  4.909e-09  1.333e+03   2.434   0.0151 *  
#dna_conc                 -2.644e-03  5.606e-04  1.429e+03  -4.716 2.64e-06 ***
#Timepoint_categoricalP28 -3.876e-02  2.764e-02  1.071e+03  -1.402   0.1611    
#Timepoint_categoricalB   -3.403e-02  3.229e-02  1.182e+03  -1.054   0.2921    
#Timepoint_categoricalM3  -4.929e-02  2.638e-02  1.056e+03  -1.869   0.0619 .  


infants <-metadata[metadata$Type=="infant",]
results_infant<-lmer(shannon ~ clean_reads_FQ_1 + dna_conc + Timepoint_categorical + (1|NEXT_ID), infants)
summary (results_infant)
#Fixed effects:
#Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)               1.416e+00  5.578e-02  2.901e+03  25.383  < 2e-16 ***
#clean_reads_FQ_1         -9.538e-09  3.633e-09  2.787e+03  -2.625 0.008711 ** 
#dna_conc                 -3.288e-03  1.275e-03  2.860e+03  -2.578 0.009973 ** 
#Timepoint_categoricalM1   1.525e-02  3.438e-02  2.499e+03   0.444 0.657390    
#Timepoint_categoricalM2   1.231e-01  3.391e-02  2.504e+03   3.632 0.000287 ***
#Timepoint_categoricalM3   1.744e-01  3.356e-02  2.551e+03   5.197 2.19e-07 ***
#Timepoint_categoricalM6   4.514e-01  3.789e-02  2.642e+03  11.913  < 2e-16 ***
#Timepoint_categoricalM9   8.407e-01  3.876e-02  2.654e+03  21.688  < 2e-16 ***
#Timepoint_categoricalM12  1.015e+00  3.880e-02  2.706e+03  26.170  < 2e-16 ***


p_values_mother <- data.frame(
  Timepoint_categorical = factor(c("P28", "B", "M3"), levels = levels(metadata$Timepoint_categorical)),
  p_value = c("ns", "ns", "ns"),
  Type = "mother"
)

p_values_infant <- data.frame(
  Timepoint_categorical = factor(c("M1", "M2", "M3", "M6", "M9", "M12"), levels = levels(metadata$Timepoint_categorical)),
  p_value = c("ns", "***", "***", "***", "***", "***"),
  Type = "infant"
)

p_values <- rbind(p_values_mother, p_values_infant)

p_values$Type<-factor(p_values$Type, levels = c("mother", "infant"))


p_value_positions <- metadata %>%
  group_by(Type, Timepoint_categorical) %>%
  summarise(y_position = max(shannon) + 1) %>%
  right_join(p_values, by = c("Type", "Timepoint_categorical"))


ggplot(metadata, aes(x = Timepoint_categorical, y = shannon, fill = Timepoint_categorical)) +
  scale_fill_manual(values = timepoint_colors) +
  scale_color_manual(values = timepoint_colors) +
  geom_violin(trim = FALSE, alpha = 0.2, aes(color = Timepoint_categorical), width = 1) +
  geom_boxplot(aes(color = Timepoint_categorical), alpha = 0.2, outlier.colour = NA, width = 0.4) +
  geom_point(aes(color = Timepoint_categorical), alpha = 0.4,
             position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0)) +
  facet_wrap(~ Type, scales = "free_x") +
  theme_bw() +
  labs(
    title = "",
    x = "Timepoint",
    y = "Shannon Diversity Index"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 20) 
  ) +
  geom_text(data = p_value_positions, aes(x = Timepoint_categorical, y = y_position, label = p_value), 
            vjust = 0.1, hjust = 2, size = 4, inherit.aes = T)

ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/figure_1/shannon_diversity_mothers_infants_all.pdf", 
       width = 6, height = 5)
# Mother gut microbiome over time 
rm (list = ls())
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
row.names(metadata)<-metadata$NG_ID
metadata$BATCH_NUMBER<-as.factor(metadata$BATCH_NUMBER)
metadata_infants<-metadata[metadata$Type=="infant", ] # n=2939
metadata_mothers<-metadata[metadata$Type=="mother", ] # n=1638
metadata_mothers <- subset(metadata_mothers, Timepoint_categorical != "M1" & Timepoint_categorical != "M2")
metadata_mothers$Timepoint_categorical <- factor(metadata_mothers$Timepoint_categorical, 
                                                 levels = c("P12", "P28", "B", "M3"))

summary (metadata_mothers$Timepoint_categorical) # Everywhere summarizing to add to data labels 
#P12 P28   B  M3 
#414 406 268 499 
taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
mother_taxa<-taxa[row.names(taxa)%in% rownames(metadata_mothers),] 
mother_taxa<-mother_taxa[match(row.names(metadata_mothers),row.names(mother_taxa)),]
mother_NEXT_ID<-metadata_mothers %>%
  select(NEXT_ID)
mother_taxa_all<-merge(mother_NEXT_ID,mother_taxa, by="row.names" )
row.names(mother_taxa_all)<-mother_taxa_all$Row.names
mother_taxa_all$Row.names=NULL
unique_counts <- sapply(mother_taxa_all, function(x) length(unique(mother_taxa_all$NEXT_ID[x >0.001]))) # Here I am counting the unique elements in the NEXT_ID column where the corresponding value in each column (i.e., x) is greater than the given cut-off. 
mother_taxa_all_filt <- mother_taxa_all[, unique_counts >= 0.3*length(unique(mother_taxa_all$NEXT_ID)) ] # Setting a 20% cut-off on prevalence (3rd august version) 
mother_taxa_all_filt$NEXT_ID=NULL

mother_taxa_SGB<-mother_taxa_all_filt[,grep("t__",colnames(mother_taxa_all_filt))]
my_pseudocount_normal=min(mother_taxa_SGB[mother_taxa_SGB!=0])/2# 

distance_mother=vegdist(mother_taxa_SGB, method = "aitchison", pseudocount=my_pseudocount_normal) 
mypcoa_CLR=cmdscale(distance_mother, k = 20, eig = T)
my_var_CLR=round(mypcoa_CLR$eig*100/sum(mypcoa_CLR$eig),2)[1:20]
barplot (my_var_CLR)
mypcoa_df_CLR=as.data.frame(mypcoa_CLR$points)
names(mypcoa_df_CLR) <- c('PC1','PC2','PC3','PC4','PC5', 'PC6','PC7','PC8','PC9','PC10',
                          'PC11','PC12','PC13','PC14','PC15', 'PC16','PC17','PC18','PC19','PC20')

time_on_pc <- mixed_models_taxa(metadata_mothers, 
                                "NG_ID", 
                                mypcoa_df_CLR, 
                                c("Timepoint_categorical"))

write.table(time_on_pc, "timepoint_versus_PC_aitchinson_mothers.txt", sep="\t", row.names = F)

all_mothers=merge(mypcoa_df_CLR,metadata_mothers, by="row.names")
row.names (all_mothers)<-all_mothers$Row.names
all_mothers$Timepoint_categorical <- factor(all_mothers$Timepoint_categorical, levels = c("P12", "P28", "B", "M3"), ordered=T)
summary (all_mothers$Timepoint_categorical)

centroids <- all_mothers %>%
  group_by(Timepoint_categorical) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

maternal_colors<- c("#f90404", "#f78310", "#fbd123", "#eca4c9")

pcoa_timepoint_mothers <- ggplot(all_mothers, aes(PC1, PC2, color = Timepoint_categorical, fill=Timepoint_categorical)) +
  geom_point(size = 2, alpha = 1) +
  stat_ellipse(aes(group = Timepoint_categorical, fill = Timepoint_categorical), 
               type = "norm", linetype = 2, geom = "polygon", alpha = 0, show.legend = F) +
  geom_point(data = centroids, aes(PC1, PC2, fill = Timepoint_categorical), shape = 23, size = 4, color = 'black') +
  xlab(paste("PCo1=", round(my_var_CLR[1], digits = 2), "%", sep = "")) +
  ylab(paste("PCo2=", round(my_var_CLR[2], digits = 2), "%", sep = "")) +
  scale_color_manual(name = NULL, 
                     breaks = c("P12", "P28", "B", "M3"),
                     labels = c("P12", "P28", "B", "M3"),
                     values = maternal_colors) +
  scale_fill_manual(name = NULL, 
                    breaks = c("P12", "P28", "B", "M3"),
                    labels = c("P12", "P28", "B", "M3"),
                    values = maternal_colors) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA, size = 16), 
        legend.key.size = unit(2, "lines"), 
        legend.text = element_text(size = 16), 
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

pcoa_timepoint_mothers <- ggExtra::ggMarginal(pcoa_timepoint_mothers, type = "boxplot", groupColour = TRUE,  fill=maternal_colors)

pcoa_timepoint_mothers


# Infant gut microbiome over time 
metadata_infants<-metadata[metadata$Type=="infant", ] # n=2939
metadata_infants$Timepoint_categorical=factor(metadata_infants$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))
infant_taxa<-taxa[row.names(taxa)%in% rownames(metadata_infants),] 
infant_taxa<-infant_taxa[match(row.names(metadata_infants),row.names(infant_taxa)),]
infant_NEXT_ID<-metadata_infants %>%
  select(NEXT_ID)
infant_taxa_all<-merge(infant_NEXT_ID,infant_taxa, by="row.names" )
row.names(infant_taxa_all)<-infant_taxa_all$Row.names
infant_taxa_all$Row.names=NULL
unique_counts <- sapply(infant_taxa_all, function(x) length(unique(infant_taxa_all$NEXT_ID[x >0.1]))) # Here I am counting the unique elements in the NEXT_ID column where the corresponding value in each column (i.e., x) is greater than the given cut-off. 
infant_taxa_all_filt <- infant_taxa_all[, unique_counts >= 0.1*length(unique(infant_taxa_all$NEXT_ID)) ] # Setting a 20% cut-off on prevalence (3rd august version) 
infant_taxa_all_filt$NEXT_ID=NULL

infant_taxa_SGB<-infant_taxa_all_filt[,grep("t__",colnames(infant_taxa_all_filt))]
my_pseudocount_normal=min(infant_taxa_SGB[infant_taxa_SGB!=0])/2# 


distance_infant=vegdist(infant_taxa_SGB, method = "aitchison", pseudocount=my_pseudocount_normal) 
mypcoa_CLR=cmdscale(distance_infant, k = 20, eig = T)
my_var_CLR=round(mypcoa_CLR$eig*100/sum(mypcoa_CLR$eig),2)[1:20]
barplot (my_var_CLR)
mypcoa_df_CLR=as.data.frame(mypcoa_CLR$points)
names(mypcoa_df_CLR) <- c('PC1','PC2','PC3','PC4','PC5', 'PC6','PC7','PC8','PC9','PC10',
                          'PC11','PC12','PC13','PC14','PC15', 'PC16','PC17','PC18','PC19','PC20')

time_on_pc <- mixed_models_taxa(metadata_infants, 
                                "NG_ID", 
                                mypcoa_df_CLR, 
                                c("Timepoint_categorical"))
write.table(time_on_pc, "timepoint_versus_PC_aitchinson_infants.txt", sep="\t", row.names = F)

all_infants=merge(mypcoa_df_CLR,metadata_infants, by="row.names")
row.names (all_infants)<-all_infants$Row.names

infant_colors<-c("#b5dd88","#41c0b4", "#4397bb", "#eca4c9", "#cb4563","#a42097", "#390962")

summary (all_infants$Timepoint_categorical)

centroids <- all_infants %>%
  group_by(Timepoint_categorical) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

pcoa_timepoint_infants <- ggplot(all_infants, aes(PC1, PC2, color = Timepoint_categorical, fill=Timepoint_categorical)) +
  geom_point(size = 2, alpha = 1) +
  stat_ellipse(aes(group = Timepoint_categorical, fill = Timepoint_categorical), 
               type = "norm", linetype = 2, geom = "polygon", alpha = 0, show.legend = F) +
  geom_point(data = centroids, aes(PC1, PC2, fill = Timepoint_categorical), shape = 23, size = 4, color = 'black') +
  xlab(paste("PCo1=", round(my_var_CLR[1], digits = 2), "%", sep = "")) +
  ylab(paste("PCo2=", round(my_var_CLR[2], digits = 2), "%", sep = "")) +
  scale_color_manual(name = NULL, 
                     breaks = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"),
                     labels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"),
                     values = infant_colors) +
  scale_fill_manual(name = NULL, 
                    breaks = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"),
                    labels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"),
                    values = infant_colors) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA, size = 16), 
        legend.key.size = unit(1, "lines"), 
        legend.key.width = unit(0.9, 'cm'),
        legend.text = element_text(size = 16), 
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

pcoa_timepoint_infants <- ggExtra::ggMarginal(pcoa_timepoint_infants, type = "boxplot", groupColour = TRUE,  fill=infant_colors)
pcoa_timepoint_infants


pcoa_both <-ggarrange(pcoa_timepoint_mothers, pcoa_timepoint_infants)

pcoa_both

ggsave(filename = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/figure_1/beta_timepoint_mothers_all.pdf", 
       width = 12, height = 5)

# Aitchisons distance mothers & infants (related, unrelated, twins, siblings and over time) 

# Figure 1D 
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/figure_1")
load("Figure_1_D.RData")


infant<-ggplot(dists_cons %>% filter(Type=="infant:infant"),
       aes(x=TimeComp,y=Distance,fill=Comp))+
  geom_violin(alpha=0.5,draw_quantiles = 0.5,position=position_dodge(width = 0.75))+
  theme_bw()+
  ylim(c(0,125))+
  theme(legend.position = "bottom",
        legend.text = element_text(size=4.5), 
        legend.title = element_text(size=5.5),
        legend.key.size = unit(0.5,"cm"))+
  scale_fill_manual(values = c("gray30","gray60","gray90"),name="Comparison type")+
  ylab("Aitchison Distance")+ xlab("Infant time point comparison")
#ggsave(filename = "results/beta_diversity_timepoint/infant2Infant.pdf",
#width = 6.5,height = 3)
infant

mother <-ggplot(dists_cons %>% filter(Type=="mother:mother"),aes(x=TimeComp,y=Distance,fill=Comp)) +
  geom_violin(alpha = 0.5, draw_quantiles = 0.5, position = position_dodge(width = 0.75)) +
  theme_bw() +
  ylim(c(0, 125)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 4.5), 
        legend.title = element_text(size = 5.5),
        legend.key.size = unit(0.5, "cm")) +
  scale_fill_manual(
    values = c("gray30", "gray60", "gray90"), 
    name = "Comparison type",
    labels = c("SameIndvSamePreg", "SameIndvDiffPreg", "Unrelated")
  ) +
  ylab("Aitchison Distance") + 
  xlab("Mother time point comparison")
#ggsave(filename = "results/beta_diversity_timepoint/mother2Mother.pdf",
#width = 3.5,height = 3)
mother

library(patchwork)
D <-   mother +infant
D

ggsave(filename = "Figure_1D.pdf",
width = 12,height = 5)


# Adding significance 
# All comparisons are FDR significant FDR =0.0000 expect for infant:infant_M6-M9 Same_vs_related FDR <0.43
# N

# infant:infant
#SameIndv Related Unrelated
#P12-P28        0       0         0
#P28-B          0       0         0
#B-M3           0       0         0
#B-W2           0       0         0
#W2-M1        242      25    153979
#M1-M2        414      47    231141
#M2-M3        442      47    274352
#M3-M6        299      37    191002
#M6-M9        293      23    116286
#M9-M12       312      25    137496


#mother:mother
#SameIndv Related Unrelated
#P12-P28      246      12    167826
#P28-B         23       8    108777
#B-M3         196      18    133518
#B-W2           0       0         0
#W2-M1          0       0         0
#M1-M2          0       0         0
#M2-M3          0       0         0
#M3-M6          0       0         0
#M6-M9          0       0         0
#M9-M12         0       0         0

filtered_dists_cons <- dists_cons %>%
  filter(Type == "mother:mother" & Comp == "Related" , TimeComp== "B-M3")

unique_ids <- filtered_dists_cons %>%
  select(NEXT_ID1) %>%
  distinct()
