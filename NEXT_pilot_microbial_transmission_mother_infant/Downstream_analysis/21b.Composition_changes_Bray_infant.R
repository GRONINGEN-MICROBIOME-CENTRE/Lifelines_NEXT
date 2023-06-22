setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the dynamics of bacteriome and virome 
# communitues within infant samples at the 
# species and vOTUs levels USING BRAY-CURTIS DISTANCE
#############################################################

##############################
# Functions
##############################
mixed_models_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list, consider_time) {
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
      
      if (consider_time=='time_as_covariate') {
        Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + infant_mode_delivery + Age_months + (1|Individual_ID)"), collapse="" )) 
      } else { # else is mainly for associating entities with time alone
        Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + (1|Individual_ID)"), collapse="" )) 
      }
      
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      
      if (consider_time=='time_as_covariate') {
        Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + infant_mode_delivery + Age_months + ",pheno2, "+ (1|Individual_ID)"), collapse="" ))
      } else { # else is mainly for associating entities with time alone
        Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + ",pheno2, "+ (1|Individual_ID)"), collapse="" ))
      }
      
      lmer(Model2, df_pheno, REML = F) -> resultmodel2
      M = "Mixed"
      as.data.frame(anova(resultmodel2, base_model))['resultmodel2','Pr(>Chisq)']->p_simp
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p <- p[! duplicated(paste0(p$Pheno, p$Bug)),]
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
  
}
##############################
# Loading libraries
##############################
library(vegan)

library(ggplot2)
library(ggExtra)
library(ggrepel)
library(ggforce)


library(dplyr)
library(tibble)
library(lme4)
library(RLRsim)
library(lmerTest)
##############################
# Input data
##############################
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)

MGS_metadata$Type <- factor(MGS_metadata$Type, levels=c('Infant', 'Mother'), ordered=T)
MGS_metadata[MGS_metadata$Type=='Mother',]$Timepoint <- 'Mother'
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels = c("M1", "M2", "M3",
                                                                    "M6", "M9", "M12", "Mother"), ordered = T)

row.names(MGS_metadata) <- MGS_metadata$Short_sample_ID_bact

microbiome <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep = '\t', header=T)
VLP_metadata$Type <- factor(VLP_metadata$Type, levels=c('Infant', 'Mother'), ordered=T)
VLP_metadata[VLP_metadata$Type=='Mother',]$Timepoint <- 'Mother'
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels = c("M1", "M2", "M3",
                                                                    "M6", "M12", "Mother"), ordered = T)

row.names(VLP_metadata) <- VLP_metadata$Short_sample_ID

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)
##############################
# ANALYSIS
##############################
#### NMDS for bacteriome
# preapring phenotypes
for_bacplot <- MGS_metadata
for_bacplot <- for_bacplot[colnames(microbiome),]
#CHOOSING TECHNICAL PHENOTYPES, SAMPLE TYPE & BACTERIAL ALPHA DIVERSITY
for_bacplot <- for_bacplot[,c("Clean_reads", "DNA_CONC", "bacterial_alpha_diversity", "Age_days", "Timepoint",
                              "metaphlan_unknown_perc")]


ord <- metaMDS(t(microbiome), distance = "bray", k=2)
en = envfit(ord, for_bacplot, permutations = 999, na.rm = TRUE)
en$factors 
en$vectors

centroids <- as.data.frame(scores(en, "factors"))
centroids$Timepoint <- c(gsub('Timepoint', '', row.names(centroids)))
centroids$Timepoint <- factor(centroids$Timepoint, levels = c("M1", "M2", "M3",
                                                              "M6", "M9", "M12", "Mother"), ordered = T)

spp.scrs <- as.data.frame(scores(en, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Type <- 'Infant'
spp.scrs$Species <- c('N clean reads', 'DNA concentration', 'Alpha diversity (bacteria)', 'Infant age', '% unknown in Metaphlan4')

### for plotting:
data.scores = as.data.frame(scores(ord, "sites"))
data.scores <- merge(data.scores, MGS_metadata, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL

pdf('./04.PLOTS/Bacterial_spp_Bray_NMDS_Timepoint_color.pdf', width=10/2.54, height=10/2.54)
ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Timepoint)) + 
  geom_point(size = 2, alpha=0.8) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = Timepoint, fill=Timepoint), linetype = 2) +
  geom_point(data=centroids, aes(fill=Timepoint),shape=23, size=4, color='black', ) + 
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "#444444") +
  geom_label_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', size = 2, alpha=0.7) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 8, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom") 
dev.off()

##### SAVE FOR PATCHING PANELS ####
write.table(data.scores, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1C_NMDS_scores_BacSp.txt", sep='\t', quote = F)
write.table(centroids, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1C_NMDS_centroids_BacSp.txt", sep='\t', quote = F)
write.table(spp.scrs, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1C_NMDS_vectors_BacSp.txt", sep='\t', quote = F)
##### SAVE FOR PATCHING PANELS ####

#### infants phenotypes

# allow only 1 dimension
ord2 <- metaMDS(t(microbiome), distance = "bray", k=1)

colnames(MGS_metadata)

infants_phenos_MGS <- mixed_models_taxa(MGS_metadata[MGS_metadata$Type=='Infant',], 
                                    "Short_sample_ID_bact", 
                                    as.data.frame(scores(ord2, "sites")), 
                                    c("infant_type_pregnancy", "infant_sex", "infant_gestational_age",
                                    "infant_birthweight", "infant_place_delivery", "infant_mode_delivery",
                                    "mother_age_years", "infant_feeding_mode_imputed_W2", "infant_feeding_mode_imputed_M1", 
                                    "infant_ever_never_breastfed", "Age_days", "perc_reads_aligned", "viral_richness_MGS", 
                                    "viral_alpha_diversity_MGS", "temperate_richness_MGS", "temperate_RA_MGS",
                                    "bacterial_alpha_diversity", "infant_ffq_feeding_mode_simple", "infant_ffq_feeding_mode_complex"), "dont_consider_time")

infants_phenos_MGS <- infants_phenos_MGS[infants_phenos_MGS$FDR<0.05,]

data.scores2 = as.data.frame(scores(ord2, "sites"))
data.scores2 <- merge(data.scores2, MGS_metadata[MGS_metadata$Type=='Infant',], by='row.names')
row.names(data.scores2) <- data.scores2$Row.names
data.scores2$Row.names <- NULL

pdf('./04.PLOTS/Bacterial_spp_Bray_NMDS1_Timepoint_Feeding_mode.pdf', width=10/2.54, height=10/2.54)
ggplot(data = data.scores2[!is.na(data.scores2$infant_ffq_feeding_mode_complex),], aes(x = Timepoint, y = NMDS1, fill=infant_ffq_feeding_mode_complex)) + 
  geom_boxplot(outlier.shape = NA,alpha=0.5) +
  geom_sina(aes(color=infant_ffq_feeding_mode_complex), size=0.6,alpha=0.7) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 8, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom") +
        labs(fill='Feeding mode', color='Feeding mode')
dev.off()

pdf('./04.PLOTS/Bacterial_spp_Bray_NMDS1_Timepoint_viral_richness_MGS.pdf', width=10/2.54, height=10/2.54)
ggplot(data = data.scores2, aes(x = viral_richness_MGS, y = NMDS1)) + 
  geom_smooth(aes(viral_richness_MGS, y = NMDS1)) +
  geom_point(aes(color=Timepoint), alpha=0.8) +
  theme_bw()+ 
  annotate("text", x = 3000, y = 1.5, label = "FDR=3.7e-27")+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 8, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom") +
  labs(x='vOTUs richness in TMs')
dev.off()

pdf('./04.PLOTS/Bacterial_spp_Bray_NMDS1_Timepoint_temperate_richness_MGS.pdf', width=10/2.54, height=10/2.54)
ggplot(data = data.scores2, aes(x = temperate_richness_MGS, y = NMDS1)) + 
  geom_smooth(aes(temperate_richness_MGS, y = NMDS1)) +
  geom_point(aes(color=Timepoint), alpha=0.8) +
  theme_bw()+ 
  annotate("text", x = 250, y = 1.5, label = "FDR=1.2e-25")+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 8, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom") +
  labs(x='Temperate phages richness in TMs')
dev.off()

### infant NMDS against bacterial species
microbiome_infant <- microbiome[,colnames(microbiome) %in% MGS_metadata[MGS_metadata$Type=='Infant',]$Short_sample_ID_bact]
microbiome_infant <- microbiome_infant[rowSums(microbiome_infant)!=0,]
microbiome_infant <- as.data.frame(t(microbiome_infant))
microbiome_infant <- merge(microbiome_infant, MGS_metadata, by='row.names')
row.names(microbiome_infant) <- microbiome_infant$Row.names
microbiome_infant$Row.names <- NULL

infant_species <- mixed_models_taxa(microbiome_infant, 
                                    "Short_sample_ID_bact", 
                                    as.data.frame(scores(ord2, "sites")), 
                                    colnames(microbiome_infant)[grep('k_',colnames(microbiome_infant))], 'dont_consider_time')
infant_species <- infant_species[infant_species$FDR<0.05,]



#### NMDS for virome
# preapring phenotypes
for_virplot <- VLP_metadata
for_virplot <- for_virplot[colnames(RPKM_counts_VLP),]
#CHOOSING TECHNICAL PHENOTYPES, SAMPLE TYPE & BACTERIAL ALPHA DIVERSITY
for_virplot <- for_virplot[,c("Clean_reads", "DNA_CONC", "viral_alpha_diversity", "Age_days", "Timepoint",
                              "bacterial_contamination_perc_reads")]


ord <- metaMDS(t(RPKM_counts_VLP), distance = "bray", k=2)
en = envfit(ord, for_virplot, permutations = 999, na.rm = TRUE)
en$factors 
en$vectors

centroids <- as.data.frame(scores(en, "factors"))
centroids$Timepoint <- c(gsub('Timepoint', '', row.names(centroids)))
centroids$Timepoint <- factor(centroids$Timepoint, levels = c("M1", "M2", "M3",
                                                              "M6", "M12", "Mother"), ordered = T)

spp.scrs <- as.data.frame(scores(en, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Type <- 'Infant'
spp.scrs <- spp.scrs[!(spp.scrs$Species %in% c('Clean_reads')),]
spp.scrs$Species <- c('DNA concentration','Alpha diversity (vOTUs)', 'Infant age', '% bacterial contamination')

### for plotting:
data.scores = as.data.frame(scores(ord, "sites"))
data.scores <- merge(data.scores, VLP_metadata, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL

pdf('./04.PLOTS/vOTUs_Bray_NMDS_Timepoint_color.pdf', width=10/2.54, height=10/2.54)
ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Timepoint)) + 
  geom_point(size = 2, alpha=0.8) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = Timepoint, fill=Timepoint), linetype = 2) +
  geom_point(data=centroids, aes(fill=Timepoint),shape=23, size=4, color='black', ) + 
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "#444444") +
  geom_label_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', size = 2, alpha=0.7) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 8, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom") 
dev.off()

##### SAVE FOR PATCHING PANELS ####
write.table(data.scores, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1B_NMDS_scores_vOTUs.txt", sep='\t', quote = F)
write.table(centroids, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1B_NMDS_centroids_vOTUs.txt", sep='\t', quote = F)
write.table(spp.scrs, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1B_NMDS_vectors_vOTUs.txt", sep='\t', quote = F)
##### SAVE FOR PATCHING PANELS ####


#### infants phenotypes

# allow only 1 dimension
ord2 <- metaMDS(t(RPKM_counts_VLP), distance = "bray", k=1)

colnames(VLP_metadata)

infants_phenos_VLP <- mixed_models_taxa(VLP_metadata[VLP_metadata$Type=='Infant',], 
                                    "Short_sample_ID", 
                                    as.data.frame(scores(ord2, "sites")), 
                                    c("infant_place_delivery", "Age_days",
                                      "bacterial_alpha_diversity", "infant_ffq_feeding_mode_complex" ), "dont_consider_time")

infants_phenos_VLP <- infants_phenos_VLP[infants_phenos_VLP$FDR<0.05,]




#ordiplot3d(ord, col = for_bacplot$Type, ax.col= "black", pch = 18)
###### OUTPUT #####
write.table(infants_phenos_MGS, '03a.RESULTS/BacSp_Bray_NMDS_phenos_mixed_models_FDR_significant.txt', sep='\t', quote=F, row.names = F)
write.table(infants_phenos_VLP, '03a.RESULTS/vOTUs_Bray_NMDS_phenos_mixed_models_FDR_significant.txt', sep='\t', quote=F, row.names = F)
