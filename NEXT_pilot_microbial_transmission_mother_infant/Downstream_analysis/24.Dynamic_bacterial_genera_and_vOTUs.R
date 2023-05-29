setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore changes in abundance of gut bacterial
# genera and aggregated vOTUs at the level of bacterial genus
# in the developing infant gut 
#############################################################
# Authors: Trishla Sinha & Sanzhima Garmaeva

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
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
  
}

##############################
# Loading libraries
##############################
library(vegan)
library(dplyr)
library(tibble)
library(lme4)
library(RLRsim)
library(lmerTest)

library(reshape2)
library(ggplot2)
library(ggforce)

##############################
# Input data
##############################
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata <- MGS_metadata[MGS_metadata$Type=='Infant',]
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels = c("M1", "M2", "M3",
                                                                    "M6", "M9", "M12"), ordered = T)
row.names(MGS_metadata) <- MGS_metadata$Short_sample_ID_bact


genera <- read.table('02.CLEAN_DATA/Microbiome_genera_unfiltred.txt', sep='\t', header=T)
genera <- genera[,MGS_metadata$Short_sample_ID_bact]
# filtering for presence in more than 10% of microbiome samples:
genera_filt <- genera[(rowSums(genera!=0) > 0.10*ncol(genera)), ]
genera_filt <- as.data.frame(t(genera_filt))
# CLR-transformation
my_pseudocount_normal=min(genera_filt[genera_filt!=0])/2
genera_filt_CLR<-decostand(genera_filt, "clr", pseudocount=my_pseudocount_normal)

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)
VLP_metadata <- VLP_metadata[VLP_metadata$Type=='Infant',]
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels = c("M1", "M2", "M3",
                                                                    "M6", "M12"), ordered = T)
row.names(VLP_metadata) <- VLP_metadata$Short_sample_ID

vOTUS_assigned_BacGen <- read.table('02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_genus.txt', sep='\t', header=T)
row.names(vOTUS_assigned_BacGen) <- vOTUS_assigned_BacGen$Host.genus
vOTUS_assigned_BacGen <- vOTUS_assigned_BacGen[,VLP_metadata$Short_sample_ID]
vOTUS_assigned_BacGen_filt <- vOTUS_assigned_BacGen[(rowSums(vOTUS_assigned_BacGen!=0) > 0.10*ncol(vOTUS_assigned_BacGen)),]
vOTUS_assigned_BacGen_filt <- as.data.frame(t(vOTUS_assigned_BacGen_filt))
#CLR-transformation
my_pseudocount_normal=min(vOTUS_assigned_BacGen_filt[vOTUS_assigned_BacGen_filt!=0])/2
vOTUS_assigned_BacGen_filt_CLR <- decostand(vOTUS_assigned_BacGen_filt, "clr", pseudocount=my_pseudocount_normal)

##############################
# ANALYSIS
##############################
genera_dynamics <- mixed_models_taxa(MGS_metadata, "Short_sample_ID_bact", genera_filt_CLR, "Age_months", "don't consider time as a covariate")

genera_phenos_UPD <- mixed_models_taxa(MGS_metadata, "Short_sample_ID_bact", genera_filt_CLR, c("infant_sex","infant_gestational_age", "infant_birthweight",
                                                                                                "infant_place_delivery", "mother_age_years",
                                                                                                "infant_ffq_feeding_mode_complex"), 'time_as_covariate')


# selection of the top 15 most dynamic genera
dynamic_genera <- genera_dynamics[genera_dynamics$FDR < 0.05 & genera_dynamics$Pheno=='Age_months',]
dynamic_genera <- dynamic_genera[order(dynamic_genera$FDR),]
dynamic_plot <- genera_filt_CLR[,dynamic_genera[1:15,]$Bug]
dynamic_plot <- dynamic_plot[,-grep('CFG', colnames(dynamic_plot))]
dynamic_plot <- merge(dynamic_plot, MGS_metadata[,c("Timepoint", "NEXT_ID")], by='row.names')
dynamic_plot_melt <- melt(dynamic_plot)
dynamic_plot_melt$NEXT_ID=factor(dynamic_plot_melt$NEXT_ID)

dynamic_plot_melt <- dynamic_plot_melt[order(dynamic_plot_melt$NEXT_ID, dynamic_plot_melt$Timepoint),]
dynamic_plot_melt <- dynamic_plot_melt %>% 
  arrange(NEXT_ID, Timepoint)

dynamic_plot_melt$variable <- gsub('.*g__','',dynamic_plot_melt$variable)
  
trial<-dynamic_plot_melt  %>% group_by(variable, Timepoint) %>% 
  summarise_all(funs(mean, n()))

colors <- c("#FF0000", "#ff0080", "#FFA500", "#0000FF", "#FFC0CB", "#008000", "#800000", "#008080", "#DC143C", "#40E0D0", "#6A5ACD","#4B0082", "#FF6347","#6A5ACD")

vars <- unique(trial$variable)
colors_subset <- colors[1:length(vars)]

pdf('./04.PLOTS/Dynamic_bacterial_genera_infants_24_05_2023.pdf', width=16/2.54, height=14/2.54)
ggplot(trial, aes(x = Timepoint, y = value_mean, color = variable, group = variable)) +
  geom_point(size=2) +
  geom_line(size=1.5, alpha=0.8) +
  scale_color_manual(values = colors_subset) +
  ggtitle("")+
  theme_bw()+
  labs(x="Timepoint", y = "Mean CLR Tranformed abundance", color='Genera')+
  theme(
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=12, angle = 60, hjust = 1), 
    legend.title = element_text(color="black", size=12, face="bold"))
dev.off()


genera_phenos_FDR <- genera_phenos_UPD[genera_phenos_UPD$FDR < 0.05,]

# selected genera dependent on feeding mode:
feeding_plot <- genera_filt_CLR[,unique(genera_phenos_FDR[genera_phenos_FDR$Pheno=='infant_ffq_feeding_mode_complex',]$Bug)]
colnames(feeding_plot) <- gsub('.*g__', '', colnames(feeding_plot))
feeding_plot <- merge(feeding_plot, MGS_metadata[,c("Timepoint", 'infant_ffq_feeding_mode_complex')], by='row.names')
feeding_plot <- feeding_plot[!is.na(feeding_plot$infant_ffq_feeding_mode_complex),]

pdf('./04.PLOTS/Genus_Acinetobacter_infants_feeding.pdf', width=13/2.54, height=9/2.54)
ggplot(feeding_plot, aes(Timepoint, Acinetobacter, fill=infant_ffq_feeding_mode_complex)) + 
  labs (y="CLR transformed abundance", x="Timepoint") + 
  geom_boxplot(outlier.shape = NA, alpha=0.5) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Feeding mode", 
                    labels=c('BF', 'FF', 'MF'),
                    values=c("#143F6B", "#F55353", "#FEB139"))
dev.off()

pdf('./04.PLOTS/Genus_Cutibacterium_infants_feeding.pdf', width=13/2.54, height=9/2.54)
ggplot(feeding_plot, aes(Timepoint, Cutibacterium, fill=infant_ffq_feeding_mode_complex)) + 
  labs (y="CLR transformed abundance", x="Timepoint") + 
  geom_boxplot(outlier.shape = NA, alpha=0.5) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Feeding mode", 
                    labels=c('BF', 'FF', 'MF'),
                    values=c("#143F6B", "#F55353", "#FEB139"))
dev.off()

# selected genus dependant on plsce of delivery:

birth_plot <- genera_filt_CLR[,unique(genera_phenos_FDR[genera_phenos_FDR$Pheno %in% c('infant_place_delivery'),]$Bug), drop=F]
colnames(birth_plot) <- gsub('.*g__', '', colnames(birth_plot))
birth_plot <- merge(birth_plot, MGS_metadata[,c("Timepoint", 'infant_place_delivery', 'infant_mode_delivery', 'Age_months', 'DNA_CONC', 'Clean_reads', 'NEXT_ID')], by='row.names')

# for additional analysis for presence of Akkermansia
Akkermansia_binary <- genera[grep('Akkermansia', row.names(genera)),] 
Akkermansia_binary[Akkermansia_binary!=0] <- 1
Akkermansia_binary <- as.data.frame(t(Akkermansia_binary))
colnames(Akkermansia_binary) <- 'Akkermansia_binary'
Akkermansia_binary$Row.names <- row.names(Akkermansia_binary)

birth_plot <- merge(birth_plot, Akkermansia_binary, by='Row.names')


place_mod0 <- lmer(Akkermansia ~ Age_months + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = birth_plot[birth_plot$infant_mode_delivery=='VG',])
place_mod1  = lmer(Akkermansia ~ Age_months + DNA_CONC + Clean_reads + infant_place_delivery + (1|NEXT_ID), REML = F, data = birth_plot[birth_plot$infant_mode_delivery=='VG',])
BIC(place_mod0, place_mod1)
summary(place_mod1) # infant_place_deliveryhospital -2.525e+00  7.358e-01  2.790e+01  -3.431  0.00189 **

place_mod2 <- lmer(Akkermansia_binary ~ Age_months + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = birth_plot[birth_plot$infant_mode_delivery=='VG',])
place_mod3  = lmer(Akkermansia_binary ~ Age_months + DNA_CONC + Clean_reads + infant_place_delivery + (1|NEXT_ID), REML = F, data = birth_plot[birth_plot$infant_mode_delivery=='VG',])
BIC(place_mod2, place_mod3)
summary(place_mod3) # infant_place_deliveryhospital -2.424e-01  7.383e-02  2.806e+01  -3.283  0.00275 **

pdf('./04.PLOTS/Genus_Akkermansia_infants_birthplace_no_CS.pdf', width=12/2.54, height=9/2.54)
ggplot(birth_plot[birth_plot$infant_mode_delivery=='VG',], aes(Timepoint, Akkermansia, fill=infant_place_delivery)) + 
  labs (y="CLR transformed abundance", x="Timepoint") + 
  geom_boxplot(outlier.shape = NA, alpha=0.5) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Place of delivery", 
                    labels=c('Home', 'Hospital'),
                    values=c("#FAE3D9", "#BBDED6"))
dev.off()

###### vOTUs-dynamic genera: 
vOTUs_aggr_gen_dynamics <- mixed_models_taxa(VLP_metadata, "Short_sample_ID", vOTUS_assigned_BacGen_filt_CLR, "Age_months", "don't consider time as a covariate")

vOTUs_gen_phenos <- mixed_models_taxa(VLP_metadata, "Short_sample_ID", vOTUS_assigned_BacGen_filt_CLR, c("infant_sex","infant_gestational_age", "infant_birthweight",
                                                                                                         "infant_place_delivery", "mother_age_years",
                                                                                                         "infant_ffq_feeding_mode_complex"), 'time_as_covariate')

dynamic_vOTUs_gen <- vOTUs_aggr_gen_dynamics[vOTUs_aggr_gen_dynamics$FDR < 0.05 & vOTUs_aggr_gen_dynamics$Pheno=='Age_months',]

dynamic_vOTUs_gen$Genus <- gsub('.*g__', '', dynamic_vOTUs_gen$Bug)

# choosing phages predicted to infect the most dynamic bacterial genera
dynamic_vOTUs_plot <- vOTUS_assigned_BacGen_filt_CLR[,c(dynamic_vOTUs_gen[dynamic_vOTUs_gen$Genus %in% unique(dynamic_plot_melt$variable),]$Bug)]

dynamic_vOTUs_plot <- merge(dynamic_vOTUs_plot, VLP_metadata[,c("Timepoint", "NEXT_ID")], by='row.names')
dynamic_vOTUs_plot_melt <- melt(dynamic_vOTUs_plot)
dynamic_vOTUs_plot_melt$NEXT_ID=factor(dynamic_vOTUs_plot_melt$NEXT_ID)

dynamic_vOTUs_plot_melt <- dynamic_vOTUs_plot_melt[order(dynamic_vOTUs_plot_melt$NEXT_ID, dynamic_vOTUs_plot_melt$Timepoint),]
dynamic_vOTUs_plot_melt <- dynamic_vOTUs_plot_melt %>% 
  arrange(NEXT_ID, Timepoint)

dynamic_vOTUs_plot_melt$variable <- gsub('.*g__','',dynamic_vOTUs_plot_melt$variable)

trial_VLP <-dynamic_vOTUs_plot_melt  %>% group_by(variable, Timepoint) %>% 
  summarise_all(funs(mean, n()))

colors_phages <- c("#FF0000", "#ff0080", "#FFA500", "#008000", "#800000", "#008080", "#40E0D0", "#4B0082")

vars <- unique(trial_VLP$variable)
colors_subset <- colors_phages[1:length(vars)]

pdf('./04.PLOTS/Dynamic_phages_aggregated_by_bacterial_genera_infants.pdf', width=16/2.54, height=14/2.54)
ggplot(trial_VLP, aes(x = Timepoint, y = value_mean, color = variable, group = variable)) +
  geom_point(size=2) +
  geom_line(size=1.5, alpha=0.8) +
  scale_color_manual(values = colors_subset) +
  ggtitle("")+
  theme_bw()+
  labs(x="Timepoint", y = "Mean CLR Tranformed abundance", color='Host Genera')+
  theme(
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=12, angle = 60, hjust = 1), 
    legend.title = element_text(color="black", size=12, face="bold"))
dev.off()

vOTUs_gen_phenos_FDR <- vOTUs_gen_phenos[vOTUs_gen_phenos$FDR < 0.05,]
###### OUTPUT #####
write.table(dynamic_genera, '03a.RESULTS/DYNAMIC_TOP_PREVALENT_RA_GENERA_FIRST_YEAR_NEXT_PILOT_24_05_2023.txt', sep='\t', row.names = F, quote=F)
write.table(genera_phenos_UPD, '03a.RESULTS/TOP_PREVALENT_RA_GENERA_phenos_mixed_models_all_results.txt_24_05_2023.txt', sep='\t', row.names = F, quote=F)
write.table(vOTUs_gen_phenos, '03a.RESULTS/TOP_PREVALENT_vOTUs_aggr_Bac_Gen_phenos_mixed_models_all_results_24_05_2023.txt', sep='\t', row.names = F, quote=F)

