setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore temperate phages of in gut microbiome and 
# virome
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

mixed_models_taxa_corr_host <- function(metadata, ID, CLR_transformed_data, pheno_list, consider_time) {
  df <- VLP_metadata[VLP_metadata$Type=="Infant",]
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(grep('vOTUs_', colnames(CLR_transformed_data), value = T))
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
        Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + infant_mode_delivery +", gsub('vOTUs_', '', Bug2)," + Age_months + (1|Individual_ID)"), collapse="" )) 
      } else { # else is mainly for associating entities with time alone
        Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + (1|Individual_ID)"), collapse="" )) 
      }
      
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      
      if (consider_time=='time_as_covariate') {
        Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + infant_mode_delivery +", gsub('vOTUs_', '', Bug2)," + Age_months + ",pheno2, "+ (1|Individual_ID)"), collapse="" ))
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

mixed_models_taxa_corr_host_and_prophage <- function(metadata, ID, CLR_transformed_data, pheno_list, consider_time) {
  df <- VLP_metadata[VLP_metadata$Type=="Infant",]
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(grep('vOTUs_', colnames(CLR_transformed_data), value = T))
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
        Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads +",paste0('prophages_',gsub('vOTUs_', '', Bug))," + infant_mode_delivery +", gsub('vOTUs_', '', Bug)," + Age_months + (1|Individual_ID)"), collapse="" )) 
      } else { # else is mainly for associating entities with time alone
        Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + (1|Individual_ID)"), collapse="" )) 
      }
      
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      
      if (consider_time=='time_as_covariate') {
        Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + ",paste0('prophages_',gsub('vOTUs_', '', Bug))," + infant_mode_delivery +", gsub('vOTUs_', '', Bug)," + Age_months + ",pheno2, "+ (1|Individual_ID)"), collapse="" ))
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

confidence_interval_bootstrap <- function(df){
  # creating the data frame with means per column and sd per column
  basic_stat <- data.frame( colMeans(df, na.rm = T), apply(df, 2, sd, na.rm=T) )
  
  # assigning colnames to the data frame
  colnames(basic_stat) <- c("mean_value", "sd_value")
  
  # Set the number of bootstrap samples
  n_bootstrap <- 1000
  
  # calculate 95% confidence interval using bootstrap
  for (i in colnames(df) ) {
    
    # Perform bootstrapping
    boot_samples <- replicate(n_bootstrap, sample(df[,i], replace = TRUE))
    
    # Calculate the mean for each bootstrap sample
    bootstrap_means <- apply(boot_samples, 2, mean, na.rm=T)
    
    # Calculate the 95% confidence interval and derive quantiles
    basic_stat[i,"q1_bootstrap"] <- unname(quantile(bootstrap_means, c(0.025), na.rm = T))
    basic_stat[i,"q2_bootstrap"] <- unname(quantile(bootstrap_means, c(0.975), na.rm = T))
    
    # Calculate bootstrapped mean
    basic_stat[i,'mean_bootstrap'] <- mean(bootstrap_means, na.rm=T)
  }
  
  return(basic_stat)
}
##############################
# Loading libraries
##############################

library(lme4)
library(RLRsim)
library(lmerTest)
library(reshape2)
library(ggplot2)
library(ggforce)
library(ggsignif)

library(dplyr)
library(tibble)

library(EnhancedVolcano)
##############################
# Input data
##############################

VLP_metadata <- read.table("02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt", sep='\t', header=T, row.names = "Short_sample_ID")
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels = c("P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M12"), ordered = T)
VLP_metadata$Type <- as.factor(VLP_metadata$Type)
VLP_metadata$Short_sample_ID <- row.names(VLP_metadata)
VLP_metadata$temperate_perc <- VLP_metadata$temperate_richness/VLP_metadata$viral_richness

MGS_metadata <- read.table("02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt", sep='\t', header=T, row.names = "Short_sample_ID")
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels = c("P3", "P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M9", "M12"), ordered = T)
MGS_metadata$Type <- as.factor(MGS_metadata$Type)
MGS_metadata$Short_sample_ID <- row.names(MGS_metadata)

host_assignment_species <- read.table('02.CLEAN_DATA/Host_prediction_to_genome_m90_refined_taxonomy_no_generalists.txt', sep='\t', header=T)

contigs_metadata <- read.table('02.CLEAN_DATA/VLP_viral_contigs_metadata.txt', sep='\t', header=T)

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)
RPKM_counts_VLP <- RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in% contigs_metadata[contigs_metadata$temperate==1,]$V1,]
RPKM_counts_VLP[RPKM_counts_VLP!=0] <- 1
temperate_binary <- RPKM_counts_VLP

temperate_binary$Virus <- row.names(temperate_binary)
temperate_binary <- merge(host_assignment_species[,c("Virus", "Species")], temperate_binary, by='Virus', all.y=T)
temperate_binary[is.na(temperate_binary$Species),]$Species <- 'Unassigned'
temperate_binary$Virus <- NULL
### aggregating temperate phages by host:
temperate_by_BacSp <- aggregate(.~Species, temperate_binary, sum)
row.names(temperate_by_BacSp) <- temperate_by_BacSp$Species
temperate_by_BacSp$Species <- NULL

#considered_infants_feed has all infants that are considered for testing the differentially prevalent temperate phages
considered_infants_feed <- VLP_metadata[VLP_metadata$Type=='Infant' & !is.na(VLP_metadata$infant_ever_never_breastfed),]$Short_sample_ID
temperate_by_BacSp_infant <- temperate_by_BacSp[,considered_infants_feed]
temperate_by_BacSp_infant <- temperate_by_BacSp_infant[rowSums(temperate_by_BacSp_infant!=0)>round(0.05*ncol(temperate_by_BacSp_infant)),]
temperate_by_BacSp_infant <- as.data.frame(t(temperate_by_BacSp_infant ))

## microbiome
microbiome <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)
microbiome <- microbiome[,MGS_metadata[MGS_metadata$Type=='Infant' & !is.na(MGS_metadata$infant_ever_never_breastfed),]$Short_sample_ID_bact]
microbiome <- microbiome[rowSums(microbiome)!=0,]
microbiome <- as.data.frame(t(microbiome))
microbiome$Short_sample_ID_bact <- row.names(microbiome)
colnames(microbiome) <- gsub('.*s__', '', colnames(microbiome))
microbiome_concurrent <- microbiome[gsub('B','',row.names(microbiome)) %in% VLP_metadata[VLP_metadata$Type=="Infant",]$Universal_fecal_ID,]

# to check the RA and richness of temperates in MGS
RPKM_counts_MGS <- read.table('02.CLEAN_DATA/RPKM_counts_MGS.txt', sep='\t', header = T)
RPKM_counts_MGS <- RPKM_counts_MGS[row.names(RPKM_counts_MGS) %in% contigs_metadata[contigs_metadata$temperate==1,]$V1,]
RPKM_counts_MGS[RPKM_counts_MGS!=0] <- 1
temperate_binary_MGS <- RPKM_counts_MGS

temperate_binary_MGS$Virus <- row.names(temperate_binary_MGS)
temperate_binary_MGS <- merge(host_assignment_species[,c("Virus", "Species")], temperate_binary_MGS, by='Virus', all.y=T)
temperate_binary_MGS[is.na(temperate_binary_MGS$Species),]$Species <- 'Unassigned'
temperate_binary_MGS$Virus <- NULL
### aggregating temperate phages by host:
temperate_by_BacSp_MGS <- aggregate(.~Species, temperate_binary_MGS, sum)
row.names(temperate_by_BacSp_MGS) <- temperate_by_BacSp_MGS$Species
temperate_by_BacSp_MGS$Species <- NULL

considered_infants_feed_MGS <- MGS_metadata[MGS_metadata$Type=='Infant' & !is.na(MGS_metadata$infant_ever_never_breastfed),]$Short_sample_ID
temperate_by_BacSp_infant_MGS <- temperate_by_BacSp_MGS[,considered_infants_feed_MGS]
temperate_by_BacSp_infant_MGS <- temperate_by_BacSp_infant_MGS[rowSums(temperate_by_BacSp_infant_MGS!=0)>round(0.05*ncol(temperate_by_BacSp_infant_MGS)),]
temperate_by_BacSp_infant_MGS <- as.data.frame(t(temperate_by_BacSp_infant_MGS ))

##############################
# ANALYSIS
##############################
# is there a difference in richness of active temperate phages between mothers and infants
type_mod0 <- lm(temperate_richness ~ Type + DNA_CONC + Clean_reads, data = VLP_metadata)
type_mod1  = lmer(temperate_richness ~ Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata)
BIC(type_mod0, type_mod1)
exactLRT(type_mod1,type_mod0)
summary(type_mod1) #Type 3.777e+02  5.980e+01  6.565e+01   6.315 2.67e-08 ***

# dependant on time in babies? (switching to the exact ages for a higher precision)
btmod0 <- lm(temperate_richness ~ Age_days + DNA_CONC + Clean_reads, data = VLP_metadata[VLP_metadata$Type=="Infant",])
btmod1  = lmer(temperate_richness ~ Age_days + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
BIC(btmod0, btmod1)
exactLRT(btmod1,btmod0)
summary(btmod0) 
summary(btmod1) #Age_days     2.573e-01  4.817e-02  7.948e+01   5.342 8.55e-07 ***

# dependant on time in mothers? NS
mtmod0 <- lm(temperate_richness ~ Timepoint_continuous+ DNA_CONC + Clean_reads, data = VLP_metadata[VLP_metadata$Type=="Mother",])
mtmod1  = lmer(temperate_richness ~ Timepoint_continuous+ DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Mother",])
BIC(mtmod0, mtmod1)
exactLRT(mtmod1,mtmod0)
summary(mtmod1) #Timepoint_continuous -3.816e-01  1.255e+01  9.017e+01  -0.030  0.97580 

# Richness of active temperate phages over time
pdf('./04.PLOTS/Dynamics_temperate_richness_infant.pdf', width=16/2.54, height=9/2.54)
ggplot(VLP_metadata, aes(Timepoint, temperate_richness, fill=Type)) + 
  labs (y="Richness of active\ntemperate phages (log10)", x="") + 
  geom_violin() + 
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(legend.position="none") + 
  facet_grid(~Type, scales="free", space = "free") +
  scale_y_continuous(trans='log10') +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type")
dev.off()


## viral richness differences
confidence_interval_bootstrap(data.frame(viral_richness_infant=VLP_metadata[VLP_metadata$Type=='Infant',]$viral_richness))
confidence_interval_bootstrap(data.frame(viral_richness_mother=VLP_metadata[VLP_metadata$Type=='Mother',]$viral_richness))

# is there a difference in richness of viruses between mothers and infants
VR_type_mod0 <- lm(viral_richness ~ Type + DNA_CONC + Clean_reads, data = VLP_metadata)
VR_type_mod1  = lmer(viral_richness ~ Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata)
BIC(VR_type_mod0, VR_type_mod1)
exactLRT(VR_type_mod1,VR_type_mod0)
summary(VR_type_mod1) 
#### FOR SUPPLEMENTARY TABLE ####
write.table(as.data.frame(summary(VR_type_mod1)$coefficients)[,1:5], '05.MANUSCRIPT/Supplementary_tables/MM_vOTUs_richness_Type.txt', sep='\t', quote=F)


pdf('./04.PLOTS/Dynamics_viral_richness_infant.pdf', width=16/2.54, height=9/2.54)
ggplot(VLP_metadata, aes(Timepoint, viral_richness, fill=Type)) + 
  labs (y="Richness of viruses (log10)", x="") + 
  geom_violin() + 
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(legend.position="none") + 
  facet_grid(~Type, scales="free", space = "free") +
  scale_y_continuous(trans='log10') +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type")
dev.off()


temperate_richnes_phenos <- mixed_models_taxa(VLP_metadata[VLP_metadata$Type=="Infant",grep('temperate_richness', colnames(VLP_metadata), invert = T)],
                                              "Short_sample_ID",
                                              VLP_metadata[VLP_metadata$Type=="Infant", "temperate_richness", drop=F], 
                                              c("infant_place_delivery", "infant_ever_never_breastfed", "infant_ffq_feeding_mode_complex"), 
                                              "time_as_covariate")
#### FOR SUPPLEMENTARY TABLE ####
write.table(temperate_richnes_phenos, '05.MANUSCRIPT/Supplementary_tables/MM_temperate_vOTUs_richness_phenos_infants.txt', sep='\t', quote=F)

temprich_inf_virdiv_corr <- lmer(temperate_richness ~ DNA_CONC + Clean_reads + infant_mode_delivery + Age_months + infant_ever_never_breastfed + viral_alpha_diversity +(1|Individual_ID),REML=F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
as.data.frame(summary(temprich_inf_virdiv_corr)$coefficients)[,1:5]
#### FOR SUPPLEMENTARY TABLE ####
write.table(as.data.frame(summary(temprich_inf_virdiv_corr)$coefficients)[,1:5], '05.MANUSCRIPT/Supplementary_tables/MM_temperate_richness_feeding_mode_corrected_viral_diversity.txt', sep='\t', quote=F)

temperate_RA_phenos <- mixed_models_taxa(VLP_metadata[VLP_metadata$Type=="Infant",grep('temperate_RA', colnames(VLP_metadata), invert = T)],
                                              "Short_sample_ID",
                                              VLP_metadata[VLP_metadata$Type=="Infant", "temperate_RA", drop=F], 
                                              c("infant_place_delivery", "infant_ever_never_breastfed", "infant_ffq_feeding_mode_complex"), 
                                              "time_as_covariate")


Richness_Feeding <- melt(VLP_metadata[VLP_metadata$Type=='Infant' & !is.na(VLP_metadata$infant_ever_never_breastfed),c("Timepoint","temperate_richness", "infant_ever_never_breastfed", "Short_sample_ID")])
Richness_Feeding$infant_ever_never_breastfed <- as.factor(Richness_Feeding$infant_ever_never_breastfed)

pdf('./04.PLOTS/Temperate_richness_infant_feeding_mode.pdf', width=13/2.54, height=9/2.54)
ggplot(Richness_Feeding, aes(Timepoint, value, fill=infant_ever_never_breastfed)) + 
  labs (y="Richness of active\ntemperate phages (log10)", x="") + 
  geom_boxplot(outlier.shape = NA) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Feeding mode", 
                    labels=c('Ever breastfed', 'Never breastfed'),
                    values=c("#FAE3D9", "#BBDED6"))
dev.off()


# is there a difference in the relative abundance of active temperate phages between mothers and infants?
type_mod0_RA <- lm(temperate_RA ~ Type + DNA_CONC + Clean_reads, data = VLP_metadata)
type_mod1_RA  = lmer(temperate_RA ~ Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata)
BIC(type_mod0_RA, type_mod1_RA)
exactLRT(type_mod1_RA,type_mod0_RA)
summary(type_mod1_RA) #Type -2.220e+01  2.923e+00  4.694e+01  -7.595 1.04e-09 ***
#### FOR SUPPLEMENTARY TABLE ####
write.table(as.data.frame(summary(type_mod1_RA)$coefficients)[,1:5], '05.MANUSCRIPT/Supplementary_tables/MM_RA_active_temperate_phages_by_Type.txt', sep='\t', quote=F)

confidence_interval_bootstrap(data.frame(temperate_RA_infant_M1_M2_M3=VLP_metadata[VLP_metadata$Type=="Infant" & VLP_metadata$Timepoint %in% c('M1', 'M2', 'M3'),]$temperate_RA))

confidence_interval_bootstrap(data.frame(temperate_RA_infant_M6=VLP_metadata[VLP_metadata$Timepoint=="M6",]$temperate_RA))

confidence_interval_bootstrap(data.frame(temperate_RA_infant_M12=VLP_metadata[VLP_metadata$Timepoint=="M12",]$temperate_RA))

confidence_interval_bootstrap(data.frame(temperate_RA_mother=VLP_metadata[VLP_metadata$Type=="Mother",]$temperate_RA))

# dependant on time in babies? (switching to the exact ages for a higher precision)
btmod0_RA <- lm(temperate_RA ~ Age_days + DNA_CONC +  Clean_reads, data = VLP_metadata[VLP_metadata$Type=="Infant",])
btmod1_RA  = lmer(temperate_RA ~ Age_days + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
BIC(btmod0_RA, btmod1_RA)
exactLRT(btmod1_RA,btmod0_RA)
summary(btmod0_RA) 
summary(btmod1_RA) #Age_days     -8.744e-02  2.310e-02  8.076e+01  -3.785 0.000294 ***
#### FOR SUPPLEMENTARY TABLE ####
write.table(as.data.frame(summary(btmod1_RA)$coefficients)[,1:5], '05.MANUSCRIPT/Supplementary_tables/MM_RA_active_temperate_phages_in_infants_over_time.txt', sep='\t', quote=F)

# dependant on time in mothers? NS
mtmod0 <- lm(temperate_RA ~ Timepoint_continuous+ DNA_CONC + Clean_reads, data = VLP_metadata[VLP_metadata$Type=="Mother",])
mtmod1  = lmer(temperate_RA ~ Timepoint_continuous+ DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Mother",])
BIC(mtmod0, mtmod1)
exactLRT(mtmod1,mtmod0)
summary(mtmod1) #Timepoint_continuous -5.295e-01  3.108e-01  9.399e+01  -1.704   0.0917 . 
#### FOR SUPPLEMENTARY TABLE ####
write.table(as.data.frame(summary(mtmod1)$coefficients)[,1:5], '05.MANUSCRIPT/Supplementary_tables/MM_RA_active_temperate_phages_in_mothers_over_time.txt', sep='\t', quote=F)

# still different at M12?
wilcox.test(VLP_metadata[VLP_metadata$Timepoint=='M12',]$temperate_RA, 
            VLP_metadata[VLP_metadata$Type=='Mother',]$temperate_RA) # W = 1731, p-value = 0.04502


# Growth of relative abundance of active temperate phages
pdf('./04.PLOTS/Dynamics_of_temperate_RelAb.pdf', width=12/2.54, height=8.5/2.54)
ggplot(VLP_metadata, aes(Timepoint, temperate_RA, fill=Type)) + 
  labs (y="Relative abundance of\nactive temperate phages", x="") + 
  #geom_violin() + 
  geom_boxplot(width=0.5,outlier.shape = NA, alpha=0.5) + 
  geom_sina(aes(color=Type), size=0.6) +
  facet_grid(~Type, scales = 'free_x') +
  theme_bw()+
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type")
dev.off()

#### richness of temperate phages depending on the feeding mode
temperate_by_BacSp_infant <- temperate_by_BacSp_infant[,colnames(temperate_by_BacSp_infant) %in% gsub('_', ' ',colnames(microbiome)),]
diff_ab_temperate <- mixed_models_taxa(VLP_metadata[VLP_metadata$Type=="Infant",], 
                                       "Short_sample_ID", 
                                       temperate_by_BacSp_infant, 
                                       c('infant_ever_never_breastfed'),
                                       "time_as_covariate")
#### FOR SUPPLEMENTARY TABLE ####
write.table(diff_ab_temperate, '05.MANUSCRIPT/Supplementary_tables/MM_diff_N_active_temperate_aggr_BacSp_feeding_mode.txt', sep='\t', quote=F)


# rerunning the association with the correction by host:
temperate_by_BacSp_infant_to_merge <- temperate_by_BacSp_infant
colnames(temperate_by_BacSp_infant_to_merge) <- gsub(' ', '_', colnames(temperate_by_BacSp_infant_to_merge))
temperate_by_BacSp_infant_to_merge <- temperate_by_BacSp_infant_to_merge[,colnames(temperate_by_BacSp_infant_to_merge) %in% colnames(microbiome_concurrent)]
row.names(temperate_by_BacSp_infant_to_merge) <- gsub('V', '', row.names(temperate_by_BacSp_infant_to_merge))

microbiome_concurrent <- microbiome_concurrent[,colnames(microbiome_concurrent) %in% colnames(temperate_by_BacSp_infant_to_merge)]
row.names(microbiome_concurrent) <- gsub('B', '', row.names(microbiome_concurrent))

colnames(temperate_by_BacSp_infant_to_merge) <- paste0('vOTUs_',colnames(temperate_by_BacSp_infant_to_merge))

temperate_by_BacSp_vs_hosts <- merge(temperate_by_BacSp_infant_to_merge, microbiome_concurrent, by='row.names')
row.names(temperate_by_BacSp_vs_hosts) <- temperate_by_BacSp_vs_hosts$Row.names
temperate_by_BacSp_vs_hosts$Row.names <- NULL

diff_ab_temperate_host_corrected <- mixed_models_taxa_corr_host(VLP_metadata[VLP_metadata$Type=="Infant",], 
                                       "Universal_fecal_ID", 
                                       temperate_by_BacSp_vs_hosts, 
                                       c('infant_ever_never_breastfed'),
                                       "time_as_covariate")
#### FOR SUPPLEMENTARY TABLE ####
write.table(diff_ab_temperate_host_corrected, '05.MANUSCRIPT/Supplementary_tables/MM_diff_N_active_temperate_aggr_BacSp_corr_host_abundance_feeding_mode.txt', sep='\t', quote=F)

diff_ab_temperate_hosts <- mixed_models_taxa(MGS_metadata[MGS_metadata$Type=="Infant",], 
                                       "Universal_fecal_ID", 
                                       microbiome_concurrent[,colSums(microbiome_concurrent)!=0], 
                                       c('infant_ever_never_breastfed'),
                                       "time_as_covariate")
#### FOR SUPPLEMENTARY TABLE ####
write.table(diff_ab_temperate_hosts, '05.MANUSCRIPT/Supplementary_tables/MM_diff_ab_BacSp_hosts_of_temperates_feeding_mode.txt', sep='\t', quote=F)

diff_ab_temperate$log2FC <- NA

# necessary to calculate log2FC in cases when a bug is not present in one of the groups:
considered_infants_metadata <- VLP_metadata[considered_infants_feed,]
ever_fed <- considered_infants_metadata[considered_infants_metadata$infant_ever_never_breastfed=='ever_BF',]$Short_sample_ID
never_fed <- considered_infants_metadata[considered_infants_metadata$infant_ever_never_breastfed=='never_BF',]$Short_sample_ID

min_ever <- min(colMeans(temperate_by_BacSp_infant[ever_fed,])[colMeans(temperate_by_BacSp_infant[ever_fed,]) > 0])
min_never <- min(colMeans(temperate_by_BacSp_infant[never_fed,])[colMeans(temperate_by_BacSp_infant[never_fed,]) > 0])

for (i in diff_ab_temperate$Bug) {
  
  A <- mean( temperate_by_BacSp_infant[ ever_fed,  i ] ) 
  B <- mean( temperate_by_BacSp_infant[ never_fed,  i ] ) 
  
  if ( abs((log2(A) - log2(B)))==Inf ) {
    if (A==0) {
      A <- min_ever
    } else {
      B <- min_never
    }
    
  }
  diff_ab_temperate[diff_ab_temperate$Bug==i,]$log2FC <-  log2(B) - log2(A) 
}

row.names(diff_ab_temperate) <- diff_ab_temperate$Bug

pdf('./04.PLOTS/Differentially_prevalent_temperate_feeding.pdf', width=20/2.54, height=15/2.54)
EnhancedVolcano(diff_ab_temperate,
                lab = rownames(diff_ab_temperate),
                x = 'log2FC',
                y = 'FDR', 
                ylim = c(0, -log10(1e-3)),
                xlim = c(-6,6),
                title = 'Differentially prevalent temperate vOTUs\naggregated by host between exclusively\nformula fed and breastfed infants',
                pCutoff = 0.05,
                #FCcutoff = 1.75,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                boxedLabels = TRUE) + 
  ggplot2::labs(y=expression(-Log["10"]~FDR))

dev.off()

# is it different when prophages are taken into account?
colnames(temperate_by_BacSp_infant_MGS) <- gsub(' ', '_', colnames(temperate_by_BacSp_infant_MGS))
temperate_by_BacSp_infant_MGS <- temperate_by_BacSp_infant_MGS[,colnames(temperate_by_BacSp_infant_MGS) %in% colnames(temperate_by_BacSp_vs_hosts)]
row.names(temperate_by_BacSp_infant_MGS) <- gsub('M','',row.names(temperate_by_BacSp_infant_MGS))

temperate_by_BacSp_infant_to_merge_to_MGS <- temperate_by_BacSp_infant_to_merge[,gsub('vOTUs_', '', colnames(temperate_by_BacSp_infant_to_merge)) %in% colnames(temperate_by_BacSp_infant_MGS)]
microbiome_concurrent_to_merge <- microbiome_concurrent[,colnames(microbiome_concurrent) %in% colnames(temperate_by_BacSp_infant_MGS)]

colnames(temperate_by_BacSp_infant_MGS) <- paste0('prophages_',colnames(temperate_by_BacSp_infant_MGS))

temperate_by_BacSp_infant_host_prophage <- merge(temperate_by_BacSp_infant_to_merge_to_MGS, microbiome_concurrent_to_merge, by='row.names')
row.names(temperate_by_BacSp_infant_host_prophage) <- temperate_by_BacSp_infant_host_prophage$Row.names
temperate_by_BacSp_infant_host_prophage$Row.names <- NULL
temperate_by_BacSp_infant_host_prophage <- merge(temperate_by_BacSp_infant_host_prophage, temperate_by_BacSp_infant_MGS, by='row.names')
row.names(temperate_by_BacSp_infant_host_prophage) <- temperate_by_BacSp_infant_host_prophage$Row.names
temperate_by_BacSp_infant_host_prophage$Row.names <- NULL


###
diff_ab_temperate_corect_host_correct_prophage <- mixed_models_taxa_corr_host_and_prophage(VLP_metadata[VLP_metadata$Type=="Infant",], 
                                                                                           "Universal_fecal_ID", 
                                                                                           temperate_by_BacSp_infant_host_prophage, 
                                                                                           c('infant_ever_never_breastfed'),
                                                                                           "time_as_covariate")
#### FOR SUPPLEMENTARY TABLE ####
write.table(diff_ab_temperate_corect_host_correct_prophage, '05.MANUSCRIPT/Supplementary_tables/MM_diff_N_temperate_phages_corrected_hosts_ab_prophages_N_feeding_mode.txt', sep='\t', quote=F)


diff_ab_temperate_MGS_host_corrected <- mixed_models_taxa_corr_host(VLP_metadata[VLP_metadata$Type=="Infant",], 
                                                                "Universal_fecal_ID", 
                                                                temperate_by_BacSp_MGS_vs_hosts, 
                                                                c('infant_ever_never_breastfed'),
                                                                "time_as_covariate")


### is RA of temperates is still different when prophages are taken into account?
btmod1_mgs_RA  = lmer(temperate_RA_MGS ~ Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = MGS_metadata)
summary(btmod1_mgs_RA) #TypeMother  -2.117e+00  1.032e+00  8.003e+01  -2.052   0.0434 *
#### FOR SUPPLEMENTARY TABLE ####
write.table(as.data.frame(summary(btmod1_mgs_RA)$coefficients)[,1:5], '05.MANUSCRIPT/Supplementary_tables/MM_RA_temperate_phages_MGS_by_Type.txt', sep='\t', quote=F)


ggplot(MGS_metadata, aes(Timepoint, temperate_RA_MGS, fill=Type)) + 
  labs (y="Relative abundance of\ntemperate phages (in MGS)", x="") + 
  #geom_violin() + 
  geom_boxplot(width=0.5,outlier.shape = NA, alpha=0.5) + 
  geom_sina(aes(color=Type), size=0.6) +
  facet_grid(~Type, scales = 'free_x') +
  theme_bw()+
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type")

### is richness of temperates is still different when prophages are taken into account?
btmod1_mgs_richness  = lmer(temperate_richness_MGS ~ Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = MGS_metadata)
summary(btmod1_mgs_richness) #TypeMother   5.324e+02  3.175e+01  6.956e+01  16.767  < 2e-16 ***

ggplot(MGS_metadata, aes(Timepoint, temperate_richness_MGS, fill=Type)) + 
  labs (y="Richness of temperate phages (in MGS)", x="") + 
  #geom_violin() + 
  geom_boxplot(width=0.5,outlier.shape = NA, alpha=0.5) + 
  geom_sina(aes(color=Type), size=0.6) +
  facet_grid(~Type, scales = 'free_x') +
  theme_bw()+
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type")


fmod1_mgs_richness  = lmer(temperate_richness_MGS ~ infant_ever_never_breastfed + bacterial_alpha_diversity  + Age_days + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = MGS_metadata[MGS_metadata$Type=="Infant" & !is.na(MGS_metadata$infant_ever_never_breastfed),])
summary(fmod1_mgs_richness) #infant_ever_never_breastfednever_BF  3.880e+01  1.583e+01  2.526e+01   2.451 0.021508 * 
#### FOR SUPPLEMENTARY TABLE ####
write.table(as.data.frame(summary(fmod1_mgs_richness)$coefficients)[,1:5], '05.MANUSCRIPT/Supplementary_tables/MM_temperate_richness_MGS_feeding_mode.txt', sep='\t', quote=F)


pdf('./04.PLOTS/Temperate_richness_MGS_infant_feeding_mode.pdf', width=13/2.54, height=9/2.54)
ggplot(MGS_metadata[considered_infants_feed_MGS,], aes(Timepoint, temperate_richness_MGS, fill=infant_ever_never_breastfed)) + 
  labs (y="Richness of temperate phages (log10), (in MGS)", x="") + 
  geom_boxplot(outlier.shape = NA, alpha=0.5) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  scale_y_continuous(trans='log10') +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Feeding mode", 
                    labels=c('Ever breastfed', 'Never breastfed'),
                    values=c("#FAE3D9", "#BBDED6"))
dev.off()


#### richness of temperate phages depending on the feeding mode in MGS
diff_ab_temperate_MGS <- mixed_models_taxa(MGS_metadata[considered_infants_feed_MGS,], "Short_sample_ID", temperate_by_BacSp_infant_MGS, c('infant_ever_never_breastfed'))

###### OUTPUT #####
write.table(diff_ab_temperate, '03a.RESULTS/Differentially_prevalent_aggregated_temperates_Mixed_Models_Feeding.txt', sep='\t', quote=F, row.names=F)
write.table(diff_ab_temperate_MGS, '03a.RESULTS/Differentially_prevalent_aggregated_temperates_MGS_Mixed_Models_Feeding.txt', sep='\t', quote=F, row.names=F)

##### RESULTS #####
write.table(temperate_richnes_phenos, '03a.RESULTS/Active_temperate_phages_richness_over_time_phenos_mixed_models.txt', sep='\t', quote=F, row.names=F)
write.table(temperate_RA_phenos, '03a.RESULTS/Active_temperate_phages_RA_over_time_phenos_mixed_models.txt', sep='\t', quote=F, row.names=F)

##### FOR VISUALIZATION #####
write.table(diff_ab_temperate, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig3C_Differentially_prevalent_aggregated_temperates_Mixed_Models_Feeding.txt', sep='\t', quote=F)
