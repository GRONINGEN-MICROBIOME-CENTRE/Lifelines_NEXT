setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the vOTUs (species-level) stability
# of infant and maternal gut 
# Apologies for a very long code, it's very simple though and
# repeated for infants and mothers.
#############################################################

##############################
# Functions
##############################
source("03.SCRIPTS/NEXT_pilot_FUP_downstream/stability_functions.R")

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
library(gplots)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(MetBrewer)

library(lme4)
library(RLRsim)
library(lmerTest)

library(ggforce)
##############################
# Input data
##############################

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels=c("P7", "B",
                                                                  "M1", "M2", "M3",
                                                                  "M6", "M12"), ordered=T)
VLP_metadata$temperate_perc <- VLP_metadata$temperate_richness/VLP_metadata$viral_richness
RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)

vOTUs_infants <- RPKM_counts_VLP[,VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID]
vOTUs_infants <- vOTUs_infants[rowSums(vOTUs_infants)!=0,]

vOTUs_mothers <- RPKM_counts_VLP[,VLP_metadata[VLP_metadata$Type=='Mother',]$Short_sample_ID]
vOTUs_mothers <- vOTUs_mothers[rowSums(vOTUs_mothers)!=0,]

Demuth <- met.brewer("Demuth")
Monet <- met.brewer('Monet')

##############################
# ANALYSIS
##############################
# vOTUs:
Infant_timepoints <- c('M1', 'M2', 'M3', 'M6', 'M12')
Mother_timepoints <- c("P7", "B", "M1", "M2", "M3")

### Infant vOTUs ###
# N bacterial species retained from M1 over different timepoints
virstability_M1 <- stability_initial(VLP_metadata[VLP_metadata$Type=='Infant',], vOTUs_infants, Infant_timepoints, 'Short_sample_ID', 'Richness')
virstability_M1_stat <- summary_stat_bootstrap(virstability_M1, Infant_timepoints)
virstability_M1_stat$Condition <- factor(virstability_M1_stat$Condition, levels = c('Retained', 'Richness', 'Not_retained'), ordered=T)
virstability_M1_stat <- virstability_M1_stat[row.names(virstability_M1_stat)!="Richness_M1",]

# what is the percentage of species detected at M1 is retained at different timepoints?
virstability_M1_perc_M1 <- virstability_M1[,c('M2', 'M3', 'M6', 'M12')]/virstability_M1[,"M1"]
virstability_M1_perc_M1_stat <- summary_stat_bootstrap(virstability_M1_perc_M1, c('M2', 'M3', 'M6', 'M12'))
colnames(virstability_M1_perc_M1) <- paste0('perc_from_M1_retained_at_', colnames(virstability_M1_perc_M1))

# what is the percentage of richness at given timepoint is occupied by species present at M1 and given timepoint?
virstability_M1_perc_richness <- virstability_M1[,c('M2', 'M3', 'M6', 'M12')]/virstability_M1[,c("Richness_M2", "Richness_M3", "Richness_M6", "Richness_M12")]
virstability_M1_perc_richness_stat <- summary_stat_bootstrap(virstability_M1_perc_richness, c('M2', 'M3', 'M6', 'M12'))
colnames(virstability_M1_perc_richness) <- paste0('M1_retained_perc_of_richness_at_', colnames(virstability_M1_perc_richness))

# what is the relative abundance of vOTUs common with M1 at the given timepoint?
virstability_M1_abundance <- stability_initial(VLP_metadata[VLP_metadata$Type=='Infant',], vOTUs_infants, Infant_timepoints, 'Short_sample_ID', 'abundance')
virstability_M1_abundance_stat <- summary_stat_bootstrap(virstability_M1_abundance, Infant_timepoints)
virstability_M1_abundance_stat <- virstability_M1_abundance_stat[row.names(virstability_M1_abundance_stat)!="Total_space_M1",]
virstability_M1_abundance_stat$Condition <- factor(virstability_M1_abundance_stat$Condition, levels=c('Retained', 'Total_space', 'Not_retained', ordered=T))

# for formal testing
virstability_M1_test <- cbind(virstability_M1_perc_M1, virstability_M1_perc_richness, virstability_M1_abundance[,c("M2", "M3", "M6", "M12")])
colnames(virstability_M1_test)[c(9:12)] <- paste0("RA_M1_retained_at_", colnames(virstability_M1_test)[c(9:12)])
virstability_M1_test <- virstability_M1_test[rowSums(is.na(virstability_M1_test))!=ncol(virstability_M1_test),]
virstability_M1_test$Individual_ID <- row.names(virstability_M1_test)
virstability_M1_test <- melt(virstability_M1_test)
virstability_M1_test$Short_sample_ID <- paste0(substr(virstability_M1_test$Individual_ID, 1,1),
                                               sprintf('%02s', gsub('.*_M', '', virstability_M1_test$variable)),
                                               substr(virstability_M1_test$Individual_ID, 2,5),
                                               'V')
virstability_M1_test$Individual_ID <- NULL
virstability_M1_test$variable <- gsub("_at_.*", "", virstability_M1_test$variable)
virstability_M1_test <- dcast(virstability_M1_test, Short_sample_ID~variable, value.var = "value")
row.names(virstability_M1_test) <- virstability_M1_test$Short_sample_ID
virstability_M1_test$Short_sample_ID <- NULL
virstability_M1_test <- virstability_M1_test[rowSums(is.na(virstability_M1_test))!=ncol(virstability_M1_test),]

# formal testing of dependency on time:
virstability_M1_time <- mixed_models_taxa(VLP_metadata, "Short_sample_ID", virstability_M1_test, "Age_months", "dont_consider_time")
virstability_M1_time$FDR <- NULL # use p-value output as these tests are independent
# formal testing of dependency on time:
virstability_M1_withtime <- mixed_models_taxa(VLP_metadata, 
                                              "Short_sample_ID", 
                                              virstability_M1_test,
                                              c("infant_place_delivery", "infant_ffq_feeding_mode_complex"), 
                                              "time_as_covariate")
# need to change how FDR is applied due to testing three independent virome traits:
virstability_M1_withtime$FDR <- unlist(lapply(colnames(virstability_M1_test), function(i) {
  p.adjust(virstability_M1_withtime[virstability_M1_withtime$Bug == i, "P"], method = "BH")
}))

### this plot contains only 11 infants that had M1 VLP sequenced
pdf('./04.PLOTS/vOTUs_retained_from_M1_over_time.pdf', width=7/2.54, height=10/2.54)
ggplot(virstability_M1_stat[virstability_M1_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_errorbar(data=virstability_M1_stat[virstability_M1_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "N detected vOTUs") + 
  labs(fill='Present vOTUs', color='', alpha='Transparency') + 
  theme_bw() + 
  theme(axis.title =element_text(size=9,face="bold"),
        legend.position = "bottom",
        legend.title = element_text(size=8,face="bold"),
        legend.text = element_text(size=7)) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from M1', 'Not present at M1')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  scale_alpha_manual(values=c(0.9,0.3)) +
  guides(color="none", 
         fill=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'), 
         alpha=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'))
dev.off()

pdf('./04.PLOTS/vOTUs_retained_from_M1_over_time_abundance.pdf', width=7/2.54, height=10/2.54)
ggplot(virstability_M1_abundance_stat[virstability_M1_abundance_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position="identity") +
  geom_errorbar(data=virstability_M1_abundance_stat[virstability_M1_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Present vOTUs', color='',alpha='Transparency') + 
  theme_bw() + 
  theme(axis.title =element_text(size=9,face="bold"),
        legend.position = "bottom",
        legend.title = element_text(size=8,face="bold"),
        legend.text = element_text(size=7)) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from M1', 'Not present at M1')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained', 'Total space')) +
  guides(color="none",
         fill=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'), 
         alpha=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'))
dev.off()


virpresence_from_M1 <- persistent_taxa(VLP_metadata[VLP_metadata$Type=='Infant',], vOTUs_infants, Infant_timepoints, "Short_sample_ID")
virpresence_from_M1_prevalent <- virpresence_from_M1[virpresence_from_M1$M1 > 0.5*length((VLP_metadata[VLP_metadata$Type=='Infant' & VLP_metadata$Timepoint=='M1',]$Short_sample_ID)),]
row.names(virpresence_from_M1_prevalent) <- gsub('.*_length','Length', row.names(virpresence_from_M1_prevalent) )
row.names(virpresence_from_M1_prevalent) <- gsub("\\...\\K\\d+", "", row.names(virpresence_from_M1_prevalent), perl = TRUE)
Okeeffe1_continuous <- met.brewer("OKeeffe1", n=9)
Okeeffe1_continuous <- Okeeffe1_continuous[9:1]

pdf('./04.PLOTS/vOTUs_persistent_over_time.pdf', width=20/2.54, height=20/2.54)
heatmap.2(as.matrix(virpresence_from_M1_prevalent), col=as.vector(Okeeffe1_continuous),
          #dendrogram ="none", 
          trace="none",
          #Rowv=FALSE, 
          Colv=FALSE,  
          margins=c(18,18)
)
dev.off()

## Personal virome fractions
# had to consider only 14 infants to exclude those that have less than 3 timepoints
p_vir_frac_infants <- personal_biome_fraction(VLP_metadata[VLP_metadata$Type=='Infant' & VLP_metadata$N_timepoints>2,], vOTUs_infants, 'Short_sample_ID')

PPB_vir_freq_infants <- as.data.frame(table(unlist(unname( p_vir_frac_infants[["PPB_bacteria"]] ))))

p_vir_frac_infants_individual <- p_vir_frac_infants[["personal_biome_fractions"]]
p_vir_frac_infants_individual$Infant <- row.names(p_vir_frac_infants_individual)
p_vir_frac_infants_individual <- melt(p_vir_frac_infants_individual[,c('PPB', 'TDB', 'Singletons', 'Infant')], id.vars = 'Infant')
p_vir_frac_infants_individual_df <- p_vir_frac_infants_individual[p_vir_frac_infants_individual$variable=='PPB',]
p_vir_frac_infants_individual_df <- p_vir_frac_infants_individual_df[order(p_vir_frac_infants_individual_df$value,decreasing = T),]
p_vir_frac_infants_individual$Infant <- factor(p_vir_frac_infants_individual$Infant, levels=p_vir_frac_infants_individual_df$Infant, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction.pdf', width=10/2.54, height=10/2.54)
ggplot(p_vir_frac_infants_individual, aes(Infant, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity") + 
  xlab(label = "Infant") + 
  ylab(label = "% of infant personal virome") + 
  labs(fill='Fraction') +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c("PPV", "TDV", "Singletons")) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.position = "bottom")
dev.off()

p_vir_frac_infants_per_sample <- p_vir_frac_infants[["personal_biome_per_sample"]]
p_vir_frac_infants_per_sample$Sample <- row.names(p_vir_frac_infants_per_sample)

p_vir_frac_infants_per_sample_ab <- melt(p_vir_frac_infants_per_sample[,c(colnames(p_vir_frac_infants_per_sample)[grep('ab_', colnames(p_vir_frac_infants_per_sample))], "Sample")], id='Sample')
p_vir_frac_infants_per_sample_ab$Infant <- VLP_metadata$Individual_ID[match(p_vir_frac_infants_per_sample_ab$Sample, VLP_metadata$Short_sample_ID)] 
p_vir_frac_infants_per_sample_ab$Timepoint <- VLP_metadata$Timepoint[match(p_vir_frac_infants_per_sample_ab$Sample, VLP_metadata$Short_sample_ID)]
p_vir_frac_infants_per_sample_ab$Infant <- factor(p_vir_frac_infants_per_sample_ab$Infant, levels=p_vir_frac_infants_individual_df$Infant, ordered = T)
p_vir_frac_infants_per_sample_ab$Timepoint <- factor(p_vir_frac_infants_per_sample_ab$Timepoint, levels=Infant_timepoints, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction_abundance.pdf', width=10/2.54, height=7/2.54)
ggplot(p_vir_frac_infants_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Infant samples") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Fraction') +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPV', 'TDV', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()


plot_vir_frac_infants_per_sample_N <- melt(p_vir_frac_infants_per_sample[,c(colnames(p_vir_frac_infants_per_sample)[grep('N_', colnames(p_vir_frac_infants_per_sample))], "Sample")], id='Sample')
plot_vir_frac_infants_per_sample_N$Infant <- VLP_metadata$Individual_ID[match(plot_vir_frac_infants_per_sample_N$Sample, VLP_metadata$Short_sample_ID)] 
plot_vir_frac_infants_per_sample_N$Timepoint <- VLP_metadata$Timepoint[match(plot_vir_frac_infants_per_sample_N$Sample, VLP_metadata$Short_sample_ID)]
plot_vir_frac_infants_per_sample_N$Infant <- factor(plot_vir_frac_infants_per_sample_N$Infant, levels=p_vir_frac_infants_individual_df$Infant, ordered = T)
plot_vir_frac_infants_per_sample_N$Timepoint <- factor(plot_vir_frac_infants_per_sample_N$Timepoint, levels=Infant_timepoints, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction_N.pdf', width=10/2.54, height=7/2.54)
ggplot(plot_vir_frac_infants_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Infant samples") + 
  ylab(label = "% of sample richness") + 
  labs(fill='Fraction') +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPV', 'TDV', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()

## FRACTION STAT: 

p_vir_frac_infants_per_sample$viral_richness <- VLP_metadata$viral_richness[match(p_vir_frac_infants_per_sample$Sample, VLP_metadata$Short_sample_ID)]
p_vir_frac_infants_per_sample <- merge(p_vir_frac_infants_per_sample, as.data.frame(colSums(RPKM_counts_VLP)), by="row.names")
row.names(p_vir_frac_infants_per_sample) <- p_vir_frac_infants_per_sample$Row.names
p_vir_frac_infants_per_sample$Row.names <- NULL
colnames(p_vir_frac_infants_per_sample)[9] <- 'Total_space'
p_vir_frac_infants_per_sample_percs <- cbind(p_vir_frac_infants_per_sample[,c("N_PPB", "N_TDB", "N_Singletons")]/p_vir_frac_infants_per_sample[,"viral_richness"],
                                             p_vir_frac_infants_per_sample[,c("ab_PPB", "ab_TDB", "ab_Singletons")]/p_vir_frac_infants_per_sample[,"Total_space"])
p_vir_frac_infants_per_sample_percs <- p_vir_frac_infants_per_sample_percs*100
p_vir_frac_infants_per_sample_percs_stat <- summary_stat_bootstrap(p_vir_frac_infants_per_sample_percs, colnames(p_vir_frac_infants_per_sample_percs))

p_vir_frac_infants_per_sample_percs$Short_sample_ID <- row.names(p_vir_frac_infants_per_sample_percs)
p_vir_frac_infants_per_sample_percs_melt <- melt(p_vir_frac_infants_per_sample_percs)
p_vir_frac_infants_per_sample_percs_melt$Individual_ID <- VLP_metadata$Individual_ID[match(p_vir_frac_infants_per_sample_percs_melt$Short_sample_ID, VLP_metadata$Short_sample_ID)]

# are the fractions different in size?
model1 <- lmer(value ~ variable + (1|Individual_ID),
               data = p_vir_frac_infants_per_sample_percs_melt[grep('N_', p_vir_frac_infants_per_sample_percs_melt$variable),], REML = F)
summary(model1)

# which is larger?
anova1 <- aov(value ~ variable,
                 data = p_vir_frac_infants_per_sample_percs_melt[grep('N_', p_vir_frac_infants_per_sample_percs_melt$variable),])
TukeyHSD(anova1)

# are the fractions different in abundance?
model2 <- lmer(value ~ variable + (1|Individual_ID),
               data = p_vir_frac_infants_per_sample_percs_melt[grep('ab_', p_vir_frac_infants_per_sample_percs_melt$variable),], REML = F)
summary(model2)
anova2 <- aov(value ~ variable,
              data = p_vir_frac_infants_per_sample_percs_melt[grep('ab_', p_vir_frac_infants_per_sample_percs_melt$variable),])
TukeyHSD(anova2)

### phenotypes: 
infant_fractions_time <- mixed_models_taxa(VLP_metadata, 
                  "Short_sample_ID", 
                  p_vir_frac_infants_per_sample_percs[,-7], 
                  c("Age_months"), 
                  "dont_consider_time")
infant_fractions_time$FDR <- unlist(lapply(c("N_", "ab_"), function(i) {
  p.adjust(infant_fractions_time[grep(i,infant_fractions_time$Bug), "P"], method = "BH")
}))

infant_fractions_withtime <- mixed_models_taxa(VLP_metadata, 
                                              "Short_sample_ID", 
                                              p_vir_frac_infants_per_sample_percs[,-7], 
                                              c("infant_place_delivery", 
                                                "infant_ffq_feeding_mode_complex"), 
                                              "time_as_covariate")
infant_fractions_withtime$FDR <- unlist(lapply(c("N_", "ab_"), function(i) {
  p.adjust(infant_fractions_withtime[grep(i,infant_fractions_withtime$Bug), "P"], method = "BH")
}))

plot_ab_PPB <- merge(p_vir_frac_infants_per_sample_percs, VLP_metadata, by="Short_sample_ID")

pdf('./04.PLOTS/Difference_in_ab_PPV_home_vs_hospital.pdf', width=5/2.54, height=7/2.54)
ggplot(plot_ab_PPB, aes(Timepoint, ab_PPB, fill=infant_place_delivery)) + 
  geom_boxplot(alpha=0.2, outlier.shape = NA) + 
  geom_sina(aes(color=infant_place_delivery), size=0.5) + 
  ylab(label = "Relative abundance of PPV") + 
  labs(fill='Place of delivery', color='Place of delivery') +
  theme_bw() + 
  scale_fill_manual(values=c("#8D9CA3","#B43E23"), labels=c('Home', 'Hospital')) +
  scale_color_manual(values=c("#8D9CA3","#B43E23"), labels=c('Home', 'Hospital')) +
  theme(axis.text=element_text(size=6), 
        axis.title =element_text(size=7,face="bold"),
        legend.position = "bottom", 
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        panel.spacing = unit(0.1, "lines")) + 
  guides(fill = guide_legend(title.position = "top"))
dev.off()

### Mother vOTUs ###
virstability_P7 <- stability_initial(VLP_metadata[VLP_metadata$Type=='Mother',], vOTUs_mothers, Mother_timepoints, 'Short_sample_ID', 'Richness')
virstability_P7_stat <- summary_stat_bootstrap(virstability_P7, Mother_timepoints)
virstability_P7_stat$Condition <- factor(virstability_P7_stat$Condition, levels = c('Retained', 'Richness', 'Not_retained'), ordered=T)
virstability_P7_stat <- virstability_P7_stat[row.names(virstability_P7_stat)!="Richness_P7",]

maternal_richness <- VLP_metadata[VLP_metadata$Type=="Mother",c("Short_sample_ID", "viral_richness")]
row.names(maternal_richness) <- maternal_richness$Short_sample_ID
maternal_richness$Short_sample_ID <- NULL
m_richness_stat <- summary_stat_bootstrap(maternal_richness, "viral_richness")

virstability_P7_perc_P7 <- virstability_P7[,c('B', 'M1', 'M2', 'M3')]/virstability_P7[,"P7"]
virstability_P7_perc_P7 <- melt(virstability_P7_perc_P7)
virstability_P7_perc_P7$variable <- NULL
virstability_P7_perc_P7_stat <- summary_stat_bootstrap(virstability_P7_perc_P7, "value")

virstability_P7_abundance <- stability_initial(VLP_metadata[VLP_metadata$Type=='Mother',], vOTUs_mothers, Mother_timepoints, 'Short_sample_ID', 'abundance')
virstability_P7_abundance_stat <- summary_stat_bootstrap(virstability_P7_abundance, Mother_timepoints)
virstability_P7_abundance_stat$Condition <- factor(virstability_P7_abundance_stat$Condition, levels = c('Retained', 'Total_space', 'Not_retained'), ordered=T)
virstability_P7_abundance_stat <- virstability_P7_abundance_stat[row.names(virstability_P7_abundance_stat)!="Total_space_P7",]

virstability_P7_abundance_melt <- melt(virstability_P7_abundance[,c('B', 'M1', 'M2', 'M3')])
virstability_P7_abundance_melt$variable <- NULL
virstability_P7_abundance_all_stat <- summary_stat_bootstrap(virstability_P7_abundance_melt, 'value')

pdf('./04.PLOTS/vOTUs_retained_from_P7_over_time.pdf', width=7/2.54, height=10/2.54)
ggplot(virstability_P7_stat[virstability_P7_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_errorbar(data=virstability_P7_stat[virstability_P7_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap,
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "N detected vOTUs") + 
  labs(fill='Present vOTUs', color='') + 
  theme_bw() + 
  theme(axis.title =element_text(size=9,face="bold"),
        legend.position = "bottom",
        legend.title = element_text(size=8,face="bold"),
        legend.text = element_text(size=7)) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from P7', 'Not present at P7')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  scale_alpha_manual(values=c(0.9,0.3)) +
  guides(color="none", 
         fill=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'), 
         alpha=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'))
dev.off()


pdf('./04.PLOTS/vOTUs_retained_from_P7_over_time_abundance.pdf', width=7/2.54, height=10/2.54)
ggplot(virstability_P7_abundance_stat[virstability_P7_abundance_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_errorbar(data=virstability_P7_abundance_stat[virstability_P7_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap,
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Present vOTUs', color='') + 
  theme_bw() + 
  theme(axis.title =element_text(size=9,face="bold"),
        legend.position = "bottom",
        legend.title = element_text(size=8,face="bold"),
        legend.text = element_text(size=7)) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from P7', 'Not present at P7')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained', 'Total space')) +
  guides(color="none", 
         fill=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'), 
         alpha=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'))
dev.off()

virpresence_from_P7 <- persistent_taxa(VLP_metadata[VLP_metadata$Type=='Mother',], vOTUs_mothers, Mother_timepoints, "Short_sample_ID")
virpresence_from_P7_prevalent <- virpresence_from_P7[virpresence_from_P7$P7 > 0.5*length((VLP_metadata[VLP_metadata$Type=='Mother' & VLP_metadata$Timepoint=='P7',]$Short_sample_ID)),]
row.names(virpresence_from_P7_prevalent) <- gsub('.*_length','Length', row.names(virpresence_from_P7_prevalent) )
row.names(virpresence_from_P7_prevalent) <- gsub("\\...\\K\\d+", "", row.names(virpresence_from_P7_prevalent), perl = TRUE)

Okeeffe1_continuous <- met.brewer("OKeeffe1", n=19)
Okeeffe1_continuous <- Okeeffe1_continuous[19:1]

pdf('./04.PLOTS/vOTUs_persistent_over_time_mothers.pdf', width=20/2.54, height=20/2.54)
heatmap.2(as.matrix(virpresence_from_P7_prevalent), col=as.vector(Okeeffe1_continuous),
          #dendrogram ="none", 
          trace="none",
          #Rowv=FALSE, 
          Colv=FALSE,  
          margins=c(18,18)
)
dev.off()



p_vir_frac_mothers <- personal_biome_fraction(VLP_metadata[VLP_metadata$Type=='Mother' & VLP_metadata$N_timepoints>2,], vOTUs_mothers, 'Short_sample_ID')

PPB_vir_freq_mothers <- as.data.frame(table(unlist(unname( p_vir_frac_mothers[["PPB_bacteria"]] ))))

p_vir_frac_mothers_individual <- p_vir_frac_mothers[["personal_biome_fractions"]]
p_vir_frac_mothers_individual$Mother <- row.names(p_vir_frac_mothers_individual)
p_vir_frac_mothers_individual <- melt(p_vir_frac_mothers_individual[,c('PPB', 'TDB', 'Singletons', 'Mother')], id.vars = 'Mother')
p_vir_frac_mothers_individual_df <- p_vir_frac_mothers_individual[p_vir_frac_mothers_individual$variable=='PPB',]
p_vir_frac_mothers_individual_df <- p_vir_frac_mothers_individual_df[order(p_vir_frac_mothers_individual_df$value,decreasing = T),]
p_vir_frac_mothers_individual$Mother <- factor(p_vir_frac_mothers_individual$Mother, levels=p_vir_frac_mothers_individual_df$Mother, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction_mothers.pdf', width=10/2.54, height=10/2.54)
ggplot(p_vir_frac_mothers_individual, aes(Mother, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity") + 
  xlab(label = "Mother") + 
  ylab(label = "% of maternal personal virome") + 
  labs(fill='Fraction') +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c("PPV", "TDV", "Singletons")) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.position = "bottom")
dev.off()

p_vir_frac_mothers_per_sample <- p_vir_frac_mothers[["personal_biome_per_sample"]]
p_vir_frac_mothers_per_sample$Sample <- row.names(p_vir_frac_mothers_per_sample)

p_vir_frac_mothers_per_sample_ab <- melt(p_vir_frac_mothers_per_sample[,c(colnames(p_vir_frac_mothers_per_sample)[grep('ab_', colnames(p_vir_frac_mothers_per_sample))], "Sample")], id='Sample')
p_vir_frac_mothers_per_sample_ab$Mother <- VLP_metadata$Individual_ID[match(p_vir_frac_mothers_per_sample_ab$Sample, VLP_metadata$Short_sample_ID)] 
p_vir_frac_mothers_per_sample_ab$Timepoint <- VLP_metadata$Timepoint[match(p_vir_frac_mothers_per_sample_ab$Sample, VLP_metadata$Short_sample_ID)]
p_vir_frac_mothers_per_sample_ab$Mother <- factor(p_vir_frac_mothers_per_sample_ab$Mother, levels=p_vir_frac_mothers_individual_df$Mother, ordered = T)
p_vir_frac_mothers_per_sample_ab$Timepoint <- factor(p_vir_frac_mothers_per_sample_ab$Timepoint, levels=Mother_timepoints, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction_abundance_mothers.pdf', width=25/2.54, height=7/2.54)
ggplot(p_vir_frac_mothers_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Maternal samples") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Fraction') +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPV', 'TDV', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()


plot_vir_frac_mothers_per_sample_N <- melt(p_vir_frac_mothers_per_sample[,c(colnames(p_vir_frac_mothers_per_sample)[grep('N_', colnames(p_vir_frac_mothers_per_sample))], "Sample")], id='Sample')
plot_vir_frac_mothers_per_sample_N$Mother <- VLP_metadata$Individual_ID[match(plot_vir_frac_mothers_per_sample_N$Sample, VLP_metadata$Short_sample_ID)] 
plot_vir_frac_mothers_per_sample_N$Timepoint <- VLP_metadata$Timepoint[match(plot_vir_frac_mothers_per_sample_N$Sample, VLP_metadata$Short_sample_ID)]
plot_vir_frac_mothers_per_sample_N$Mother <- factor(plot_vir_frac_mothers_per_sample_N$Mother, levels=p_vir_frac_mothers_individual_df$Mother, ordered = T)
plot_vir_frac_mothers_per_sample_N$Timepoint <- factor(plot_vir_frac_mothers_per_sample_N$Timepoint, levels=Mother_timepoints, ordered=T)

pdf('./04.PLOTS/vOTUs_personal_fraction_N_mothers.pdf', width=25/2.54, height=7/2.54)
ggplot(plot_vir_frac_mothers_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Maternal samples") + 
  ylab(label = "% of sample richness") + 
  labs(fill='Fraction') +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPV', 'TDV', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()

## FRACTION STAT: 
p_vir_frac_mothers_per_sample$viral_richness <- VLP_metadata$viral_richness[match(p_vir_frac_mothers_per_sample$Sample, VLP_metadata$Short_sample_ID)]

p_vir_frac_mothers_per_sample <- merge(p_vir_frac_mothers_per_sample, as.data.frame(colSums(RPKM_counts_VLP)), by="row.names")
p_vir_frac_mothers_per_sample$Row.names <- NULL
colnames(p_vir_frac_mothers_per_sample)[9] <- 'Total_space'
row.names(p_vir_frac_mothers_per_sample) <- p_vir_frac_mothers_per_sample$Sample
p_vir_frac_mothers_per_sample_percs <- cbind(p_vir_frac_mothers_per_sample[,c("N_PPB", "N_TDB", "N_Singletons")]/p_vir_frac_mothers_per_sample[,"viral_richness"],
                                             p_vir_frac_mothers_per_sample[,c("ab_PPB", "ab_TDB", "ab_Singletons")]/p_vir_frac_mothers_per_sample[,"Total_space"])

p_vir_frac_mothers_per_sample_percs <- p_vir_frac_mothers_per_sample_percs*100
p_vir_frac_mothers_per_sample_percs_stat <- summary_stat_bootstrap(p_vir_frac_mothers_per_sample_percs, colnames(p_vir_frac_mothers_per_sample_percs))

p_vir_frac_mothers_per_sample_percs$Sample <- row.names(p_vir_frac_mothers_per_sample_percs)
p_vir_frac_mothers_per_sample_percs_melt <- melt(p_vir_frac_mothers_per_sample_percs)
p_vir_frac_mothers_per_sample_percs_melt$Individual_ID <- VLP_metadata$Individual_ID[match(p_vir_frac_mothers_per_sample_percs$Sample, VLP_metadata$Short_sample_ID)]

model3 <- lmer(value ~ variable + (1|Individual_ID),
               data = p_vir_frac_mothers_per_sample_percs_melt[grep('N_', p_vir_frac_mothers_per_sample_percs_melt$variable),], REML = F)
summary(model3)

anova3 <- aov(value ~ variable,
              data = p_vir_frac_mothers_per_sample_percs_melt[grep('N_', p_vir_frac_mothers_per_sample_percs_melt$variable),])
TukeyHSD(anova3)
boxplot(p_vir_frac_mothers_per_sample_percs[,c("N_PPB", "N_TDB", "N_Singletons")])

model4 <- lmer(value ~ variable + (1|Individual_ID),
               data = p_vir_frac_mothers_per_sample_percs_melt[grep('ab_', p_vir_frac_mothers_per_sample_percs_melt$variable),], REML = F)
summary(model4)
anova4 <- aov(value ~ variable,
              data = p_vir_frac_mothers_per_sample_percs_melt[grep('ab_', p_vir_frac_mothers_per_sample_percs_melt$variable),])
TukeyHSD(anova4)
boxplot(p_vir_frac_mothers_per_sample_percs[,c("ab_PPB", "ab_TDB", "ab_Singletons")])
##### OUTPUT #####
write.table(virstability_M1, '03a.RESULTS/N_retained_vOTUs_from_M1_infants.txt', sep='\t', quote=F)
write.table(virstability_M1_abundance, '03a.RESULTS/Abundance_retained_vOTUs_from_M1_infants.txt', sep='\t', quote=F)
write.table(virpresence_from_M1, '03a.RESULTS/vOTUs_retained_from_M1_over_time_infants.txt', sep='\t', quote=F)
write.table(PPB_vir_freq_infants, '03a.RESULTS/vOTUs_members_of_PPV_infants_freq.txt', sep='\t', quote=F)
write.table(p_vir_frac_infants_per_sample, '03a.RESULTS/Fractions_of_personal_infant_virome_N_and_abundance.txt', sep='\t', quote=F)

write.table(virstability_P7, '03a.RESULTS/N_retained_vOTUs_from_P7_mothers.txt', sep='\t', quote=F)
write.table(virstability_P7_abundance, '03a.RESULTS/Abundance_retained_vOTUs_from_P7_mothers.txt', sep='\t', quote=F)
write.table(virpresence_from_P7, '03a.RESULTS/vOTUs_retained_from_P7_over_time_mothers.txt', sep='\t', quote=F)
write.table(PPB_vir_freq_mothers, '03a.RESULTS/vOTUs_members_of_PPV_mothers_freq.txt', sep='\t', quote=F)
write.table(p_vir_frac_mothers_per_sample, '03a.RESULTS/Fractions_of_personal_maternal_virome_N_and_abundance.txt', sep='\t', quote=F)
