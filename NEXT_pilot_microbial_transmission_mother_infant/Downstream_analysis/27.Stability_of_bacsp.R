setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the bacterial species-level stability
# of infant and maternal gut 


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

MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels=c("P3", "P7", "B",
                                                                  "M1", "M2", "M3",
                                                                  "M6", "M9", "M12"), ordered=T)

species <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)

species_infants <- species[,MGS_metadata[MGS_metadata$Type=='Infant',]$Short_sample_ID_bact]
species_infants <- species_infants[rowSums(species_infants)!=0,]

species_mothers <- species[,MGS_metadata[MGS_metadata$Type=='Mother',]$Short_sample_ID_bact]
species_mothers <- species_mothers[rowSums(species_mothers)!=0,]

Demuth <- met.brewer("Demuth")
Monet <- met.brewer('Monet')

##############################
# ANALYSIS
##############################
# bacteria:
Infant_timepoints <- c('M1', 'M2', 'M3', 'M6', 'M9', 'M12')
Mother_timepoints <- c("P3", "P7", "B", "M1", "M2", "M3")

### Infant Bacteria ###
# N bacterial species retained from M1 over different timepoints
bacstability_M1 <- stability_initial(MGS_metadata[MGS_metadata$Type=='Infant',], species_infants, Infant_timepoints, 'Short_sample_ID_bact', 'Richness')
bacstability_M1_stat <- summary_stat_bootstrap(bacstability_M1, Infant_timepoints)
bacstability_M1_stat$Condition <- factor(bacstability_M1_stat$Condition, levels = c('Retained', 'Richness', 'Not_retained'), ordered=T)
bacstability_M1_stat <- bacstability_M1_stat[row.names(bacstability_M1_stat)!="Richness_M1",]

# what is the percentage of species detected at M1 is retained at different timepoints?
bacstability_M1_perc_M1 <- bacstability_M1[,c('M2', 'M3', 'M6', 'M9', 'M12')]/bacstability_M1[,"M1"]
bacstability_M1_perc_M1_stat <- summary_stat_bootstrap(bacstability_M1_perc_M1, c('M2', 'M3', 'M6', 'M9','M12'))
colnames(bacstability_M1_perc_M1) <- paste0('perc_from_M1_retained_at_', colnames(bacstability_M1_perc_M1))

# what is the percentage of richness at given timepoint is occupied by species present at M1 and given timepoint?
bacstability_M1_perc_richness <- bacstability_M1[,c('M2', 'M3', 'M6', 'M9', 'M12')]/bacstability_M1[,c("Richness_M2", "Richness_M3", "Richness_M6", "Richness_M9", "Richness_M12")]
bacstability_M1_perc_richness_stat <- summary_stat_bootstrap(bacstability_M1_perc_richness, c('M2', 'M3', 'M6', 'M9','M12'))
colnames(bacstability_M1_perc_richness) <- paste0('M1_retained_perc_of_richness_at_', colnames(bacstability_M1_perc_richness))

# what is the relative abundance of vOTUs common with M1 at the given timepoint?
bacstability_M1_abundance <- stability_initial(MGS_metadata[MGS_metadata$Type=='Infant',], species_infants, Infant_timepoints, 'Short_sample_ID_bact', 'abundance')
bacstability_M1_abundance_stat <- summary_stat_bootstrap(bacstability_M1_abundance, Infant_timepoints)
bacstability_M1_abundance_stat$Condition <- factor(bacstability_M1_abundance_stat$Condition, levels = c('Retained', 'Total_space', 'Not_retained'), ordered=T)
bacstability_M1_abundance_stat <- bacstability_M1_abundance_stat[row.names(bacstability_M1_abundance_stat)!="Total_space_M1",]

# for formal testing
bacstability_M1_test <- cbind(bacstability_M1_perc_M1, bacstability_M1_perc_richness, bacstability_M1_abundance[,c("M2", "M3", "M6","M9", "M12")])
colnames(bacstability_M1_test)[c(11:15)] <- paste0("RA_M1_retained_at_", colnames(bacstability_M1_test)[c(11:15)])
bacstability_M1_test <- bacstability_M1_test[rowSums(is.na(bacstability_M1_test))!=ncol(bacstability_M1_test),]
bacstability_M1_test$Individual_ID <- row.names(bacstability_M1_test)
bacstability_M1_test <- melt(bacstability_M1_test)
bacstability_M1_test$Short_sample_ID_bact <- paste0(substr(bacstability_M1_test$Individual_ID, 1,1),
                                               sprintf('%02s', gsub('.*_M', '', bacstability_M1_test$variable)),
                                               substr(bacstability_M1_test$Individual_ID, 2,5),
                                               'B')
bacstability_M1_test$Individual_ID <- NULL
bacstability_M1_test$variable <- gsub("_at_.*", "", bacstability_M1_test$variable)
bacstability_M1_test <- dcast(bacstability_M1_test, Short_sample_ID_bact~variable, value.var = "value")
row.names(bacstability_M1_test) <- bacstability_M1_test$Short_sample_ID_bact
bacstability_M1_test$Short_sample_ID_bact <- NULL
bacstability_M1_test <- bacstability_M1_test[rowSums(is.na(bacstability_M1_test))!=ncol(bacstability_M1_test),]

# formal testing of dependency on time:
bacstability_M1_time <- mixed_models_taxa(MGS_metadata, "Short_sample_ID_bact", bacstability_M1_test, "Age_months", "dont_consider_time")
bacstability_M1_time$FDR <- NULL # use p-value output as these tests are independent
# formal testing of dependency on time:
bacstability_M1_withtime <- mixed_models_taxa(MGS_metadata, 
                                              "Short_sample_ID_bact", 
                                              bacstability_M1_test,
                                              c("infant_place_delivery", "infant_ffq_feeding_mode_complex"), 
                                              "time_as_covariate")

# need to change how FDR is applied due to testing three independent virome traits:
bacstability_M1_withtime$FDR <- unlist(lapply(colnames(bacstability_M1_test), function(i) {
  p.adjust(bacstability_M1_withtime[bacstability_M1_withtime$Bug == i, "P"], method = "BH")
}))


pdf('./04.PLOTS/BacSp_retained_from_M1_over_time.pdf', width=8/2.54, height=10/2.54)
ggplot(bacstability_M1_stat[bacstability_M1_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_errorbar(data=bacstability_M1_stat[bacstability_M1_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "N detected species") + 
  labs(fill='Present species', color='', alpha='Transparency') + 
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


pdf('./04.PLOTS/BacSp_retained_from_M1_over_time_abundance.pdf', width=8/2.54, height=10/2.54)
ggplot(bacstability_M1_abundance_stat[bacstability_M1_abundance_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_errorbar(data=bacstability_M1_abundance_stat[bacstability_M1_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Present species', color='', alpha='Transparency') + 
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


######## will not be used probably ########
stability <- as.data.frame(matrix(NA, nrow=length( unique(MGS_metadata$Individual_ID) ), ncol=6))
row.names(stability) <- unique(MGS_metadata$Individual_ID)
colnames(stability) <- c(paste0('Richness_', Infant_timepoints))

timepoint_pairs <- Map(function(x, y) c(x, y), Infant_timepoints[-length(Infant_timepoints)], Infant_timepoints[-1])

for (i in row.names(stability)) {
 
  for (timepoint in Infant_timepoints) {
    if ( length(MGS_metadata[MGS_metadata$Individual_ID == i & MGS_metadata$Timepoint==timepoint,]$Short_sample_ID_bact) != 0 ) {
      stability[i,paste0('Richness_', timepoint)] <- sum(species_infants[, MGS_metadata[MGS_metadata$Individual_ID == i & MGS_metadata$Timepoint==timepoint,]$Short_sample_ID_bact ] != 0)
    }
  }
  
  
  for (tp_pair in timepoint_pairs) {
    A <- MGS_metadata[MGS_metadata$Individual_ID == i & MGS_metadata$Timepoint == tp_pair[1], ]$Short_sample_ID_bact
    B <- MGS_metadata[MGS_metadata$Individual_ID == i & MGS_metadata$Timepoint == tp_pair[2], ]$Short_sample_ID_bact
    
    if (length(A) != 0 && length(B) != 0) {
      columns <- c(A, B)
      presence <- rowSums(species_infants[, columns] != 0) == 2
      stability[i, paste0(tp_pair[1], "_to_", tp_pair[2])] <- sum(presence)
      
    }
  }
}

# relocating M1_to_M2 for better representation:
stability <- stability[,c(  1:grep('M12', colnames(stability))[1], 
                            grep('M1_to_M2', colnames(stability)),  
                            (grep('M12', colnames(stability))[1] + 1):(grep('M1_to_M2', colnames(stability))-1)) ]

# for visualization:
for (i in Infant_timepoints) {
  
  
  if (i!='M1') {
    stability[,paste0('Not_retained_', i)] <- stability_M1[,paste0('Richness_', i)] - stability[,colnames(stability)[grep(i, colnames(stability))[2] ] ]
  } else {
    stability[,paste0('Not_retained_', i)] <- 0 
  }
  
  
}

stability_stat <- data.frame( colMeans(stability, na.rm = T), apply(stability, 2, sd, na.rm=T) )

colnames(stability_stat) <- c("mean_value", "sd_value")
stability_stat$Timepoint <- gsub(".*_", "", row.names(stability_stat))
stability_stat$Timepoint <- factor(stability_stat$Timepoint, levels=Infant_timepoints, ordered = T)
stability_stat["M1_to_M1",] <- stability_stat[grep('Richness_M1', row.names(stability_stat))[1],]
stability_stat$Condition <- gsub("_M.*", "", row.names(stability_stat))
stability_stat[grep('to', stability_stat$Condition),]$Condition <- 'Retained'

pdf('./04.PLOTS/BacSp_retained_from_previous_over_time_N.pdf', width=7/2.54, height=9/2.54)
ggplot(stability_stat[stability_stat$Condition!='Richness',], aes(Timepoint, mean_value, fill=Condition)) + 
  geom_bar(stat = "identity", position="stack") +
  geom_errorbar(data=stability_stat[stability_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=mean_value - sd_value, 
                    ymax=mean_value + sd_value, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "N detected species") + 
  labs(fill='Present species', color='') + 
  theme_bw() + 
  theme(axis.title =element_text(size=10,face="bold"),
        legend.position = 'bottom', 
        legend.text = element_text(size=5),
        legend.title = element_text(size=6)) + 
  scale_fill_manual(values = c(Demuth[6:5], "#FFFFFF"), labels=c('Not retained', 'Retained from previous', '')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  guides(color="none")
dev.off()

######## will not be used probably ########

bacpresence_from_M1 <- persistent_taxa(MGS_metadata[MGS_metadata$Type=='Infant',], species_infants, Infant_timepoints, "Short_sample_ID_bact")
bacpresence_from_M1_prevalent <- bacpresence_from_M1[bacpresence_from_M1$M1 > 0.5*length((MGS_metadata[MGS_metadata$Type=='Infant' & MGS_metadata$Timepoint=='M1',]$Short_sample_ID_bact)),]
row.names(bacpresence_from_M1_prevalent) <- gsub('.*s__','', row.names(bacpresence_from_M1_prevalent) )
row.names(bacpresence_from_M1_prevalent) <- gsub('_', ' ', row.names(bacpresence_from_M1_prevalent))
Okeeffe1_continuous <- met.brewer("OKeeffe1", n=26)
Okeeffe1_continuous <- Okeeffe1_continuous[26:1]

pdf('./04.PLOTS/BacSp_persistent_over_time.pdf', width=20/2.54, height=20/2.54)
heatmap.2(as.matrix(bacpresence_from_M1_prevalent), col=as.vector(Okeeffe1_continuous),
          #dendrogram ="none", 
          trace="none",
          #Rowv=FALSE, 
          Colv=FALSE,  
          margins=c(18,18)
          )
dev.off()


## Personal microbiome fractions

p_bac_frac_infants <- personal_biome_fraction(MGS_metadata[MGS_metadata$Type=='Infant',], species_infants, 'Short_sample_ID_bact')

PPB_bac_freq_infants <- as.data.frame(table(unlist(unname( p_bac_frac_infants[["PPB_bacteria"]] ))))

p_bac_frac_infants_individual <- p_bac_frac_infants[["personal_biome_fractions"]]
p_bac_frac_infants_individual$Infant <- row.names(p_bac_frac_infants_individual)
p_bac_frac_infants_individual <- melt(p_bac_frac_infants_individual[,c('PPB', 'TDB', 'Singletons', 'Infant')], id.vars = 'Infant')
p_bac_frac_infants_individual_df <- p_bac_frac_infants_individual[p_bac_frac_infants_individual$variable=='PPB',]
p_bac_frac_infants_individual_df <- p_bac_frac_infants_individual_df[order(p_bac_frac_infants_individual_df$value,decreasing = T),]
p_bac_frac_infants_individual$Infant <- factor(p_bac_frac_infants_individual$Infant, levels=p_bac_frac_infants_individual_df$Infant, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction.pdf', width=10/2.54, height=10/2.54)
ggplot(p_bac_frac_infants_individual, aes(Infant, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity") + 
  xlab(label = "Infant") + 
  ylab(label = "% of infant personal bacteriome") + 
  labs(fill='Fraction') +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5]) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.position = "bottom")
dev.off()

p_bac_frac_infants_per_sample <- p_bac_frac_infants[["personal_biome_per_sample"]]
p_bac_frac_infants_per_sample$Sample <- row.names(p_bac_frac_infants_per_sample)

p_bac_frac_infants_per_sample_ab <- melt(p_bac_frac_infants_per_sample[,c(colnames(p_bac_frac_infants_per_sample)[grep('ab_', colnames(p_bac_frac_infants_per_sample))], "Sample")], id='Sample')
p_bac_frac_infants_per_sample_ab$Infant <- MGS_metadata$Individual_ID[match(p_bac_frac_infants_per_sample_ab$Sample, MGS_metadata$Short_sample_ID_bact)] 
p_bac_frac_infants_per_sample_ab$Timepoint <- MGS_metadata$Timepoint[match(p_bac_frac_infants_per_sample_ab$Sample, MGS_metadata$Short_sample_ID_bact)]
p_bac_frac_infants_per_sample_ab$Infant <- factor(p_bac_frac_infants_per_sample_ab$Infant, levels=p_bac_frac_infants_individual_df$Infant, ordered = T)
p_bac_frac_infants_per_sample_ab$Timepoint <- factor(p_bac_frac_infants_per_sample_ab$Timepoint, levels=Infant_timepoints, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction_abundance.pdf', width=25/2.54, height=7/2.54)
ggplot(p_bac_frac_infants_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Infant samples") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Fraction') +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPB', 'TDB', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()


plot_bac_frac_infants_per_sample_N <- melt(p_bac_frac_infants_per_sample[,c(colnames(p_bac_frac_infants_per_sample)[grep('N_', colnames(p_bac_frac_infants_per_sample))], "Sample")], id='Sample')
plot_bac_frac_infants_per_sample_N$Infant <- MGS_metadata$Individual_ID[match(plot_bac_frac_infants_per_sample_N$Sample, MGS_metadata$Short_sample_ID_bact)] 
plot_bac_frac_infants_per_sample_N$Timepoint <- MGS_metadata$Timepoint[match(plot_bac_frac_infants_per_sample_N$Sample, MGS_metadata$Short_sample_ID_bact)]
plot_bac_frac_infants_per_sample_N$Infant <- factor(plot_bac_frac_infants_per_sample_N$Infant, levels=p_bac_frac_infants_individual_df$Infant, ordered = T)
plot_bac_frac_infants_per_sample_N$Timepoint <- factor(plot_bac_frac_infants_per_sample_N$Timepoint, levels=Infant_timepoints, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction_N.pdf', width=25/2.54, height=7/2.54)
ggplot(plot_bac_frac_infants_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Infant samples") + 
  ylab(label = "% of sample richness") + 
  labs(fill='Fraction') +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPB', 'TDB', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()

## FRACTION STAT: 
p_bac_frac_infants_per_sample <- merge(p_bac_frac_infants_per_sample, as.data.frame(colSums(species_infants > 0)), by="row.names")
row.names(p_bac_frac_infants_per_sample) <- p_bac_frac_infants_per_sample$Row.names
p_bac_frac_infants_per_sample <- merge(p_bac_frac_infants_per_sample, as.data.frame(colSums(species_infants)), by="row.names")
p_bac_frac_infants_per_sample$Row.names <- NULL
row.names(p_bac_frac_infants_per_sample) <- p_bac_frac_infants_per_sample$Row.names
p_bac_frac_infants_per_sample$Row.names <- NULL
colnames(p_bac_frac_infants_per_sample)[c(8,9)] <- c('bacterial_richness', 'Total_space')
p_bac_frac_infants_per_sample_percs <- cbind(p_bac_frac_infants_per_sample[,c("N_PPB", "N_TDB", "N_Singletons")]/p_bac_frac_infants_per_sample[,"bacterial_richness"],
                                             p_bac_frac_infants_per_sample[,c("ab_PPB", "ab_TDB", "ab_Singletons")]/p_bac_frac_infants_per_sample[,"Total_space"])

p_bac_frac_infants_per_sample_percs <- p_bac_frac_infants_per_sample_percs*100
p_bac_frac_infants_per_sample_percs_stat <- summary_stat_bootstrap(p_bac_frac_infants_per_sample_percs, colnames(p_bac_frac_infants_per_sample_percs))
p_bac_frac_infants_per_sample_percs$Short_sample_ID_bact <- row.names(p_bac_frac_infants_per_sample_percs)
p_bac_frac_infants_per_sample_percs_melt <- melt(p_bac_frac_infants_per_sample_percs)
p_bac_frac_infants_per_sample_percs_melt$Individual_ID <- MGS_metadata$Individual_ID[match(p_bac_frac_infants_per_sample_percs_melt$Short_sample_ID_bact, MGS_metadata$Short_sample_ID_bact)]

# are the fractions different in size?
model1 <- lmer(value ~ variable + (1|Individual_ID),
               data = p_bac_frac_infants_per_sample_percs_melt[grep('N_', p_bac_frac_infants_per_sample_percs_melt$variable),], REML = F)
summary(model1)

# which is larger?
anova1 <- aov(value ~ variable,
              data = p_bac_frac_infants_per_sample_percs_melt[grep('N_', p_bac_frac_infants_per_sample_percs_melt$variable),])
TukeyHSD(anova1)
boxplot(data=p_bac_frac_infants_per_sample_percs_melt[grep('N_', p_bac_frac_infants_per_sample_percs_melt$variable),], value ~ variable)
# are the fractions different in abundance?
model2 <- lmer(value ~ variable + (1|Individual_ID),
               data = p_bac_frac_infants_per_sample_percs_melt[grep('ab_', p_bac_frac_infants_per_sample_percs_melt$variable),], REML = F)
summary(model2)
anova2 <- aov(value ~ variable,
              data = p_bac_frac_infants_per_sample_percs_melt[grep('ab_', p_bac_frac_infants_per_sample_percs_melt$variable),])
TukeyHSD(anova2)
boxplot(data=p_bac_frac_infants_per_sample_percs_melt[grep('ab_', p_bac_frac_infants_per_sample_percs_melt$variable),], value ~ variable)
### phenotypes: 
infant_fractions_time <- mixed_models_taxa(MGS_metadata, 
                                           "Short_sample_ID_bact", 
                                           p_bac_frac_infants_per_sample_percs[,-7], 
                                           c("Age_months"), 
                                           "dont_consider_time")
infant_fractions_time$FDR <- unlist(lapply(c("N_", "ab_"), function(i) {
  p.adjust(infant_fractions_time[grep(i,infant_fractions_time$Bug), "P"], method = "BH")
}))

infant_fractions_withtime <- mixed_models_taxa(MGS_metadata, 
                                               "Short_sample_ID_bact", 
                                               p_bac_frac_infants_per_sample_percs[,-7], 
                                               c("infant_place_delivery", 
                                                 "infant_ffq_feeding_mode_complex"), 
                                               "time_as_covariate")
infant_fractions_withtime$FDR <- unlist(lapply(c("N_", "ab_"), function(i) {
  p.adjust(infant_fractions_withtime[grep(i,infant_fractions_withtime$Bug), "P"], method = "BH")
})) # no significant associations with phenotypes

### Mother Bacteria ###
bacstability_P3 <- stability_initial(MGS_metadata[MGS_metadata$Type=='Mother',], species_mothers, Mother_timepoints, 'Short_sample_ID_bact', 'Richness')
bacstability_P3_stat <- summary_stat_bootstrap(bacstability_P3, Mother_timepoints)
bacstability_P3_stat$Condition <- factor(bacstability_P3_stat$Condition, levels = c('Retained', 'Richness', 'Not_retained'), ordered=T)
bacstability_P3_stat <- bacstability_P3_stat[row.names(bacstability_P3_stat)!="Richness_P3",]


bacstability_P3_abundance <- stability_initial(MGS_metadata[MGS_metadata$Type=='Mother',], species_mothers, Mother_timepoints, 'Short_sample_ID_bact', 'abundance')
bacstability_P3_abundance_stat <- summary_stat_bootstrap(bacstability_P3_abundance, Mother_timepoints)
bacstability_P3_abundance_stat$Condition <- factor(bacstability_P3_abundance_stat$Condition, levels = c('Retained', 'Total_space', 'Not_retained'), ordered=T)
bacstability_P3_abundance_stat <- bacstability_P3_abundance_stat[row.names(bacstability_P3_abundance_stat)!="Total_space_P3",]


pdf('./04.PLOTS/BacSp_retained_from_P3_over_time.pdf', width=8/2.54, height=10/2.54)
ggplot(bacstability_P3_stat[bacstability_P3_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_errorbar(data=bacstability_P3_stat[bacstability_P3_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "N detected species") + 
  labs(fill='Present species', color='', alpha='Transparency') + 
  theme_bw() + 
  theme(axis.title =element_text(size=9,face="bold"),
        legend.position = "bottom",
        legend.title = element_text(size=8,face="bold"),
        legend.text = element_text(size=7)) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from P3', 'Not present at P3')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  scale_alpha_manual(values=c(0.9,0.3)) +
  guides(color="none",
         fill=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'), 
         alpha=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'))
dev.off()

pdf('./04.PLOTS/BacSp_retained_from_P3_over_time_abundance.pdf', width=8/2.54, height=10/2.54)
ggplot(bacstability_P3_abundance_stat[bacstability_P3_abundance_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_errorbar(data=bacstability_P3_abundance_stat[bacstability_P3_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  xlab(label = "Timepoints") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Present species', color='', alpha='Transparency') + 
  theme_bw() + 
  theme(axis.title =element_text(size=9,face="bold"),
        legend.position = "bottom",
        legend.title = element_text(size=8,face="bold"),
        legend.text = element_text(size=7)) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from P3', 'Not present at P3')) + 
  scale_color_manual(values = c(Demuth[c(3,8)])) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained', 'Total space')) +
  guides(color="none",
         fill=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'), 
         alpha=guide_legend(nrow=2,byrow=TRUE, title.position = 'top'))
dev.off()

bacpresence_from_P3 <- persistent_taxa(MGS_metadata[MGS_metadata$Type=='Mother',], species_mothers, Mother_timepoints, "Short_sample_ID_bact")
bacpresence_from_P3_prevalent <- bacpresence_from_P3[bacpresence_from_P3$P3 > 0.5*length((MGS_metadata[MGS_metadata$Type=='Mother' & MGS_metadata$Timepoint=='P3',]$Short_sample_ID_bact)),]
row.names(bacpresence_from_P3_prevalent) <- gsub('.*s__','', row.names(bacpresence_from_P3_prevalent) )
row.names(bacpresence_from_P3_prevalent) <- gsub('_', ' ', row.names(bacpresence_from_P3_prevalent))
Okeeffe1_continuous <- met.brewer("OKeeffe1", n=21)
Okeeffe1_continuous <- Okeeffe1_continuous[21:1]

pdf('./04.PLOTS/BacSp_persistent_over_time_mothers.pdf', width=20/2.54, height=20/2.54)
heatmap.2(as.matrix(bacpresence_from_P3_prevalent), col=as.vector(Okeeffe1_continuous),
          #dendrogram ="none", 
          trace="none",
          #Rowv=FALSE, 
          Colv=FALSE,  
          margins=c(18,18)
)
dev.off()

p_bac_frac_mothers <- personal_biome_fraction(MGS_metadata[MGS_metadata$Type=='Mother',], species_mothers, 'Short_sample_ID_bact')

PPB_bac_freq_mothers <- as.data.frame(table(unlist(unname( p_bac_frac_mothers[["PPB_bacteria"]] ))))

p_bac_frac_mothers_individual <- p_bac_frac_mothers[["personal_biome_fractions"]]
# since it is impossible to judge about the stability in 2 timepoints
p_bac_frac_mothers_individual <- p_bac_frac_mothers_individual[p_bac_frac_mothers_individual$N_timepoints > 2,]
p_bac_frac_mothers_individual$Mother <- row.names(p_bac_frac_mothers_individual)
p_bac_frac_mothers_individual <- melt(p_bac_frac_mothers_individual[,c('PPB', 'TDB', 'Singletons', 'Mother')], id.vars = 'Mother')
p_bac_frac_mothers_individual_df <- p_bac_frac_mothers_individual[p_bac_frac_mothers_individual$variable=='PPB',]
p_bac_frac_mothers_individual_df <- p_bac_frac_mothers_individual_df[order(p_bac_frac_mothers_individual_df$value,decreasing = T),]
p_bac_frac_mothers_individual$Mother <- factor(p_bac_frac_mothers_individual$Mother, levels=p_bac_frac_mothers_individual_df$Mother, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction_mothers.pdf', width=10/2.54, height=10/2.54)
ggplot(p_bac_frac_mothers_individual, aes(Mother, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity") + 
  xlab(label = "Mother") + 
  ylab(label = "% of maternal personal bacteriome") + 
  labs(fill='Fraction') +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5]) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.position = "bottom")
dev.off()

p_bac_frac_mothers_per_sample <- p_bac_frac_mothers[["personal_biome_per_sample"]]
p_bac_frac_mothers_per_sample$Sample <- row.names(p_bac_frac_mothers_per_sample)

p_bac_frac_mothers_per_sample_ab <- melt(p_bac_frac_mothers_per_sample[,c(colnames(p_bac_frac_mothers_per_sample)[grep('ab_', colnames(p_bac_frac_mothers_per_sample))], "Sample")], id='Sample')
p_bac_frac_mothers_per_sample_ab$Mother <- MGS_metadata$Individual_ID[match(p_bac_frac_mothers_per_sample_ab$Sample, MGS_metadata$Short_sample_ID_bact)] 
# excluding those that have only 2 timepoints: 
p_bac_frac_mothers_per_sample_ab <- p_bac_frac_mothers_per_sample_ab[p_bac_frac_mothers_per_sample_ab$Mother %in% p_bac_frac_mothers_individual$Mother,]
p_bac_frac_mothers_per_sample_ab$Timepoint <- MGS_metadata$Timepoint[match(p_bac_frac_mothers_per_sample_ab$Sample, MGS_metadata$Short_sample_ID_bact)]
p_bac_frac_mothers_per_sample_ab$Mother <- factor(p_bac_frac_mothers_per_sample_ab$Mother, levels=p_bac_frac_mothers_individual_df$Mother, ordered = T)
p_bac_frac_mothers_per_sample_ab$Timepoint <- factor(p_bac_frac_mothers_per_sample_ab$Timepoint, levels=Mother_timepoints, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction_abundance_mothers.pdf', width=25/2.54, height=7/2.54)
ggplot(p_bac_frac_mothers_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Maternal samples") + 
  ylab(label = "Relative abundance") + 
  labs(fill='Fraction') +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPB', 'TDB', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()


plot_bac_frac_mothers_per_sample_N <- melt(p_bac_frac_mothers_per_sample[,c(colnames(p_bac_frac_mothers_per_sample)[grep('N_', colnames(p_bac_frac_mothers_per_sample))], "Sample")], id='Sample')
plot_bac_frac_mothers_per_sample_N$Mother <- MGS_metadata$Individual_ID[match(plot_bac_frac_mothers_per_sample_N$Sample, MGS_metadata$Short_sample_ID_bact)] 
# excluding those that have only 2 timepoints: 
plot_bac_frac_mothers_per_sample_N <- plot_bac_frac_mothers_per_sample_N[plot_bac_frac_mothers_per_sample_N$Mother %in% p_bac_frac_mothers_individual$Mother,]
plot_bac_frac_mothers_per_sample_N$Timepoint <- MGS_metadata$Timepoint[match(plot_bac_frac_mothers_per_sample_N$Sample, MGS_metadata$Short_sample_ID_bact)]
plot_bac_frac_mothers_per_sample_N$Mother <- factor(plot_bac_frac_mothers_per_sample_N$Mother, levels=p_bac_frac_mothers_individual_df$Mother, ordered = T)
plot_bac_frac_mothers_per_sample_N$Timepoint <- factor(plot_bac_frac_mothers_per_sample_N$Timepoint, levels=Mother_timepoints, ordered=T)

pdf('./04.PLOTS/BacSp_personal_fraction_N_mothers.pdf', width=25/2.54, height=7/2.54)
ggplot(plot_bac_frac_mothers_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  xlab(label = "Maternal samples") + 
  ylab(label = "% of sample richness") + 
  labs(fill='Fraction') +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Monet[3:5], labels=c('PPB', 'TDB', 'Singletons')) +
  theme(axis.text.y=element_text(size=12), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=14,face="bold"),
        strip.text.x = element_blank(),
        legend.position = "bottom", 
        panel.spacing = unit(0.1, "lines"))
dev.off()

###### OUTPUT #####
write.table(bacstability_M1, '03a.RESULTS/N_retained_BacSp_from_M1_infants.txt', sep='\t', quote=F)
write.table(bacstability_M1_abundance, '03a.RESULTS/Abundance_retained_BacSp_from_M1_infants.txt', sep='\t', quote=F)
write.table(bacpresence_from_M1, '03a.RESULTS/BacSp_retained_from_M1_over_time_infants.txt', sep='\t', quote=F)
write.table(PPB_bac_freq_infants, '03a.RESULTS/BacSp_members_of_PPB_infants_freq.txt', sep='\t', quote=F)
write.table(p_bac_frac_infants_per_sample, '03a.RESULTS/Fractions_of_personal_infant_bacteriome_N_and_abundance.txt', sep='\t', quote=F)

write.table(bacstability_P3, '03a.RESULTS/N_retained_BacSp_from_P3_mothers.txt', sep='\t', quote=F)
write.table(bacstability_P3_abundance, '03a.RESULTS/Abundance_retained_BacSp_from_P3_mothers.txt', sep='\t', quote=F)
write.table(bacpresence_from_P3, '03a.RESULTS/BacSp_retained_from_P3_over_time_mothers.txt', sep='\t', quote=F)
write.table(PPB_bac_freq_mothers, '03a.RESULTS/BacSp_members_of_PPB_mothers_freq.txt', sep='\t', quote=F)
write.table(p_bac_frac_mothers_per_sample, '03a.RESULTS/Fractions_of_personal_maternal_bacteriome_N_and_abundance.txt', sep='\t', quote=F)
