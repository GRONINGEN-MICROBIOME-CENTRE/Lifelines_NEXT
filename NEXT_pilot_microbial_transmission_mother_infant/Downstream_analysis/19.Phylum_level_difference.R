setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the bacterial phyla prevalence in 
# infant and adult gut
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################
library(vegan)


library(dplyr)
library(lme4)
library(lmerTest)
library(tibble)

library(reshape2)
library(ggplot2)
library(ggforce)
library(ggforestplot)

library(data.table)
library(stringr)

##############################
# Input data
##############################

MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_with_phenos.txt', sep='\t', header=T)
row.names(MGS_metadata) <- MGS_metadata$Short_sample_ID_bact

microbiome_phyla <- read.table('02.CLEAN_DATA/Microbiome_phyla_unfiltred.txt', sep='\t', header=T)
# filtering for presence in more than 5% of microbiome samples:
microbiome_phyla_filt <- microbiome_phyla[(rowSums(microbiome_phyla!=0) > 0.05*ncol(microbiome_phyla)),  ]
microbiome_phyla_filt <- as.data.frame(t(microbiome_phyla_filt))
# CLR-transformation
my_pseudocount_normal=min(microbiome_phyla_filt[microbiome_phyla_filt!=0])/2
microbiome_phyla_filt_CLR<-decostand(microbiome_phyla_filt, "clr", pseudocount=my_pseudocount_normal)


##############################
# ANALYSIS
##############################

# find differential abundant ones
# see which are differential abundant in viruses
# correlate the phyla and or genera in bacteria and phages

df<-merge(MGS_metadata, microbiome_phyla_filt_CLR, by="row.names")
row.names(df)<-df$Row.names
df$Row.names=NULL

Prevalent= names (microbiome_phyla_filt_CLR)
names (MGS_metadata)

pheno_list= "Type"

Overall_result =tibble() 

for (Bug in Prevalent){
  if (! Bug %in% colnames(df)){ next }
  #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
  # print (c(Bug, Prevalence))
  Bug2 = paste(c("`",Bug, "`"), collapse="")
  for ( pheno in pheno_list){
    pheno2 = paste(c("`",pheno, "`"), collapse="")
    df$NG_ID[is.na(df[colnames(df) == pheno]) == F] -> To_keep
    df_pheno = filter(df, NG_ID %in% To_keep )
    Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + (1|NEXT_ID)"), collapse="" )) 
    lmer(Model0, df_pheno) -> resultmodel0
    base_model=resultmodel0
    Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads +",pheno2, "+ (1|NEXT_ID)"), collapse="" ))
    lmer(Model2, df_pheno, REML = F) -> resultmodel2
    M = "Mixed"
    as.data.frame(anova(resultmodel2, base_model))[2,8]->p_simp
    as.data.frame(summary(resultmodel2)$coefficients)[3,1:5] -> Summ_simple
    Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
    rbind(Overall_result, temp_output) -> Overall_result
  }
}

p =as.data.frame(Overall_result)
p$FDR<-p.adjust(p$P, method = "BH")
p$Phylum <- gsub('.*p__', '', p$Bug)
p$significance_level <- NA
p[p$FDR > 0.05,]$significance_level <- 'ns'
p[p$FDR <= 0.05,]$significance_level <- '*'
p[p$FDR <= 0.01,]$significance_level <- '**'
p[p$FDR <= 0.001,]$significance_level <- '***'

#### ALL DETECTED
microbiome_phyla_plot <- as.data.frame(t(microbiome_phyla_filt_CLR))
microbiome_phyla_plot$Phylum <- row.names(microbiome_phyla_plot)
microbiome_phyla_plot_melt <- melt(microbiome_phyla_plot)
microbiome_phyla_plot_melt$Type <- MGS_metadata$Type[match(microbiome_phyla_plot_melt$variable, MGS_metadata$Short_sample_ID_bact)]
microbiome_phyla_plot_melt$Kingdom <- gsub('.*k__','', sapply(strsplit(microbiome_phyla_plot_melt$Phylum, '\\|'), "[", 1))
microbiome_phyla_plot_melt$Phylum <- gsub('.*p__','',microbiome_phyla_plot_melt$Phylum)
microbiome_phyla_plot_melt <- microbiome_phyla_plot_melt[order(microbiome_phyla_plot_melt$Phylum),]
microbiome_phyla_plot_melt$significance_level <- p$significance_level[match(microbiome_phyla_plot_melt$Phylum, p$Phylum )]

pdf('./04.PLOTS/Phyla_difference.pdf', width=10/2.54, height=20/2.54)
ggplot(microbiome_phyla_plot_melt, aes(value, reorder(Phylum, desc(Phylum)), fill=Type)) +
  geom_boxplot(outlier.shape = NA,alpha=0.5) +
  geom_sina(aes(color=Type), size=0.6,alpha=0.5) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(size = 10),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  geom_stripes(odd = "#33333333", even = "#00000000") + 
  facet_grid(rows = vars(Kingdom), scales = "free_y", space='free_y') + 
  geom_text(aes(label = significance_level, x = 14, y = Phylum), size = 4, angle=270) + 
  labs (y="Phylum", x="CLR abundance")
dev.off()

#### ALL SIGNIFICANT
microbiome_phyla_plot <- as.data.frame(t(microbiome_phyla_filt_CLR[,(colnames(microbiome_phyla_filt_CLR) %in% p[p$significance_level!='ns',]$Bug)]))
microbiome_phyla_plot$Phylum <- row.names(microbiome_phyla_plot)
microbiome_phyla_plot_melt <- melt(microbiome_phyla_plot)
microbiome_phyla_plot_melt$Type <- MGS_metadata$Type[match(microbiome_phyla_plot_melt$variable, MGS_metadata$Short_sample_ID_bact)]
microbiome_phyla_plot_melt$Kingdom <- gsub('.*k__','', sapply(strsplit(microbiome_phyla_plot_melt$Phylum, '\\|'), "[", 1))
microbiome_phyla_plot_melt$Phylum <- gsub('.*p__','',microbiome_phyla_plot_melt$Phylum)
microbiome_phyla_plot_melt <- microbiome_phyla_plot_melt[order(microbiome_phyla_plot_melt$Phylum),]
microbiome_phyla_plot_melt$significance_level <- p$significance_level[match(microbiome_phyla_plot_melt$Phylum, p$Phylum )]

pdf('./04.PLOTS/Phyla_difference_significant.pdf', width=14/2.54, height=11/2.54)
ggplot(microbiome_phyla_plot_melt, aes(value, reorder(Phylum, desc(Phylum)), fill=Type)) +
  geom_boxplot(outlier.shape = NA,alpha=0.5) +
  geom_sina(aes(color=Type), size=0.6,alpha=0.5) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(size = 10),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  geom_stripes(odd = "#33333333", even = "#00000000") + 
  facet_grid(rows = vars(Kingdom), scales = "free_y", space='free_y') + 
  geom_text(aes(label = significance_level, x = 14, y = Phylum), size = 4, angle=270) + 
  labs (y="Phylum", x="CLR abundance")
dev.off()

##### OUTPUT #####
write.table(p, '03a.RESULTS/Phylum_vs_Type_Mixed_Models.txt', sep='\t', quote=F, row.names=F)
write.table(microbiome_phyla_filt_CLR, '02.CLEAN_DATA/Microbiome_phyla_filtred_CLR_trans.txt', sep='\t', quote=F)

