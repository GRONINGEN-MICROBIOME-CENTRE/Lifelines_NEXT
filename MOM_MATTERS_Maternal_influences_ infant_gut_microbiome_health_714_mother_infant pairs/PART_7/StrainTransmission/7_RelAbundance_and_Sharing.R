library(tidyverse)
library(lmerTest)
source("Common_function.R")


#Get RDS file from command line
args <- commandArgs(trailingOnly = TRUE)
file = args[1]
sp4 = read_rds(file)
#Get SGB name
SGB_name = str_remove(basename(file), "\\.rds") 

#Abundance in early/late timepoints with transmitted species
Abundance = read.table("metadata/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt") %>% rownames_to_column("NG_ID") %>% as_tibble()

grepl(SGB_name , colnames(Abundance) ) -> C_int
if (sum(C_int) == 0) { print("SGB not in abundance table") ; q() }
colnames(Abundance)[ C_int ] -> N


print(paste0(SGB_name," ", N))
Abundance %>% select(c("NG_ID", N)) -> Abundance_merge


Mother_infant_formatting(sp4) -> motherinfant
motherinfant %>% filter(Type_2 == "infant" )  -> Reversed
Remove_col =  colnames(Reversed)[grepl("_1", colnames(Reversed))]
Reversed %>% select(-Remove_col) -> Reversed
colnames(Reversed) = str_replace(colnames(Reversed), "_2", "")

motherinfant %>% filter(Type_1 == "infant" )  -> Forward
Remove_col =  colnames(Forward)[grepl("_2", colnames(Forward))]
Forward %>% select(-Remove_col) -> Forward
colnames(Forward) = str_replace(colnames(Forward), "_1", "")

rbind(Forward %>% rename(NEXT_ID=NEXT_ID_long ) , Reversed) -> motherinfant


left_join(motherinfant, Abundance_merge, by=c("Sample_ID" = "NG_ID") ) %>% filter(! is.na(!!sym(N)) ) -> For_analysis

#Check if enough samples per level
For_analysis %>% group_by(infant_early_late_timepoint, Strain_sharing) %>% summarise(N = n()) -> I
if ( dim(I)[1] < 4){ print("Not enough samples per level") ; print(I) }
I %>% filter(N < 5) -> I2
if ( dim(I2)[1] != 0 ){ print("Not enough samples per level") ; print(I) }


Model = as.formula( paste0( N, " ~ Strain_sharing+infant_early_late_timepoint+ Strain_sharing:infant_early_late_timepoint + (1|NEXT_ID)" ) )
For_analysis[,N] = scale(For_analysis[,N])
lmer(Model, For_analysis) -> Res
summary(Res)$coefficients %>% as.data.frame() %>% rownames_to_column("Feature") %>% mutate(SGB = SGB_name) %>% as_tibble() %>% mutate( N_Early_noshare=filter(I, infant_early_late_timepoint=="early", Strain_sharing=="no")$N  , N_Early_share=filter(I, infant_early_late_timepoint=="early", Strain_sharing=="yes")$N  , N_Late_noshare=filter(I, infant_early_late_timepoint=="late", Strain_sharing=="no")$N  , N_Late_share=filter(I, infant_early_late_timepoint=="late", Strain_sharing=="yes")$N  ) -> Res

print(Res)

Out = paste0("Results/AbundanceAssociation/", SGB_name, ".tsv")
write_tsv(Res, Out)

print("Saving")
For_analysis %>% ggplot(aes_string(x="Strain_sharing", y=N)) + geom_boxplot(outlier.shape = NA, aes(fill=Strain_sharing)) + ggforce::geom_sina(alpha=0.2) + facet_wrap(~infant_early_late_timepoint) + theme_bw() + scale_fill_manual(values= c("#3d67d1", "#cb343f")) + ylab(paste0( str_replace(N, "_", " "), " (CLR)") ) + xlab("Mother-Infant strain sharing") + guides(fill = FALSE)  ->Plot



Out_plot = paste0("plots_associations/AbundanceInfant_", SGB_name, ".pdf")

ggsave(Out_plot, Plot,  height=5, width=5)

