library(tidyverse)
library(lmerTest)
library(ggforce)
source("Common_function.R")


#Get RDS file from command line
args <- commandArgs(trailingOnly = TRUE)
file = args[1]
sp4 = read_rds(file)
#Get SGB name
SGB_name = str_remove(basename(file), "\\.rds") 

Abundance_mother = read.table("metadata/NEXT_metaphlan_4_CLR_transformed_fil_SGB_mothers_03_08_2023.txt") %>% rownames_to_column("NG_ID") %>% as_tibble()


grepl(SGB_name , colnames(Abundance_mother) ) -> C_int
if (sum(C_int) == 0) { print("SGB not in abundance table") ; q() }
colnames(Abundance_mother)[ C_int ] -> N


print(paste0(SGB_name," ", N))
Abundance_mother %>% select(c("NG_ID", N)) -> Abundance_mother_merge

Mother_infant_formatting(sp4,select_mother = F ) -> motherinfant
motherinfant %>% filter(infant_early_late_timepoint == "early") -> motherinfant

motherinfant %>% filter(Type_2 == "mother" )  -> Reversed
Remove_col =  colnames(Reversed)[grepl("_1", colnames(Reversed))]
Reversed %>% select(-Remove_col) -> Reversed
colnames(Reversed) = str_replace(colnames(Reversed), "_2", "")

motherinfant %>% filter(Type_1 == "mother" )  -> Forward
Remove_col =  colnames(Forward)[grepl("_2", colnames(Forward))]
Forward %>% select(-Remove_col) -> Forward
colnames(Forward) = str_replace(colnames(Forward), "_1", "")

rbind(Forward %>% rename(NEXT_ID=NEXT_ID_long ) , Reversed) -> motherinfant


left_join(motherinfant, Abundance_mother_merge, by=c("Sample_ID" = "NG_ID") ) %>% filter(! is.na(!!sym(N)) ) -> For_analysis

#Check only transmission in the first two timepoints
#For_analysis %>% filter(infant_timepoint_categorical %in% c("M1", "W2" ) ) -> For_analysis


#Check if enough samples per level
For_analysis %>% group_by(Strain_sharing) %>% summarise(N = n()) -> I
if ( dim(I)[1] < 2){ print("Not enough samples per level") ; print(I) }
I %>% filter(N < 5) -> I2
if ( dim(I2)[1] != 0 ){ print("Not enough samples per level") ; print(I) }

left_join(motherinfant, Abundance_mother_merge, by=c("Sample_ID" = "NG_ID") ) %>% filter(! is.na(N) ) -> For_analysis
For_analysis[,N] = scale(For_analysis[,N])
For_analysis %>% mutate(Strain_sharing =  factor(Strain_sharing, levels=c("no", "yes") )) -> For_analysis
Model = as.formula( paste0( "Strain_sharing ~ infant_timepoint_categorical+mother_timepoint_categorical  +  (1|NEXT_ID) +", N ) )

Res = Fit_glmer(Model, For_analysis)
if (Res[[2]] == F) { Res = summary(Res[[1]])$coefficients %>% as.data.frame()
} else { q() }

Res %>% rownames_to_column("Feature") %>% filter(Feature == N) %>% as_tibble() %>% mutate(SGB = SGB_name) %>% mutate( Sharing_families = filter(I, Strain_sharing == "yes")$N, NoSharing_families = filter(I, Strain_sharing == "no")$N) -> Res

print(Res)

Out = paste0("Results/AbundanceAssociationMother/", SGB_name, ".tsv")
write_tsv(Res, Out)

For_analysis %>% ggplot(aes_string(x="Strain_sharing", y=N)) + geom_boxplot(outlier.shape = NA, aes(fill=Strain_sharing) ) + ggforce::geom_sina(alpha=0.2) + theme_bw() + scale_fill_manual(values= c("#3d67d1", "#cb343f")) + ylab(paste0( str_replace(N, "_", " "), " (CLR)") ) + xlab("Mother-Infant strain sharing") + guides(fill = FALSE) ->Plot

Out_plot = paste0("plots_associations/AbundanceMother_", SGB_name, ".pdf")
print(Out_plot)
ggsave(Out_plot, Plot, height=5, width=2.5)


For_analysis %>% ggplot(aes_string(x= "log10(dist)", y= N )) + geom_point(aes(col=Strain_sharing)) + theme_bw() + geom_smooth(method='lm', formula= y~x) -> Plot2
Out_plot2 = paste0("plots_associations/AbundanceMother_", SGB_name, "_scatter.pdf")
ggsave(Out_plot2, Plot2)


