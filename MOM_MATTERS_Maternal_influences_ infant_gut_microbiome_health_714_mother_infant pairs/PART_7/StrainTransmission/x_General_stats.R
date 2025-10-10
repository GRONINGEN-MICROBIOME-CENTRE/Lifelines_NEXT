library(tidyverse)

read_tsv("SGB_babies.tsv") -> To_check
To_check %>% filter(N >= 20) -> To_check

read_rds("tmp/Complete_long_sharing.rds") -> AllSharing

AllSharing %>% filter(SGB %in% To_check$SGB) -> AllSharing
print("Number sharing events")
AllSharing %>% filter(SGB %in% To_check$SGB) %>% group_by(Strain_sharing) %>% summarise(n())
print("Number of SGBs with at least 1 sharing event")
AllSharing %>% filter(Strain_sharing==1) -> AtLeast1
AtLeast1 %>% group_by(SGB) %>% summarise(n())
print("Number of SGBs")
AllSharing %>% group_by(SGB) %>% summarise(n())

print("Whats the range of sharing rate per species")
AllSharing %>% print()
AllSharing %>% group_by(SGB, infant_timepoint_categorical) %>% summarise(Rate= mean(Strain_sharing), N = n() ) -> SharingRateSpTime
AllSharing %>% group_by(SGB) %>% summarise(Rate= mean(Strain_sharing), N = n() ) -> SharingRateSp

SharingRateSp$Rate %>% summary()
SharingRateSp %>% arrange(desc(Rate)) %>% print()
SharingRateSp %>% arrange(Rate) %>% print()

#Save as table
SharingRateSpTime %>% ungroup() %>% rbind(SharingRateSp %>% mutate(infant_timepoint_categorical="overall") ) %>% pivot_wider(names_from = infant_timepoint_categorical, values_from = c(Rate, N)) %>% write_tsv("Results/SGB_sharingRate_table.tsv")

q()

print("Low sharing species")
SharingRateSp %>% filter(SGB %in% c("t__SGB8007_group", "t__SGB8002") ) %>% print()

SharingRateSpTime %>% filter(SGB == "t__SGB1888") %>% print()



print("How many species change sharing through time")
read_tsv("Results/LMERTEST_merged_output.tsv") -> SummaryModelPheno
SummaryModelPheno %>% filter(FDR_all<0.05, Phenotype == "Age_days_infant") -> Res
Res %>% group_by(Estimate <0) %>% summarise(n())
print(Res)


print("Which species have more abundance if they are shared")

read_tsv("Results/AbundanceAssociation/Abundance_vs_sharing_merged_output.tsv") -> Abundanceifshared
Abundanceifshared %>% filter(!Feature == "infant_early_late_timepointlate") %>% filter(N_Early_share>5)  %>% mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr") ) -> Abundanceifshared
Abundanceifshared %>% filter(FDR < 0.05) -> Abundanceifshared_sign
Abundanceifshared_sign %>% filter(Feature == "Strain_sharingyes") %>% print()
Abundanceifshared_sign %>% filter(Feature != "Strain_sharingyes") %>% select(-c(N_Early_noshare, N_Early_share, N_Late_noshare, N_Late_share,FDR_all) ) %>% print()

print("Which species are more likely to be shared in mother has more abundance")

read_tsv("Results/AbundanceAssociationMother/AbundanceMother_vs_sharing_merged_output.tsv") ->  Share_and_Abundance_mother
Share_and_Abundance_mother %>% dim() %>% print()
Share_and_Abundance_mother %>% filter(FDR_all<0.05) -> Share_and_Abundance_mother_sign
print(Share_and_Abundance_mother_sign)


print("Whcih phenotypes are assocaited with individidual species")
SummaryModelPheno %>% filter(Phenotype != "Age_days_infant") -> SummaryModelPheno_other
SummaryModelPheno_other %>% group_by(SGB) %>% summarise(n())
SummaryModelPheno_other %>% filter(FDR_all<0.05) ->  SummaryModelPheno_other_sign
SummaryModelPheno_other_sign %>% group_by(SGB) %>% summarise(n())
SummaryModelPheno_other_sign %>% group_by(Phenotype) %>% summarise(n())
SummaryModelPheno_other_sign %>% arrange(Phenotype) %>% print()

#Mode of delivery plot
#t__SGB17248       birth_deliverybirthcard_mode_binaryVG          3.29        0.737      4.47 0.00000799 NA      birth_deliverybirthcard_mode_binary        no             0.000370 Bifidobacterium_lo…
#__SGB1861        birth_deliverybirthcard_mode_binaryVG



#Supplementary table threhsolds
library(purrr)

# List all files matching the pattern *_threshold.tsv in the Distance_tables directory
file_list <- list.files("Distance_tables", pattern = "_threshold.tsv", full.names = TRUE)
# Read all files and merge them using rbind
merged_data <- map_dfr(file_list, read_tsv)
write_tsv(merged_data, "Results/SupTableThresholds.tsv")



