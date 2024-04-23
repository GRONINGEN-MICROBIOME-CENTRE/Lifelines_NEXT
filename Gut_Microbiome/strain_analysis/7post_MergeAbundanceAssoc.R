library(tidyverse)

folder ="Results/AbundanceAssociation/"
Output = paste0(folder, "Abundance_vs_sharing_merged_output.tsv")

if (! file.exists(Output)){
        file_list = list.files(path=folder,pattern="*.tsv$") # Names of all the species you have distance matrices for 
        Summary_statistics = tibble()
        for (i in file_list){
                if ( grepl("merged", i) ){ next }
                read_tsv(paste0(folder, "/", i)) -> d
                if (dim(d)[1] == 1 ){ 
                        if (is.na(d$Estimate)) { next }
                 }  
                d %>% rbind(Summary_statistics, .) -> Summary_statistics        
        }
write_tsv(Summary_statistics, Output)
} else {
 Summary_statistics = read_tsv(Output)
}

Summary_statistics %>% filter(! Feature == "(Intercept)") %>% mutate(FDR_all = p.adjust(`Pr(>|t|)`, "fdr")) -> Summary_statistics
Summary_statistics %>% arrange(`Pr(>|t|)`) -> Summary_statistics

Info = read_tsv("metadata/link_SGB_species_names.txt")
Summary_statistics %>% left_join(Info, by="SGB") -> Summary_statistics

write_tsv(Summary_statistics, Output) 

print(Summary_statistics)



