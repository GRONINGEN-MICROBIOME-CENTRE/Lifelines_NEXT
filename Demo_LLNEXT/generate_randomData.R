################### Generate random data for testing purpose ##########################
# Author: A. Kurilshikov
# Last updated: 21th Nov, 2024

## Choosing simulation parameters:  

Nfamilies = 500
Nbac = 1000
set.seed(1234)

Nsamples.babies = sapply(1:Nfamilies,\(x) sample(1:5,1))
Nsamples.mothers = sapply(1:Nfamilies,\(x) sample(1:3,1))
sample_babies = sapply(Nsamples.babies,\(x) sample(c("W2","M1","M3","M6","M12"),x))
sample_mothers = sapply(Nsamples.mothers,\(x) sample(c("B","P12","P28"),x))

metadata_Rand = rbind(
  data.frame(Type = "Infant",do.call(rbind,
                              lapply(1:length(sample_babies),
                                     \(x) data.frame(
                                       Family = paste0("FAM",sprintf("%03d",x)),
                                       Timepoint = sample_babies[[x]]
                                       )
                                     )
                              )
             ),
  data.frame(Type = "Mother",do.call(rbind,
                              lapply(1:length(sample_babies),
                                     \(x) data.frame(
                                       Family = paste0("FAM",sprintf("%03d",x)),
                                       Timepoint = sample_mothers[[x]]
                                     )
                              )
  )
  )
)
rownames(metadata_Rand) = paste0("Sample",sprintf("%04d",1:nrow(metadata_Rand)))


#Simulating metadata
metadata_Rand = data.frame(metadata_Rand,
                           Covariate1= sample(c("Covar1.1","Covar1.2","Covar1.3"),
                                              nrow(metadata_Rand),replace = T),
                           Covariate2= rnorm(nrow(metadata_Rand)),
                           Binary_phenotype_dynamic = sample(c("BinaryDynamic1","BinaryDynamic2"),
                                                             nrow(metadata_Rand),replace = T),
                           Factor_phenotype_dynamic = sample(c("LevelDynamic1","LevelDynamic2","LevelDynamic3"),
                                                             nrow(metadata_Rand),replace = T),
                           Quant_phenotype_dynamic = rnorm(nrow(metadata_Rand)),
                           Binary_phenotype_static = NA,
                           Factor_phenotype_static = NA,
                           Quant_phenotype_static = NA
)

## generatic static phenotypes
static_phenos = data.frame(Family = unique(metadata_Rand$Family),
                           Binary_phenotype_static = sample(c("BinaryStatic1","BinaryStatic2"),
                                                            length(unique(metadata_Rand$Family)),replace = T),
                           Factor_phenotype_static = sample(c("LevelStatic1","LevelStatic2","LevelStatic3"),
                                                            length(unique(metadata_Rand$Family)),replace = T),
                           Quant_phenotype_static = rnorm(length(unique(metadata_Rand$Family))))
metadata_Rand$Binary_phenotype_static = static_phenos[match(metadata_Rand$Family,static_phenos$Family),"Binary_phenotype_static"]
metadata_Rand$Factor_phenotype_static = static_phenos[match(metadata_Rand$Family,static_phenos$Family),"Factor_phenotype_static"]
metadata_Rand$Quant_phenotype_static = static_phenos[match(metadata_Rand$Family,static_phenos$Family),"Quant_phenotype_static"]



#Simulating metaphlan
if(!require("VGAM")) {
  install.packages("VGAM")
  library("VGAM")}
  
metaphlan_Rand = do.call(cbind,
                         lapply(1:Nbac,
                                \(x) rzipois(nrow(metadata_Rand), 
                                             lambda = rnorm(1,runif(1)* 100,7), 
                                             pstr0 = runif(1))))

rownames(metaphlan_Rand) = rownames(metadata_Rand)
colnames(metaphlan_Rand) = paste0("Bacterium",sprintf("%04d",1:Nbac))
metaphlan_Rand = metaphlan_Rand[,!is.na(colSums(metaphlan_Rand))]
metaphlan_Rand = sweep(metaphlan_Rand,1,rowSums(metaphlan_Rand),"/")

# Saving results
write.table(metadata_Rand,file = "./simulated_metadata.txt",sep="\t")                           
write.table(metaphlan_Rand,file = "./simulated_metaphlan.txt",sep="\t")                           
                           
                           
                           
                           
                           
                           
                           
                           