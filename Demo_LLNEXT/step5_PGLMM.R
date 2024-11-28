############################ PGLMM ANALYSIS  ###########
# Authors: Alex Kurilshikov
# Last update: 26st of Nov, 2024 
# PGLMM analysis of association of phenotypes to strain phylogeny

library(ape)
library(phangorn)
library(MCMCglmm)

## Functions 

run_MCMC_singlePair = function(tree,covars,data,s_nitt = 60000, s_burnin = 30000,s_thin = 50) {
  invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))} 
  pheno_name = colnames(data)
  dataset = data.frame(covars,phenotype=data[,1])
  # filtering
  dataset = dataset[!is.na(dataset$phenotype),,drop = F]
  name_tree = names(tree)
  #tree = tree[[1]]
  tree_filtered = keep.tip(tree, intersect(tree$tip.label,dataset$NG_ID))
  dataset = dataset[match(tree_filtered$tip.label,dataset$NG_ID),,drop =F]
  dataset$Timepoint.fac = as.integer(factor(dataset$Timepoint,levels = c("W2","M1","M3","M6","M12")))
  dataset = do.call(rbind,lapply(unique(dataset$Family),\(x) {
    subset = dataset[dataset$Family == x,]
    subset[which.min(subset$Timepoint.fac),]
  }))
  tree_filtered = try(keep.tip(tree_filtered,dataset$NG_ID))
  
  
  if(class(tree_filtered)[1] == "phylo") {
    if (length(table(dataset$Family))<4) {
      NA
    } else {
      Ainv.1 = try(inverseA(tree_filtered,nodes = "TIPS",scale = T)$Ainv)
      if(class(Ainv.1)[1]=="try-error") Ainv.1 = try(inverseA(tree_filtered,nodes = "TIPS",scale = T)$Ainv)
      if(class(Ainv.1)[1]=="try-error") Ainv.1 = inverseA(force.ultrametric(remove.zero.brlen(tree_filtered),nodes = "TIPS",scale = T))$Ainv
      
      family.mcmc = ifelse(class(dataset$phenotype)=="numeric"|class(dataset$phenotype)=="integer",
                           "gaussian",
                           "categorical")
      if(family.mcmc=="gaussian") {
        dataset$phenotype = invrank(dataset$phenotype)
        p.var = var(dataset$phenotype,na.rm = T)
        prior1 <- list(
          G=list(
            G1=list(V=1,nu=0.002)),
          R=list(V=1,nu=0.002))
        Nlevels <- 2
        
      } else {dataset$phenotype = factor(dataset$phenotype)} 
      
      if(family.mcmc=="categorical") {
        Nlevels = length(table(dataset$phenotype))
        I <- diag(Nlevels-1)
        J <- matrix(rep(1, (Nlevels-1)^2), c(Nlevels-1, Nlevels-1))
        RV = (I+J)/Nlevels
        prior1 <- list(G=list(
          G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)),
          R=list(V=RV,fix=1))
      }
      
      if(family.mcmc=="categorical" & Nlevels >2) {
        modelMCMC <- try(MCMCglmm(phenotype ~ Covariate1 + Covariate1,
                                  random =~ NG_ID,
                                  ginverse = list(NG_ID = Ainv.1),
                                  rcov = ~idh(trait):units,
                                  data = dataset,
                                  prior = prior1,
                                  thin=s_thin,
                                  family = family.mcmc,nitt = s_nitt,burnin=s_burnin))
        modelMCMC.time <- try(MCMCglmm(phenotype ~ Covariate1 + Covariate1 + Timepoint,
                                       random =~ NG_ID,
                                       ginverse = list(NG_ID = Ainv.1),
                                       rcov = ~idh(trait):units,
                                       data = dataset,
                                       prior = prior1,
                                       thin=s_thin,
                                       family = family.mcmc,nitt = s_nitt,burnin=s_burnin))
      } else {
        modelMCMC <- try(MCMCglmm(phenotype ~ Covariate1 + Covariate1,
                                  random =~ NG_ID,
                                  ginverse = list(NG_ID = Ainv.1),
                                  prior = prior1,
                                  data = dataset,
                                  thin=s_thin,
                                  family = family.mcmc,nitt = s_nitt,burnin=s_burnin))
        modelMCMC.time <- try(MCMCglmm(phenotype ~ Covariate1 + Covariate1 + Timepoint,
                                       random =~ NG_ID,
                                       ginverse = list(NG_ID = Ainv.1),
                                       prior = prior1,
                                       data = dataset,
                                       thin=s_thin,
                                       family = family.mcmc,nitt = s_nitt,burnin=s_burnin))
      }
      #modelMCMC$VCV = modelMCMC$VCV / rowSums(modelMCMC$VCV)
      list(MCMC=modelMCMC,
           MCMC.time = modelMCMC.time,
           N = nrow(dataset),
           Nbabies = length(unique(dataset$Family)),
           pheno_name = pheno_name,
           bacteria = name_tree)
    }
  } else {
    list(MCMC=NA,
         MCMC.time = NA,
         N = 0,
         Nbabies = 0,
         pheno_name = pheno_name,
         bacteria = name_tree)
  }
}

parse_mcmc_single = function(mcmc_object){
  #load(mcmc_name)
  #mcmc_object = result
  if(length(mcmc_object)<6) {
    mcmc_object = list(NA,NA,NA,NA,NA,NA)
    class(mcmc_object[[1]]) = "try-error"
    class(mcmc_object[[2]]) = "try-error"
    
  }
  if (
    class(mcmc_object[[1]])[1] != "MCMCglmm") {
    outline = data.frame(
      family = NA,
      #bacterium = mcmc_object[[6]],
      phenotype = mcmc_object[[5]],
      N=mcmc_object[[3]],
      Nbaby = mcmc_object[[4]],
      N.iter = NA, 
      time.N.iter = NA,   
      R2.random.lower = NA,
      R2.random.median = NA,
      R2.random.mode = NA,
      R2.random.higher = NA,
      R2.total.lower = NA,
      R2.total.median = NA,
      R2.total.mode = NA,
      signal2noise.mode = NA,
      time.R2.random.lower = NA,
      time.R2.random.median = NA,
      time.R2.random.mode = NA,
      time.R2.random.higher = NA,
      time.R2.total.lower = NA,
      time.R2.total.median = NA,
      time.R2.total.mode = NA,
      time.signal2noise.mode = NA
    )
  }else {
    family = mcmc_object[[1]]$family[1]
    if (family =="gaussian") {
      VCV1 = mcmc_object[[1]]$VCV / rowSums(mcmc_object[[1]]$VCV)
      VCV2 = try(mcmc_object[[2]]$VCV / rowSums(mcmc_object[[2]]$VCV))
      VCV1.total = mcmc_object[[1]]$VCV / (
        apply(mcmc_object[[1]]$Sol %*% t(as.matrix(mcmc_object[[1]]$X)), 1,var) +
          rowSums(mcmc_object[[1]]$VCV)
      )
      VCV2.total = try(mcmc_object[[2]]$VCV / (
        apply(mcmc_object[[2]]$Sol %*% t(as.matrix(mcmc_object[[2]]$X)), 1,var) +
          rowSums(mcmc_object[[2]]$VCV)))
    } else {
      VCV1 = mcmc_object[[1]]$VCV / (rowSums(mcmc_object[[1]]$VCV) + (pi^2)/3)
      VCV2 = try(mcmc_object[[2]]$VCV / (rowSums(mcmc_object[[2]]$VCV) + (pi^2)/3))
      VCV1.total = mcmc_object[[1]]$VCV / (
        apply(mcmc_object[[1]]$Sol %*% t(mcmc_object[[1]]$X), 1,var) +
          rowSums(mcmc_object[[1]]$VCV) + (pi^2)/3
      )
      VCV2.total = try(mcmc_object[[2]]$VCV / (
        apply(mcmc_object[[2]]$Sol %*% t(mcmc_object[[2]]$X), 1,var) +
          rowSums(mcmc_object[[2]]$VCV)+ (pi^2)/3
      )
      )
    }
    outline = data.frame(family = family,
                         #bacterium = mcmc_object[[6]],
                         phenotype = mcmc_object[[5]],
                         N=mcmc_object[[3]],
                         Nbaby = mcmc_object[[4]],
                         N.iter = effectiveSize(VCV1)["NG_ID"], 
                         time.N.iter = ifelse(class(VCV2)[1] == "try-error",NA,effectiveSize(VCV2)["NG_ID"]),  
                         R2.random.lower = HPDinterval(VCV1)[1,1],
                         R2.random.median = median(VCV1[,"NG_ID"]),
                         R2.random.mode = posterior.mode(VCV1)["NG_ID"],
                         R2.random.higher = HPDinterval(VCV1)[1,2],
                         R2.total.lower = HPDinterval(VCV1.total)[1,1],
                         R2.total.median = median(VCV1.total[,"NG_ID"]),
                         R2.total.mode = posterior.mode(VCV1.total)["NG_ID"],
                         signal2noise.mode = posterior.mode(mcmc_object[[1]]$VCV[,1] / rowSums(mcmc_object[[1]]$VCV[,2,drop = F])),
                         time.R2.random.lower = ifelse(class(VCV2)[1] == "try-error",NA,HPDinterval(VCV2)[1,1]),
                         time.R2.random.median = ifelse(class(VCV2)[1] == "try-error",NA,median(VCV2[,"NG_ID"])),
                         time.R2.random.mode = ifelse(class(VCV2)[1] == "try-error",NA,posterior.mode(VCV2.total)["NG_ID"]),
                         time.R2.random.higher = ifelse(class(VCV2)[1] == "try-error",NA,HPDinterval(VCV2.total)[1,2]),
                         time.R2.total.lower = ifelse(class(VCV2.total)[1] == "try-error",NA,HPDinterval(VCV2.total)[1,1]),
                         time.R2.total.median = ifelse(class(VCV2.total)[1] == "try-error",NA,median(VCV2.total[,"NG_ID"])),
                         time.R2.total.mode = ifelse(class(VCV2.total)[1] == "try-error",NA,posterior.mode(VCV2.total)["NG_ID"]),
                         time.signal2noise.mode =  ifelse(class(VCV2)[1] == "try-error",NA,posterior.mode(mcmc_object[[2]]$VCV[,1] / rowSums(mcmc_object[[2]]$VCV[,2,drop = F])))
    )
  }
  outline  
}


## Data loading

metadata = read.table("simulated_metadata.txt",header=T)
align = read.dna("alignment_example.aln",format = "fasta")

#generating ultrametric GTR tree
tree = pml_bb(align,model = "GTR",method = "ultrametric",
              rearrangement="NNI")$tree


# preparing dataset for loading. 
metadata_infant = metadata[metadata$Type=="Infant",]
metadata_infant$NG_ID = rownames(metadata_infant)
covars = metadata_infant[,c("Covariate1","Covariate2","NG_ID","Timepoint","Family")]

## running PGLMM analysis 
results_mcmc = lapply(c("Binary_phenotype_dynamic","Quant_phenotype_dynamic","Factor_phenotype_dynamic"),
       \(x) run_MCMC_singlePair(tree,
                    covars,
                    metadata_infant[,x,drop = F]))

## parsing PGLMM results
do.call(rbind,lapply(results_mcmc, parse_mcmc_single))
