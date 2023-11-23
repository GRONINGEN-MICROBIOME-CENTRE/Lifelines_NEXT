########### Infant eczema prediction using CODACORE ##############################
# Author: Trishla Sinha 
# Date: 23rd October 
# Last update: 23rd  November

setwd("~/Desktop/LLNEXT/Analysis/prediction/")

####Load libraries####

library(tidyverse)
library(vegan)
library(ggrepel)
library(tidyverse)
library(codacore)
library(reticulate)
use_condaenv("for_coda_core", required = TRUE)
library(tensorflow)


Pseudocount = function(Taxa){
  return(min(Taxa[!Taxa==0])/2)
}

ggplot <- function(...) ggplot2::ggplot(...) + scale_color_manual(values=met.brewer("Manet",5)[c(1,5, 3)] ) +
  scale_fill_manual(values=met.brewer("Manet", 5)[c(1,5)] )

metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_03_08_2023.txt")
taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
phenos <-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_cross_sectional_2023_11_15.txt")
eczema<-phenos %>%
  select(next_id_mother, infant_health_eczema_diagnosis_strict)
names (eczema)[1]<-"NEXT_ID"
eczema<-eczema[!duplicated(eczema$NEXT_ID), ]
metadata<-left_join(metadata, eczema)
row.names(metadata)<-metadata$NG_ID


# Filtering of Mother taxa to most abundant and prevalent 
row.names(metadata)<-metadata$NG_ID
mother_metadata<-metadata[metadata$Type=="mother",]
mother_taxa<-taxa[row.names(taxa)%in% rownames(mother_metadata),] 
mother_taxa<-mother_taxa[match(row.names(mother_taxa),row.names(mother_metadata)),]
mother_NEXT_ID<-mother_metadata %>%
  select(NEXT_ID)
mother_taxa_all<-merge(mother_NEXT_ID,mother_taxa, by="row.names" )
row.names(mother_taxa_all)<-mother_taxa_all$Row.names
mother_taxa_all$Row.names=NULL
unique_counts <- sapply(mother_taxa_all, function(x) length(unique(mother_taxa_all$NEXT_ID[x >0.001]))) # Here I am counting the unique elements in the NEXT_ID column where the corresponding value in each column (i.e., x) is greater than the given cut-off. 
mother_taxa_all_filt <- mother_taxa_all[, unique_counts >= 0.3*length(unique(mother_taxa_all$NEXT_ID)) ] # Setting a cut-off on prevalence 
mother_taxa_all_filt$NEXT_ID=NULL
mother_taxa_all_filt$UNCLASSIFIED=NULL

# Only SGB level 
mother_taxa_SGB_filt=mother_taxa_all_filt[,grep("t__",colnames(mother_taxa_all_filt))]

#Simplify names
colnames(mother_taxa_SGB_filt)=gsub(".*s__","",colnames(mother_taxa_SGB_filt))

# Only birth samples 
mother_metadata_B<-mother_metadata[mother_metadata$Timepoint_categorical=="B",]
mother_metadata_B<- mother_metadata_B%>% drop_na(infant_health_eczema_diagnosis_strict)
mother_taxa_SGB_filt_B<-mother_taxa_SGB_filt[row.names(mother_taxa_SGB_filt)%in% rownames(mother_metadata_B),] 
mother_taxa_SGB_filt_B<-mother_taxa_SGB_filt_B[match(row.names(mother_taxa_SGB_filt_B),row.names(mother_metadata_B)),]


# Merge the taxa data with the metadata
eczema<-mother_metadata_B %>%
  select(infant_health_eczema_diagnosis_strict)
final_data <- merge(eczema, mother_taxa_SGB_filt_B, by = "row.names")
rownames(final_data) <- final_data$Row.names
final_data$Row.names <- NULL

# Prepare the data for codacore
Y <- final_data$infant_health_eczema_diagnosis_strict
X_old <- as.matrix(final_data[, colnames(final_data) != "infant_health_eczema_diagnosis_strict"])

# Add a pseudocount 
X <- X_old + Pseudocount(mother_taxa_SGB_filt_B) 

# Split data into training and test sets
set.seed(123)  # for reproducibility
train_indices <- sample(1:nrow(X), size = floor(0.8 * nrow(X)))
X_train <- X[train_indices, ]
Y_train <- Y[train_indices]
X_test <- X[-train_indices, ]
Y_test <- Y[-train_indices]

# Fit the codacore model
optimizer = tensorflow::tf$keras$optimizers$legacy$SGD()
tflow <- tensorflow::tf$compat$v1
tflow$enable_eager_execution()
print(tflow$executing_eagerly())
model <- codacore(
  X_train,
  Y_train,
  logRatioType = 'balances', # can also use 'amalgamations'
  lambda = 1                 # regularization parameter (1 corresponds to "1SE rule")
)

print(model)
plot(model)
plotROC(model)

N1 = colnames(X_train)[getNumeratorParts(model, 1)]
N2 = colnames(X_train)[getNumeratorParts(model, 2)]
D1 = colnames(X_train)[getDenominatorParts(model, 1)]
D2 = colnames(X_train)[getDenominatorParts(model, 2)]

final_data %>% select(-infant_health_eczema_diagnosis_strict) %>% as_tibble() %>% dplyr::select(N1) %>% apply(1, function(x){ x= x+ Pseudocount(mother_taxa_SGB_filt_B)  ;exp(mean(log(x)))  }  ) -> Numerator1
final_data %>% select(-infant_health_eczema_diagnosis_strict)  %>% as_tibble() %>% dplyr::select(N2) %>% apply(1, function(x){ x= x+ Pseudocount(mother_taxa_SGB_filt_B)  ;exp(mean(log(x)))  }  ) -> Numerator2
final_data %>% select(-infant_health_eczema_diagnosis_strict)  %>% as_tibble() %>% dplyr::select(D1) %>% apply(1, function(x){ x= x+ Pseudocount(mother_taxa_SGB_filt_B)  ;exp(mean(log(x)))  }  ) -> Denominator1
final_data %>% select(-infant_health_eczema_diagnosis_strict) %>% as_tibble() %>% dplyr::select(D2) %>% apply(1, function(x){ x= x+ Pseudocount(mother_taxa_SGB_filt_B)  ;exp(mean(log(x)))  }  ) -> Denominator2
final_data %>% select(-infant_health_eczema_diagnosis_strict)  %>% as.data.frame() %>% rownames_to_column("ID")  %>% as_tibble() %>% mutate(Ratio1 = log10(Numerator1/Denominator1), Ratio2=log10(Numerator2/Denominator2)  ) %>% select(ID ,Ratio1, Ratio2) -> for_logistic
mother_metadata_B %>% rownames_to_column("ID") %>% select(ID, shannon ) %>% as_tibble() %>% left_join(for_logistic, . ) %>%
  mutate(Y = as.factor(final_data$infant_health_eczema_diagnosis_strict) ) -> For_predicition

For_predicition %>% filter(ID %in% rownames(X_train) ) -> Train_for_model
For_predicition %>% filter(! ID %in% rownames(X_train) ) -> Test_for_model


Train_for_model %>% glm(Y ~ shannon + Ratio1 + Ratio2, . , family="binomial" ) -> model1
Train_for_model %>% glm(Y ~ shannon , . , family="binomial" ) -> model2
Train_for_model %>% glm(Y ~ Ratio1 + Ratio2 , . , family="binomial" ) -> model3

predicted_probs <- predict(model1, newdata = Test_for_model, type = "response")
roc_curve <- roc(Test_for_model$Y, predicted_probs)
auc_value <- auc(roc_curve)

predicted_probs <- predict(model2, newdata = Test_for_model, type = "response")
roc_curve <- roc(Test_for_model$Y, predicted_probs)
auc_value <- auc(roc_curve)

predicted_probs <- predict(model3, newdata = Test_for_model, type = "response")
roc_curve <- roc(Test_for_model$Y, predicted_probs)
auc_value <- auc(roc_curve)



# Predict on the test set
predictions <- predict(model, newx = X_test, logits=F)

# Evaluate the model performance
yHat <- predict(model,X_test, logits = F)
cat("Test set AUC =", pROC::auc(pROC::roc(Y_test,
                                          yHat, quiet = T)))

yHat <- predict(model, X_test, logits = F,
                numLogRatios = 1)
cat("Test set AUC =", pROC::auc(pROC::roc(Y_test,
                                          yHat, quiet = T)))

