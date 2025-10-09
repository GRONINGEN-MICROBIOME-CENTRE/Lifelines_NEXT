##Coda4Microbiome helper functions

##########################
###Longitudinal Functions#
##########################
#coda_glmnet_longitudinal2 : this one is the one usually caled
#coda_glmnet_longitudinal0_2 : this one is a short version that is usually called in the null distribtuion function
#coda_glmnet_longitudinal_null_2


coda_glmnet_longitudinal2 = function (x, y, x_time, subject_id, ini_time, end_time, covar = NULL, lambda = "lambda.1se", nvar = NULL, alpha = 0.9, nfolds = 5, Do_Leave_One_Out = F, showPlots = TRUE, coef_threshold = 0) {
  #Includes a few changes from the original version
  #1. Annotations
  #2. Leave-one-out CV scheme for AUC assessment
  #3. Produces folds for glmnet that keep proportions of y levels, instead of using the random split default in cv.glmnet
  #4. Returns the glmnet model, and the LOO results on top of the previous ones
  
  yini <- y
  #0 imputation - half of the minimum non-0 value
  if (sum(x == 0) > 0) {
    x <- impute_zeros(x)
  }
  #log of the resulting pseudocounted counts
  logX1 = log(x)
  nsubjects = length(unique(subject_id))
  indexUser = seq_along(subject_id)[!duplicated(subject_id)]
  y_unique <- y[indexUser]
  if (!is.null(covar)) {
    covar = covar[indexUser, ]
  }
  y.binary <- ifelse(dim(table(y_unique)) == 2, TRUE, FALSE)
  if (y.binary == TRUE) {
    y <- factor(y)
    y_unique <- as.factor(y_unique)
  }
  if (is.factor(y)) {
    labelsy <- levels(y)
    y_unique <- factor(y_unique, labels = labelsy)
  }
  alpha0 <- alpha
  kselect <- ncol(x)
  taxaselect <- (1:ncol(x))
  k <- ncol(x)
  intLogX <- NULL
  #Calculate inegrals per each of the taxa through time
  #Output: matrix of N_samples x N_species, but the value of each species in the inegral
  for (ki in (1:(ncol(logX1)))) {
    #Iterate through bugs
    print(paste("Taxa=", ki))
    yy = as.numeric(logX1[, ki])
    #We calculate the integral per each subject (no sample!); so: only one value per individual
    integrals = coda4microbiome:::integralFun(x_time, yy, subject_id, a = ini_time, 
                                              b = end_time)
    intLogX <- cbind(intLogX, matrix(integrals))
  }
  print(dim(intLogX))
  m <- length(y_unique)
  lrcolnames <- NULL
  lrX <- matrix(0, m, k * (k - 1)/2)
  idlrX <- matrix(0, k * (k - 1)/2, 2)
  nameslrX <- matrix(0, k * (k - 1)/2, 2)
  colnamesx <- colnames(x)
  lloc <- 0
  #Per each of the integrals, get the difference between each pair of taxa. Differences of integrals of logs represent ratios
  #Final size N_samples x (N_taxa*(N_taxa-1)/2)
  for (i in (1:(k - 1))) {
    for (j in ((i + 1):k)) {
      lloc = lloc + 1
      idlrX[lloc, ] <- c(i, j)
      nameslrX[lloc, ] <- c(colnamesx[i], colnamesx[j])
      lrX[, lloc] <- intLogX[, i] - intLogX[, j]
      lrcolnames <- c(lrcolnames, paste(paste("lr", i, 
                                              sep = ""), j, sep = "."))
    }
  }
  colnames(lrX) <- lrcolnames
  idlrXsub <- idlrX
  lrXsub <- lrX
  y.binary <- ifelse(dim(table(y)) == 2, TRUE, FALSE)
  #And now, run glmnet!
  if (y.binary == TRUE) {
    Folds = generate_stratified_folds_balanced(y_unique, nfolds)
    summary_table <- table(Folds, y_unique)
    print(summary_table)
    #For binary data
    if (is.null(covar)) {
      lassocv <- cv.glmnet2(lrXsub, y_unique, family = "binomial", foldid=Folds,
                                   alpha = alpha0, type.measure = "auc", nfolds = nfolds, 
                                   keep = TRUE)
    }
    else {
      df0 <- data.frame(y_unique, covar)
      model0 <- glm(y_unique ~ ., family = "binomial", 
                    data = df0)
      x0 <- predict(model0)
      lassocv <- cv.glmnet2(lrXsub, y_unique, family = "binomial", 
                                   offset = x0, alpha = alpha0, type.measure = "auc", 
                                   nfolds = nfolds, keep = TRUE)
    }
    if (Do_Leave_One_Out == T){
      probs = LOO_glment(lrXsub, y_unique, alpha0, nfolds, lambda, foldlist=Folds)
      Results_loo = Get_scores_model(probs, truth = y_unique )
    } else { Results_loo = list() }
    
  } else {  #For continuous data
    #No need for fold stratification keeping proportions
    if (is.null(covar)) {
      lassocv <- cv.glmnet2(lrXsub, y_unique, alpha = alpha0, 
                                   type.measure = "deviance", nfolds = nfolds, 
                                   keep = TRUE)
    }
    else {
      df0 <- data.frame(y_unique, covar)
      model0 <- lm(y_unique ~ ., data = df0)
      x0 <- predict(model0)
      lassocv <- cv.glmnet2(lrXsub, y_unique, offset = x0, 
                                   alpha = alpha0, type.measure = "deviance", nfolds = nfolds, 
                                   keep = TRUE)
    }
    #No implementation of LOO-CV for continuous data atm (tbd)
  }
  if (showPlots == TRUE) {
    plot(lassocv)
  }
  if (!is.null(nvar)) {
    rowlasso <- max(which(lassocv$glmnet.fit$df <= nvar))
    lambda <- lassocv$glmnet.fit$lambda[rowlasso]
  }
  lambdavalue <- lambda
  if (is.character(lambda)) {
    if (lambda == "lambda.1se") 
      lambdavalue <- lassocv$lambda.1se
    if (lambda == "lambda.min") 
      lambdavalue <- lassocv$lambda.min
  }
  idrow <- max(which(lassocv$glmnet.fit$lambda >= lambdavalue))
  coeflr <- as.vector(coef(lassocv, s = lambda))[-1]
  lrselect <- which(coeflr != 0)
  idlrXsub[lrselect, ]
  coeflogcontrast <- rep(0, ncol(x))
  for (i in (1:length(coeflr))) {
    coeflogcontrast[idlrXsub[i, 1]] <- coeflogcontrast[idlrXsub[i, 1]] + coeflr[i]
    coeflogcontrast[idlrXsub[i, 2]] <- coeflogcontrast[idlrXsub[i, 2]] - coeflr[i]
  }
  varlogcontrast <- which(abs(coeflogcontrast) > coef_threshold)
  coeflogcontrast <- coeflogcontrast[varlogcontrast]
  (names.select <- colnames(x)[varlogcontrast])
  (positive <- ifelse(coeflogcontrast > 0, 1, 0))
  positive <- factor(positive, levels = c(0, 1), labels = c("negative",  "positive"))
  logcontrast = as.matrix(lrXsub[, lrselect]) %*% coeflr[lrselect]
  if (is.null(covar)) {
    predictions <- as.numeric(predict(lassocv, lrX, s = lambdavalue))
  }
  else {
    predictions <- as.numeric(predict(lassocv, lrX, s = lambdavalue, 
                                      newoffset = x0))
  }
  coeflogcontrast <- 2 * coeflogcontrast/sum(abs(coeflogcontrast))
  if (y.binary == TRUE) {
    AUC_signature <- pROC::auc(pROC::roc(y_unique, as.numeric(predictions), quiet = TRUE))[[1]]
    if (length(varlogcontrast) == 0) 
      AUC_signature <- 0.5
    mcvAUC <- lassocv$cvm[idrow]
    sdcvAUC <- lassocv$cvsd[idrow]
  }
  else {
    mcvDev <- lassocv$cvm[idrow]
    sdcvDev <- lassocv$cvsd[idrow]
    Rsq <- 0
    if (length(varlogcontrast) > 0) {
      Rsq <- as.numeric(cor(predictions, y)^2)
    }
  }
  y <- y_unique
  plot1 <- NULL
  plot2 <- NULL
  if (length(lrselect > 0)) {
    plot1 <- plot_prediction(predictions, y, showPlots = showPlots)
    plot2 <- plot_signature(names.select, coeflogcontrast, showPlots = showPlots)
  }
  else {
    print("No variables are selected. The prediction and the signature plots are not displayed.")
  }
  plot3 <- NULL
  if (showPlots == TRUE) {
    plot3 <- plot_signature_curves(varlogcontrast, coeflogcontrast,  x = x, y = yini, x_time, subject_id, ini_time, end_time)
  }
  if (y.binary == TRUE) {
    results <- list(taxa.num = varlogcontrast, taxa.name = names.select, 
                    `log-contrast coefficients` = coeflogcontrast, predictions = predictions, 
                    `apparent AUC` = AUC_signature, `mean cv-AUC` = mcvAUC, 
                    `sd cv-AUC` = sdcvAUC, `predictions plot` = plot1, 
                    `signature plot` = plot2, `trajectories plot` = plot3, Model=lassocv, LOO_results=Results_loo, Y = y_unique  , Y_hat = probs  )
  }
  else {
    results <- list(taxa.num = varlogcontrast, taxa.name = names.select, 
                    `log-contrast coefficients` = coeflogcontrast, predictions = predictions, 
                    `apparent Rsq` = Rsq, `mean cv-Deviance` = mcvDev, 
                    `sd cv-Deviance` = sdcvDev, `predictions plot` = plot1, 
                    `signature plot` = plot2, `trajectories plot` = plot3 )
  }
  return(results)
}

coda_glmnet_longitudinal0_2 = function (x, lrX, idlrX, nameslrX, y, x_time, subject_id, ini_time, end_time, Folds, covar = NULL, ktop = NULL, lambda = "lambda.1se", alpha = 0.9, nfolds = 10,Do_Leave_One_Out = F) {
 #Shorter version of  coda_glmnet_longitudinal in the coda4microbiome package. It directly accept the table with the difference of the integrals of each pair of taxa
 #I have similarly included a LOO-CV statment, and the definition of folds within the function
  y.binary <- ifelse(dim(table(y)) == 2, TRUE, FALSE)
  if (sum(x == 0) > 0) {
    x <- impute_zeros(x)
  }
  logX1 = log(x)
  nsubjects = length(unique(subject_id))
  indexUser = seq_along(subject_id)[!duplicated(subject_id)]
  if (!is.null(covar)) {
    covar = covar[indexUser, ]
  }
  y_unique <- y[indexUser]
  if (y.binary == TRUE) {
    y <- factor(y)
    y_unique <- as.factor(y_unique)
  }
  if (is.factor(y)) {
    labelsy <- levels(y)
    y_unique <- factor(y_unique, labels = labelsy)
  }
  alpha0 <- alpha
  kselect <- ncol(x)
  idlrXsub <- idlrX
  lrXsub <- lrX
  if (y.binary == TRUE) {
    if (is.null(covar)) {
      lassocv <- cv.glmnet2(lrXsub, y_unique, family = "binomial", 
                                   alpha = alpha0, type.measure = "auc", nfolds = nfolds, foldid = Folds,
                                   keep = TRUE)
    }
    else {
      df0 <- data.frame(y_unique, covar)
      model0 <- glm(y_unique ~ ., family = "binomial", 
                    data = df0)
      x0 <- predict(model0)
      lassocv <- cv.glmnet2(lrXsub, y_unique, family = "binomial", foldid = Folds,
                                   offset = x0, alpha = alpha0, type.measure = "auc", 
                                   nfolds = nfolds, keep = TRUE)
    }
    if (Do_Leave_One_Out == T){
      probs = LOO_glment(lrXsub, y_unique, alpha0, nfolds, lambda, foldlist=Folds)
      Results_loo = Get_scores_model(probs, truth = y_unique )
    } else { Results_loo = list() }
  }
  else {
    if (is.null(covar)) {
      lassocv <- cv.glmnet2(lrXsub, y_unique, alpha = alpha0, 
                                   type.measure = "deviance", nfolds = nfolds, keep = TRUE)
    }
    else {
      df0 <- data.frame(y_unique, covar)
      model0 <- lm(y_unique ~ ., data = df0)
      x0 <- predict(model0)
      lassocv <- cv.glmnet2(lrXsub, y_unique, offset = x0, 
                                   alpha = alpha0, type.measure = "deviance", nfolds = nfolds, 
                                   keep = TRUE)
    }
  }
  lambdavalue <- lambda
  if (is.character(lambda)) {
    if (lambda == "lambda.1se") 
      lambdavalue <- lassocv$lambda.1se
    if (lambda == "lambda.min") 
      lambdavalue <- lassocv$lambda.min
  }
  idrow <- max(which(lassocv$glmnet.fit$lambda >= lambdavalue))
  coeflr <- as.vector(coef(lassocv, s = lambda))[-1]
  lrselect <- which(coeflr != 0)
  idlrXsub[lrselect, ]
  coeflogcontrast <- rep(0, ncol(x))
  for (i in (1:length(coeflr))) {
    coeflogcontrast[idlrXsub[i, 1]] <- coeflogcontrast[idlrXsub[i, 
                                                                1]] + coeflr[i]
    coeflogcontrast[idlrXsub[i, 2]] <- coeflogcontrast[idlrXsub[i, 
                                                                2]] - coeflr[i]
  }
  varlogcontrast <- which(abs(coeflogcontrast) > 0)
  coeflogcontrast <- coeflogcontrast[varlogcontrast]
  (names.select <- colnames(x)[varlogcontrast])
  (positive <- ifelse(coeflogcontrast > 0, 1, 0))
  positive <- factor(positive, levels = c(0, 1), labels = c("negative", 
                                                            "positive"))
  logcontrast = as.matrix(lrXsub[, lrselect]) %*% coeflr[lrselect]
  if (is.null(covar)) {
    predictions <- as.numeric(predict(lassocv, lrX, s = lambdavalue))
  }
  else {
    predictions <- as.numeric(predict(lassocv, lrX, s = lambdavalue, 
                                      newoffset = x0))
  }
  coeflogcontrast <- 2 * coeflogcontrast
  if (y.binary == TRUE) {
    AUC_signature <- pROC::auc(pROC::roc(y_unique, as.numeric(predictions), 
                                         quiet = TRUE))[[1]]
    if (length(varlogcontrast) == 0) 
      AUC_signature <- 0.5
    mcvAUC <- lassocv$cvm[idrow]
    sdcvAUC <- lassocv$cvsd[idrow]
  }
  else {
    mcvDev <- lassocv$cvm[idrow]
    sdcvDev <- lassocv$cvsd[idrow]
    Rsq <- 0
    if (length(varlogcontrast) > 0) {
      Rsq <- as.numeric(cor(predictions, y)^2)
    }
  }
  y <- y_unique
  if (y.binary == TRUE) {
    results <- list(taxa.num = varlogcontrast, taxa.name = names.select, 
                    `log-contrast coefficients` = coeflogcontrast, predictions = predictions, 
                    `apparent AUC` = AUC_signature, `mean cv-AUC` = mcvAUC, 
                    `sd cv-AUC` = sdcvAUC, LOO_results=Results_loo)
  }
  else {
    results <- list(taxa.num = varlogcontrast, taxa.name = names.select, 
                    `log-contrast coefficients` = coeflogcontrast, predictions = predictions, 
                    `apparent Rsq` = Rsq, `mean cv-Deviance` = mcvDev, 
                    `sd cv-Deviance` = sdcvDev)
  }
  return(results)
}  
 
coda_glmnet_longitudinal_null_2 = function (x, y, x_time, subject_id, ini_time, end_time, niter = 100, covar = NULL, alpha = 0.9, lambda = "lambda.1se", nfolds = 10, sig = 0.05, Do_Leave_One_Out = F) {
  #Only change, call coda_glmnet_longitudinal0_2, and extract LOO results if done
  #ADD: Now I also shuffle the X and not the Y, because I define the folds here
  y.binary <- ifelse(dim(table(y)) == 2, TRUE, FALSE)
  if (sum(x == 0) > 0) {
    x <- impute_zeros(x)
  }
  logX1 = log(x)
  nsubjects = length(unique(subject_id))
  indexUser = seq_along(subject_id)[!duplicated(subject_id)]
  if (!is.null(covar)) {
    covar = covar[indexUser, ]
  }
  y_unique <- y[indexUser]
  if (y.binary == TRUE) {
    y <- factor(y)
    y_unique <- as.factor(y_unique)
  }
  if (is.factor(y)) {
    labelsy <- levels(y)
    y_unique <- factor(y_unique, labels = labelsy)
  }
  k <- ncol(x)
  intLogX <- NULL
  for (ki in (1:(ncol(logX1)))) {
    print(paste("ind=", ki))
    yy = as.numeric(logX1[, ki])
    integrals = coda4microbiome:::integralFun(x_time, yy, subject_id, a = ini_time, 
                            b = end_time)
    intLogX <- cbind(intLogX, matrix(integrals))
  }
  print(dim(intLogX))
  m <- length(y_unique)
  lrcolnames <- NULL
  lrX <- matrix(0, m, k * (k - 1)/2)
  idlrX <- matrix(0, k * (k - 1)/2, 2)
  nameslrX <- matrix(0, k * (k - 1)/2, 2)
  colnamesx <- colnames(x)
  lloc <- 0
  for (i in (1:(k - 1))) {
    for (j in ((i + 1):k)) {
      lloc = lloc + 1
      idlrX[lloc, ] <- c(i, j)
      nameslrX[lloc, ] <- c(colnamesx[i], colnamesx[j])
      lrX[, lloc] <- intLogX[, i] - intLogX[, j]
      lrcolnames <- c(lrcolnames, paste(paste("lr", i, 
                                              sep = ""), j, sep = "."))
    }
  }
  colnames(lrX) <- lrcolnames
  y1 <- y
  accuracy <- rep(0, niter)
  accuracy_loo <- rep(0, niter)
  accuracy_balanced_loo <- rep(0, niter)

  fv = generate_stratified_folds_balanced(y_unique, nfolds)
  summary_table <- table(fv, y_unique)
  print(summary_table)
  lrX2 =lrX
  for (i in (1:niter)) {
    #y1 <- sample(y1)
    lrX2 = lrX2[sample(nrow(lrX2)), ]
    lr <- coda_glmnet_longitudinal0_2(x, lrX2, idlrX, nameslrX, 
                                    y = y1, x_time, subject_id, ini_time, end_time, alpha = alpha, 
                                    lambda = lambda, covar = covar, nfolds = nfolds,Do_Leave_One_Out =Do_Leave_One_Out, Folds = fv)
    if (y.binary == TRUE) {
      res <- lr$`mean cv-AUC`
    }
    else {
      res <- lr$`mean cv-MSE`
    }
    accuracy[i] <- res
    if (Do_Leave_One_Out == T){
      accuracy_loo[i] = lr$LOO_results[['AUC']]
    }
    print(paste("iter", i))
  }
  results <- list(accuracy = accuracy, `confidence interval` = quantile(accuracy,  c(sig/1, 1 - (sig/2))))
  return(results)
}


#############################
###Cross Sectional Functions#
#############################

#coda_glmnet : this one is the one usually caled
#coda_glmnet_longitudinal0_2 : this one is a short version that is usually called in the null distribtuion function
#coda_glmnet_longitudinal_null_2

coda_glmnet_2 = function(x, y, covar = NULL, lambda = "lambda.1se", nvar = NULL,  alpha = 0.9, nfolds = 10, showPlots = TRUE, coef_threshold = 0, Do_Leave_One_Out = F){
  x <- impute_zeros(x)
  kselect <- ncol(x)
  taxaselect <- (1:ncol(x))
  lrmatrix <- logratios_matrix(x)
  lrX <- lrmatrix[[1]]
  idlrX <- lrmatrix[[2]]
  nameslrX <- lrmatrix[[3]]
  y.binary <- ifelse(dim(table(y)) == 2, TRUE, FALSE)
  alpha0 <- alpha
  if (y.binary == TRUE) {
    #Added predefinition of folds
    Folds = generate_stratified_folds_balanced(y, nfolds)
    summary_table <- table(Folds, y)
    #

    if (is.null(covar)) {
      #added foldid argument
      lassocv <- cv.glmnet2(lrX, y, family = "binomial", foldid=Folds,
        alpha = alpha, type.measure = "auc", nfolds = nfolds, 
        keep = TRUE)
    }
    else {
      df0 <- data.frame(y, covar)
      model0 <- glm(y ~ ., family = "binomial", data = df0)
      x0 <- predict(model0)
      #added foldid argument
      lassocv <- cv.glmnet2(lrX, y, family = "binomial", foldid=Folds,
        offset = x0, alpha = alpha, type.measure = "auc", 
        nfolds = nfolds, keep = TRUE)
    }
    #Add LOO
    if (Do_Leave_One_Out == T){
      probs = LOO_glment(lrX, y, alpha0, nfolds, lambda, foldlist=Folds)
      Results_loo = Get_scores_model(probs, truth = y )
    } else { Results_loo = list() }

  }
  else {
    #Continuous data: No need for fold stratification keeping proportions
    #No implementation of LOO so far
    if (is.null(covar)) {
      lassocv <- glmnet::cv.glmnet(lrX, y, alpha = alpha, 
        type.measure = "deviance", nfolds = nfolds, 
        keep = TRUE)
    }
    else {
      df0 <- data.frame(y, covar)
      model0 <- lm(y ~ ., data = df0)
      x0 <- predict(model0)
      lassocv <- glmnet::cv.glmnet(lrX, y, offset = x0, 
        alpha = alpha, type.measure = "deviance", nfolds = nfolds, 
        keep = TRUE)
    }
  }
  if (showPlots == TRUE) {
    plot(lassocv)
  }
  if (!is.null(nvar)) {
    rowlasso <- max(which(lassocv$glmnet.fit$df <= nvar))
    lambda <- lassocv$glmnet.fit$lambda[rowlasso]
  }
  lambdavalue <- lambda
  if (is.character(lambda)) {
    if (lambda == "lambda.1se") 
      lambdavalue <- lassocv$lambda.1se
    if (lambda == "lambda.min") 
      lambdavalue <- lassocv$lambda.min
  }
  idrow <- max(which(lassocv$glmnet.fit$lambda >= lambdavalue))
  coeflr <- as.vector(coef(lassocv, s = lambda))[-1]
  lrselect <- which(coeflr != 0)
  coeflogcontrast <- rep(0, ncol(x))
  for (i in (1:length(coeflr))) {
    coeflogcontrast[idlrX[i, 1]] <- coeflogcontrast[idlrX[i, 
      1]] + coeflr[i]
    coeflogcontrast[idlrX[i, 2]] <- coeflogcontrast[idlrX[i, 
      2]] - coeflr[i]
  }
  varlogcontrast <- which(abs(coeflogcontrast) > coef_threshold)
  coeflogcontrast <- coeflogcontrast[varlogcontrast]
  (names.select <- colnames(x)[varlogcontrast])
  (positive <- ifelse(coeflogcontrast > 0, 1, 0))
  positive <- factor(positive, levels = c(0, 1), labels = c("negative", 
    "positive"))
  logcontrast = as.matrix(log(x)[, varlogcontrast]) %*% coeflogcontrast
  if (is.null(covar)) {
    predictions <- as.numeric(predict(lassocv, lrX, s = lambdavalue))
  }
  else {
    predictions <- as.numeric(predict(lassocv, lrX, s = lambdavalue, 
      newoffset = x0))
  }
  coeflogcontrast <- 2 * coeflogcontrast/sum(abs(coeflogcontrast))
  if (y.binary == TRUE) {
    AUC_signature <- pROC::auc(pROC::roc(y, as.numeric(predictions), 
      quiet = TRUE))[[1]]
    if (length(varlogcontrast) == 0) 
      AUC_signature <- 0.5
    mcvAUC <- lassocv$cvm[idrow]
    sdcvAUC <- lassocv$cvsd[idrow]
  }
  else {
    mcvMSE <- lassocv$cvm[idrow]
    sdcvMSE <- lassocv$cvsd[idrow]
    Rsq <- 0
    if (length(varlogcontrast) > 0) {
      Rsq <- as.numeric(cor(predictions, y)^2)
    }
  }
  plot1 <- NULL
  plot2 <- NULL
  if (length(lrselect > 0)) {
    plot1 <- plot_prediction(predictions, y, showPlots = showPlots)
    plot2 <- plot_signature(names.select, coeflogcontrast, 
      showPlots = showPlots)
  }
  else {
    print("No variables are selected. The prediction and the signature plots are not displayed.")
  }
  if (y.binary == TRUE) {
    #Add model and LOO result
    results <- list(taxa.num = varlogcontrast, taxa.name = names.select, 
      `log-contrast coefficients` = coeflogcontrast, predictions = predictions, 
      `apparent AUC` = AUC_signature, `mean cv-AUC` = mcvAUC, 
      `sd cv-AUC` = sdcvAUC, `predictions plot` = plot1, 
      `signature plot` = plot2, Model=lassocv, LOO_results=Results_loo,  Y = y  , Y_hat = probs )
  }
  else {
    results <- list(taxa.num = varlogcontrast, taxa.name = names.select, 
      `log-contrast coefficients` = coeflogcontrast, predictions = predictions, 
      `apparent Rsq` = Rsq, `mean cv-MSE` = mcvMSE, `sd cv-MSE` = sdcvMSE, 
      `predictions plot` = plot1, `signature plot` = plot2)
  }
  return(results)

}



coda_glmnet0_2 = function (x, lrX, idlrX, nameslrX, y, Folds, covar = NULL, lambda = "lambda.1se", alpha = 0.9, Do_Leave_One_Out = F) 
{
  if (sum(x == 0) > 0) {
    x <- impute_zeros(x)
  }
  kselect <- ncol(x)
  idlrXsub <- idlrX
  lrXsub <- lrX
  y.binary <- ifelse(dim(table(y)) == 2, TRUE, FALSE)
  alpha0 <- alpha
  if (y.binary == TRUE) {
    if (is.null(covar)) {
      lassocv <- cv.glmnet2(lrXsub, y, family = "binomial", foldid = Folds,
        alpha = alpha, type.measure = "auc", keep = TRUE)
    }
    else {
      df0 <- data.frame(y, covar)
      model0 <- glm(y ~ ., family = "binomial", data = df0)
      x0 <- predict(model0)
      lassocv <- cv.glmnet2(lrXsub, y, family = "binomial", foldid = Folds,
        offset = x0, alpha = alpha, type.measure = "auc", 
        keep = TRUE)
    }
    if (Do_Leave_One_Out == T){
      probs = LOO_glment(lrXsub, y, alpha0, nfolds, lambda, foldlist=Folds)
      Results_loo = Get_scores_model(probs, truth = y )
    } else { Results_loo = list() }
  }
  else {
    if (is.null(covar)) {
      lassocv <- glmnet::cv.glmnet(lrXsub, y, alpha = alpha, 
        type.measure = "deviance", keep = TRUE)
    }
    else {
      df0 <- data.frame(y, covar)
      model0 <- lm(y ~ ., data = df0)
      x0 <- predict(model0)
      lassocv <- glmnet::cv.glmnet(lrXsub, y, offset = x0, 
        alpha = alpha, type.measure = "deviance", keep = TRUE)
    }
  }
  lambdavalue <- lambda
  if (is.character(lambda)) {
    if (lambda == "lambda.1se") 
      lambdavalue <- lassocv$lambda.1se
    if (lambda == "lambda.min") 
      lambdavalue <- lassocv$lambda.min
  }
  idrow <- max(which(lassocv$glmnet.fit$lambda >= lambdavalue))
  coeflr <- as.vector(coef(lassocv, s = lambda))[-1]
  lrselect <- which(coeflr != 0)
  idlrXsub[lrselect, ]
  coeflogcontrast <- rep(0, ncol(x))
  for (i in (1:length(coeflr))) {
    coeflogcontrast[idlrXsub[i, 1]] <- coeflogcontrast[idlrXsub[i, 
      1]] + coeflr[i]
    coeflogcontrast[idlrXsub[i, 2]] <- coeflogcontrast[idlrXsub[i, 
      2]] - coeflr[i]
  }
  varlogcontrast <- which(abs(coeflogcontrast) > 0)
  coeflogcontrast <- coeflogcontrast[varlogcontrast]
  (names.select <- colnames(x)[varlogcontrast])
  (positive <- ifelse(coeflogcontrast > 0, 1, 0))
  positive <- factor(positive, levels = c(0, 1), labels = c("negative", 
    "positive"))
  logcontrast = as.matrix(log(x)[, varlogcontrast]) %*% coeflogcontrast
  if (is.null(covar)) {
    predictions <- as.numeric(predict(lassocv, lrX, s = lambdavalue))
  }
  else {
    predictions <- as.numeric(predict(lassocv, lrX, s = lambdavalue, 
      newoffset = x0))
  }
  coeflogcontrast <- 2 * coeflogcontrast/sum(abs(coeflogcontrast))
  if (y.binary == TRUE) {
    AUC_signature <- pROC::auc(pROC::roc(y, as.numeric(predictions), 
      quiet = TRUE))[[1]]
    if (length(varlogcontrast) == 0) 
      AUC_signature <- 0.5
    mcvAUC <- lassocv$cvm[idrow]
    sdcvAUC <- lassocv$cvsd[idrow]
  }
  else {
    mcvMSE <- lassocv$cvm[idrow]
    sdcvMSE <- lassocv$cvsd[idrow]
    Rsq <- 0
    if (length(varlogcontrast) > 0) {
      Rsq <- as.numeric(cor(predictions, y)^2)
    }
  }
  if (y.binary == TRUE) {
    results <- list(taxa.num = varlogcontrast, taxa.name = names.select, 
      `log-contrast coefficients` = coeflogcontrast, predictions = predictions, 
      `apparent AUC` = AUC_signature, `mean cv-AUC` = mcvAUC, 
      `sd cv-AUC` = sdcvAUC, LOO_results=Results_loo)
  }
  else {
    results <- list(taxa.num = varlogcontrast, taxa.name = names.select, 
      `log-contrast coefficients` = coeflogcontrast, predictions = predictions, 
      `apparent Rsq` = Rsq, `mean cv-MSE` = mcvMSE, `sd cv-MSE` = sdcvMSE)
  }
  return(results)
}


coda_glmnet_null_2 = function (x, y, niter = 100, covar = NULL, lambda = "lambda.1se",  alpha = 0.9, sig = 0.05, Do_Leave_One_Out = F)  {
  alpha0 <- alpha
  lambda0 <- lambda
  covar0 <- covar
  y.binary <- ifelse(dim(table(y)) == 2, TRUE, FALSE)
  y1 <- y
  lrmatrix <- logratios_matrix(x)
  lrX <- lrmatrix[[1]]
  idlrX <- lrmatrix[[2]]
  nameslrX <- lrmatrix[[3]]
  accuracy <- rep(0, niter)
  #add funcitons for loo
  accuracy_loo <- rep(0, niter)
  accuracy_balanced_loo <- rep(0, niter)

  fv = generate_stratified_folds_balanced(y, nfolds)
  summary_table <- table(fv, y)


  lrX1 = lrX
  for (i in (1:niter)) {
    #y1 <- sample(y1)
    lrX1 = lrX1[sample(nrow(lrX1)), ]
    lr <- coda_glmnet0_2(x = x, lrX = lrX, idlrX = idlrX, 
      nameslrX = nameslrX, y = y, lambda = lambda0, covar = covar0, 
      alpha = alpha0, Do_Leave_One_Out =Do_Leave_One_Out, Folds = fv)
    if (y.binary == TRUE) {
      res <- lr$`mean cv-AUC`
    }
    else {
      res <- lr$`mean cv-MSE`
    }
    accuracy[i] <- res
    if (Do_Leave_One_Out == T){
      accuracy_loo[i] = lr$LOO_results[['AUC']]
      accuracy_balanced_loo[i] = lr$LOO_results$`Balanced Accuracy`
    }
    print(paste("iter", i))
  }
  results <- list(accuracy = accuracy, accuracy_loo = accuracy_loo, balanced_accuracy_loo=accuracy_balanced_loo, `confidence interval` = quantile(accuracy, 
    c((sig/2), (1 - (sig/2)))))
  return(results)
}




######################
####COMMON FUNCTIONS##
######################
#generate_stratified_folds: Folds keepin the proportion of each level constant
#generate_stratified_folds_balanced
#LOO_glment
#cv.glmnet2
#Get_scores_model

generate_stratified_folds <- function(y, nfolds = 3) {
  #Generates a vector of fold belonging for each y value, that keeps the proportions of each level 
  #this vector can be given as an input in cv.glmnet in the argument 'foldid'
  
  set.seed(11233)
  # Get unique classes in y
  y_classes <- unique(y)
  
  # Initialize foldid vector
  foldid <- numeric(length(y))
  
  # For each class, randomly assign samples to folds
  for (class in y_classes) {
    class_indices <- which(y == class)  # Indices of samples of this class
    class_folds <- sample(rep(seq(nfolds), length = length(class_indices)))  # Randomly assign to folds
    foldid[class_indices] <- class_folds  # Assign fold numbers to corresponding class samples
  }
  
  return(foldid)
}

generate_stratified_folds_balanced <- function(y, nfolds = 3, max_attempts = 100) {
  y_classes <- unique(y)
  
  for (attempt in seq_len(max_attempts)) {
    # Initialize foldid vector
    foldid <- numeric(length(y))
    
    # For each class, assign samples to folds in a balanced way
    for (class in y_classes) {
      class_indices <- which(y == class)
      shuffled_indices <- sample(class_indices)
      
      # Round-robin assignment of samples to folds
      fold_assignment <- rep(seq(nfolds), length.out = length(class_indices))
      foldid[shuffled_indices] <- fold_assignment
    }
    
    # Check if any fold has zero samples from any class
    summary_table <- table(foldid, y)
    zero_counts <- any(summary_table == 0)
    
    # If no folds have zero samples from any class, return the foldid
    if (!zero_counts) {
      return(foldid)
    } 
  }
  
  stop("Failed to generate balanced folds after maximum attempts.")
}


LOO_glment = function(lrXsub, y_unique,alpha0, nfolds,lambda, foldlist){
  #Estimates the predictions using a LOO-CV scheme (traina model in all data -1 , test in the 1, and do this for all samples)
  #Returns a probability vector obtained in which each prediction is the output of one iteration
  set.seed(123)
  Log_odds = c()
  for (Subject in seq(1,nrow(lrXsub)) ){
    print( paste0("Running in Subject: ", Subject) )
    #Each subject will be predicted individually, in a model trained with all others 
    Test = lrXsub[Subject, ,  drop=F] 
    Test_y = y_unique[Subject]
    Train = lrXsub[-Subject,]
    Train_y =  y_unique[-Subject]
    fold_vector = foldlist[-Subject]
    lasso_loo = cv.glmnet2( x=Train, y=Train_y, family = "binomial", alpha = alpha0, type.measure = "auc", nfolds = nfolds, keep=F, foldid = fold_vector )
    #Get the model with the best lambda
    lambdavalue <- lambda
    if (is.character(lambda)) {
      if (lambda == "lambda.1se") 
        lambdavalue <- lasso_loo$lambda.1se
      if (lambda == "lambda.min") 
        lambdavalue <- lasso_loo$lambda.min
    }
    sample_pred =  predict(lasso_loo, Test, s=lambdavalue, type="response")
    Log_odds =  c( Log_odds, sample_pred )
  }
  return(Log_odds)
}

Get_scores_model = function(probs, truth=y_unique){
  #Given a probability vector and a ground truth, calculates different metrics of model performance
  #roc_obj = pROC::roc(truth,probs, levels = c(FALSE, TRUE), direction = ">" )
  roc_obj = pROC::roc(controls = probs[truth==F], cases=probs[truth==T]  )
  auc_value = pROC::auc(roc_obj)

  #Other scores
  truth = factor(truth, levels=c(F, T ) )
  predicted <- ifelse(probs > 0.5, levels(truth)[2], levels(truth)[1])
  # Confusion matrix
  #conf_matrix <- caret::confusionMatrix(as.factor(predicted), as.factor(truth))
  conf_matrix <- caret::confusionMatrix(factor(predicted, levels = c(T, F) ), factor(truth,levels = c(T, F)), positive = "TRUE")
  # Extract key scores
  accuracy <- conf_matrix$overall['Accuracy']
  sensitivity <- conf_matrix$byClass['Sensitivity']  # Also known as Recall
  specificity <- conf_matrix$byClass['Specificity']
  precision <- conf_matrix$byClass['Precision']      # Precision (PPV)
  f1_score <- conf_matrix$byClass['F1']
  balanced_accuracy <- conf_matrix$byClass['Balanced Accuracy']
  return(list(
    'AUC' = auc_value,
    'Accuracy' = accuracy,
    'Sensitivity' = sensitivity,
    'Specificity' = specificity,
    'Precision' = precision,
    'F1 Score' = f1_score,
    'Balanced Accuracy' = balanced_accuracy,
    'Confusion Matrix' = conf_matrix,
    'ROC Plot' = plot(roc_obj)
  ))
}


cv.glmnet2 = function (x, y, weights = NULL, offset = NULL, lambda = NULL, 
          type.measure = c("default", "mse", "deviance", "class", "auc", "mae", "C"), nfolds = 10, foldid = NULL, alignment = c("lambda", "fraction"), grouped = TRUE, keep = FALSE, parallel = FALSE, 
          gamma = c(0, 0.25, 0.5, 0.75, 1), relax = FALSE, trace.it = 0, ...) {
  #Allow 2 fold CV
  
  type.measure = match.arg(type.measure)
  alignment = match.arg(alignment)
  if (!is.null(lambda) && length(lambda) < 2) 
    stop("Need more than one value of lambda for cv.glmnet")
  if (!is.null(lambda) && alignment == "fraction") {
    warning("fraction of path alignment not available if lambda given as argument; switched to alignment=`lambda`")
    alignment = "lambda"
  }
  N = nrow(x)
  if (is.null(weights)) 
    weights = rep(1, N)
  else weights = as.double(weights)
  y = drop(y)
  cv.call = glmnet.call = match.call(expand.dots = TRUE)
  which = match(c("type.measure", "nfolds", "foldid", "grouped", 
                  "keep"), names(glmnet.call), FALSE)
  if (any(which)) 
    glmnet.call = glmnet.call[-which]
  glmnet.call[[1]] = as.name("glmnet")
  if (glmnet:::glmnet.control()$itrace) 
    trace.it = 1
  else {
    if (trace.it) {
      glmnet:::glmnet.control(itrace = 1)
      on.exit(glmnet.control(itrace = 0))
    }
  }
  if (is.null(foldid)) 
    foldid = sample(rep(seq(nfolds), length = N))
  else nfolds = max(foldid)
  #if (nfolds < 3) 
  #  stop("nfolds must be bigger than 3; nfolds=10 recommended")
  if (relax) 
    glmnet:::cv.relaxed.raw(x, y, weights, offset, lambda, type.measure, 
                   nfolds, foldid, alignment, grouped, keep, parallel, 
                   trace.it, glmnet.call, cv.call, gamma, ...)
  else glmnet:::cv.glmnet.raw(x, y, weights, offset, lambda, type.measure, 
                     nfolds, foldid, alignment, grouped, keep, parallel, trace.it, 
                     glmnet.call, cv.call, ...)
}

Fig_xgboost_LOO = function(X, Y){
  X %>%  coda4microbiome::logratios_matrix() -> lrmatrix
  lrX <- lrmatrix[[1]]
  idlrX <- lrmatrix[[2]]
  nameslrX <- lrmatrix[[3]]

  X_matrix <- as.matrix(lrX)
  Y_vector <- as.numeric(Y)
  loo_predictions <- numeric(nrow(X_matrix))

  for (i in 1:nrow(X_matrix)) {
    print(i)
    model <- xgboost(
      data = X_matrix[-i, ],
      label = Y_vector[-i],
      nrounds = 100,
      objective = "binary:logistic",
      verbose = 0
    )
    loo_predictions[i] <- predict(model, X_matrix[i, , drop = FALSE])
  }
  return( loo_predictions )
}


Fig_xgboost_LOO_multiclass <- function(X, Y){
  X %>% coda4microbiome::logratios_matrix() -> lrmatrix
  lrX <- lrmatrix[[1]]
  X_matrix <- as.matrix(lrX)
  
  Y_factor <- factor(Y)
  Y_vector <- as.numeric(Y_factor) - 1
  num_classes <- length(levels(Y_factor))
  
  # Matrix: rows = samples, cols = classes
  loo_prob_matrix <- matrix(NA, nrow = nrow(X_matrix), ncol = num_classes)
  
  for (i in 1:nrow(X_matrix)) {
    cat("Sample", i, "\n")
    model <- xgboost(
      data = X_matrix[-i, ],
      label = Y_vector[-i],
      nrounds = 100,
      objective = "multi:softprob",
      num_class = num_classes,
      verbose = 0
    )
    
    # Predict returns a vector of probabilities (one per class)
    prob_i <- predict(model, X_matrix[i, , drop = FALSE])
    loo_prob_matrix[i, ] <- prob_i
  }
  
  colnames(loo_prob_matrix) <- levels(Y_factor)
  return(list(probabilities = loo_prob_matrix, true_labels = Y_factor))
}
