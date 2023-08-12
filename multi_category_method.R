library(doSNOW)
library(dplyr) 
library(energy) # dcov.test
library(randomForest)
library(WeightIt)
library(nnet) # multinomial regression for weights estimation
library(glmnet)
library(rlist)


##### Modified AD/SABD Learning for effective ITR ------------------------------

## Modified AD/SABD Learning ---------------------------------------------------

simplex <- function(K = 4){

  ## Create vertex simplex.  
  # K : number of treatments.
  A <- -(1 + sqrt(K+1)) / (K^1.5)
  B <- sqrt((K+1) / K)
  vertices <- matrix(0, K, K+1)
  vertices[, 1] <- K ^(-0.5) * matrix(1, K)
  vertices[, 2:(K+1)] <- A * matrix(1, K, K) + B * diag(1, K)
  return(t(vertices))
  
}

modified_adlearn <- function(X, A, Y, w, K = 4, lasso = TRUE) {
  
  ## Perform modified AD-Learning.
  # X : matrix of covariates.
  # A : vector of treatments.
  # Y : vector of observed outcomes.
  # w : vector of weights.
  # K : number of treatments.
  # lasso : a logical indicator of whether to use Lasso. if FALSE, use unpenalized linear model.
  
  set.seed(1)
  n <- nrow(X)
  p <- ncol(X)
  vertices <- simplex(K-1)
  treatment_specific_vertices <- vertices[A, ]
  Y_modified <- (K / (K-1)) * Y * treatment_specific_vertices
  if (lasso == TRUE) {
    fit <- cv.glmnet(X, Y_modified, alpha = 1, standardize.response = T, 
                     weights = w, family = "mgaussian")    
  } else {
    fit <- lm(Y_modified ~ X, weights = w)
  }
  return(fit)
}

modified_sabdlearn <- function(X, A, Y, w, modified, K = 4, lasso = TRUE) {
  
  ## Perform modified SABD-Learning.
  # X : matrix of covariates.
  # A : vector of treatments.
  # Y : vector of observed outcomes.
  # w : vector of weights.
  # modified : a logical indicator of whether to use modified method.
  # K : number of treatments.
  # lasso : a logical indicator of whether to use Lasso. if FALSE, use unpenalized linear model.
  
  set.seed(1)
  
  # 1. Initial AD-Learning Step
  vertices <- simplex(K-1)
  treatment_specific_vertices <- vertices[A, ]
  fit <- modified_adlearn(X, A, Y, w, K, lasso)  
  
  if (lasso == TRUE){
    B <- as.matrix(list.cbind(coef(fit, "lambda.min")))
  } else {
    B <- as.matrix(coef(fit))
  }
  XB <- cbind(1, X) %*% B
  pred_all <- rowSums(treatment_specific_vertices * XB)
  resid_squared <- ((K/(K-1))*Y - pred_all)^2
  AX = cbind(X, A = as.factor(A))
  Y_modified <- (K / (K-1)) * Y * treatment_specific_vertices
  
  # 2. Residual modeling
  resid_fit <- randomForest(resid_squared ~ ., data = AX)
  resid_preds <- predict(resid_fit)
    
  # 3. Modified weights and SABD-Learning fitting step
  if (modified == TRUE){ # proposed method
    w_new <- w / resid_preds  
  } else { # original method
    w_new <- 1 / resid_preds  
  }
  if (lasso == TRUE){
    fit2 <- cv.glmnet(X, Y_modified, alpha = 1, standardize.response = T, 
                      weights = w_new, family = "mgaussian")  
  } else {
    fit2 <- lm(Y_modified ~ X, weights = w_new)
  }
  return(fit2)
}

get_ITR <- function(fit, X, K = 4, lasso = TRUE) {
  
  ## Obtain optimal ITR.
  # X : matrix of covariates in the test set. 
  # K : number of treatments.
  # lasso : a logical indicator of whether Lasso was used.
  
  vertices <- simplex(K-1)
  if (lasso == TRUE) {
    B <- as.matrix(list.cbind(coef(fit, "lambda.min")))  
  } else {
    B <- as.matrix(coef(fit))
  }
  XB <- cbind(1, X) %*% B
  cate_est <- XB %*% t(vertices) 
  itr <- apply(cate_est, 1, function(x) which.max(x))
  return(itr)  
}

## Compute covariate balancing weights -----------------------------------------

get_IPW_weights <- function(data, method = "RF"){
  
  ## Compute inverse probability weighting.
  # data : dataframe with covariates and treatment A.
  # method : method to compute weights, using either random forest or multinomial regression.
  
  set.seed(1)
  
  if (method == "RF"){
    model_rf <- randomForest(A ~ ., data)
    newdat <- as.data.frame(data[,-1]) 
    colnames(newdat) <- colnames(data)[-1] 
    ps.multi <- predict(model_rf, newdat, 'prob')
  } 
  
  if (method == "multinomial"){
    mod_fit <- multinom(A ~ ., data)
    ps.multi <- drop(predict(mod_fit, type = "probs"))
  }
  
  ipw <- rep(0, nrow(ps.multi))
  for (i in 1:nrow(ps.multi)){
    ipw[i] <- 1/as.numeric(ps.multi[i, data$A[i]])
  }
  return(ipw)
}

normalized_weights <- function(w, trt){
  
  ## Obtain normalized weights so that the sum of weights equals to sample size.
  # w : vector of weights.
  # trt : vector of treatments.
  
  trt_list <- unique(trt)
  for (i in 1:length(trt_list)){
    w[trt == trt_list[i]] <- length(trt) * 
      (w[trt == trt_list[i]])/ sum(w[trt == trt_list[i]])
  }
  return(w)
}

finding_weights <- function(y_training, trt_training, x_training, new_x_training, 
                            K = 4, random_set = TRUE){
  
  ## Obtain a n x 4 dataframe of weights with IPW and energy balancing weights for original/screened coviarates.
  # y_training : vector of outcomes in the training set.
  # trt_training : vector of treatments in the training set.
  # x_training : vector of covariates in the training set.
  # new_x_training : vector of screened covariates in the training set.
  # random_set : a logical indicator of whether it is a randomized trial.
  
  set.seed(1)
  AX_training <- data.frame(A = factor(trt_training), x_training)
  AX_training2 <- data.frame(A = factor(trt_training), new_x_training)
  
  energy_w <- tryCatch(weightit(A ~ ., data = AX_training, method = "energy", 
                                estimand = "ATE")$weights,
                       error=function(e) rep(0, nrow(AX_training)))
  energy_w2 <- tryCatch(weightit(A ~ ., data = AX_training2, method = "energy",
                                 estimand = "ATE")$weights,
                        error=function(e) rep(0, nrow(AX_training2)))
  
  wts <- tibble(
    IPW = 1/rep(1/K, length(trt_training)),
    energy = normalized_weights(energy_w, trt_training),
    IPW_s = 1/rep(1/K, length(trt_training)),
    energy_s = normalized_weights(energy_w2, trt_training),
  )
  
  if (random_set == FALSE){
    ipw <- get_IPW_weights(AX_training)
    ipw_s <- get_IPW_weights(AX_training2)
    
    estimated_IPW <- normalized_weights(ipw, trt_training)
    estimated_IPW_s <- normalized_weights(ipw_s, trt_training)
    
    wts$IPW <- estimated_IPW
    wts$IPW_s <- estimated_IPW_s
  }
  return(wts)
}


## Conduct variable screening --------------------------------------------------

distance_cov_test <- function(Y, A, X, K = 4){
  
  ## Perform variable screening based on distance covariance test.
  # Y : vector of observed outcomes.
  # A : vector of treatments.
  # X : matrix of covariates.
  # K : number of treatments.

  set.seed(1)
  pvalue <- numeric(ncol(X))
  for (i in 1:ncol(X)) {
    temp_p <- numeric(K)
    for (j in seq_len(K)) {
      temp_p[j] <- dcov.test(X[A==j,i], Y[A==j], R = 99)$p
    }
    pvalue[i] <- min(temp_p)
  }
  return(which(pvalue<0.05))
}


## Conduct outcome augmentation ------------------------------------------------

augmented_function <- function(Y, X, split_aug = TRUE){
  
  ## Obtain augmented otucome.
  # Y : vector of observed outcomes.
  # X : matrix of covariates.
  # split_aug : a logical indicator of whether to split data.
  
  set.seed(1)
  YX <- data.frame(Y = Y, X)
  n <- nrow(YX)
  
  if (split_aug == TRUE){
    split_idx <- sample(c(rep(TRUE, n/2), rep(FALSE, n/2)))
    YX1 <- YX[split_idx,]
    YX2 <- YX[!split_idx,]
    main_eff_fit1 <- randomForest(Y ~ ., data = YX1)
    main_eff_fit2 <- randomForest(Y ~ ., data = YX2)
    pred1 <- predict(main_eff_fit2, YX1[,-1])
    pred2 <- predict(main_eff_fit1, YX2[,-1])
    res1 <- YX1$Y - pred1
    res2 <- YX2$Y - pred2
    res <- c(res1, res2)
  } else {
    main_eff_fit1 <- randomForest(Y ~ ., data = YX)
    pred1 <- predict(main_eff_fit1, YX[,-1])
    res <- YX$Y - pred1
  }
  return(res[order(as.numeric(names(res)))])
}


## Compute empirical value -----------------------------------------------------

pred_model <- function(Y, X_training, A, X_test, d, K = 4){
  
  ## Estimate treatment-free effect.
  # Y : vector of observed outcomes.
  # X_training : matrix of covariates in the training set.
  # A : vector of treatments.
  # X_test : matrix of covariates in the test set.
  # d : vector of optimal ITR.
  # K : number of treatments.
  
  set.seed(1)
  YXA_training <- data.frame(Y = Y, X_training, 
                             trt = factor(A, levels = c(1:K)))
  cond_mean <- randomForest(Y ~ ., data = YXA_training)
  YXA_pred <- data.frame(X_test, trt = factor(d, levels = c(1:K)))
  pred1 <- predict(cond_mean, YXA_pred)
  return(pred1)
  
}

calculate_value <- function(pred_trt, potential_outcome, augmented, pred1 = NA){
  
  ## Compute empirical value based on potential outcome.
  # pred_trt : vector of predicted optimal treatment.
  # potential outcome : matrix of potential outcomes.
  # augmented : a logical indicator of whether it is for augmented outcome.
  # pred1 : vector of predicted treatment-free effect, needed when augmented = TRUE.
  
  n <- nrow(potential_outcome)
  
  if (augmented == F) { # original outcome
    val <- mean(sapply(1:n, function(i) potential_outcome[i, pred_trt[i]]))
  } else { # augmented outcome
    augmented_y <- potential_outcome - pred1
    val <- mean(sapply(1:n, function(i) augmented_y[i, pred_trt[i]])) + mean(pred1)
  }
  return(val)
}


## Method Evaluation -----------------------------------------------------------

evaluation <- function(wts, x_training, trt_training, y_training, x_test, K = 4,
                       new_x_training, new_x_test, d, 
                       augmented, y_test, trt_test, y_training_old = NA,
                       potential_outcome_test){

  ## Evaluate methods based on misclassification and empirical value.
  # wts : n x 4 dataframe of weights with IPW and energy balancing weights for original/screened coviarates.
  # x_training : vector of covariates in the training set.
  # trt_training : vector of treatments in the training set.
  # y_training : vector of outcomes in the training set.
  # x_test : vector of covariates in the test set. 
  # K : number of treatments.
  # new_x_training : vector of screened covariates in the training set.
  # new_x_test : vector of screened covariates in the test set.
  # d : vector of optimal ITR.
  # augmented : a logical indicator of whether it is for augmented outcome.
  # y_test : vector of outcomes in the test set.
  # y_training_old : vector of original outcomes in the training set, needed when augmented = TRUE.
  # potential_outcome_test : matrix of potential outcomes in the test set.
  
  set.seed(1)

  # AD-Learning for original covariates
  itrs1 <-  lapply(wts[,1:2], function(w) tryCatch(get_ITR(modified_adlearn(
    x_training, trt_training, y_training, w = w, K = K), x_test, K = K), 
    error = function(e) rep(NA, nrow(x_test))))  
  
  # AD-Learning for screened covariates
  itrs2 <-  lapply(wts[,3:4], function(w) tryCatch(get_ITR(modified_adlearn(
    new_x_training, trt_training, y_training, w = w, K = K), new_x_test, K = K), 
    error = function(e) rep(NA, nrow(new_x_test))))
  
  # SABD-Learning with IPW
  itrs3 <- tryCatch(get_ITR(modified_sabdlearn(
   x_training, trt_training, y_training, w = wts$IPW, modified = FALSE, K = K), 
   x_test, K = K), error=function(e) rep(NA, nrow(x_test)))  
  
  # SABD-Learning with energy balancing weights
  itrs4 <- tryCatch(get_ITR(modified_sabdlearn(
    x_training, trt_training, y_training, w = wts$energy, modified = TRUE, K = K), 
    x_test, K = K), error=function(e) rep(NA, nrow(x_test)))  
  
  # SABD-Learning for screened covariates with IPW 
  itrs5 <- tryCatch(get_ITR(modified_sabdlearn(
    new_x_training, trt_training, y_training, w = wts$IPW_s, modified = FALSE, K = K), 
    new_x_test, K = K), error=function(e) rep(NA, nrow(new_x_test)))  

  # SABD-Learning for screened covariates with energy balancing weights
  itrs6 <- tryCatch(get_ITR(modified_sabdlearn(
    new_x_training, trt_training, y_training, w = wts$energy_s, modified = TRUE, K = K), 
    new_x_test, K = K), error=function(e) rep(NA, nrow(new_x_test)))

  itrs_AD <- c(itrs1, itrs2)
  itrs_SABD <- list(itrs3, itrs4, itrs5, itrs6)
  names(itrs_SABD) <- paste0("SABD_", names(itrs_AD))
  names(itrs_AD) <- paste0("AD_", names(itrs_AD))
  itrs <- c(itrs_AD, itrs_SABD)
  
  # misclassification rate
  misclassification <- unlist(lapply(itrs, function(x) tryCatch(1-sum(x==d)/length(d),
                                                          error = function(e) NA)))

  # empirical value
  if (augmented == TRUE){
    pred1 <- pred_model(y_training_old, x_training, trt_training, x_test, d, K = K)
    pred2 <- pred_model(y_training_old, new_x_training, trt_training, new_x_test, d, K = K)
  } else {
    pred1 = pred2 = NA
  }

  # empirical value for original covariates
  emp_value1 <- lapply(itrs[c(1,2,5,6)], function(w) tryCatch(calculate_value(
    w, potential_outcome_test, augmented, pred1), error = function(e) rep(NA, 4)))
 
  # empirical value for screened covariates
  emp_value2 <- lapply(itrs[c(3,4,7,8)], function(w) tryCatch(calculate_value(
    w, potential_outcome_test, augmented, pred2), error = function(e) rep(NA, 4)))
  
  emp_value <- unlist(c(emp_value1[1:2], emp_value2[1:2], emp_value1[3:4], emp_value2[3:4]))

  # bench mark empirical value
  bench_mark_set <- tibble(
    `emp_value_opt` = d,
    `emp_value_1` = rep(1, 10000),
    `emp_value_2` = rep(2, 10000),
    `emp_value_3` = rep(3, 10000),
    `emp_value_4` = rep(4, 10000),
  )
    
  ideal_emp_value <- unlist(lapply(bench_mark_set, function(w) tryCatch(calculate_value(
    w, potential_outcome_test, augmented, pred1), error = function(e) rep(NA, 5))))
    
  return(c(misclassification, emp_value, ideal_emp_value))
}

## Summarize results -----------------------------------------------------------

summary_result <- function(result, p0 = c(20,40,60)){
  
  ## Calculate averages and standard deviations from the results.
  # result : dataframe of results obtained from main function.
  
  result0 <- data.frame(t(round(apply(result[[1]], 2, function(x) mean(x, na.rm=T)), 3)))
  result1 <- data.frame(t(round(apply(result[[1]], 2, function(x) sd(x, na.rm=T)), 3)))
  for (i in 2:length(p0)){
    temp <- t(round(apply(result[[i]], 2, function(x) mean(x, na.rm=T)), 3))
    temp1 <- t(round(apply(result[[i]], 2, function(x) sd(x, na.rm=T)), 3))
    result0 <- rbind(result0, temp)
    result1 <- rbind(result1, temp1)
  }
  rownames(result0) = rownames(result1) = p0
  return(list(mean = result0, sd = result1))
}

## Run simulation --------------------------------------------------------------

main <- function(case = 1, n_train = 200,
                 augmented = TRUE, split_aug = FALSE,
                 case_desc = "nonlinear_", 
                 p0 = c(20, 40, 60), K = 4, n_itr = 100){
  
  ## Run simulation.
  # case : a numerical value of scenario (details in 'itr_simu_data.R')
  # n_train : number of observations in the training set.
  # augmented : a logical indicator of whether to use outcome augmentation.
  # split_aug : a logical indicator of whether to split data for outcome augmentation.
  # case_desc : description of case.
  # p0 : a set of covariate dimensions.
  # K : number of treatments.
  # n_itr : number of iterations.
  
  cl <- makeCluster(4);
  registerDoSNOW(cl)
  
  result0_list = list()
  
  for (j in 1:length(p0)){
    p <- p0[j]
    
    pb <- txtProgressBar(max = n_itr, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    result0 <- foreach(i = 1:n_itr, .combine = rbind, .options.snow = opts) %dopar% {
      
      .GlobalEnv$case <- case
      
      source("multi_category_method.R")
      source("itr_simu_data.R")
      
      set.seed(i)
      
      dat <- make_data(n = (10000+n_train), p = p, case = case, sigma = 1)
      
      x <- dat$x
      trt <- dat$trt
      y <- dat$y
      opt_trt <- dat$opt_trt
      potential_outcome <- dat$potential_outcome
      true_var_set <- dat$true_var_set
      random_set <- ifelse(case %in% c(1, 3, 5, 7), TRUE, FALSE)
      
      in_test <- sample(c(rep(TRUE, 10000), rep(FALSE, n_train)))
      x_test <- x[in_test,]
      trt_test <- trt[in_test]
      y_test <- y[in_test]
      d <- opt_trt[in_test]
      potential_outcome_test <- potential_outcome[in_test,]
      
      x_training <- x[!in_test,]
      trt_training <- trt[!in_test]
      y_training <- y[!in_test]
      
      # screening
      selected_cov <- distance_cov_test(y_training, trt_training, x_training, K = 4)
      if (length(selected_cov)!=0){
        new_x_training <- x_training[,selected_cov]
        new_x_test <- x_test[,selected_cov]
      }
      
      # augmented for main effects
      if (augmented == TRUE){
        y_training_old <- y_training
        y_training <- as.numeric(augmented_function(y_training, new_x_training, split_aug = split_aug))
      } else {
        y_training_old <- NA
      } 
      
      # weights
      wts <- finding_weights(y_training, trt_training, x_training,
                             new_x_training, K = 4, random_set)
      
      # Evaluation
      evaluation(wts, x_training, trt_training, y_training, x_test, K = 4,
                 new_x_training, new_x_test, d, augmented, y_test, trt_test, 
                 y_training_old, potential_outcome_test)
      
      
    }  
    result0 <- data.frame(result0)
    method_names <- c("AD", "AD_e", "AD_s", "AD_e_s", 
                      "SABD", "SABD_e", "SABD_s", "SABD_e_s")
    colnames(result0)[1:16] <- c(paste0("misclassification_", method_names),
                                 paste0("emp_val_", method_names))
    colnames(result0)[17:21] <- paste0("emp_val_", colnames(result0)[17:21])
    result0_list[[j]] <- result0
    
  }      
  return(result0_list)
  stopCluster(cl)
}

## Fit modified AD/SABD-Learning -----------------------------------------------

ITR_Learning <- function(x_training, trt_training, y_training, type = "SABD", K,
                         use_ebw = TRUE, screened = TRUE, augmented = TRUE){
  
  ## Fit modified AD/SABD-Learning.
  
  ## input
  # x_training : vector of covariates in the training set.
  # trt_training : vector of treatments in the training set.
  # y_training : vector of outcomes in the training set.
  # type : a character string either "AD" or "SABD".
  # K : number of treatments.
  # use_ebw : a logical indicator of whether to use distributional covariate balancing weight, i.e. energy balancing weight.
  # screened : a logical indicator of whether to use variable screenning.
  # augmented : a logical indicator of whether to use outcome augmentation.
  
  ## output
  # 1. AD/SABD-Learning fit
  # 2. index of covariates used for model fitting
  
  set.seed(1)
  
  # screening
  if (screened == TRUE) {
    selected_cov <- distance_cov_test(y_training, trt_training, x_training, K = K)
    if (length(selected_cov)!=0){
      x_training <- x_training[,selected_cov]
    }  
  } else {
    p <- ncol(x_training)
    selected_cov <- 1:p
  }
  
  # augmented for main effects
  if (augmented == TRUE){
    y_training <- as.numeric(augmented_function(y_training, x_training, split_aug = FALSE))  
  }
  
  # weights
  AX_training <- data.frame(A = factor(trt_training), x_training)
  if (use_ebw == TRUE) {
    energy_w <- tryCatch(weightit(A ~ ., data = AX_training, method = "energy", 
                                  estimand = "ATE")$weights,
                         error=function(e) rep(0, nrow(AX_training)))
    w <- normalized_weights(energy_w, trt_training)
  } else {
    ipw <- get_IPW_weights(AX_training)
    w <- normalized_weights(ipw, trt_training)
  }
  
  # ITR-Learning
  if (type == "AD"){
    fit0 = modified_adlearn(x_training, trt_training, y_training, w = w, K = K)
    return(list(fit0, "used_covariate" = selected_cov))
  } 
  
  if (type == "SABD"){
    fit0 = modified_sabdlearn(x_training, trt_training, y_training, w = w, modified = FALSE, K = K)
    return(list(fit0, "used_covariate" = selected_cov))
  }
}



predict_ITR <- function(fit, x_test, K){
  
  ## Provide predicted individual treatment rule (ITR).
  
  ## input
  # fit : the output of ITR-Learning.
  # x_test : vector of covariates in the test set.
  # K : number of treatments.
  
  fit0 <- fit[[1]]
  selected_cov <- fit[[2]]
  x_test <- x_test[,selected_cov]
  itr <- get_ITR(fit0, x_test, K = K) 
  return(itr)
}
