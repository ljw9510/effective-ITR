
library(doSNOW)
library(dplyr) 
library(energy) # dcov.test
library(randomForest)
library(WeightIt)
library(nnet) # multinomial regression for weights estimation
library(glmnet)
library(rlist)

# for policy
library(policytree)
library(causalDML)
library(lmtest)
library(sandwich)

# for AD binary
library(grpreg)

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
  # lasso : a logical indicator of whether to use Lasso. 
  #         if FALSE, use unpenalized linear model.
  
  set.seed(1)
  n <- nrow(X)
  p <- ncol(X)
  vertices <- simplex(K-1)
  treatment_specific_vertices <- vertices[A, ]
  Y_modified <- (K / (K-1)) * Y * treatment_specific_vertices
  if ((lasso == TRUE) & (p > 1)) {
    fit <- cv.glmnet(X, Y_modified, alpha = 1, standardize.response = T, 
                     weights = w, family = "mgaussian")    
  } else {
    fit <- lm(Y_modified ~ data.matrix(X), weights = w)
  }
  return(fit)
}

modified_sabdlearn <- function(fit, X, A, Y, w, modified, K = 4, lasso = TRUE) {
  
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
  
  B <- fit
  if (sum(is.na(B))>0){
    XB <- cbind(1, X)[,-which(is.na(B))] %*% B[-which(is.na(B)),]  
  } else {
    XB <- data.matrix(cbind(1, X)) %*% B
  }
  
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
  if ((lasso == TRUE) & (ncol(X) > 1)) {
    fit2 <- cv.glmnet(X, Y_modified, alpha = 1, standardize.response = T, 
                      weights = w_new, family = "mgaussian")  
  } else {
    fit2 <- lm(Y_modified ~ data.matrix(X), weights = w_new)
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
    XB <- cbind(1, X) %*% B
  } else {
    #B <- as.matrix(coef(fit))
    B <- fit
    if (sum(is.na(B))>0){
      XB <- cbind(1, X)[,-which(is.na(B))] %*% B[-which(is.na(B)),]
    } else {
      XB <- data.matrix(cbind(1, X)) %*% B
    }
  }
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

augmented_function <- function(Y, X){
  
  ## Obtain augmented otucome.
  # Y : vector of observed outcomes.
  # X : matrix of covariates.
  
  set.seed(1)
  YX <- data.frame(Y = Y, X)
  n <- nrow(YX)
  
  main_eff_fit1 <- randomForest(Y ~ ., data = YX)
  pred1 <- predict(main_eff_fit1, YX[,-1])
  res <- Y - pred1
  return(res)
}

## Generate simulated dataset --------------------------------------------------

generate_data <- function(case = 1, n_train = 200, p = 60, no.seed = 1,
                          binary = FALSE){
  
  ## Generate simulated data with 4 treatments.
  # case: (integer) specified the simulation scenario. It takes values from 1 to 8.
  #       odd numbers correspond to randomized trials.
  #       even numbers correspond to observational studies.
  # n_train: (integer) specifies the size of the training set.
  # p: (integer) specifies the dimension of covariates.
  # no.seed: (integer) seed number used for data generation.
  # binary: (logical) type of outcome variable. The default is continuous (TRUE).
  
  set.seed(no.seed)
  
  # Data generation
  if (binary == TRUE){
    dat <- make_data_binary(n = (10000+n_train), p = p, sigma = 1, case = case)
  } else {
    dat <- make_data(n = (10000+n_train), p = p, sigma = 1, case = case)  
  }
  
  x <- dat$x
  trt <- dat$trt
  y <- dat$y
  opt_trt <- dat$opt_trt
  
  in_test <- sample(c(rep(TRUE, 10000), rep(FALSE, n_train)))
  x_test <- x[in_test,]
  trt_test <- trt[in_test]
  y_test <- y[in_test]
  d <- opt_trt[in_test]
  
  x_training <- x[!in_test,]
  trt_training <- trt[!in_test]
  y_training <- y[!in_test]
  
  output <- list(x_training, trt_training, y_training, 
                 x_test, trt_test, y_test, d)
  names(output) <- c("x_training", "trt_training", "y_training", 
                     "x_test", "trt_test", "y_test", "true_opt_trt_test")
  return(output)
  
}

## ITR-Learning with other methods for simulated dataset -----------------------

all_evaluation <- function(x_training, trt_training, y_training,
                           x_test, trt_test, y_test, true_opt_trt_test = NA,
                           randomized_trial = FALSE, binary = FALSE){
  
  ## Evaluate all methods including benchmark approaches based on accuracy.
  # (approaches) AD-Learning, SABD-Learning, policytree, causalDML
  
  # x_training: (data.matrix) data matrix of covariates in the training set.
  # trt_training: (vector) vector of treatments in the training set.
  # y_training: (vector) vector of outcomes in the training set.
  # x_test: (data.matrix) data matrix of covariates in the test set.
  # trt_test: (vector) vector of treatments in the test set.
  # y_test: (vector) vector of outcomes in the test set.
  # true_opt_trt_test: (vector) vector of true optimal treatments in the test set.
  #                    The default is NA for real data analysis.
  #                    If true_opt_trt_test is given, accuracy is provided.
  # randomized_trial: (logical) either randomized trial or observational study (default). 
  # binary: (logical) type of outcome variable. The default is continuous (TRUE).
  
  # Specify other factors
  K <- length(unique(trt_training))
  n_train <- nrow(x_training)
  
  # Finding ITR
  set.seed(1)
  
  if (binary == FALSE) { # continuous outcome
    # Variable screening
    selected_cov <- distance_cov_test(y_training, trt_training, x_training, K)
    new_x_training <- data.matrix(x_training[,selected_cov])
    new_x_test <- data.matrix(x_test[,selected_cov])
    
    # Outcome augmentation for variance reduction
    y_training_aug_s <- as.numeric(augmented_function(y_training, new_x_training))
    
    # Estimate weights
    ## IPW
    if (randomized_trial  == TRUE) {
      ipw <- 1/rep(0.25, n_train)
    } else {
      AX_training <- data.frame(A = factor(trt_training), x_training)
      ipw <- get_IPW_weights(AX_training)
      ipw <- normalized_weights(ipw, trt_training)
    }
    ## Distributional covariate balancing weights (DCBW) with energy balancing weights
    AX_training2 <- data.frame(A = factor(trt_training), new_x_training)
    energy_w <- tryCatch(weightit(A ~ ., data = AX_training2, method = "energy", 
                                  estimand = "ATE")$weights,
                         error=function(e) rep(0, nrow(AX_training2)))
    energy_w <- normalized_weights(energy_w, trt_training)
    
    # ITR estimation
    
    ## AD-Learning
    m_AD <- as.matrix(list.cbind(coef(
      modified_adlearn(x_training, trt_training, y_training, w = ipw, K = K, lasso = T),
      "lambda.min")))
    estimated_AD <- get_ITR(m_AD, x_test, K, lasso = F)
    
    ## SABD-Learning
    m_SABD <- as.matrix(list.cbind(coef(
      modified_sabdlearn(m_AD, x_training, trt_training, y_training, w = ipw, 
                         modified = FALSE, K = K, lasso = T),
      "lambda.min")))
    estimated_SABD <- get_ITR(m_SABD, x_test, K = K, lasso = F) 
    
    ## Our proposed method with penalized optimization using glmnet
    m_pre_ours <- as.matrix(list.cbind(coef(
      modified_adlearn(new_x_training, trt_training, y_training_aug_s, w = energy_w, K = K, lasso = T),
      "lambda.min")))
    m_ours <- as.matrix(list.cbind(coef(
      modified_sabdlearn(m_pre_ours, new_x_training, trt_training, y_training_aug_s, 
                         w = energy_w, modified = TRUE, K = K, lasso = T),
      "lambda.min")))
    estimated_ours <- get_ITR(m_ours, new_x_test, K = K, lasso = F)  
    
    ## Our proposed method with constrained optimization (need to optimize)
    #print("")
    
    ## policytree
    estimated_DR <- fit_DR(x_training, y_training, trt_training, x_test, depth = 2)
    
    ## causalDML
    estimated_DML <- fit_DML(x_training, y_training, trt_training, x_test, depth = 2)
    
    # Calculate accuracy
    if (sum(is.na(true_opt_trt_test)) != TRUE){
      AD <- mean(estimated_AD == true_opt_trt_test)
      SABD <- mean(estimated_SABD == true_opt_trt_test)
      ours <- mean(estimated_ours == true_opt_trt_test)
      DR <- mean(estimated_DR == true_opt_trt_test)
      DML <- mean(estimated_DML == true_opt_trt_test, na.rm = T)
      
      estimated_ITR <- list(estimated_AD, estimated_SABD, estimated_ours,
                            estimated_DR, estimated_DML)
      accuracy <- c(AD, SABD, ours, DR, DML)

      names(estimated_ITR) <- c("AD", "SABD", "Proposed", "policytree", "causalDML")
      names(accuracy) <- c("AD", "SABD", "Proposed", "policytree", "causalDML")
      result <- list("accuracy" = accuracy, "estimated_ITR" = estimated_ITR)
      return(result)
      
    } else {
      
      estimated_ITR <- list(estimated_AD, estimated_SABD, estimated_ours,
                            estimated_DR, estimated_DML)
      names(estimated_ITR) <- c("AD", "SABD", "Proposed", "policytree", "causalDML")
      return(estimated_ITR)
    }

  } else { # binary outcome
    print("")
  }
}

## Fit effective ITR-Learning -----------------------------------------------

ITR_Learning <- function(x_training, trt_training, y_training,
                         binary = FALSE){
  
  ## Fit the effective and robust ITR model.
  
  ## input
  # x_training: (data.matrix) data matrix of covariates in the training set.
  # trt_training: (vector) vector of treatments in the training set.
  # y_training: (vector) vector of outcomes in the training set.
  # binary: (logical) Type of outcome variable. The default is continuous (TRUE).
  
  ## output
  # 1. Model fit (B matrix)
  # 2. index of covariates used for model fitting
  
  set.seed(1)
  
  # Specify other factors
  K <- length(unique(trt_training))
  
  # Finding ITR
  if (binary == FALSE) { # continuous outcome
    # Variable screening
    selected_cov <- distance_cov_test(y_training, trt_training, x_training, K)
    new_x_training <- data.matrix(x_training[,selected_cov])
    
    # Outcome augmentation for variance reduction
    y_training_aug_s <- as.numeric(augmented_function(y_training, new_x_training))
    
    # Estimate Distributional covariate balancing weights (DCBW) with energy balancing weights
    AX_training2 <- data.frame(A = factor(trt_training), new_x_training)
    energy_w <- tryCatch(weightit(A ~ ., data = AX_training2, method = "energy", 
                                  estimand = "ATE")$weights,
                         error=function(e) rep(0, nrow(AX_training2)))
    energy_w <- normalized_weights(energy_w, trt_training)
    
    m_pre_ours <- as.matrix(list.cbind(coef(
      modified_adlearn(new_x_training, trt_training, y_training_aug_s, w = energy_w, K = K, lasso = T),
      "lambda.min")))
    m_ours <- as.matrix(list.cbind(coef(
      modified_sabdlearn(m_pre_ours, new_x_training, trt_training, y_training_aug_s, 
                         w = energy_w, modified = TRUE, K = K, lasso = T),
      "lambda.min")))
    
  } else { # binary outcome
    print("Not yet!")
  }
  
  return(list("estimated_B" = m_ours, 
              "used_covariate" = selected_cov,
              "no.treatments" = K))
  
}

predict_ITR <- function(fit, x_test){
  
  ## Provide predicted individual treatment rule (ITR).
  
  ## input
  # fit : the output of ITR-Learning.
  # x_test : vector of covariates in the test set.
  
  ## output
  # estimated optimal treatments
  
  fitted_model <- fit$estimated_B
  selected_cov <- fit$used_covariate
  K <- fit$no.treatments
  
  new_x_test <- data.matrix(x_test[,selected_cov])
  estimated_optimal_trt <- get_ITR(fitted_model, new_x_test, K = K, lasso = F)  
  return(estimated_optimal_trt)
}

### policy learning with double machine learning (DML) -------------------------

fit_DML <- function(x_training, y_training, trt_training, x_test, depth = 2){
  
  # Create the ML method to be used for nuisance parameter prediction
  forest <- create_method("forest_grf",args=list(tune.parameters = "all",seed=1234))
  
  # Run the main function that outputs nuisance parameters, APO and ATE
  cDML <- tryCatch(DML_aipw(y_training,as.factor(trt_training),x_training,
                            ml_w=list(forest),ml_y=list(forest)),
                   error = function(e) NA) # need enough sample size to fit cDML
  
  if (sum(is.na(cDML))==T) {
    return (rep(NA, nrow(x_test)))
    
  } else {
    opt.tree <- policy_tree(x_training, cDML$APO$gamma, depth)  
    itr_policy_dml <- predict(opt.tree, x_test) # Predict treatment
    return (itr_policy_dml)
  }
  
}

### policy learning with doubly robust score (DR) ------------------------------

fit_DR <- function(x_training, y_training, trt_training, x_test, depth = 2){
  
  multi.forest <- grf::multi_arm_causal_forest(x_training, y_training, factor(trt_training))
  Gamma.matrix <- double_robust_scores(multi.forest)   # Compute doubly robust reward estimates.
  opt.tree_dr <- policy_tree(x_training, Gamma.matrix, depth) # Fit a depth 2 tree on a random training subset.
  itr_policy_dr <- predict(opt.tree_dr, x_test)
  
  return (itr_policy_dr)
  
}
