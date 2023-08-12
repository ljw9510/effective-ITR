
##### Scenario -----------------------------------------------------------------

## Case 1
# 1 : linear interaction, simple treatment-free, randomized trial
# 2 : linear interaction, simple treatment-free, observational study

## Case 2
# 3 : linear interaction, complicated treatment-free, randomized trial
# 4 : linear interaction, complicated treatment-free, observational study

## Case 3
# 5 : non-linear interaction, simple treatment-free, randomized trial
# 6 : non-linear interaction, simple treatment-free, observational study

## Case 4
# 7 : non-linear interaction, complicated treatment-free, randomized trial
# 8 : non-linear interaction, complicated treatment-free, observational study

#-------------------------------------------------------------------------------


library(truncnorm)

eff_main_fn <- function(x, case) { # main (treatment-free) effect function
  main <- switch(
    case,
    `1` =  1.0 + 2.0 * x[,1] + 2.0 * x[,2],
    `2` =  1.0 + 1.0 * x[,5] + 3.0 * x[,6] + 2.0 * x[,1] * x[,2],
    `3` =  1.0 + 1.0 * x[,5] + 1.0 * x[,5]^2 + 2.0 * exp(-x[,1]*x[,2]) + sin(x[,3]),
    `4` =  1.0 + 1.0 * x[,5] + 1.0 * x[,5]^2 + 2.0 * exp(-x[,1]*x[,2]) + sin(x[,3]), 
    `5` =  1.0 + 2.0 * x[,1] + 2.0 * x[,2] + 2.0 * x[,4] - 2.0 * x[,4]^2 + 2.0 * x[,1]*x[,2],
    `6` =  1.0 + 2.0 * x[,1] + 2.0 * x[,2] + 2.0 * x[,4] - 2.0 * x[,4]^2 + 2.0 * x[,1]*x[,2], 
    `7` =  1.0 + 2.0 * x[,1] + 2.0 * x[,2] + 2.0 * x[,4] - 2.0 * x[,4]^2 + 2.0 * x[,1]*x[,2] + 2.0 * exp(-x[,1]*x[,2]) + sin(x[,3]), 
    `8` =  1.0 + 2.0 * x[,1] + 2.0 * x[,2] + 2.0 * x[,4] - 2.0 * x[,4]^2 + 2.0 * x[,1]*x[,2] + 2.0 * exp(-x[,1]*x[,2]) + sin(x[,3]), 
  )
  return(main)
}

eff_interaction_fn <- function(x, a, case) { # interaction function
  interaction <- switch(
    case,
    `1` = case_when(a == 1 ~ 0.75 + 1.5 * x[,1] + 1.5 * x[,2] + 1.5 * x[,3] + 1.5 * x[,4],
                    a == 2 ~ 0.75 + 1.5 * x[,1] - 1.5 * x[,2] - 1.5 * x[,3] + 1.5 * x[,4],
                    a == 3 ~ 0.75 + 1.5 * x[,1] - 1.5 * x[,2] + 1.5 * x[,3] - 1.5 * x[,4],
                    a == 4 ~ 0.75 - 1.5 * x[,1] + 1.5 * x[,2] - 1.5 * x[,3] - 1.5 * x[,4]),
    `2` = case_when(a == 1 ~ 0.50 + 2.0 * x[,1] + 1.0 * x[,2] + 1.0 * x[,3],
                    a == 2 ~ 1.00 + 1.0 * x[,1] - 1.0 * x[,2] - 1.0 * x[,3],
                    a == 3 ~ 1.50 + 3.0 * x[,1] - 1.0 * x[,2] + 1.0 * x[,3],
                    a == 4 ~ 1.00 - 1.0 * x[,1] - 1.0 * x[,2] + 1.0 * x[,3]),
    `3` = case_when(a == 1 ~ 0.75 + 1.5 * x[,1] + 1.5 * x[,2] + 1.5 * x[,3] + 1.5 * x[,4],
                    a == 2 ~ 0.75 + 1.5 * x[,1] - 1.5 * x[,2] - 1.5 * x[,3] + 1.5 * x[,4],
                    a == 3 ~ 0.75 + 1.5 * x[,1] - 1.5 * x[,2] + 1.5 * x[,3] - 1.5 * x[,4],
                    a == 4 ~ 0.75 - 1.5 * x[,1] + 1.5 * x[,2] - 1.5 * x[,3] - 1.5 * x[,4]), 
    `4` = case_when(a == 1 ~ 0.50 + 2.0 * x[,1] + 1.0 * x[,2] + 1.0 * x[,3],
                    a == 2 ~ 1.00 + 1.0 * x[,1] - 1.0 * x[,2] - 1.0 * x[,3],
                    a == 3 ~ 1.50 + 3.0 * x[,1] - 1.0 * x[,2] + 1.0 * x[,3],
                    a == 4 ~ 1.00 - 1.0 * x[,1] - 1.0 * x[,2] + 1.0 * x[,3]),
    `5` = case_when(a == 1 ~ 0.50 + 1.0 * x[,1] - 2.0 * x[,4] + 0.5 * x[,4]^2 ,
                    a == 2 ~ 1.00 + 1.0 * x[,1] + 1.0 * x[,4] - 1.0 * x[,4]^2 ,
                    a == 3 ~ 1.50 + 2.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2,
                    a == 4 ~ 1.00 - 1.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2),    
    `6` = case_when(a == 1 ~ 0.50 + 1.0 * x[,1] - 2.0 * x[,4] + 0.5 * x[,4]^2 ,
                    a == 2 ~ 1.00 + 1.0 * x[,1] + 1.0 * x[,4] - 1.0 * x[,4]^2 ,
                    a == 3 ~ 1.50 + 2.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2,
                    a == 4 ~ 1.00 - 1.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2),    
    `7` = case_when(a == 1 ~ 0.50 + 1.0 * x[,1] - 2.0 * x[,4] + 0.5 * x[,4]^2 ,
                    a == 2 ~ 1.00 + 1.0 * x[,1] + 1.0 * x[,4] - 1.0 * x[,4]^2 ,
                    a == 3 ~ 1.50 + 2.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2,
                    a == 4 ~ 1.00 - 1.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2),    
    `8` = case_when(a == 1 ~ 0.50 + 1.0 * x[,1] - 2.0 * x[,4] + 0.5 * x[,4]^2 ,
                    a == 2 ~ 1.00 + 1.0 * x[,1] + 1.0 * x[,4] - 1.0 * x[,4]^2 ,
                    a == 3 ~ 1.50 + 2.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2,
                    a == 4 ~ 1.00 - 1.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2),    
    
  )
  return(interaction)
}

counterfactuals_fn <- function(x, case) { # to get optimal ITR
  cate <- switch(
    case,
    `1` = data.frame(
      A1 = 0.75 + 1.5 * x[,1] + 1.5 * x[,2] + 1.5 * x[,3] + 1.5 * x[,4],
      A2 = 0.75 + 1.5 * x[,1] - 1.5 * x[,2] - 1.5 * x[,3] + 1.5 * x[,4],
      A3 = 0.75 + 1.5 * x[,1] - 1.5 * x[,2] + 1.5 * x[,3] - 1.5 * x[,4],
      A4 = 0.75 - 1.5 * x[,1] + 1.5 * x[,2] - 1.5 * x[,3] - 1.5 * x[,4]), 
    `2` = data.frame(
      A1 = 0.5 + 2.0 * x[,1] + 1.0 * x[,2] + 1.0 * x[,3],
      A2 = 1.0 + 1.0 * x[,1] - 1.0 * x[,2] - 1.0 * x[,3],
      A3 = 1.5 + 3.0 * x[,1] - 1.0 * x[,2] + 1.0 * x[,3],
      A4 = 1.0 - 1.0 * x[,1] - 1.0 * x[,2] + 1.0 * x[,3]),
    `3` = data.frame(
      A1 = 0.75 + 1.5 * x[,1] + 1.5 * x[,2] + 1.5 * x[,3] + 1.5 * x[,4],
      A2 = 0.75 + 1.5 * x[,1] - 1.5 * x[,2] - 1.5 * x[,3] + 1.5 * x[,4],
      A3 = 0.75 + 1.5 * x[,1] - 1.5 * x[,2] + 1.5 * x[,3] - 1.5 * x[,4],
      A4 = 0.75 - 1.5 * x[,1] + 1.5 * x[,2] - 1.5 * x[,3] - 1.5 * x[,4]), 
    `4` = data.frame(
      A1 = 0.5 + 2.0 * x[,1] + 1.0 * x[,2] + 1.0 * x[,3],
      A2 = 1.0 + 1.0 * x[,1] - 1.0 * x[,2] - 1.0 * x[,3],
      A3 = 1.5 + 3.0 * x[,1] - 1.0 * x[,2] + 1.0 * x[,3],
      A4 = 1.0 - 1.0 * x[,1] - 1.0 * x[,2] + 1.0 * x[,3]),
    `5` = data.frame( 
      A1 =  0.5 + 1.0 * x[,1] - 2.0 * x[,4] + 0.5 * x[,4]^2 ,
      A2 =  1.0 + 1.0 * x[,1] + 1.0 * x[,4] - 1.0 * x[,4]^2 ,
      A3 =  1.5 + 2.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2,
      A4 =  1.0 - 1.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2),    
    `6` = data.frame( 
      A1 =  0.5 + 1.0 * x[,1] - 2.0 * x[,4] + 0.5 * x[,4]^2 ,
      A2 =  1.0 + 1.0 * x[,1] + 1.0 * x[,4] - 1.0 * x[,4]^2 ,
      A3 =  1.5 + 2.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2,
      A4 =  1.0 - 1.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2),    
    `7` = data.frame( 
      A1 =  0.5 + 1.0 * x[,1] - 2.0 * x[,4] + 0.5 * x[,4]^2 ,
      A2 =  1.0 + 1.0 * x[,1] + 1.0 * x[,4] - 1.0 * x[,4]^2 ,
      A3 =  1.5 + 2.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2,
      A4 =  1.0 - 1.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2),    
    `8` = data.frame(
      A1 =  0.5 + 1.0 * x[,1] - 2.0 * x[,4] + 0.5 * x[,4]^2 ,
      A2 =  1.0 + 1.0 * x[,1] + 1.0 * x[,4] - 1.0 * x[,4]^2 ,
      A3 =  1.5 + 2.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2,
      A4 =  1.0 - 1.0 * x[,1] - 1.0 * x[,4] - 1.0 * x[,4]^2),    
  )
  return(cate)
}

err_fn <- function(x, a, case){ # heterogeneous variance, function of covariate and treatment
  err <- switch(
    case,
    `1` = 0.25 + 0.2 * (1.5 - x[,2])^2,
    `2` = 0.25 + 2*x[,2]*(x[,2]>0) + x[,3]*(x[,3]>0)*(a==1) + x[,4]*(x[,4]>0)*(a==2),
    `3` = 0.25 + 0.2 * (1.5 - x[,2])^2,
    `4` = 0.25 + 2*x[,2]*(x[,2]>0) + x[,3]*(x[,3]>0)*(a==1) + x[,4]*(x[,4]>0)*(a==2),
    `5` = 0.25 + 0.2 * (1.5 - x[,2])^2,
    `6` = 0.25 + 2*x[,2]*(x[,2]>0) + x[,3]*(x[,3]>0)*(a==1) + x[,4]*(x[,4]>0)*(a==2),
    `7` = 0.25 + 0.2 * (1.5 - x[,2])^2,
    `8` = 0.25 + 2*x[,2]*(x[,2]>0) + x[,3]*(x[,3]>0)*(a==1) + x[,4]*(x[,4]>0)*(a==2),
  )
  return(err)
}

make_data <- function(n, p = 20, case = 1, sigma = 1) {
  
  ## covariates from truncated Normal(0, 1, -3, 3)
  x <- matrix(rtruncnorm(n*p, a=-3, b=3, mean = 0, sd = 1), nrow = n, ncol = p)
  
  ## treatment assignment
  prob_trt <- rep(0.25, n) # meaningful only for randomized trial
  
  trt <- rep(0, n)
  if (case %in% c(1, 3, 5, 7)) {
    trt <- sample(c(1, 2, 3, 4), n, replace = TRUE)
  } else {
    for (i in 1:n) {
      if (x[i,1] < 0) {a <- sample(x = c(1,2,3,4), 1, replace = T, prob = c(0.25, 0.25, 0.25, 0.25))}
      if (x[i,1] > 0) {a <- sample(x = c(1,2,3,4), 1, replace = T, prob = c(0.4, 0.2, 0.2, 0.2))}
      trt[i] <- a
    } 
  }
  
  ## outcome
  eff_main <- eff_main_fn(x, case)
  eff_interaction <- eff_interaction_fn(x, trt, case) 
  sigma_x <- err_fn(x, trt, case)
  err <- rnorm(n, 0, sd = sigma) * sqrt(sigma_x)
  y <- eff_main + eff_interaction + err

  eff_cate <- counterfactuals_fn(x, case) 
  opt_trt <- max.col(eff_cate)
  potential_outcome <- eff_main + eff_cate + err
  
  # need only for variable screening evaluation
  if (case %in% c(1, 3)) {
    true_var_set <- 1:4
  } else if (case %in% c(2, 4)) {
    true_var_set <- 1:3
  } else {
    true_var_set <- c(1, 4)
  }  
  
  list(x   = x,
       trt = as.integer(trt),
       y   = y,
       opt_trt = opt_trt,
       true_var_set = true_var_set,
       potential_outcome = potential_outcome)
  
}



