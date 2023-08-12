
##### Main code for simulation -------------------------------------------------

## input
# case : a numerical value of scenario (details in 'itr_simu_data.R')
# n_train : number of observations in the training set.
# augmented : a logical indicator of whether to use outcome augmentation.
# split_aug : a logical indicator of whether to split data for outcome augmentation.
# case_desc : description of case.
# p0 : a set of covariate dimensions.
# K : number of treatments.
# n_itr : number of iterations.

source("itr_simu_data.R")
source("multi_category_method.R")

case = 8; 
n_train = 200; 
augmented = TRUE; 
split_aug = FALSE; 
case_desc = ""

system.time({results <- main(case = case, 
                             n_train = n_train,
                             augmented = augmented, 
                             split_aug = split_aug,
                             case_desc = case_desc, 
                             p0 = c(20, 40, 60), 
                             K = 4, 
                             n_itr = 20)})

summary_result(results)[[1]][,1:8] # misclassification rate
summary_result(results)[[1]][,-c(1:8)] # empirical values

random_set <- ifelse(case %in% c(1, 3, 5, 7), TRUE, FALSE)
save.image(paste(case_desc, "case_", case, ifelse(random_set == TRUE, "_random_", "_observational_"), 
                 "n_", n_train, "_", ifelse(augmented, "augmented", "original"), ".RData", sep=""))

