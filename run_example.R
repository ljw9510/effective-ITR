
##### Run example for ITR estimation -------------------------------------------

source("multi_category_method.R")
source("itr_simu_data.R")

set.seed(1)
dat <- make_data(n = 400, p = 10, case = 1, sigma = 1)
x <- dat$x; trt <- dat$trt; y <- dat$y
K <- 4 

x_test <- x[1:10,]
x_training <- x[-c(1:10),]
trt_training <- trt[-c(1:10)]
y_training <- y[-c(1:10)]

fit1 = ITR_Learning(x_training, trt_training, y_training, type = "SABD", K,
                    use_ebw = TRUE, screened = TRUE, augmented = TRUE)
predict_ITR(fit1, x_test, K)

