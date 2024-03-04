# #######################################################################################################################

# Author of code: Glen P. Martin.

# This is code runs the worked examples in a manuscript entitled: 
# Statistical Primer: Sample Size Requirements for Developing and Validating Risk Prediction Models

# Authors:
#   Glen P. Martin
#   Richard D. Riley
#   Joie Ensor
#   Stuart W. Grant

# #######################################################################################################################

###-----------------------------------------------------------------------------
### Worked Example 1: model development - binary outcome
###-----------------------------------------------------------------------------

## Pre-specify the required information for the sample size calculations:
Y_prev <- 3/100 #the anticipated binary outcome prevalence
P <- 25 #the number of predictor parameters
anticipated_R2 <- 0.023 #the anticipated level of performance
SF <- 0.9 #the level of overfitting that is targeted 

## Now run the sample size calculation:
# install.packages("pmsampsize")
library(pmsampsize)
pmsampsize(type = "b",
           csrsquared = anticipated_R2,
           parameters = P,
           shrinkage = SF,
           prevalence = Y_prev)

## An alternative approach is to pre-specifying performance based on an
## anticipated c-statistic:
anticipated_c_stat <- 0.74
pmsampsize(type = "b",
           cstatistic = anticipated_c_stat,
           parameters = P,
           shrinkage = 0.9,
           prevalence = Y_prev)

## In the above, the required sample size is too big. Option (i) is to accept a
## higher level of overfitting. For example, dropping the shrinkage from 0.9 to
## 0.8 gives the following required sample size:
pmsampsize(type = "b",
           csrsquared = anticipated_R2,
           parameters = P,
           shrinkage = 0.8,
           prevalence = Y_prev)

## Option (ii) is to consider reducing the number of predictor parameters:
P <- 20
pmsampsize(type = "b",
           csrsquared = anticipated_R2,
           parameters = P,
           shrinkage = SF,
           prevalence = Y_prev)
P <- 18
pmsampsize(type = "b",
           csrsquared = anticipated_R2,
           parameters = P,
           shrinkage = SF,
           prevalence = Y_prev)
P <- 15
pmsampsize(type = "b",
           csrsquared = anticipated_R2,
             parameters = P,
             shrinkage = SF,
             prevalence = Y_prev)



###-----------------------------------------------------------------------------
### Worked Example 2: model validation - simulation-based approach
# See Box 2 of Snell et al. doi: 10.1016/j.jclinepi.2021.02.011
###-----------------------------------------------------------------------------
set.seed(1234)

y_prev <- 0.031 #set the outcome proportion

#specify distribution of LP:
mu_LP <- -4.02 #mean
sd_LP <- 1.17 #standard deviation

#Specify a sequence of sample sizes to try for the size of the validation study:
N <- seq(from = 2500, to = 13000, by = 100) 

n_iter <- 500 #number of simulation repeats per sample size, N

OE_width <- matrix(NA, ncol = length(N), nrow = n_iter)
cal_slope_width <- matrix(NA, ncol = length(N), nrow = n_iter)
c_stat_width <- matrix(NA, ncol = length(N), nrow = n_iter)
for(n in 1:length(N)){
  for(i in 1:n_iter) {
    # generate Y from assumed LP distribution
    LP <- rnorm(N[n], mu_LP, sd_LP)
    PR <- exp(LP) / (1 + exp(LP))
    Y <- rbinom(N[n], size = 1, prob = PR)
    
    #calculate O:E ratio:
    O_E <- mean(Y) / mean(PR)
    SE_OE <- sqrt((1 - y_prev) / (N[n]*y_prev))
    OE_width[i, n] <- 2 * qnorm(0.025, lower = F) * SE_OE
    
    #calculate calibration slope:
    cal_model <- glm(Y ~ LP, family = binomial(link = "logit"))
    cal_slope <- as.numeric(coef(cal_model)[2])
    cal_slope_SE <- as.numeric(sqrt(diag(vcov(cal_model)))[2])
    cal_slope_width[i, n] <- 2 * qnorm(0.025, lower = F) * cal_slope_SE
    
    #calculate the c-statistic
    roc_estimate <- pROC::roc(response = Y,
                        predictor = PR,
                        direction = "<",
                        levels = c(0,1),
                        ci = TRUE)
    c_stat <- as.numeric(pROC::auc(roc_estimate))
    c_stat_width[i, n] <- roc_estimate$ci[2] - roc_estimate$ci[1]
    
  }
}
simulation_sample_size_results <- data.frame("Sample Size" = N,
                                             "Width of O:E Ratio 95% CI" = colMeans(OE_width),
                                             "Width of Calibration Slope 95% CI" = colMeans(cal_slope_width),
                                             "Width of C-Stat 95% CI" = colMeans(c_stat_width))

#step i
min(simulation_sample_size_results$Sample.Size[which(
  simulation_sample_size_results$Width.of.O.E.Ratio.95..CI < 0.2
)])

#step ii
min(simulation_sample_size_results$Sample.Size[which(
  simulation_sample_size_results$Width.of.Calibration.Slope.95..CI < 0.2
)])

#step iii
min(simulation_sample_size_results$Sample.Size[which(
  simulation_sample_size_results$Width.of.C.Stat.95..CI < 0.1
)])

### final sample size
min(simulation_sample_size_results$Sample.Size[which(
  simulation_sample_size_results$Width.of.O.E.Ratio.95..CI < 0.2 &
    simulation_sample_size_results$Width.of.Calibration.Slope.95..CI < 0.2 &
    simulation_sample_size_results$Width.of.C.Stat.95..CI < 0.1
)])

###-----------------------------------------------------------------------------
### Worked Example 2: model validation - closed-form solution
###-----------------------------------------------------------------------------

# install.packages("pmvalsampsize")
library(pmvalsampsize)
pmvalsampsize(type = "b", 
              prevalence = 0.031, 
              cstatistic = 0.74, 
              lpnormal = c(-4.02,1.17), 
              oeciwidth = 0.2, 
              csciwidth = 0.2, 
              cstatciwidth = 0.1)
