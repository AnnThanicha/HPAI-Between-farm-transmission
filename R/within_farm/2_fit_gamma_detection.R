
# To fit gamma distribution to infectious duration  ------------------------------

#include libraries ####
packages <- c("ggplot2","deSolve","tidyverse","bbmle", "fitdistrplus","dplyr")


## Load all package
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)



# Function to get negative log-likelihood of gamma distribution------------------
# Negative log-likelihood function with gamma
neg_gamma_ll_gamma <- function(intercept,coef, inf_t, size, var) {
  
  mean <- intercept + (coef*size)
  var <- var # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  
    # Negative log-likelihood
  log_likelihood <- sum(dgamma(inf_t, shape = shape, rate = rate, log = TRUE)) 
  return(-log_likelihood) 
}

# second gamma distribution with polynomial ------------------------------------
neg_gamma_ll_gamma2 <- function(intercept,coef1,coef2, inf_t, size, var) {
  
  mean <- intercept + coef1*size + coef2*size^2
  var <- var # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  
  # Negative log-likelihood
  log_likelihood <- sum(dgamma(inf_t, shape = shape, rate = rate, log = TRUE)) 
  return(-log_likelihood) 
}

# Negative log-likelihood with normal distribution------------------------------
neg_normal_ll_norm <- function(intercept_mean,coef_mean, inf_t, size, sd) {
  # Calculate mu as a linear function of the predictor (x)
  mean <- intercept_mean + coef_mean * size
  
  # Calculate the negative log-likelihood using dnorm
  log_likelihood <- sum(dnorm(inf_t, mean = mean, sd = sd, log = TRUE))
  
  return(-log_likelihood)  # Return negative log-likelihood
}


# mle is used to find the minimum of the negative log-likelihood, the function should be a negative log-likelihood
# The log-likelihood is often used because it transforms products of probabilities into sums, 
# which are computationally more stable and easier to work with.
# The log-likelihood provides a more intuitive understanding of the goodness-of-fit. 
# Higher values of the log-likelihood (closer to zero) indicate better model fit. 
# By minimizing the negative log-likelihood, you seek the parameters that lead to the highest likelihood of observing the given data.



# Import output from within-farm model------------------------------------------
output_chicken_fullday <- readRDS("./results/within_farm/output_chicken_fullday.rds")
output_turkey_fullday <- readRDS("./results/within_farm/output_turkey_fullday.rds")
output_duck_fullday <- readRDS("./results/within_farm/output_duck_fullday.rds")



# Fit distribution to chicken --------------------------------------------------

data <- output_chicken_fullday
#Summary detection time 
detect_time <- data %>% group_by(farm_size,run) %>% summarise(detect_time = max(time),
                                                                     farm_size = max(farm_size))  
# detect_time$detect_time <- detect_time$detect_time +0.1 # add 0.1 to make it full day
detect_time$log_farm_size <- log(detect_time$farm_size)

# Visualize data
plot(detect_time$farm_size,detect_time$detect_time )
plot(detect_time$log_farm_size,detect_time$detect_time )

# It follows the linear,after log(4) (approximately farm size higher than 60), 
# so let's exclude the farm size lower than 60 and fit the model again
detect_time <- detect_time[detect_time$farm_size> 60,]



# Fit distributions
fit_normal <- fitdist(detect_time$detect_time, "norm")
fit_exponential <- fitdist(detect_time$detect_time, "exp")
fit_gamma <- fitdist(detect_time$detect_time, "gamma")

# Compare fits visually
plot(fit_normal)
plot(fit_exponential)
plot(fit_gamma)

# Goodness-of-fit statistics
gofstat(list(fit_normal, fit_exponential, fit_gamma))

# from the goodness of fit test (rad how to interpret in One note), normal distribution is the best one, but we might have a problem with negative value if we use it, 
# so I would play safe and select gamma distribution // Noted dicuss with Egil on the realistic detecetion, and run this again

# fit gamma distribution
fit <- fitdist(detect_time$detect_time, "gamma")
summary(fit) # use shape and rate 
summary(rgamma(2000,shape = 34.14, rate = 3.03 )) # this is estimate value from fitdist


  
summary(detect_time$detect_time)
var(detect_time$detect_time)
mean(detect_time$detect_time)


# Use bbmle to minimize the negative log-likelihood with log_farm_size
fit1 <-  mle2(neg_gamma_ll_gamma,start = list(coef_mean = 1, var = 1), 
             data = list(inf_t = detect_time$detect_time, size = detect_time$log_farm_size ),
             method = "L-BFGS-B",                                      # Use constrained optimization
             lower = list( coef_mean = 0.1, var = 0) )
summary(fit1)
warnings()

fit2 <-  mle2(neg_gamma_ll_gamma2 ,start = list(intercept = 1, coef1 =0.1, coef2 = 0.1,  var = 1), 
              data = list(inf_t = detect_time$detect_time, size = detect_time$log_farm_size ),
              method = "L-BFGS-B",                                      # Use constrained optimization
              lower = list( intercept = 0.1, coef1 =0, coef2 = 0,  var = 0.1) )

summary(fit2)



# function to simulate t_inf from fit value
t_inf_samp1 <- function(size) {
  mean <- 1.05856673*size
  var <-  0.47 # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# function to simulate t_inf from fit value with polynomial
t_inf_samp2 <- function(size) {
  mean <-1.0516315 + 0.3599556*size + 0.0602508*size^2
  var <-  0.1954552 # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# simulate with different size from parameters that got from fit 
detect_time$t_inf <- sapply(detect_time$log_farm_size, FUN= t_inf_samp  )
detect_time$t_inf2 <- sapply(detect_time$log_farm_size, FUN= t_inf_samp2  )

plot(detect_time$farm_size, detect_time$t_inf)
plot(detect_time$farm_size, detect_time$t_inf2,col = 3)
points( detect_time$farm_size,detect_time$detect_time,col = 2)

# check RMSE between the two fitting
sqrt(sum((detect_time$detect_time-detect_time$t_inf)^2)/nrow(detect_time))
sqrt(sum((detect_time$detect_time-detect_time$t_inf2)^2)/nrow(detect_time))

# the second model with polynomial fit much better


# Fit data with turkey-----------------------------------------------------------
data <- output_turkey_fullday
#Summary detection time 
detect_time <- data %>% group_by(farm_size,run) %>% summarise(detect_time = max(time),
                                                              farm_size = max(farm_size))  
# detect_time$detect_time <- detect_time$detect_time +0.1 # add 0.1 to make it full day
detect_time$log_farm_size <- log(detect_time$farm_size)

# Visualize data
plot(detect_time$farm_size,detect_time$detect_time )
plot(detect_time$log_farm_size,detect_time$detect_time )

# It follows the linear,after log(4) (approximately farm size higher than 60), 
# so let's exclude the farm size lower than 60 and fit the model again
detect_time <- detect_time[detect_time$farm_size> 60,]

# fit with gamma and polynomial 
fit2 <-  mle2(neg_gamma_ll_gamma2 ,start = list(intercept = 1, coef1 =0.1, coef2 = 0.1,  var = 1), 
              data = list(inf_t = detect_time$detect_time, size = detect_time$log_farm_size ),
              method = "L-BFGS-B",                                      # Use constrained optimization
              lower = list( intercept = 0.1, coef1 =0, coef2 = 0,  var = 0.1) )

summary(fit2)

# function to simulate t_inf from fit value with polynomial
t_inf_samp2 <- function(size) {
  mean <-1.1439855 + 0.8948005*size + 0.0365499*size^2
  var <-  0.2061337 # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# simulate with different size from parameters that got from fit 
detect_time$t_inf2 <- sapply(detect_time$log_farm_size, FUN= t_inf_samp2  )

plot(detect_time$farm_size, detect_time$t_inf2,col = 3)
points( detect_time$farm_size,detect_time$detect_time,col = 2)

# check RMSE between the two fitting
sqrt(sum((detect_time$detect_time-detect_time$t_inf2)^2)/nrow(detect_time))


# Fit data with duck------------------------------------------------------------

neg_gamma_ll_gamma3 <- function(intercept,coef, inf_t, size, var) {
  
  mean <- intercept + (coef*size)
  var <- var*size # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  
  # Negative log-likelihood
  log_likelihood <- sum(dgamma(inf_t, shape = shape, rate = rate, log = TRUE)) 
  return(-log_likelihood) 
}


data <- output_duck_fullday
#Summary detection time 
detect_time <- data %>% group_by(farm_size,run) %>% summarise(detect_time = max(time),
                                                              farm_size = max(farm_size))  
# detect_time$detect_time <- detect_time$detect_time +0.1 # add 0.1 to make it full day
detect_time$log_farm_size <- log(detect_time$farm_size)

# Visualize data
plot(detect_time$farm_size,detect_time$detect_time )

plot(detect_time$log_farm_size,detect_time$detect_time )

# It follows the linear,after log(4) (approximately farm size higher than 60), 
# so let's exclude the farm size lower than 60 and fit the model again
detect_time <- detect_time[detect_time$farm_size> 60,]

# Fit distributions
fit_normal <- fitdist(detect_time$detect_time, "norm")
fit_exponential <- fitdist(detect_time$detect_time, "exp")
fit_gamma <- fitdist(detect_time$detect_time, "gamma")

# Goodness-of-fit statistics
gofstat(list(fit_normal, fit_exponential, fit_gamma))


# fit with gamma distribution
fit1 <-  mle2(neg_gamma_ll_gamma,start = list(intercept = 1 ,coef = 1, var = 0.1), 
              data = list(inf_t = detect_time$detect_time, size = detect_time$log_farm_size ),
              method = "L-BFGS-B",                                      # Use constrained optimization
              lower = list( intercept = 0.1 ,coef = 0.1, var = 0.01) )

summary(fit1)

# fit with gamma distribution
fit3 <-  mle2(neg_gamma_ll_gamma3,start = list(intercept = 1 ,coef = 1, var = 0.1), 
              data = list(inf_t = detect_time$detect_time, size = detect_time$log_farm_size ),
              method = "L-BFGS-B",                                      # Use constrained optimization
              lower = list( intercept = 0.1 ,coef = 0.1, var = 0.001) )

summary(fit3)


# function to simulate t_inf from fit value
t_inf_samp1 <- function(size) {
  mean <- 2.1781620 + 0.47698461*size
  var <-  0.057 # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}


# function to simulate t_inf from fit value with polynomial
t_inf_samp3 <- function(size) {
  mean <-2.93114591 + 0.46725553*size 
  var <-  0.01038639*size # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# simulate with different size from parameters that got from fit 
detect_time$t_inf <- sapply(detect_time$log_farm_size, FUN= t_inf_samp  )
detect_time$t_inf3 <- sapply(detect_time$log_farm_size, FUN= t_inf_samp3  )

plot(detect_time$farm_size, detect_time$t_inf)
plot(detect_time$farm_size, detect_time$t_inf3,col = 3)
points( detect_time$farm_size,detect_time$detect_time,col = 2)

# check RMSE between the two fitting
sqrt(sum((detect_time$detect_time-detect_time$t_inf)^2)/nrow(detect_time)) 
sqrt(sum((detect_time$detect_time-detect_time$t_inf3)^2)/nrow(detect_time))

# gamma 3 with variance increase with the size perform best
