# To fit gamma distribution to human exposure ----------------------------------
# Human exposure is calculated from cumulative sick-animal days

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

# Import output from within-farm model------------------------------------------
output_chicken_fullday <- readRDS("./results/within_farm/output_chicken_fullday.rds")
output_turkey_fullday <- readRDS("./results/within_farm/output_turkey_fullday.rds")
output_duck_fullday <- readRDS("./results/within_farm/output_duck_fullday.rds")



# For chicken -------------------------------------
# Calculate cumulative infectious-day animals
data <- output_chicken_fullday

# Explore the relationship----------------
human_exposure <- data %>% dplyr::group_by(farm_size,run) %>% summarise(sick_day = sum(I)*0.1,# *0.1 from deltat
                                                              suscep_day = max(farm_size) * max(time),
                                                              farm_size = max(farm_size),
                                                              detect_time = max(time))

plot(human_exposure$sick_day)
plot(human_exposure$farm_size, human_exposure$sick_day)
plot(human_exposure$detect_time ,human_exposure$sick_day)
plot(human_exposure$suscep_day ,human_exposure$sick_day)

human_exposure$prop_day <- human_exposure$sick_day/human_exposure$suscep_day
plot(human_exposure$prop_day )

# the proportion of sick animal day and suscep animal day is high in small farms 
# because the number of index case is high (10) compared to farm size
# Let exclude the small farm size
human_exposure <- human_exposure[human_exposure$farm_size >= 100,]

# Try take log farm size and sick animal day to see relationship
human_exposure$log_farm_size <- log(human_exposure$farm_size)
human_exposure$log_detect_time <- log(human_exposure$detect_time)
human_exposure$log_suscep_day <- log(human_exposure$suscep_day)

plot(human_exposure$log_farm_size,human_exposure$sickl_day)
plot(human_exposure$log_detect_time ,human_exposure$sick_day)
plot(human_exposure$log_suscep_day,human_exposure$sick_day)

# fit gamma distribution -------------------------
gamma_ll_1 <- function(intercept1,intercept2,coef1,coef2,  suscep_day, sick_day) {
  
  mean <- intercept1 + coef1*suscep_day 
  var <- intercept2 + coef2*suscep_day  # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  
  # Negative log-likelihood
  log_likelihood <- sum(dgamma(sick_day, shape = shape, rate = rate, log = TRUE)) 
  return(-log_likelihood) 
}



# Fit distribution
fit1 <-  mle2(gamma_ll_1, start = list(intercept1=1 ,intercept2 = 1 ,coef1 = 3,coef2 = 1), 
              data = list(suscep_day = human_exposure$suscep_day, sick_day = human_exposure$sick_day ),
              method = "L-BFGS-B", lower = c(intercept1 = 0.001, intercept2 = 0.001, coef1 = 0, coef2 = 0))
summary(fit1)



# Check fit

# function to simulate prediction from fit value
pred_gamma_ll_1 <- function(suscep_day) {
  mean <-  2.2742e+02 + 3.1337e-02*suscep_day
  var <-  8.6092e-01 + 1.4177e+01*suscep_day # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# simulate with different size from parameters that got from fit 
pred_fit1 <- sapply(human_exposure$suscep_day, FUN= pred_gamma_ll_1 )
plot(pred_fit1)

human_exposure$pred_fit1 <- pred_fit1 
hist(pred_fit1 )
summary(pred_fit1)
hist(human_exposure$sick_day)
summary(human_exposure$sick_day)
plot(human_exposure$sick_day)
# check RMSE between the two fitting
sqrt(sum((human_exposure$sick_day-pred_fit1)^2)/nrow(human_exposure))

plot(human_exposure$suscep_day,human_exposure$sick_day,col = 3)
points( human_exposure$suscep_day,human_exposure$pred_fit1,col = 2)

plot( human_exposure$suscep_day,human_exposure$pred_fit1,col = 2)
points(human_exposure$suscep_day,human_exposure$sick_day,col = 3)


# Try with different coefficients for 
gamma_ll_2 <- function(intercept1,intercept2,coef1,coef2, coef3, coef4, suscep_day, detect_time,farm_size,sick_day) {
  
  mean <- intercept1 + coef1*farm_size +coef2*detect_time 
  var <- intercept2 + coef3*farm_size +coef4*detect_time   # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  
  # Negative log-likelihood
  log_likelihood <- sum(dgamma(sick_day, shape = shape, rate = rate, log = TRUE)) 
  return(-log_likelihood) 
}

fit1.1 <-  mle2(gamma_ll_2, start = list(intercept1=1 ,intercept2 = 1 ,coef1 = 3,coef2 = 1, coef3 =1, coef4 = 1), 
                data = list(suscep_day = human_exposure$suscep_day, sick_day = human_exposure$sick_day, detect_time =  human_exposure$detect_time,farm_size =  human_exposure$farm_size),
                method = "L-BFGS-B", lower = c(intercept1 = 0.001, intercept2 = 0.001, coef1 = 0, coef2 = 0, coef3 =0, coef4 = 0))
summary(fit1.1)

# function to simulate prediction from fit value
pred_gamma_ll_1.1 <- function(suscep_day, detect_time, farm_size) {
  mean <-  8.3863e+00 + 3.9156e-01*farm_size +4.4606e+01*detect_time 
  var <-  1.0001e+00 + 2.1993e+02*farm_size +1.0291e+00*detect_time # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

  # simulate with different size from parameters that got from fit 
pred_fit1.1 <- mapply(pred_gamma_ll_1.1, human_exposure$suscep_day, human_exposure$detect_time,human_exposure$farm_size )
  plot(pred_fit1.1)
  human_exposure$pred_fit1.1 <- pred_fit1.1 
  hist(pred_fit1.1 )
  summary(pred_fit1.1)
  hist(human_exposure$sick_day)
  summary(human_exposure$sick_day)
  plot(human_exposure$sick_day)
  # check RMSE between the two fitting
  sqrt(sum((human_exposure$sick_day-pred_fit1.1)^2)/nrow(human_exposure))
  
  plot(human_exposure$sick_day,col = 3)
  points(human_exposure$pred_fit1.1,col = 2)
  
  plot( human_exposure$pred_fit1.1,col = 2)
  points(human_exposure$sick_day,col = 3)

  
  

# For Turkey --------------------------------------------------------------------

data <-  output_turkey_fullday

# Explore the relationship
human_exposure <- data %>% dplyr::group_by(farm_size,run) %>% summarise(sick_day = sum(I)*0.1,# *0.1 from deltat
                                                                        suscep_day = max(farm_size) * max(time),
                                                                        farm_size = max(farm_size),
                                                                        detect_time = max(time))

plot(human_exposure$sick_day)
plot(human_exposure$farm_size, human_exposure$sick_day)
plot(human_exposure$detect_time ,human_exposure$sick_day)
plot(human_exposure$suscep_day ,human_exposure$sick_day)

human_exposure$prop_day <- human_exposure$sick_day/human_exposure$suscep_day
plot(human_exposure$prop_day )

# the proportion of sick animal day and suscep animal day is high in small farms 
# because the number of index case is high (10) compared to farm size
# Let exclude the small farm size
human_exposure <- human_exposure[human_exposure$farm_size >= 100,]

# Try take log farm size and sick animal day to see relationship
human_exposure$log_farm_size <- log(human_exposure$farm_size)
human_exposure$log_detect_time <- log(human_exposure$detect_time)
human_exposure$log_suscep_day <- log(human_exposure$suscep_day)

plot(human_exposure$log_farm_size,human_exposure$sickanimal_day)
plot(human_exposure$log_detect_time ,human_exposure$sickanimal_day)
plot(human_exposure$log_suscep_day,human_exposure$sickanimal_day)

summary(human_exposure)
# fit gamma distribution -------------------------
gamma_ll_1 <- function(intercept1,intercept2,coef1,coef2,  suscep_day, sick_day) {
  
  mean <- intercept1 + coef1*suscep_day 
  var <- intercept2 + coef2*suscep_day  # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  
  # Negative log-likelihood
  log_likelihood <- sum(dgamma(sick_day, shape = shape, rate = rate, log = TRUE)) 
  return(-log_likelihood) 
}



# Fit distribution
fit2 <-  mle2(gamma_ll_1, start = list(intercept1 = 1, intercept2 = 1, coef1 = 1, coef2 = 1), 
              data = list(suscep_day = human_exposure$suscep_day, sick_day = human_exposure$sick_day ),
              method = "L-BFGS-B",lower = c(intercept1 = 0.001, intercept2 = 0.001, coef1 = 0, coef2 = 0))
summary(fit2)

# Check fit
# function to simulate t_inf from fit value
pred_gamma_ll_1 <- function(suscep_day) {
  mean <-  1.0758e+00 + 1.3799e-01*suscep_day
  var <-  1.0005e+00 + 8.1179e+01*suscep_day # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# simulate with different size from parameters that got from fit 
pred_fit1 <- sapply(human_exposure$suscep_day, FUN= pred_gamma_ll_1 )
plot(pred_fit1)
human_exposure$pred_fit1 <- pred_fit1 
hist(pred_fit1 )
summary(pred_fit1)
hist(human_exposure$sick_day)
summary(human_exposure$sick_day)
# check RMSE between the two fitting
sqrt(sum((human_exposure$sick_day-pred_fit1)^2)/nrow(human_exposure))

plot(human_exposure$suscep_day,human_exposure$sick_day,col = 3)
points( human_exposure$suscep_day,human_exposure$pred_fit1,col = 2)

plot( human_exposure$suscep_day,human_exposure$pred_fit1,col = 2)
points(human_exposure$suscep_day,human_exposure$sick_day,col = 3)


# Try with different coefficients for infectious duration and farm size
gamma_ll_2 <- function(intercept1,intercept2,coef1,coef2, coef3, coef4, suscep_day, detect_time,farm_size,sick_day) {
  
  mean <- intercept1 + coef1*farm_size +coef2*detect_time 
  var <- intercept2 + coef3*farm_size +coef4*detect_time   # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  
  # Negative log-likelihood
  log_likelihood <- sum(dgamma(sick_day, shape = shape, rate = rate, log = TRUE)) 
  return(-log_likelihood) 
}

fit1.1 <-  mle2(gamma_ll_2, start = list(intercept1=1 ,intercept2 = 1 ,coef1 = 3,coef2 = 1, coef3 =1, coef4 = 1), 
                data = list(suscep_day = human_exposure$suscep_day, sick_day = human_exposure$sick_day, detect_time =  human_exposure$detect_time,farm_size =  human_exposure$farm_size),
                method = "L-BFGS-B", lower = c(intercept1 = 0.001, intercept2 = 0.001, coef1 = 0, coef2 = 0, coef3 =0, coef4 = 0))
summary(fit1.1)

# function to simulate prediction from fit value
pred_gamma_ll_1.1 <- function(suscep_day, detect_time, farm_size) {
  mean <-  1.3074e+00 + 1.8588e+00*farm_size +6.9339e+01*detect_time 
  var <-  1.0541e+00 + 1.0685e+03*farm_size +1.7156e+00*detect_time # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# simulate with different size from parameters that got from fit 
pred_fit1.1 <- mapply(pred_gamma_ll_1.1, human_exposure$suscep_day, human_exposure$detect_time,human_exposure$farm_size )
plot(pred_fit1.1)
human_exposure$pred_fit1.1 <- pred_fit1.1 
hist(pred_fit1.1 )
summary(pred_fit1.1)
hist(human_exposure$sick_day)
summary(human_exposure$sick_day)
# check RMSE between the two fitting
sqrt(sum((human_exposure$sick_day-pred_fit1.1)^2)/nrow(human_exposure))

plot(human_exposure$sick_day,col = 3)
points( human_exposure$pred_fit1.1,col = 2)

plot( human_exposure$pred_fit1.1,col = 2)
points(human_exposure$sick_day,col = 3)



# For Duck ---------------------------------------------------------------------
# Calculate cumulative infectious-day animals
data <-  output_duck_fullday

# Explore the relationship
human_exposure <- data %>% dplyr::group_by(farm_size,run) %>% summarise(sick_day = sum(I)*0.1,# *0.1 from deltat
                                                                        suscep_day = max(farm_size) * max(time),
                                                                        farm_size = max(farm_size),
                                                                        detect_time = max(time))

plot(human_exposure$sick_day)
plot(human_exposure$farm_size, human_exposure$sick_day)
plot(human_exposure$detect_time ,human_exposure$sick_day)
plot(human_exposure$suscep_day ,human_exposure$sick_day)

human_exposure$prop_day <- human_exposure$sick_day/human_exposure$suscep_day
plot(human_exposure$prop_day )

# the proportion of sick animal day and suscep animal day is high in small farms 
# because the number of index case is high (10) compared to farm size
# Let exclude the small farm size
human_exposure <- human_exposure[human_exposure$farm_size >= 100,]

# Try take log farm size and sick animal day to see relationship
human_exposure$log_farm_size <- log(human_exposure$farm_size)
human_exposure$log_detect_time <- log(human_exposure$detect_time)
human_exposure$log_suscep_day <- log(human_exposure$suscep_day)

plot(human_exposure$log_farm_size,human_exposure$sickanimal_day)
plot(human_exposure$log_detect_time ,human_exposure$sickanimal_day)
plot(human_exposure$log_suscep_day,human_exposure$sickanimal_day)

summary(human_exposure)
# fit gamma distribution -------------------------
gamma_ll_1 <- function(intercept1,intercept2,coef1,coef2,  suscep_day, sick_day) {
  
  mean <- intercept1 + coef1*suscep_day 
  var <- intercept2 + coef2*suscep_day  # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  
  # Negative log-likelihood
  log_likelihood <- sum(dgamma(sick_day, shape = shape, rate = rate, log = TRUE)) 
  return(-log_likelihood) 
}



# Fit distribution
fit3 <-  mle2(gamma_ll_1, start = list(intercept1 = 1, intercept2 = 1, coef1 = 1, coef2 = 1), 
              data = list(suscep_day = human_exposure$suscep_day, sick_day = human_exposure$sick_day ),
              method = "L-BFGS-B",lower = c(intercept1 = 0.001, intercept2 = 0.001, coef1 = 0, coef2 = 0))
summary(fit3)

# Check fit
# function to simulate t_inf from fit value
pred_gamma_ll_1 <- function(suscep_day) {
  mean <-  2.2198e+03 + 4.8221e-01*suscep_day
  var <-  7.4070e-01 + 9.1822e+01*suscep_day # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# simulate with different size from parameters that got from fit 
pred_fit1 <- sapply(human_exposure$suscep_day, FUN= pred_gamma_ll_1 )
plot(pred_fit1)
human_exposure$pred_fit1 <- pred_fit1 
hist(pred_fit1 )
summary(pred_fit1)
hist(human_exposure$sick_day)
summary(human_exposure$sick_day)
# check RMSE between the two fitting
sqrt(sum((human_exposure$sick_day-pred_fit1)^2)/nrow(human_exposure))

plot(human_exposure$sick_day,col = 3)
points( human_exposure$pred_fit1,col = 2)

plot(human_exposure$pred_fit1,col = 2)
points(human_exposure$sick_day,col = 3)




