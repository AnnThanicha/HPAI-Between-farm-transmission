# Prepare Q_init_matrix and T_inf_matrix----------------------------------------
# Create Q_init_matrix and T_inf_matrix to use for outbreak  simulation
# Read https://doi.org/10.1371/journal.pcbi.1008009

# import farm data
farm_2024 <- readRDS("./data/farm_gen.rds")

# Prepare tinf and Qinf
# Random Q_inf
numsim <-1000 # no. of iteration
totpoints <- nrow(farm_2024)
set.seed(001)
Q_init_matrix <- matrix(0,nrow=numsim,ncol=totpoints)
for (ii in 1:numsim){
  Q_init_matrix[ii,] <- rexp(totpoints, rate = 1) # random threshold to infection
}

saveRDS(Q_init_matrix,"./data/Q_init_matrix.rds")


## Create T_inf from fitted model of within farm model 
# function for T_inf of chicken 
t_inf_chick <- function(log_size) {
  mean <-1.0516315 + 0.3599556*log_size + 0.0602508*log_size^2
  var <-  0.1954552 # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# function for T_inf of turkey
t_inf_turkey <- function(log_size) {
  mean <-1.1439855 + 0.8948005*log_size + 0.0365499*log_size^2
  var <-  0.2061337 # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# function for T_inf of duck
t_inf_duck <- function(log_size) {
  mean <-2.93114591 + 0.46725553*log_size 
  var <-  0.01038639*log_size # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# Indicate types
chick_id <- which(farm_2024$type%in% c("LAYER","LAYER_BREEDER" ,"BROILER","BROILER_BREEDER"))
turkey_id <- which(farm_2024$type%in% c("TURKEY"))
duck_id <- which(farm_2024$type%in% c("DUCK"))
# Create matrix of T_inf
T_inf_matrix <- matrix(0,nrow=numsim,ncol=totpoints)
for (ii in 1:numsim){
  T_inf_matrix[ii,chick_id] <- sapply(log(farm_2024$size[chick_id]), FUN= t_inf_chick) 
  T_inf_matrix[ii,turkey_id] <- sapply(log(farm_2024$size[turkey_id]), FUN= t_inf_turkey) 
  T_inf_matrix[ii,duck_id] <- sapply(log(farm_2024$size[duck_id]), FUN= t_inf_duck) 
}


saveRDS(T_inf_matrix,"./data/T_inf_matrix.rds") 
