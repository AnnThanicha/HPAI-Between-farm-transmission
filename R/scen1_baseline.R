# Simulate Scen 1 b -----------------------

# First specify the packages of interest
packages = c("ggplot2", "sf", "sp","raster", "readxl",
             "dplyr", "RandomFields","RColorBrewer","ggridges", "cowplot")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Prepare data -----------------------------------------------------------------
## Import data 
farm_2024 <- readRDS("./data/farm_gen.rds")
Q_init_matrix <- readRDS("./data/Q_init_matrix.rds")
T_inf_matrix <- readRDS("./data/T_inf_matrix.rds")
distancematrix <- readRDS("./data/distancematrix.rds")

## Function transmission kernel 
# Calibrated parameters
h0 <- 0.0005741044 # calibrated H0
r0  <- 2.5
alpha  <- 2.2
h_kernel <- function(r){h0/(1 + (r/r0)^alpha)} ; # transmission kernel as a function of r
beta<-1
theta = 7490
numsim <- 1000 # number of iteration

# Run simulation for baseline with high density index case ---------------------
farm_df <- farm_2024

# Matrix of susceptibility and infectivity
het_matrix <-  as.matrix(abs(outer(farm_df$suscep,farm_df$infect,"*")))

# create an hazard matrix evaluating for each host j the chance to be infected by host i as a function of distance 
hazardmatrix_before  <- as.matrix(apply(distancematrix,MARGIN=c(1,2),FUN=h_kernel))

# Create matrix of farm size
farm_df$mod_size <- 1-exp(-1*farm_df$size/theta)
size_matrix <- as.matrix(abs(outer(farm_df$mod_size,farm_df$mod_size,"*")))
diag(size_matrix) <- 0



result <- list()
source("./r/event.R")

# Start for-loop of outbreak simulation
for(i in 1:numsim) {
  
  # Multiplied
  hazardmatrix<- hazardmatrix_before * het_matrix * size_matrix
  diag(hazardmatrix) <- 0 # because the chance of infecting itself is 0
  
  # take T_inf and Q_inf
  T_inf <- T_inf_matrix[i,]   # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
  Q_init <- Q_init_matrix[i,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1
  # Define index case 
  K <- sample(which(farm_df$index=="high"), size = 1) # random index case within high density area
  
  source("./r/InitSim.R") # initialization and setting of the first infected
  while(nrow(Queue)!=0){
    source("./r/Simloop.R")
    result[[i]] <- History
    
  }
  print(i)
  
}

# Summarise result  ------------------------------------------------------------
results <- bind_rows(result, .id = "iter")
results$per_od_layer <- 47.5 # %outdoor layer farms
results$per_od_broiler <- 6.5  # %outdoor broiler farms
results$scen <- "Baseline" # Scenario name
results$cluster <- "Baseline"
results$index <- "high"


# number of infected farms
n_infected <- results  %>% filter(Type_event ==2) %>% group_by(iter) %>% summarise(n=n())
hist(n_infected$n)
summary(n_infected$n)
# outbreak duration
duration <- results %>% group_by(iter) %>% summarise(duration=max(Event_time))
hist(duration$duration)
summary(duration$duration)

