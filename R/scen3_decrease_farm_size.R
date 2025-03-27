# Simulate Scen 3 tfarm size reduction --------------

# First specify the packages of interest
packages = c("ggplot2", "sf", "sp","raster", "readxl","cluster","factoextra",
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


## Function transmission kernel 
# Calibrated parameters
h0 <- 0.0005741044 # calibrated H0
r0  <- 2.5
alpha  <- 2.2
h_kernel <- function(r){h0/(1 + (r/r0)^alpha)} ; # transmission kernel as a function of r
beta<-1
theta = 7490
numsim <- 1000 # number of iteration
totpoints <- nrow(farm_2024)



# ## Simulate model ------------------------------------------------------------
# We will create the situation when all farms are outdoor 
#change all farm to outdoor
farm_df <- farm_2024


# Express the coordinates as complex numbers for fast calculation of the euclidean distance
matrix_points <- farm_df[,c("X","Y")]
totpoints <- nrow(matrix_points) # total number of points
colnames(matrix_points) <- c("xcoord","ycoord")
# add a column for the index
Index_points <- c(1:totpoints)
Coord <- (complex(length.out=2,real=matrix_points$xcoord,imaginary=matrix_points$ycoord));
distancematrix <- as.matrix(abs(outer(Coord,Coord,"-")))
distancematrix <- distancematrix/1000 # change to km

# create an hazard matrix evaluating for each host j the chance to be infected by host i as a function of distance 
hazardmatrix_before  <- as.matrix(apply(distancematrix,MARGIN=c(1,2),FUN=h_kernel))


# To make it comparable with scen 2.1, we use the random farms that we selected 
# from scen 2.1
files <- list.files(path ="./data/scen2.1")

for(j in 1:length(files)){
  result <- list()
  source("./R/event.R")
  suscep <- readRDS(paste0("./data/scen2.1/",files[j]))
  
  
  for(i in 1:numsim) {
    
    # create matrix for susceptibility
    het_matrix <-  as.matrix(abs(outer( suscep[[i]], farm_df$infect,"*")))
    
    # Create reduce farm size
    size <-  farm_df$size
    size[farm_df$type2 =="LAYER" & suscep[[i]] == 6.3] <- ceiling(size[farm_df$type2 =="LAYER" & suscep[[i]] == 6.3] * 0.65)
    size[farm_df$type2 =="BROILER" & suscep[[i]] == 0.8442] <- ceiling(size[farm_df$type2 =="BROILER" & suscep[[i]] == 0.8442] * 0.5)
    mod_size <- 1-exp(-1*size/theta) # size matrix calculate from reduce size
    size_matrix <- as.matrix(abs(outer(mod_size,mod_size,"*")))
    diag(size_matrix) <- 0
    
    
    hazardmatrix <- hazardmatrix_before * size_matrix * het_matrix
    diag(hazardmatrix) <- 0 # because the chance of infecting itself is 0
    
    # take T_inf and Q_inf
    T_inf <- T_inf_matrix[i,]   # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
    Q_init <- Q_init_matrix[i,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1
    # Define index case 
    K <- sample(which(farm_df$index=="high"), size = 1) # random index case, but not choosing thin out farms
    
    source("./R/InitSim.R") # initialization and setting of the first infected
    while(nrow(Queue)!=0){
      source("./R/Simloop.R")
      result[[i]] <- History
      
    }
    print(i)
  }
  
  # Summarise outcomes
  results <- bind_rows(result, .id = "iter")
  per_od_layer <- as.numeric(sub(".*_layer([0-9]+\\.?[0-9]*)_.*", "\\1", files[j]))
  per_od_broiler  <- as.numeric(sub(".*_broiler([0-9]+\\.?[0-9]*).*", "\\1", files[j]))
  results$per_od_layer <- per_od_layer
  results$per_od_broiler <- per_od_broiler 
  results$scen <- "3"
  results$cluster <- "Random"
  results$index <- "high"
  saveRDS(results, paste0("./results/scen3/","sim3","_layer",per_od_layer,"_broiler",per_od_broiler,"_high",".rds" )) 
  gc()
}
