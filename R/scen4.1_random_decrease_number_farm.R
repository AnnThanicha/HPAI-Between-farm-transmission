# Simulate Scen 4.1 randomly decrease number of farms--------------------------------------

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



# Randomly select farms fro removal --------------------------------------------
farm_df <- farm_2024

# Change all layer and broiler farms to outdoor
farm_df$suscep[farm_df$type2 == "BROILER"] <- "0.8442"
farm_df$suscep[farm_df$type2 == "LAYER"] <- "6.3"
summary(as.factor(farm_df$suscep))

# We will thin out farms, by setting infectivity of the farms to 0
# First let's randomly select farms
per_thin_layer <- c(0,25,50,75,100) # percent layer farms removal
per_thin_broiler <-c(0,50,100) # percent broiler farms removal
per_thin_out <- expand.grid(per_thin_layer, per_thin_broiler) # combination of these two % transition


for(j in 1:nrow(per_thin_out)){
  
  rm_farm <- list()
  per_layer <- per_thin_out[j,1]
  per_broiler <- per_thin_out[j,2]
  
  for(i in 1:numsim){
    infect <- farm_df$infect
    suscep <- farm_df$suscep
    rm_layer <- sample(which(suscep == 6.3), size = (sum(suscep == 6.3)*per_layer)/100, replace = FALSE )
    rm_broiler <- sample(which(suscep == 0.8442), size = (sum(suscep == 0.8442)*per_broiler)/100, replace = FALSE )
    infect[c(rm_layer,rm_broiler)] <- 0
    rm_farm[[i]] <- infect
  }
  
  saveRDS(rm_farm, paste0("./data/scen4.1/",'rm6.1','_layer',per_thin_out[j,1],'_broiler',per_thin_out[j,2],'.rds'))
}



# Run outbreak simulation ------------------------------------------------------

farm_df <- farm_2024

# Express the coordinates as complex numbers for fast calculation of the euclidean distance
matrix_points <- farm_df[,c("x","y")]
totpoints <- nrow(matrix_points) # total number of points
colnames(matrix_points) <- c("xcoord","ycoord")
# add a column for the index
Index_points <- c(1:totpoints)
Coord <- (complex(length.out=2,real=matrix_points$xcoord,imaginary=matrix_points$ycoord));
distancematrix <- as.matrix(abs(outer(Coord,Coord,"-")))
distancematrix <- distancematrix/1000 # change to km

# create an hazard matrix evaluating for each host j the chance to be infected by host i as a function of distance 
hazardmatrix_before  <- as.matrix(apply(distancematrix,MARGIN=c(1,2),FUN=h_kernel))


# Create farm size matrix
farm_df$mod_size <- 1-exp(-1*farm_df$size/theta)
size_matrix <- as.matrix(abs(outer(farm_df$mod_size,farm_df$mod_size,"*")))
diag(size_matrix) <- 0

# Change all layer and broiler farms to outdoor
farm_df$suscep[farm_df$type2 == "BROILER"] <- "0.8442"
farm_df$suscep[farm_df$type2 == "LAYER"] <- "6.3"
suscep_outdoor <- as.numeric(as.character(farm_df$suscep))
het_matrix <-  as.matrix(abs(outer(suscep_outdoor,farm_df$infect,"*")))

# Multiplied hazard matrix
hazardmatrix_before_thin <- hazardmatrix_before * het_matrix * size_matrix
diag(hazardmatrix_before_thin) <- 0

# list files of randomly remove farms that we generated before
files <- list.files("./data/scen4.1")

for(j in 1:length(files)){
  result <- list()
  source("./R/event.R")
  infect <- readRDS(paste0("./data/scen4.1/",files[j]))
  
  
  for(i in 1:numsim) {
    # thin out matrix
    thin_martrix <- outer(infect[[i]], infect[[i]], "*")
    
    # With censor matrix
    hazardmatrix <- hazardmatrix_before_thin * thin_martrix
    diag(hazardmatrix) <- 0
    
    diag(hazardmatrix) <- 0 # because the chance of infecting itself is 0
    
    # take T_inf and Q_inf
    T_inf <- T_inf_matrix[i,]   # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
    Q_init <- Q_init_matrix[i,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1
    # Define index case 
    K <- sample(which(farm_df$index=="high" & infect[[i]] !=0), size = 1) # random index case, but not choosing thin out farms
    
    source("./R/InitSim.R") # initialization and setting of the first infected
    while(nrow(Queue)!=0){
      source("./R/Simloop.R")
      result[[i]] <- History
      
    }
    print(i)
  }
  # Summarise outcomes -------
  # bind list
  results <- bind_rows(result, .id = "iter")
  per_thin_layer <- as.numeric(sub(".*_layer([0-9]+)_.*", "\\1", files[j]))
  per_thin_broiler  <- as.numeric(sub(".*_broiler([0-9]+)\\.rds", "\\1", files[j]))
  results$per_thin_layer <- per_thin_layer
  results$per_thin_broiler <- per_thin_broiler
  results$scen <- "scen4.1"
  results$cluster <- "random"
  results$index <- "high"
  saveRDS(results, paste0("./results/scen4.1/","Sim4.1_layer",per_thin_layer,"_broiler", per_thin_broiler,"_high",".rds" )) 
  gc()
}

# Simulate outbreak high index-------------------------

farm_df <- farm_024

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

# Create matrix for farm size
theta = 7490
# Create farm size matrix
farm_df$mod_size <- 1-exp(-1*farm_df$size/theta)
size_matrix <- as.matrix(abs(outer(farm_df$mod_size,farm_df$mod_size,"*")))
diag(size_matrix) <- 0

numsim <- 1000

# Now 47.5% of layer farms is outdoor and 6.5% of broiler is outdoor, this will be baseline
# if we calculate from the (outdoor broilers+layers)/(total broilers+layers) = (46+421)/1598 = 0.29



# Change all layer and broiler farms to outdoor
farm_df$suscep[farm_df$type2 == "BROILER"] <- "0.8442"
farm_df$suscep[farm_df$type2 == "LAYER"] <- "6.3"
suscep_outdoor <- as.numeric(as.character(farm_df$suscep))
het_matrix <-  as.matrix(abs(outer(suscep_outdoor,farm_df$infect,"*")))

# Multiplied hazard matrix
hazardmatrix_before_thin <- hazardmatrix_before * het_matrix * size_matrix
diag(hazardmatrix_before_thin) <- 0

# Run model --------------
for(j in 1:length(Scen6.2_intersect_farm)){
  result <- list()
  source("event.R")
  infect <- farm_df$infect
  infect[Scen6.2_intersect_farm[[j]]$host_id] <- 0
  
  # thin out matrix
  thin_martrix <- outer(infect, infect, "*")
  
  # With censor matrix
  hazardmatrix <- hazardmatrix_before_thin * thin_martrix
  diag(hazardmatrix) <- 0 # because the chance of infecting itself is 0
  
  for(i in 1:numsim) {
    
    # take T_inf and Q_inf
    T_inf <- T_inf_matrix[i,]   # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
    Q_init <- Q_init_matrix[i,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1
    # Define index case 
    K <- sample(which(farm_df$index=="high" & infect !=0), size = 1) # random index case, but not choosing thin out farms
    
    source("InitSim.R") # initialization and setting of the first infected
    while(nrow(Queue)!=0){
      source("Simloop.R")
      result[[i]] <- History
      
    }
    print(i)
  }
  # Summarise outcomes -------
  # bind list
  results <- bind_rows(result, .id = "iter")
  per_thin_layer <- ceiling(sum(farm_df$host_id %in% c(Scen6.2_intersect_farm[[j]]$host_id) & farm_df$type2 =="LAYER")/sum(farm_df$type2 =="LAYER")*100)
  results$per_thin_layer  <- per_thin_layer 
  per_thin_broiler <- ceiling(sum(farm_df$host_id %in% c(Scen6.2_intersect_farm[[j]]$host_id) & farm_df$type2 =="BROILER")/sum(farm_df$type2 =="BROILER")*100)
  results$per_thin_broiler <- per_thin_broiler
  results$scen <- "scen6.2"
  results$cluster <- "buffer"
  results$dist <- buff_dist[j]
  results$index <- "high"
  saveRDS(results, paste0("Sim6.2_layer",per_thin_layer,"_broiler", per_thin_broiler,"_high",".rds" )) 
  gc()
}
