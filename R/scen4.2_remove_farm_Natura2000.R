# Simulate Scen 4.2 remove farms around Natura2000 --------------

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
# create sf object from farm 2024
farm_2024_sf <- st_as_sf(farm_2024, coords = c("x", "y"), crs = 28992)

# import natura2000
natura2000 <- read_sf("./data/Natura2000/Natura2000_NL.shp")

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




# Select farm removal around Natura 2000 ---------------------------------------
buff_dist <- c( 0, 1000) # set distance area from Natura 2000 in meter 

# at 21 km all farms are within the buffer
farm_in_buffer <- data.frame (buff_dist = buff_dist/1000, no_farm = rep(NA,length(buff_dist)) )
intersect_farm <- list()

for( i in 1:length(buff_dist)){
  nature_buffer <- sf::st_buffer(natura2000, dist = buff_dist[i])
  # select farms within buffer zone
  intersect_farm[[i]] <- st_intersection(farm_2024_sf, nature_buffer)
 }




# Simulate outbreak ------------------------------------------------------------

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

# Run model --------------
for(j in 1:length(intersect_farm)){
  result <- list()
  source("./R/event.R")
  infect <- farm_df$infect
  infect[intersect_farm[[j]]$host_id] <- 0
  
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
  # percent remove layer farms
  per_thin_layer <- ceiling(sum(farm_df$host_id %in% c(intersect_farm[[j]]$host_id) & farm_df$type =="LAYER")/sum(farm_df$type =="LAYER")*100)
  results$per_thin_layer  <- per_thin_layer 
  # percent remove broiler farms
  per_thin_broiler <- ceiling(sum(farm_df$host_id %in% c(intersect_farm[[j]]$host_id) & farm_df$type =="BROILER")/sum(farm_df$type =="BROILER")*100)
  results$per_thin_broiler <- per_thin_broiler
  results$scen <- "scen4.2"
  results$cluster <- "buffer"
  results$dist <- buff_dist[j]
  results$index <- "high"
  saveRDS(results, paste0("./results/scen4.2/","Sim4.2_layer",per_thin_layer,"_broiler", per_thin_broiler,"_high",".rds" )) 
  gc()
}

