# Simulate Scen 2.3 transition to outdoor farms around Natura2000 --------------

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


# Define farms within the certain distance of Natura 2000-----------------------
buff_dist <- c( 250 ,500 ,1000, 2500, 5000, 10000, 15000) # set buffer area in meter 
# at 21 km all farms are within the buffer
farm_in_buffer <- data.frame (buff_dist = buff_dist/1000, no_farm = rep(NA,7) )
intersect_farm <- list()

for( i in 1:length(buff_dist)){
  nature_buffer <- sf::st_buffer(natura2000, dist = buff_dist[i])
  intersect_farm[[i]] <- st_intersection(farm_2024_sf, nature_buffer)
}

saveRDS(intersect_farm, "./data/scen2.3/scen2.3_intersect_farm.rds")



# Run model --------------------------------------------------------------------
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




for(j in 1:length(intersect_farm)){
  
  farm_df <- farm_2024
  transit_farm <- intersect_farm[[j]] # intersect_farm
  farm_df$suscep[farm_df$host_id %in% transit_farm$host_id & farm_df$type=="LAYER" ] <- 6.3
  farm_df$system[farm_df$host_id %in% transit_farm$host_id & farm_df$type=="LAYER" ] <- "Outdoor"
  farm_df$suscep[farm_df$host_id %in% transit_farm$host_id & farm_df$type=="BROILER" ] <- 0.8442
  farm_df$system[farm_df$host_id %in% transit_farm$host_id & farm_df$type=="BROILER" ] <- "Outdoor"
  
  suscep <- farm_df$suscep
  infect <- farm_df$infect
  het_matrix <-  as.matrix(abs(outer(suscep,infect,"*")))
  diag(het_matrix) <- 0
  
  result <- list()
  source("./R/event.R")
  
  for(i in 1:numsim) {
    
    
    # Multiplied
    hazardmatrix<- hazardmatrix_before * het_matrix * size_matrix
    diag(hazardmatrix) <- 0 # because the chance of infecting itself is 0
    
    # take T_inf and Q_inf
    T_inf <- T_inf_matrix[i,]   
    Q_init <- Q_init_matrix[i,] 
    # Define index case 
    K <- sample(which(farm_df$index=="high"), size = 1) # index case from high density
    
    source("./R/InitSim.R") # initialization and setting of the first infected
    while(nrow(Queue)!=0){
      source("./R/Simloop.R")
      result[[i]] <- History
      
    }
    print(i)
  }
  # Summarise outcomes
  results <- bind_rows(result, .id = "iter")
  per_od_layer <- round(sum(farm_df$type=="LAYER"&farm_df$system=="Outdoor")/sum(farm_df$type=="LAYER")*100)
  per_od_broiler  <-  round(sum(farm_df$type=="BROILER"&farm_df$system=="Outdoor")/sum(farm_df$type=="BROILER")*100)
  results$per_od_layer <- per_od_layer
  results$per_od_broiler <- per_od_broiler 
  results$scen <- "2.3"
  results$cluster <- "buffer"
  results$index <- "high"
  results$buff_dist <- buff_dist[j]
  saveRDS(results, paste0("./results/scen2.3/","sim2.3_buffer",buff_dist[j],"m_high",".rds" )) 
  gc()
}
