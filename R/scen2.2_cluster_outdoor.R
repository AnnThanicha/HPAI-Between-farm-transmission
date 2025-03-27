# Simulate Scen 2.2 clustered transition to outdoor farms -----------------------

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


# Define the number of clusters ------------------------------------------------
# Find the cluster of outdoor farms both layer and broiler
# when the percent outdoor decrease, the outdoor farms that are the furthest from medoids change first
# when the percent outdoor increase, the outdoor farms that are the closet from medoids change first

cluster_farm <- data.frame(host_id = 1:nrow(farm_2024), X= farm_2024$x, Y = farm_2024$y, type = farm_2024$type, size = farm_2024$size, system = farm_2024$system)
cluster_farm_outdoor <- cluster_farm [cluster_farm$system == "Outdoor",] # only select layer farms

# choose the number of clusters
fviz_nbclust(cluster_farm_outdoor[, c("X","Y")], pam, method = "wss") # the elbow points at  4 clusters
fviz_nbclust(cluster_farm_outdoor[, c("X","Y")], pam, method = "silhouette") #  k means silhouette score maximum at 4
#calculate gap statistic based on number of clusters
gap_stat <- clusGap(cluster_farm_outdoor[, c("X","Y")],
                    FUN = pam,
                    K.max = 8, #max clusters to consider
                    B = 50) #total bootstrapped iterations

# plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat) # gap statistic is highest at 4 before plateau

# Perform k-medoids and calculate distance to medoids---------------------------
# Here, we used number of cluster (k) = 4, 
# but you should decide based on the gap statistics and shoulder test above
set.seed(1)
#perform k-medoids clustering with k = 4 clusters
kmed <- pam(cluster_farm_outdoor[, c("X","Y")], k = 4)

# plot clusters
fviz_cluster(kmed, data = cluster_farm_outdoor[, c("X","Y")])

# get the medoids to sf
medoids_sf <- st_as_sf(data.frame(kmed$medoids), coords = c("X", "Y"), crs = 28992)

# Calculate the distance from farms to these four medoids --> the closet medoids is the cluster of that farms
# change to sf object to calculate distance
farm_2024_sf <- st_as_sf(farm_2024, coords = c("x", "y"), crs = 28992)

farm_2024$med1 <- as.numeric(st_distance(farm_2024_sf, medoids_sf[1,]))
farm_2024$med2 <- as.numeric(st_distance(farm_2024_sf, medoids_sf[2,]))
farm_2024$med3 <- as.numeric(st_distance(farm_2024_sf, medoids_sf[3,]))
farm_2024$med4 <- as.numeric(st_distance(farm_2024_sf, medoids_sf[4,]))

# Define cluster based on shortest distance to a medoid
farm_2024$group <- apply(farm_2024[, c("med1","med2","med3","med4")], 1, function(row) {which.min(row)})
farm_2024$dist_medoid <- apply(farm_2024[, c("med1","med2","med3","med4")], 1, function(row) {min(row)})




# Create susceptibility with different % of outdoor system----------------------
farm_df <- farm_2024

# Now 47.5% of layer farms is outdoor and 6.5% of broiler is outdoor, this will be baseline

# Get the id of layer and broiler farms beforehand to save time for sampling
layer_id <- farm_df$host_id[farm_df$type =="LAYER"]
layer_id_outdoor <- farm_df$host_id[farm_df$type =="LAYER" & farm_df$system =="Outdoor"]
layer_id_indoor <- farm_df$host_id[farm_df$type =="LAYER" & farm_df$system =="Indoor"]
broiler_id <- farm_df$host_id[farm_df$type =="BROILER"]
broiler_id_outdoor <- farm_df$host_id[farm_df$type =="BROILER" & farm_df$system =="Outdoor"]
broiler_id_indoor <- farm_df$host_id[farm_df$type =="BROILER" & farm_df$system =="Indoor"]

# Create clustered transition
per_od_layer <- c(0,25,50,75,100) # percent outdoor layer farms
per_od_broiler <- c(6.5,100) # percent outdoor broiler farms
per_transit <- expand.grid(per_od_layer,per_od_broiler) # combination of these two % transition


for(j in 1:nrow(per_transit)){
  
  transit_suscep <- list()
  per_transit_layer <- per_transit[j,1]
  per_transit_broiler <- per_transit[j,2]
  
  
  
  suscep <- farm_df$suscep
  # For layer farms 
  # in case 0% outdoor layer farms, all layer farms change to indoor
  if (per_transit_layer==0){
    id_change <- farm_df$host_id[farm_df$type =="LAYER" & farm_df$system =="Outdoor" ]
    suscep[id_change] <- 1
  }
  
  
  # in case 25% outdoor layer farms, reduce outdoor farm that furthest from the cluster
  if (per_transit_layer> 0 & per_transit_layer < 47.5){
    size <- length(layer_id_outdoor)- round(per_transit_layer *length(layer_id)/100)
    dist_threshold <- min(sort(farm_df$dist_medoid[layer_id_outdoor], decreasing = TRUE)[1:size]) # the outdoor farm above this threshold will be change to indoor
    id_change <- farm_df$host_id[farm_df$host_id %in% layer_id_outdoor & farm_df$dist_medoid>= dist_threshold]
    suscep[id_change] <- 1
  }
  
  # if per_transit_layer = 47.5 (baseline) do nothing
  if (per_transit_layer == 47.5) {}
  
  # in case 50%,75% outdoor layer farms, the shortest distance indoor layer farms to be outdoor
  if (per_transit_layer > 47.5){
    size <- round(per_transit_layer *length(layer_id)/100)- length(layer_id_outdoor)
    dist_threshold <- max(sort(farm_df$dist_medoid[layer_id_indoor], decreasing = FALSE)[1:size]) # the indoor farm below this threshold will be change to outdoor
    id_change <- farm_df$host_id[farm_df$host_id %in% layer_id_indoor & farm_df$dist_medoid<= dist_threshold]
    suscep[id_change] <- 6.3
  }
  
  # in case 100% outdoor layer farms, random indoor layer farms to be outdoor
  if (per_transit_layer == 100){
    id_change <- farm_df$host_id[farm_df$type2 =="LAYER" & farm_df$system =="Indoor" ]
    suscep[id_change] <- 6.3
  }
  
  #For broiler farms
  # in case 0%, all broiler farms change to indoor
  if (per_transit_broiler==0){
    id_change <- farm_df$host_id[farm_df$type =="BROILER" & farm_df$system =="Outdoor" ]
    suscep[id_change] <-  0.134
  }
  
  # if per_transit_layer = 47.5 (baseline) do nothing
  if (per_transit_broiler == 6.5) {}
  
  # in case 25%,50% and 75% , choose the closet indoor broiler farms to be outdoor
  if (per_transit_broiler > 6.5 & per_transit_broiler < 100){
    size <- round(per_transit_broiler *length(broiler_id)/100)- length(broiler_id_outdoor)
    dist_threshold <- max(sort(farm_df$dist_medoid[broiler_id_indoor], decreasing = FALSE)[1:size]) # the outdoor farm below this threshold will be change to indoor
    id_change <- farm_df$host_id[farm_df$host_id %in% broiler_id_indoor & farm_df$dist_medoid<= dist_threshold]
    suscep[id_change] <- 0.8442
  }
  
  
  # in case 100% layer farms, all indoor layer farms to be outdoor
  if (per_transit_broiler == 100){
    id_change <- farm_df$host_id[farm_df$type =="BROILER" & farm_df$system =="Indoor" ]
    suscep[id_change] <- 0.8442
  }
  
  
  saveRDS(suscep, paste0('./data/scen2.2/','scen2.2','_layer',per_transit[j,1],'_broiler',per_transit[j,2],'.rds'))
  
}




# Run simulation ---------------------------------------------------------------
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



# Now 47.5% of layer farms is outdoor and 6.5% of broiler is outdoor, this will be baseline
# if we calculate from the (outdoor broilers+layers)/(total broilers+layers) = (46+421)/1598 = 0.29

files <- list.files(path="./data/scen2.2")



for (j in 1:length(files)){
  
  infect <- farm_df$infect
  suscep <- readRDS(paste0("./data/scen2.2/",files[j]))
  het_matrix <-  as.matrix(abs(outer(suscep,farm_2024$infect,"*")))
  
  # Multiplied
  hazardmatrix<- hazardmatrix_before * het_matrix * size_matrix
  diag(hazardmatrix) <- 0 # because the chance of infecting itself is 0
  
  result <- list()
  source("./R/event.R")
  
  for(i in 1:numsim) {
    
    # take T_inf and Q_inf
    T_inf <- T_inf_matrix[i,]   
    Q_init <- Q_init_matrix[i,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1
    # Define index case 
    K <- sample(which(farm_df$index=="high"), size = 1) # random index case 
    
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
  per_od_layer <- as.numeric(sub(".*_layer([0-9]+\\.?[0-9]*)_.*", "\\1", files[j]))
  per_od_broiler  <- as.numeric(sub(".*_broiler([0-9]+\\.?[0-9]*).*", "\\1",  files[j]))
  results$per_od_layer <- per_od_layer
  results$per_od_broiler <- per_od_broiler 
  results$scen <- "2.2"
  results$cluster <- "cluster"
  results$index <- "high"
  saveRDS(results, paste0("./results/scen2.2/","sim2.2_layer",per_od_layer,"_broiler",per_od_broiler,"_high",".rds" )) 
  gc()
}

# read results
results <- readRDS("./results/scen2.2/sim2.2_layer0_broiler6.5_high.rds")
# number of infected farms
n_infected <- results %>% filter(Type_event ==2) %>% group_by(iter) %>% summarise(n=n())
hist(n_infected$n)
summary(n_infected$n)
# outbreak duration
duration <- results %>% group_by(iter) %>% summarise(duration=max(Event_time))
hist(duration$duration)
summary(duration$duration)
