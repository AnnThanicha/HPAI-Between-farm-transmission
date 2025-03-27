# Simulate Scen 2.1 randomly transition to outdoor farms -----------------------

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


## Function transmission kernel 
# Calibrated parameters
h0 <- 0.0005741044 # calibrated H0
r0  <- 2.5
alpha  <- 2.2
h_kernel <- function(r){h0/(1 + (r/r0)^alpha)} ; # transmission kernel as a function of r
beta<-1
theta = 7490
numsim <- 1000 # number of iteration




# Randomly select the farms transit to outdoor----------------------------------
# Get the id of layer and broiler farms beforehand to save time for sampling
layer_id <- farm_df$host_id[farm_df$type =="LAYER"]
layer_id_outdoor <- farm_df$host_id[farm_df$type =="LAYER" & farm_df$system =="Outdoor"]
broiler_id <- farm_df$host_id[farm_df$type =="BROILER"]
broiler_id_outdoor <- farm_df$host_id[farm_df$type =="BROILER" & farm_df$system =="Outdoor"]

# Set percentge of outdoor transition ------------------
per_od_layer <- c(0,20,40, 47.5,60,80,100) # percent outdoor  layer farms
per_od_broiler <-c(6.5,100) # percent outdoor broiler farms
per_transit <- expand.grid(per_od_layer,per_od_broiler) # combination of these two % transition

set.seed(1)
for(j in 1:nrow(per_transit)){
  
  transit_suscep <- list()
  per_transit_layer <- per_transit[j,1]
  per_transit_broiler <- per_transit[j,2]
  
  for(i in 1:numsim){
    suscep <- farm_df$suscep
    # For layer farms 
    # in case 0% outdoor layer farms, all layer farms change to indoor
    if (per_transit_layer==0){
      id_change <- farm_df$host_id[farm_df$type =="LAYER" & farm_df$system =="Outdoor" ]
      suscep[id_change] <- 1
    }
    
    # in case per_transit_layer>0 & per_transit_layer < 47.5
    if (per_transit_layer>0 & per_transit_layer < 47.5 ){
      id_change <- sample(farm_df$host_id[farm_df$type =="LAYER" & farm_df$system =="Outdoor" ], size = length(layer_id_outdoor)- round(0.25 *length(layer_id)), replace =FALSE)
      suscep[id_change] <- 1
    }
    
    # in case per_transit_layer == 47.5 (baseline)
    if (per_transit_layer== 47.5 ){}
    
    # in case per_transit_layer > 47.5
    if (per_transit_layer > 47.5){
      id_change <- sample(farm_df$host_id[farm_df$type2 =="LAYER" & farm_df$system =="Indoor" ], size = round(per_transit_layer *length(layer_id)/100)- length(layer_id_outdoor), replace =FALSE)
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
      id_change <- farm_df$host_id[farm_df$type2 =="BROILER" & farm_df$system =="Outdoor" ]
      suscep[id_change] <- 0.134
    }
    
    if (per_transit_broiler == 6.5 ){}
    
    # in case 25% , random half of outdoor layer farms to be indoor
    if (per_transit_broiler > 6.5 & per_transit_broiler < 100){
      id_change <- sample(farm_df$host_id[farm_df$type2 =="BROILER" & farm_df$system =="Indoor" ], size = round(per_transit_broiler *length(broiler_id)/100)- length(broiler_id_outdoor), replace =FALSE)
      suscep[id_change] <- 0.8442
    }
    
    
    # in case 100% layer farms, random indoor layer farms to be outdoor
    if (per_transit_broiler == 100){
      id_change <- farm_df$host_id[farm_df$type2 =="BROILER" & farm_df$system =="Indoor" ]
      suscep[id_change] <- 0.8442
    }
    
    transit_suscep[[i]] <- suscep
  }
  # save the randomly outdoor farms
  saveRDS(transit_suscep, paste0('./data/scen2.1/','scen2.1','_layer',per_transit[j,1],'_broiler',per_transit[j,2],'.rds'))
}



# Outbreak simulation ----------------------------------------------------------
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


# Run model --------------
files <- list.files(path ="./data/scen1")

source("./R/event.R")

for(j in 1:length(files)){
  result <- list()
  infect <- farm_df$infect
  suscep <- readRDS(paste0("./data/scen1/",files[j]))
  
  
  for(i in 1:numsim) {
    # read suscep from randomly transit matrix
    het_matrix <-  as.matrix(abs(outer(suscep[[i]],farm_df$infect,"*")))
    
    # Multiplied
    hazardmatrix<- hazardmatrix_before * het_matrix * size_matrix
    diag(hazardmatrix) <- 0 # because the chance of infecting itself is 0
    
    # take T_inf and Q_inf
    T_inf <- T_inf_matrix[i,]   # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
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
  # Summarise outcomes 
  # bind list
  results <- bind_rows(result, .id = "iter")
  per_od_layer <- as.numeric(sub(".*_layer([0-9]+\\.?[0-9]*)_.*", "\\1", files[j]))
  per_od_broiler  <- as.numeric(sub(".*_broiler([0-9]+\\.?[0-9]*).*", "\\1", files[j]))
  results$per_od_layer <- per_od_layer
  results$per_od_broiler <- per_od_broiler 
  results$scen <- "1.1"
  results$cluster <- "random"
  results$index <- "high"
  saveRDS(results, paste0("./results/scen2.1/","sim2.1_layer",per_od_layer,"_broiler",per_od_broiler,"_high",".rds" )) 
  gc()
}


# read results
results <- readRDS("./results/scen2.1/Sim2.1_layer0_broiler6.5_high.rds")
# number of infected farms
n_infected <- results %>% filter(Type_event ==2) %>% group_by(iter) %>% summarise(n=n())
hist(n_infected$n)
summary(n_infected$n)
# outbreak duration
duration <- results %>% group_by(iter) %>% summarise(duration=max(Event_time))
hist(duration$duration)
summary(duration$duration)


