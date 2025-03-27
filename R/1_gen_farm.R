# Generate farm using log-gaussion cox process and CBS data --------------------

# Due to anonymity, we could not share the real poultry farm data to you.
# But we will generate poultry farm data using log-gaussion cox process and CBS data.
# This generated farm data can be used for outbreak simulation.

# The code for log-gaussion cox process is adapt from: DOI: 10.1371/journal.pcbi.1008009
# The CBS data contains the number of poultry farms in each agricultural areas can be downloaded from: https://opendata.cbs.nl/statline/#/CBS/en/dataset/80783eng/table?ts=1721205355440
 


# Preparation ------------------------------------------------------------------
# First specify the packages of interest
packages = c("ggplot2", "sf", "sp","lwgeom", "spatstat.random", "raster", 
             "RandomFields","maptools", "spatstat.geom", "reshape", "dplyr", 
             "spatstat.random", "spatstat", "readxl", "tidyr")

# Load library
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


# Download data
NL <- read_sf("./data/clean_NL3/clean_NL3.shp") 
# NL agricultural area map excluding water, city and forest areas. This is the area that the farms could potentially be located.

farm_census_wide <- read_excel("./data/farm_census_wide.xlsx") # farm data from CBS website
farm_census <- pivot_longer(farm_census_wide, cols = 2:66)# change to long data format



# Generate farm location using log-gaussion cox process-------------------------
# Function to calculate farm location using log-gaussion cox process
# The number of generated farm follows th number of farms in each agricultural area.
# The distribution of farm location can be customized via var, scale and alpha parameters.
# read DOI: 10.1371/journal.pcbi.1008009 for explanation

# Parameters: 
#   var = variance of the random field
#   scale = correlation distance. Note: the effective correlation distance depends on the model
#   alpha = 'smoothness' parameter (0, 2]. alpha = 1 is exponential decay

gen_LGC <- function(scalepar,varpar, farm_census, NL) {
  # to save points output
  point_list <- list() 
  
  # for-loop to generate point -----
  for (i in 1:length(NL$statcode)) {
    
    # if the number of farm = 0 in that province we skip that province
    if (farm_census$value[farm_census$name == NL$statcode[i] & farm_census$Type ==  "Total"]==0) next
    
    polygon <- as.owin(as_Spatial(NL[NL$statcode == NL$statcode[i], ]))
    # Define grids
    grid_size_x <- 100  # Number of grid points in the x direction
    grid_size_y <- 100  # Number of grid points in the y direction
    
    desired_points <- farm_census$value[farm_census$name == NL$statcode[i] & farm_census$Type ==  "Total"] # Desired number of points
    
    
    # Create a grid of points over the bounding box of the polygon
    x <- seq(polygon$xrange[1], polygon$xrange[2], length.out = grid_size_x)
    y <- seq(polygon$yrange[1], polygon$yrange[2], length.out = grid_size_y)
    grid_points <- expand.grid(x = x, y = y)
    
    # Filter points that lie inside the polygon
    inside <- inside.owin(x = grid_points$x, y = grid_points$y, w = polygon)
    grid_points <- grid_points[inside, ]
    
    n.cells <- nrow(grid_points)
    n.points <- desired_points # number of farms
    
    
    # Set spConform to FALSE to enable faster simulation
    # Set seed to NA to make seed arbitrary
    # See help(RFoptions)
    RFoptions(spConform = FALSE, seed = NA)
    
    
    # Covariance model. Whittle-Matern Covariance Model following Bennica et al.
    cov.model <- RMwhittle(var = scalepar, scale = varpar , nu = 1)
    # Calculate the mean intensity 
    # This is the mean number of points per grid cell: n.points/n.cells
    # However, because of the exponent in the rf term, we should divide by the smearing factor: exp(0.5*var)
    mean.intensity <- n.points/n.cells/exp(0.5*cov.model@par.general$var)
    
    result.data  <- grid_points
    result.data <- within(result.data, {
      # Generate a realisation of the random field
      rf <- RFsimulate(model = cov.model, x = x, y = y)
      # Calculate log intensity
      log.intensity <- log(mean.intensity) + rf
      # Intensity (needs to be positive)
      intensity <- exp(log.intensity)
    })
    
    # The sum of all intensities is approximately equal to n.points
    # We take care of this later
    sum(result.data$intensity)
    
    # Simulate inhomogeneous Poisson process------
    # For each grid cell, draw a number of points from a Multinomial distribution with a given probability vector
    # We use the link between the Poisson distribution and the Multinomial distribution to condition on the total n.points
    # The probability = intensity/sum(intensity).
    # Within the grid cell, the location of the points is uniform in space
    result.data <- within(result.data, {
      # The number of points within each cell
      n <- rmultinom(n = 1, size = n.points, prob = intensity/sum(intensity))
    })
    
    # What is the size of cells?
    x.dim <- diff(x)[1]
    y.dim <- diff(y)[1]
    
    # Draw points
    points.list <- with(result.data, mapply(
      x = x, y = y, n = n,
      FUN = function(x, y, n) {
        # Only draw points if n > 0
        if (n > 0) {
          return(cbind(
            x = runif(n = n, min = x - x.dim/2, max = x + x.dim/2),
            y = runif(n = n, min = y - y.dim/2, max = y + y.dim/2)))
        }
      }))
    
    points.mat <- do.call(what = "rbind", args = points.list)
    
    point_list[[i]] <- points.mat
    print(i)
    
  }
  
  point_list <- point_list[lengths(point_list) != 0]
  point_df <- as.data.frame( do.call("rbind", point_list))
  point_df$area <- rep(NL$statcode, farm_census$value[farm_census$Type ==  "Total"])
  
  return(point_df)
}

# Here, we will generated farm location using scale = 16 and var = 16, 
# but  you can change as you want.

farm_gen <- gen_LGC(scalepar = 0.1,varpar = 1, farm_census, NL)

# Add farm type
# Order the data by area first
farm_gen <- farm_gen[order(farm_gen$area), ]
farm_census <- farm_census[order(farm_census$name), ]
farm_census <- farm_census[farm_census$Type!= "Total",]
farm_gen$type <- rep(farm_census$Type, farm_census$value)

# Add farm system (indoor, outdoor)
# From the poultry farm data 47.5% of layer farms are outdoor and 6.5% of broiler farms are outdoor
farm_gen$system <- "Indoor"
farm_gen$system[sample(which(farm_gen$type =="LAYER"), 0.475*sum(farm_gen$type =="LAYER"))] <- "Outdoor"
farm_gen$system[sample(which(farm_gen$type =="BROILER"), 0.065*sum(farm_gen$type =="BROILER"))] <- "Outdoor"
table(farm_gen$type  ,farm_gen$system) # Check data

# Add farm size 
# We random farm size based on size of poultry farm 2024 data 
farm_gen$size[farm_gen$type == "LAYER"] <- round(runif(min = 1, max =500000,n=sum(farm_gen$type == "LAYER")), 0)
farm_gen$size[farm_gen$type == "BROILER"] <- round(runif(min = 1, max =490000,n = sum(farm_gen$type == "BROILER")), 0)                                             
farm_gen$size[farm_gen$type == "BROILER_BREEDER"] <- round(runif(min = 1, max = 115017,n=sum(farm_gen$type == "BROILER_BREEDER")), 0) 
farm_gen$size[farm_gen$type == "LAYER_BREEDER"] <- round(runif(min = 1, max = 133473 ,n=sum(farm_gen$type == "LAYER_BREEDER")), 0) 
farm_gen$size[farm_gen$type == "DUCK"] <- round(runif(min = 1, max = 101500 ,n=sum(farm_gen$type == "DUCK")), 0) 
farm_gen$size[farm_gen$type == "TURKEY"] <- round(runif(min = 1, max =49265  ,n=sum(farm_gen$type == "TURKEY")), 0) 

# Add susceptibility
# Susceptibility follows the parameters for the table in the paper
farm_gen <- farm_gen %>% mutate (suscep = case_when( type == "LAYER" & system == "Outdoor" ~ 6.3,
                                                     type == "LAYER" & system == "Indoor" ~ 1,
                                                     type == "BROILER" & system == "Outdoor" ~ 0.84,
                                                     type == "BROILER" & system == "Indoor" ~ 0.134,
                                                     type == "DUCK" & system == "Indoor" ~ 0.377,
                                                     type == "TURKEY" & system == "Indoor" ~ 3,
                                                     type == "LAYER_BREEDER" & system == "Indoor" ~ 1,
                                                     type == "BROILER_BREEDER" & system == "Indoor" ~ 0.134))

# Add infectivity
# Now we assume infectivity of all farm is 1.
farm_gen$infect <- 1

# Define index case. There are two types of index case: 
# high populated density area and low populated density area.
# Define index case ----------
# Index farms were selected randomly in a high density poultry area (HDPA) or
# a low density poultry area (LDPA). DPPA farms were random selected from those farms with the 5% percentile highest point density 
# within 1 km from a farm. SPPA were randomly selected among the 75% percentile lowest point density within 1 km.

# So first, we need to start by finding how many farm with radius 1 km of a farm
farm_df <- farm_gen

# Express the coordinates as complex numbers for fast calculation of the euclidean distance
matrix_points <- farm_df[,c("x","y")]
totpoints <- nrow(matrix_points) # total number of points
colnames(matrix_points) <- c("xcoord","ycoord")
# add a column for the index
Index_points <- c(1:totpoints)
Coord <- (complex(length.out=2,real=matrix_points$xcoord,imaginary=matrix_points$ycoord));
distancematrix <- as.matrix(abs(outer(Coord,Coord,"-")))
distancematrix <- distancematrix/1000 # change to km

# From the distance matrix count the distance between farm that below 1 km for each column
farm_df$point_1km <- colSums(distancematrix  <= 1)
summary(farm_df$point_1km)

# Define if it is in high or low density area
farm_df$index <- NA
farm_df$index[farm_df$point_1km >= quantile(farm_df$point_1km, 0.95)] <- "high" # the farms above 5% are in high density areas
farm_df$index[farm_df$point_1km <= quantile(farm_df$point_1km, 0.75)] <- "low" # the farms below 75% are in low  density areas
# add index case to farm data
farm_gen$index <-farm_df$index


# Add farm_id
farm_gen$host_id <- 1:nrow(farm_gen)

# Save generated farm data for outbreak simulation later
saveRDS(farm_gen, "./data/farm_gen.rds")
# Save distance matrix for outbreak simulation later
saveRDS(distancematrix, "./data/distancematrix.rds")