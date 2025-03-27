# Calculate human exposure------------------------------------------------------
# We calculated human exposure from output of between-farm model and 
# fitted gamma distribution from within-farm model
# Read supplementary in the report for fitted gamma distribution

# Specify the packages of interest
packages = c("ggplot2", "sf", "sp","raster", "readxl","ggpubr",
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

# import fam data
farm_2024 <- readRDS("./data/farm_gen.rds")

# Function for calculating human exposure ---------------------------------------
hex_cal <- function (dat) {
 
  a <- dat %>% group_by(iter, per_od_layer, per_od_broiler, cluster, host_id, index) %>% summarise (inf_t = diff(Event_time))
  
   # get the farm size
  a <- a  %>% left_join( farm_2024[c("size", "host_id", "type")], by=c('host_id'='host_id'))
  
  
  # Calculate mean for gamma distribution
  a <- a %>%  mutate(
    hex_mean = case_when(
      type %in% c("LAYER","LAYER_BREEDER" ,"BROILER","BROILER_BREEDER") ~ 8.3863e+00 + 3.9156e-01*size +4.4606e+01*inf_t,
      type%in% c("MEAT_TURKEY") ~ 1.3074e+00 + 1.8588e+00*size +6.9339e+01*inf_t,
      type%in% c("MEAT_DUCK", "DUCK_BREEDER") ~ 2.2198e+03 + 4.8221e-01*inf_t*size))   
 
  
  
  # Calculate mean and variance for gamma distribution
  a <- a %>%  mutate(
    hex_var = case_when(
      type %in% c("LAYER","LAYER_BREEDER" ,"BROILER","BROILER_BREEDER") ~ 1.0001e+00 + 2.1993e+02*size +1.0291e+00*inf_t,
      type%in% c("MEAT_TURKEY") ~ 1.0541e+00 + 1.0685e+03*size +1.7156e+00*inf_t,
      type%in% c("MEAT_DUCK", "DUCK_BREEDER") ~ 7.4070e-01 + 9.1822e+01*inf_t*size))   
  
  # If farm size < 100 the human exposure = sick-animal days
  a <- a %>%  mutate(hex = if_else(size < 100, 
                                   inf_t*size, 
                                   rgamma(1, shape = hex_mean^2/hex_var,  rate = hex_mean/hex_var)))
  
  
  a <- a %>%  mutate(hex = if_else(inf_t*size < hex, 
                                   inf_t*size, 
                                   hex))
  
  return(a)
  
}




# Calculate human exposure ------------------------------------------------------

# Import output of between-farm model from scenario that you want to calculate
result <- readRDS("./results/scen2.1/Sim2.1_layer0_broiler6.5_high.rds")
result <- dplyr::bind_rows(result)


# Calculate human exposure 

hex <- hex_cal(dat = result)
# 
mhex <- hex %>%  group_by(iter) %>% summarise (sum_hex = sum(hex)) 

# column sum)hex is cumulative sick animal-days for wach iteration