# HPAI within-farm model -------------------------------------------------------                        

#include libraries ####
packages <- c("ggplot2","deSolve","tidyverse","bbmle", "fitdistrplus","readxl")


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

# parameters--------------------------------------------------------------------

# Parameters for chicken
param.list.chicken <- list (time = 0,   
                    L0 = 0,
                    R0 = 0, 
                    DS0 = 0, 
                    DL0= 0, 
                    DI0 = 0, 
                    DR0 = 0,
                    mortRate = 0.0005,
                    deltat = 0.1, #time step#
                    pdie = 0.99,
                    beta = 1.31,
                    latency.period = 1,#choose very short latency period
                    k.latency =2,##k parameter Erlang distribution
                    infectious.period = 3,#Duration infectious period 
                    k.infectious = 20, #k parameter Erlang distribution
                    seed = 1)


# Parameters for turkeys
param.list.turkey <- list (time = 0,   
                            L0 = 0,
                            R0 = 0, 
                            DS0 = 0, 
                            DL0= 0, 
                            DI0 = 0, 
                            DR0 = 0,
                            mortRate = 0.0005,
                            deltat = 0.1, #time step#
                            pdie = 0.9,
                            beta = 1.26,
                            latency.period = 1,#choose very short latency period
                            k.latency =2,##k parameter Erlang distribution
                            infectious.period = 6.2,#Duration infectious period 
                            k.infectious = 20, #k parameter Erlang distribution
                            seed = 1)


# Parameters for ducks
param.list.duck <- list (time = 0,   
                            L0 = 0,
                            R0 = 0, 
                            DS0 = 0, 
                            DL0= 0, 
                            DI0 = 0, 
                            DR0 = 0,
                            mortRate = 0.0005,
                            deltat = 0.1, #time step#
                            pdie = 0.7,
                            beta = 4.1,
                            latency.period = 0.17,
                            k.latency =2,##k parameter Erlang distribution
                            infectious.period = 4.3,#Duration infectious period 
                            k.infectious = 20, #k parameter Erlang distribution
                            seed = 1)





# Stochastic SEIR model ------------------------------
# Load necessary package
#####simulation model#######
# for this function the detection was conduct at a complete-full day
 SEIR_HPAI_stochastic <- function(param.list,runs){
   with(c(param.list),{
     #calculate the transition rates in the latency and infectious states
     lLat <- k.latency/latency.period
     lInf <- k.infectious/infectious.period
     currun <- 1
     N0 <- N0
     ifelse (N0 <= 10, I0 <- N0, I0 <-10) 
     S0  <- N0-I0
     
     
     while(currun <= runs){
       
       #initialize the state of the system with k.latency and k.infectious states
       #set state of the system with initial values
       state <- NULL
       state$time <- 1 # time will start at 1 and add 1 at a time, we will convert it to deltat scale later
       state$dt <-0
       state$S <- S0
       state$L <- matrix(0, 1,k.latency)
       state$L[,1] <- L0
       state$I<- matrix(0, 1,k.infectious)
       state$I[,1] <- I0
       state$R <- R0
       state$DS <- DS0
       state$DL <- DL0
       state$DI <- DI0
       state$DR <- DR0
       state$N <- N0
       state$day_mortal <- 0 # dead animal start at 0
       state$passive_detect <- 0
       detect_threshold <- ceiling(state$N*0.005)
       
       # Run model until the run setting and detect criteria is not fulfilled
       while(rowSums(state$I)+rowSums(state$L)>0 & state$passive_detect < 2){
         #Create list to hold the simulated data
         #record output
         if(!exists("output")) output <- NULL;
         output$time <- rbind(output$time, state$time);
         output$run <- rbind(output$run, c(currun));
         output$S <- rbind(output$S, c(state$S));
         output$L <- rbind(output$L, rowSums(state$L));
         output$I <- rbind(output$I, rowSums(state$I));
         output$DS <- rbind(output$DS, c(state$DS));
         output$DL<- rbind(output$DL, c(state$DL));
         output$DI<- rbind(output$DI, c(state$DI));
         output$R <- rbind(output$R, c(state$R));
         output$DR<- rbind(output$DR, c(state$DR));
         output$N <- rbind(output$N, c(state$N));
         
         
         # calculate force of infection
         foi <- (beta * rowSums(state$I)/state$N) 
         
         
         # calculate detect threshold this threshold will be updated every completed 1-day time step
         output$detect_threshold <- rbind(output$detect_threshold, c(detect_threshold) ) 
         if(state$time %% 10 == 1){detect_threshold <- ceiling(state$N * 0.005)}
         
         # reset cumulative mortality every completed 1-day time step
         output$day_mortal <- rbind(output$day_mortal, c(state$day_mortal));
         if(state$time %% 10 == 1) {state$day_mortal <- 0}
         
         # add detect status
         # update detect when complete 1-day time step, the detect time will stop at completed 1 day
         output$passive_detect <-  rbind(output$passive_detect, c(state$passive_detect));
         
         #determine transitions to next state or death
         transS <- rbinom(1, state$S,1-exp(-(foi+mortRate)*deltat))
         transL <- matrix(rbinom(rep(1,k.latency), state$L, 1-exp(-(lLat+mortRate)*deltat)), nrow = 1)
         transI <- matrix(rbinom(rep(1,k.infectious), state$I, 1-exp(-(lInf+mortRate)*deltat)), nrow = 1)
         mortR <- rbinom(1, state$R, 1-exp(-mortRate*deltat))
         
         #divide between death and transition
         mortS <- round(transS*mortRate/(foi+10^-100 + mortRate)) #prevent dividing by 0
         StoL <- transS - mortS
         mortL <- round(transL*mortRate/(lLat + mortRate))
         LtoNext <- cbind(StoL, transL - mortL)
         mortI <- round(transI*mortRate/(lInf + mortRate)) 
         ItoNext <- cbind(LtoNext[,ncol(LtoNext)], transI - mortI)
         
         
         
         #Store the change
         state$S <- state$S - transS
         state$DS <-state$DS + mortS
         state$L <- state$L - transL + LtoNext[,1:(dim(LtoNext)[2]-1)]
         state$DL <- state$DL + rowSums(mortL)
         state$I <- state$I - transI + ItoNext[,1:(dim(ItoNext)[2]-1)] 
         
         #those dying at recovery
         mortIdisease <- rbinom(1,ItoNext[,dim(ItoNext)[2]], pdie)
         state$R <- state$R + c(ItoNext[,dim(ItoNext)[2]]) - mortIdisease - mortR
         state$DI<- state$DI + rowSums(mortI) + mortIdisease
         state$DR <- state$DR + mortR
         # calculate alive animals
         state$N <- N0 - state$DS -state$DL -state$DI - state$DR
         
         # this mortality will keep update until a completes 1-daytime step, 
         # then reset to get a new cumulative death on the next day
         state$day_mortal <- state$day_mortal + mortS + rowSums(mortL) + rowSums(mortI) + mortR + mortIdisease
         
         if(state$time %% 10 == 0 & state$day_mortal >= detect_threshold) {state$passive_detect <- state$passive_detect+1}
         if(state$time %% 10 == 0 & state$day_mortal < detect_threshold) {state$passive_detect <- 0}
         
         #update time
         state$time <- state$time + 1
         
         
       }
       currun <- currun + 1
     }
     output$time <- output$time*0.1-0.1
     return(output)
   })
   
   
 }



# Check farm size for each poultry types---------------------------------------

# Import poultry data
fix_poultry_data_prop <- read_excel("./data/poultry_data_anonymous.xlsx")

# farm size for chicken
summary(fix_poultry_data_prop$farm_size[fix_poultry_data_prop$Animal_type%in% c("BROILER", "LAYER", "BROILER_BREEDER", "LAYER_BREEDER")])
boxplot(fix_poultry_data_prop$farm_size[fix_poultry_data_prop$Animal_type%in% c("BROILER", "LAYER", "BROILER_BREEDER", "LAYER_BREEDER")])



# farm size for duck
summary(fix_poultry_data_prop$farm_size[fix_poultry_data_prop$Animal_type%in% c("MEAT_DUCK", "DUCK_BREEDER")])
boxplot(fix_poultry_data_prop$farm_size[fix_poultry_data_prop$Animal_type%in% c("MEAT_DUCK", "DUCK_BREEDER")])


# farm size for turkey
summary(fix_poultry_data_prop$farm_size[fix_poultry_data_prop$Animal_type%in% c("MEAT_TURKEY")])
boxplot(fix_poultry_data_prop$farm_size[fix_poultry_data_prop$Animal_type%in% c("MEAT_TURKEY")])



# set farm size based on data above
# instead of using fix bin we will set the sample of farm size based on percentile
# to make sure that the within-farm model covered all range of  farm size 
farm_size_chick <- quantile(fix_poultry_data_prop$farm_size[fix_poultry_data_prop$Animal_type%in% c("BROILER", "LAYER", "BROILER_BREEDER", "LAYER_BREEDER")], seq(0,1,0.005))
farm_size_chick <- unique(as.numeric(round(farm_size_chick))) # round number and remove duplicated numbers

farm_size_duck <- quantile(fix_poultry_data_prop$farm_size[fix_poultry_data_prop$Animal_type%in% c("MEAT_DUCK", "DUCK_BREEDER")], seq(0,1,0.005))
farm_size_duck <-  unique(as.numeric(round(farm_size_duck))) # round number and remove duplicated numbers

farm_size_turkey <- quantile(fix_poultry_data_prop$farm_size[fix_poultry_data_prop$Animal_type%in% c("MEAT_TURKEY")], seq(0,1,0.005)) 
farm_size_turkey  <-  unique(as.numeric(round(farm_size_turkey ))) # round number and remove duplicated numbers





## Within farm model for chicken -----------------------------------------------
output_size <- list()
for (i in 1:length(farm_size_chick)) {
  N0 <- farm_size_chick[[i]]
  output_list <- SEIR_HPAI_stochastic(param.list = param.list.chicken , runs = 100)
  output_df <- as.data.frame(do.call(cbind, output_list))
  colnames(output_df) <- names(output_list)
  output_df$farm_size <- N0 
  output_df$species <- "chicken"
  output_size[[i]] <- output_df
  print(i)
}

# Merge output
output_chicken <- as.data.frame(do.call(rbind, output_size))
saveRDS(output_chicken,"./results/within_farm/output_chicken_fullday.rds")

## Summary detection time 
detect_time <- output_size_df %>% group_by(farm_size,run) %>% summarise(detect_time = max(time),
                                                      farm_size = max(farm_size))  


plot(detect_time$farm_size,detect_time$detect_time)

# prepare data for plotting graph  
sum_detect_time <- detect_time %>% group_by(farm_size) %>% summarise(max_time = max(detect_time),
                                                                          median_time = median(detect_time),
                                                                         min_time = min(detect_time))  
sum_detect_time$farm_size <-as.numeric(sum_detect_time$farm_size)
summary(sum_detect_time)
# plot to show detect time and farm size
ggplot(sum_detect_time, aes(x = farm_size)) +
  geom_line(aes(y = median_time, color = "Median"), linetype = "solid") + 
  geom_line(aes(y = min_time, color = "Min"), linetype = "twodash") +
  geom_line(aes(y = max_time, color = "Max"), linetype = "twodash") +
  scale_color_manual(name = "", values = c("Median" = "black", "Min" = "steelblue", "Max" = "red")) +
  scale_y_continuous(breaks = seq(0, 16, 1)) +
  theme_bw()+ 
  ggtitle("Plot of detection time and farm size") +
  xlab("farm size") + ylab("detection time (day)")





# Calculate human exposure from within farm model 
# Calculate human exposure as cumulative day of sick birds
human_exoposure_DI <- output_size_df %>% group_by(farm_size,run) %>% summarise(human_exposure_DI = sum(I)*0.1) # divide 10 because deltat = 0.1
sum_human_exoposure_DI <- human_exoposure_DI %>% group_by(farm_size)%>% summarise(max_human_exposure_DI = max(human_exposure_DI),
                                                                            median_human_exposure_DI = median(human_exposure_DI),
                                                                            min_human_exposure_DI = min(human_exposure_DI))

# plot human exposure as day-infectious
ggplot(sum_human_exoposure_DI , aes(x = farm_size)) +
  geom_line(aes(y = median_human_exposure_DI, color = "Median"), linetype = "solid") + 
  geom_line(aes(y = min_human_exposure_DI, color = "Min"), linetype = "twodash") +
  geom_line(aes(y = max_human_exposure_DI, color = "Max"), linetype = "twodash") +
  scale_color_manual(name = "", values = c("Median" = "black", "Min" = "steelblue", "Max" = "red")) +
  theme_bw()+
  ggtitle("Plot of human exposure") +
  xlab("farm size") + ylab("Cumulative day-infectious")




## Within farm model for Turkeys -----------------------------------------------
# Use a for loop to apply the function to each farm_size
output_size <- list()
for (i in 1:length(farm_size_turkey)) {
  N0 <- farm_size_turkey[[i]]
  output_list <- SEIR_HPAI_stochastic(param.list = param.list.turkey , runs = 100)
  output_df <- as.data.frame(do.call(cbind, output_list))
  colnames(output_df) <- names(output_list)
  output_df$farm_size <- N0 
  output_df$species <- "turkey"
  output_size[[i]] <- output_df
  print(i)
}

# Merge output
output_turkey <- as.data.frame(do.call(rbind, output_size))
saveRDS(output_turkey,"./results/within_farm/output_turkey_fullday.rds")


output_size_df <- output_turkey

# Summary detection time 
detect_time <- output_size_df %>% group_by(farm_size,run) %>% summarise(detect_time = max(time),
                                                                        farm_size = max(farm_size))  
a <- output_size_df[output_size_df$farm_size== 12001,]

# prepare data for plotting graph  
sum_detect_time <- detect_time %>% group_by(farm_size) %>% summarise(max_time = max(detect_time),
                                                                     median_time = median(detect_time),
                                                                     min_time = min(detect_time))  
sum_detect_time$farm_size <-as.numeric(sum_detect_time$farm_size)

# plot to show detect time and farm size
ggplot(sum_detect_time, aes(x = farm_size)) +
  geom_line(aes(y = median_time, color = "Median"), linetype = "solid") + 
  geom_line(aes(y = min_time, color = "Min"), linetype = "twodash") +
  geom_line(aes(y = max_time, color = "Max"), linetype = "twodash") +
  scale_color_manual(name = "", values = c("Median" = "black", "Min" = "steelblue", "Max" = "red")) +
  scale_y_continuous(breaks = seq(0, 16, 1)) +
  theme_bw()+
  ggtitle("Plot of detection time and farm size") +
  xlab("farm size") + ylab("detection time (day)")


# Calculate human exposure from within farm model
# human exposure = sum(infectious*beta*dt)
human_exoposure <- output_size_df %>% group_by(farm_size,run) %>% mutate(sum_foi = I * param.list.chicken$beta * param.list.chicken$deltat ) %>% 
  summarise(human_exposure = sum(sum_foi)) 
sum_human_exoposure <- human_exoposure %>% group_by(farm_size)%>% summarise(max_human_exposure = max(human_exposure),
                                                                            median_human_exposure = median(human_exposure),
                                                                            min_human_exposure = min(human_exposure))

# plot human exposure
ggplot(sum_human_exoposure , aes(x = farm_size)) +
  geom_line(aes(y = median_human_exposure, color = "Median"), linetype = "solid") + 
  geom_line(aes(y = min_human_exposure, color = "Min"), linetype = "twodash") +
  geom_line(aes(y = max_human_exposure, color = "Max"), linetype = "twodash") +
  scale_color_manual(name = "", values = c("Median" = "black", "Min" = "steelblue", "Max" = "red")) +
  theme_bw()+
  ggtitle("Plot of human exposure") +
  xlab("farm size") + ylab("Cumulative force of infection")


# Calculate human exposure as cumulative day of sick birds
human_exoposure_DI <- output_size_df %>% group_by(farm_size,run) %>% summarise(human_exposure_DI = sum(I)*0.1) # divide 10 because deltat = 0.1
sum_human_exoposure_DI <- human_exoposure_DI %>% group_by(farm_size)%>% summarise(max_human_exposure_DI = max(human_exposure_DI),
                                                                                  median_human_exposure_DI = median(human_exposure_DI),
                                                                                  min_human_exposure_DI = min(human_exposure_DI))

# plot human exposure as day-infectious
ggplot(sum_human_exoposure_DI , aes(x = farm_size)) +
  geom_line(aes(y = median_human_exposure_DI, color = "Median"), linetype = "solid") + 
  geom_line(aes(y = min_human_exposure_DI, color = "Min"), linetype = "twodash") +
  geom_line(aes(y = max_human_exposure_DI, color = "Max"), linetype = "twodash") +
  scale_color_manual(name = "", values = c("Median" = "black", "Min" = "steelblue", "Max" = "red")) +
  theme_bw()+
  ggtitle("Plot of human exposure") +
  xlab("farm size") + ylab("Cumulative day-infectious")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~--------------------------------------------------

## Within farm model for Ducks -------------------------------------------------
# Use a for loop to apply the function to each farm_size
output_size <- list()
for (i in 1:length(farm_size_duck)) {
  N0 <- farm_size_duck[[i]]
  output_list <- SEIR_HPAI_stochastic1(param.list = param.list.duck , runs = 100)
  output_df <- as.data.frame(do.call(cbind, output_list))
  colnames(output_df) <- names(output_list)
  output_df$farm_size <- N0 
  output_df$species <- "duck"
  output_size[[i]] <- output_df
  print(i)
}

# Merge output
output_duck <- as.data.frame(do.call(rbind, output_size))
saveRDS(output_duck,"./results/within_farm/output_duck_fullday.rds")


output_size_df <- output_duck

# Summary detection time 
detect_time <- output_size_df %>% group_by(farm_size,run) %>% summarise(detect_time = max(time),
                                                                        farm_size = max(farm_size))  

# prepare data for plotting graph  
sum_detect_time <- detect_time %>% group_by(farm_size) %>% summarise(max_time = max(detect_time),
                                                                     median_time = median(detect_time),
                                                                     min_time = min(detect_time))  
sum_detect_time$farm_size <-as.numeric(sum_detect_time$farm_size)

# plot to show detect time and farm size
ggplot(sum_detect_time, aes(x = farm_size)) +
  geom_line(aes(y = median_time, color = "Median"), linetype = "solid") + 
  geom_line(aes(y = min_time, color = "Min"), linetype = "twodash") +
  geom_line(aes(y = max_time, color = "Max"), linetype = "twodash") +
  scale_color_manual(name = "", values = c("Median" = "black", "Min" = "steelblue", "Max" = "red")) +
  scale_y_continuous(breaks = seq(0, 9, 1)) +
  theme_bw()+
  ggtitle("Plot of detection time and farm size") +
  xlab("farm size") + ylab("detection time (day)")



# Calculate human exposure from within farm model 
# human exposure = sum(infectious*beta*dt)
human_exoposure <- output_size_df %>% group_by(farm_size,run) %>% mutate(sum_foi = I * param.list.chicken$beta * param.list.chicken$deltat ) %>% 
  summarise(human_exposure = sum(sum_foi)) 
sum_human_exoposure <- human_exoposure %>% group_by(farm_size)%>% summarise(max_human_exposure = max(human_exposure),
                                                                            median_human_exposure = median(human_exposure),
                                                                            min_human_exposure = min(human_exposure))

# plot human exposure
ggplot(sum_human_exoposure , aes(x = farm_size)) +
  geom_line(aes(y = median_human_exposure, color = "Median"), linetype = "solid") + 
  geom_line(aes(y = min_human_exposure, color = "Min"), linetype = "twodash") +
  geom_line(aes(y = max_human_exposure, color = "Max"), linetype = "twodash") +
  scale_color_manual(name = "", values = c("Median" = "black", "Min" = "steelblue", "Max" = "red")) +
  scale_y_continuous(breaks = seq(0, 9, 1)) +
  theme_bw()+
  ggtitle("Plot of human exposure") +
  xlab("farm size") + ylab("Cumulative force of infection")


# Calculate human exposure as cumulative day of sick birds
human_exoposure_DI <- output_size_df %>% group_by(farm_size,run) %>% summarise(human_exposure_DI = sum(I)*0.1) # divide 10 because deltat = 0.1
sum_human_exoposure_DI <- human_exoposure_DI %>% group_by(farm_size)%>% summarise(max_human_exposure_DI = max(human_exposure_DI),
                                                                                  median_human_exposure_DI = median(human_exposure_DI),
                                                                                  min_human_exposure_DI = min(human_exposure_DI))

# plot human exposure as day-infectious
ggplot(sum_human_exoposure_DI , aes(x = farm_size)) +
  geom_line(aes(y = median_human_exposure_DI, color = "Median"), linetype = "solid") + 
  geom_line(aes(y = min_human_exposure_DI, color = "Min"), linetype = "twodash") +
  geom_line(aes(y = max_human_exposure_DI, color = "Max"), linetype = "twodash") +
  scale_color_manual(name = "", values = c("Median" = "black", "Min" = "steelblue", "Max" = "red")) +
  theme_bw()+
  ggtitle("Plot of human exposure") +
  xlab("farm size") + ylab("Cumulative day-infectious")

