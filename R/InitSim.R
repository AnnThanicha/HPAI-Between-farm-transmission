# initialize vector status with 1 (all susceptibles)
Status <- matrix(1,nrow=totpoints)      # status is a vector recording the state of each host.
colnames(Status) <- "status"
#1= susceptible, 2=infectious, 3 culled
Queue <- {} # queue has elements of length 4 (time, event type,host id, extra slot), last slot contains infector in case of transmission event (=event 2)
History <- {} # The {} denotes an empty list in R
Current <-{}
infected_over_time<-{}
time_vector<-{}  # Initialize the Cumulative Force of Infection at the begininng of the epidemic
CFI <-matrix(0,nrow=totpoints)
CFI_matrix <- matrix(0,ncol=totpoints,nrow=10000)
index_new_event <- 0
counter<-0
# create vectors for the infected and the susceptibles, to use to compute the cumulative force of infection
# These are index used for the calculation of the cumulative force of infection.
# They are used to keep track which are infected or susceptible.
# I could have also directly used the vector status, but this is more clear to understand
List_to_remove <-{}
List_to_infect <- {}
indexI <- matrix(0,nrow=totpoints)#indexI==0, not yet infected, indexI==1 infected, indexI==3 culled
indexS <- matrix(0,nrow=totpoints)   # indexS==0, not yet infected, indexS==1 infected, indexS=3 culled
t_infection<-matrix(0,nrow=totpoints)
# coefficient of increase
bb <- matrix(0,nrow=totpoints)
tt<-0 # to keep track of time
next_infection_time <- 0
next_infection_host <- 0

#######  # initialize with the first one to be infected
firstone<-K
Status[firstone] <- 2 # now Status= 2 (infected)
indexI[firstone] <- 1 # infected
indexS[firstone]  <-1 # not susceptible anymore
#calculate coefficient of increase
bb[which(indexS==0)] <- apply(matrix(hazardmatrix[which(indexI==1),which(indexS==0)],nrow=length(which(indexI==1)),ncol=length(which(indexS==0))), MARGIN=2,FUN=sum) # In the rows (i) the infected. in the j the supectibles
# Ann's note: Q_init is the threshold for infection. 
# Susceptible hosts are infected as soon as the cumulative infection pressure exceeds a stochastic threshold to infection Qj. Qj ~ Exp(1)
t_infection[which(indexS==0)] <- (Q_init[which(indexS==0)]- CFI[which(indexS==0)])/bb[which(indexS==0)]

t_infection[which(indexS==1)] <- 10000000   # I set an extremely high number, because it cannot infect itself
next_infection_time <- min(t_infection)
next_infection_host <- which.min(t_infection)
# add this infected host to the history vector
# explain on Re and Im function https://en.wikipedia.org/wiki/Complex_number (complex number)
History <<-rbind(History,data.frame(Event_time=tt,Type_event=2,host_id=firstone,x_coord=Re(Coord[as.numeric(firstone)]),y_coord=Im(Coord[as.numeric(firstone)]))) # record it in the history vector. In the history vector add the coord
# update the list_to_infect and the list_to remove
List_to_infect <- rbind(List_to_infect,data.frame(Event_time= tt +next_infection_time,Type_event=2,id_host= next_infection_host))
List_to_infect <- List_to_infect[order(List_to_infect[,1]),]
List_to_remove <- rbind(List_to_remove,data.frame(Event_time=tt +T_inf[firstone],Type_event=3,id_host=firstone))
List_to_remove <- List_to_remove[order(List_to_remove[,1]),]
next_events<-rbind(List_to_infect[1,],List_to_remove[1,])
index_next_event<- which.min(next_events[,1])
Queue <- rbind(Queue,next_events[index_next_event,])
# now remove this event from the list_to_infect or list_to_remove
# Ann's note: Check List_to_infect and List_to_remove with Queue --> it is equal to Queue --> remove it from list
if(all.equal(cbind(List_to_infect[1,1],List_to_infect[1,2],List_to_infect[1,3]),cbind(Queue[1,1],Queue[1,2],Queue[1,3]))==TRUE){
  List_to_infect <-List_to_infect[-c(1),]
} else if(all.equal(cbind(List_to_remove[1,1],List_to_remove[1,2],List_to_remove[1,3]),cbind(Queue[1,1],Queue[1,2],Queue[1,3]))==TRUE){
  List_to_remove <-List_to_remove[-c(1),]
}
# I add the current CFI (which is 0 at the beginning to the CFI_matrix)
index_new_event <- index_new_event+1
CFI_matrix[index_new_event,] <- CFI
timevector<-rbind(time_vector,tt)
infected_over_time<-rbind(infected_over_time,c(0,length(indexI[indexI==1])))
