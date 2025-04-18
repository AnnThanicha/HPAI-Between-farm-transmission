### Define the function handling the events in the spatial transmission model
# Event function has four entries: event_time= time at which an event occurs, eventtype= type of event, statustype=status of host,id_= ID of host
# eventtype:  2=infection; 3=removal;
# statustype: 1= susceptible; 2= infectious; 3= culled;
event <- function(time_event,eventype,statustype,id_){
  if (eventype==2 & statustype==1){
    # calculate the CFI up to that moment
    # I let the CFI grow also for the infected hosts, because they are infected. Of course, the already infected hosts will not be considered in the calculation of the next infection time, because they are already infected 
    CFI <<- CFI+ apply(matrix(hazardmatrix[which(indexI==1),],nrow=length(which(indexI==1)),totpoints),   MARGIN=2,FUN=sum)*(time_event-tt)     # In the rows (i) the infected, in the columns (j) the susceptibles
    # save the CFI in the matrix of CFI
    index_new_event <<- index_new_event + 1
    CFI_matrix[index_new_event,] <<- CFI
    timevector <<- rbind(timevector,time_event)
    # update the status vector and the indices vectors for S and for I
    Status[id_] <<- 2 # now Status= 2 (infected) 
    indexI[id_] <<- 1 # infected
    indexS[id_] <<- 1 # not susceptible anymore
    # save the number of infected over time in a list
    infected_over_time <<- rbind(infected_over_time,c(time_event,length(indexI[indexI==1])))
    # calculate the slope of the force of infection from this moment onwards
    bb[which(indexS==0)] <<- apply(matrix(hazardmatrix[which(indexI==1),which(indexS==0)],nrow=length(which(indexI==1)),ncol=length(which(indexS==0))), MARGIN=2,FUN=sum)       
    # in the rows (i) the infected, in thecolumns (j) the susceptibles
    # calculate the next infection events
    t_infection[which(indexS==0)] <<- (Q_init[which(indexS==0)]-CFI[which(indexS==0)])/bb[which(indexS==0)]
    t_infection[which(indexS==1|indexS==3)] <<- 10000000 # I set a very high number which is not going to happen
    # update the list of points to infect
    next_infection_host <<- which.min(t_infection)  
    next_infection_time <<- t_infection[which.min(t_infection)]
    #update the time
    tt <<- time_event
    #update the list of hosts to infect
    List_to_infect <<- rbind(List_to_infect,data.frame(Event_time = tt+next_infection_time,   Type_event = rep(2,length(next_infection_host)), id_host = next_infection_host))
    List_to_infect <<- List_to_infect[order(List_to_infect[,1]),] 
    #update the list of hosts to remove
    List_to_remove <<- rbind(List_to_remove,data.frame(Event_time=tt+T_inf[id_],Type_event=3,id_host=id_))
    List_to_remove <<- List_to_remove[order(List_to_remove[,1]),] 
    return(1)} else if (eventype==2 & statustype==2) {# already infected
      return(0)} else if (eventype==2 & statustype==3) {# if it has been culled it cannot be infected
        return(0)} else if (eventype==3 & statustype==1) {# it does not occur
          return(0)}else if (eventype==3 & statustype==2){
            # calculate the CFI up to that moment
            # I let the CFI grow also for the infected hosts, because they are infected. 
            # Of course they will not be considered in the calculation of the next infection time, because they are already infected 
            CFI <<- CFI + apply(matrix(hazardmatrix[which(indexI==1),],nrow=length(which(indexI==1)),totpoints),MARGIN=2,FUN=sum)*(time_event-tt) 
            index_new_event <<- index_new_event + 1
            CFI_matrix[index_new_event,] <<- CFI
            # track the CFI over time
            timevector <<- rbind(timevector,time_event)  # track the time
            # update the status vector and the indices vectors for S and for I
            Current[,3] <<- id_
            Status[id_] <<- 3   #culled 
            indexI[id_] <<- 3   #culled, it will not contribute to the infectious matrix anymore
            indexS[id_] <<- 3   #culled, it is not susceptible anymore
            # save the number of infected over time in a list
            infected_over_time <<- rbind(infected_over_time,c(time_event,length(indexI[indexI==1])))
            if(length(which(indexI==1))!=0){ # if there are still individual infected
              # update the slope of the force of infection 
              bb[which(indexS==0)] <<- apply(matrix(hazardmatrix[which(indexI==1),which(indexS==0)],nrow=length(which(indexI==1)), ncol=length(which(indexS==0))),MARGIN=2,FUN=sum) 
              # calculate for each susceptible the next possible infection event
              t_infection[which(indexS==0)] <<- (Q_init[which(indexS==0)]-CFI[which(indexS==0)])/bb[which(indexS==0)]
              t_infection[which(indexS==1|indexS==3)]  <<- 10000000
              # update the time
              tt <<- time_event
              next_infection_host <<- which.min(t_infection)  
              next_infection_time <<- t_infection[which.min(t_infection)]
              # if the next infection host is already in the queue, you need to remove the old one and add the new one with the     updated time of infection
              if ((next_infection_host%in%List_to_infect$id_host)==TRUE){
                List_to_infect <<- List_to_infect[-(which((List_to_infect$id_host%in%next_infection_host)==TRUE)),]
              }
              List_to_infect <<- rbind(List_to_infect,data.frame(Event_time=tt+next_infection_time,Type_event=rep(2,length(next_infection_time)),id_host=next_infection_host))
              List_to_infect <<- List_to_infect[order(List_to_infect[,1]),] 
            } else if (length(which(indexI==1)) == 0){# if there are no infection events anymore the cumulative force infection of the susceptible should be set to 0 
              #update the slope of the force of infection 
              bb[which(indexS==0)] <<- 0           
              List_to_infect <<- {}
            } 
            return(1)}
}
