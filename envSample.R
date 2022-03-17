
# Function for environmental filtering from Varela et al. (2014)

envSample<- function (coord, filters, res, do.plot=TRUE){
  
  n<- length (filters)
  pot_points<- list ()
  for (i in 1:n){
    k<- filters [[i]] [!is.na(filters[[i]])]
    ext1<- range (k)
    ext1 [1]<- ext1[1]- 1
    x<- seq(ext1[1],ext1[2], by=res[[i]])
    pot_points[[i]]<- x
  }
  pot_p<- expand.grid(pot_points)
  
  ends<- NULL
  for (i in 1:n){
    fin<- pot_p [,i] + res[[i]]
    ends<- cbind (ends, fin)
  }
  
  pot_pp<- data.frame (pot_p, ends)
  pot_pp<- data.frame (pot_pp, groupID=c(1:nrow (pot_pp)))
  rows<- length (filters[[1]])
  filter<- data.frame(matrix(unlist(filters), nrow=rows))
  real_p<- data.frame (coord, filter)
  
  names_real<- c("lon", "lat")
  names_pot_st<- NULL
  names_pot_end<- NULL
  sql1<- NULL
  for (i in 1:n){
    names_real<- c(names_real, paste ("filter", i, sep=""))
    names_pot_st<- c(names_pot_st, paste ("start_f", i, sep=""))
    names_pot_end<- c(names_pot_end, paste ("end_f", i, sep=""))
    sql1<- paste (sql1, paste ("real_p.filter", i, sep=""), sep=", ")   
  }
  
  names (real_p)<- names_real
  names (pot_pp)<- c(names_pot_st, names_pot_end, "groupID")
  
  conditions<- paste ("(real_p.filter", 1, "<= pot_pp.end_f", 1,") and (real_p.filter", 1, "> pot_pp.start_f", 1, ")", sep="")
  for (i in 2:n){
    conditions<- paste (conditions, 
                        paste ("(real_p.filter", i, "<= pot_pp.end_f", i,") and (real_p.filter", i, "> pot_pp.start_f", i, ")", sep=""), 
                        sep="and")
  }
  
  selection_NA<- sqldf(paste ("select real_p.lon, real_p.lat, pot_pp.groupID",   
                              sql1, "from pot_pp left join real_p on", conditions, sep=" "))
  
  selection<- selection_NA [complete.cases(selection_NA),]
  
  final_points<- selection[!duplicated(selection$groupID), ]
  coord_filter<- data.frame (final_points$lon, final_points$lat) 
  names (coord_filter)<- c("lon", "lat")
  
  if (do.plot==TRUE){
    par (mfrow=c(1,2), mar=c(4,4,0,0.5))
    plot (filters[[1]], filters[[2]], pch=19, 
          col="grey50", xlab="Filter 1", ylab="Filter 2")
    points (final_points$filter1, final_points$filter2, 
            pch=19, col="black", cex=1)
    plot (coord, pch=19, col="grey50")
    map(add=T)
    points (coord_filter, pch=19, col="black", cex=1)
    
  }
  coord_filter
}
