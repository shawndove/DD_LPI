######################
### Author: Shawn Dove
######################

# Function to calculate population and species growth rates from time series

growth_rate_calc_fn <- function(merged.matrix, c, model=FALSE) {
  
  # create empty vectors to store max, min, and mean growth rates
  maxr2 <- vector()
  geomeanr2 <- vector()
  armeanr2 <- vector()
  minr2 <- vector()
  sd.mr2 <- vector()
  specgr <- vector()
  specsd <- vector()
  gratelist2 <- list()
  popgeogrlist2 <- list()
  popargrlist2 <- list()
  popsdlist2 <- list()
  
  counter <- 1
  
  # loop to get max and mean growth rates
  for (i in unique(merged.matrix$SpecID)) {
    
    # create empty vectors to store max, min, and mean growth rates
    maxr <- vector()
    geomeanr <- vector()
    armeanr <- vector()
    minr <- vector()
    sd.mr <- vector()
    gratelist <- list()
    
    tempa <- merged.matrix[merged.matrix$SpecID==i,1:c] # get non-na values
    
    for (j in 1:nrow(tempa)) {
      
      temp <- tempa[j,]
      
      if (model==TRUE) {
        
        temp <- temp[!is.na(temp)] # exclude NAs
        
        temp <- temp[temp > 0] # exclude values of 0 or less
        
      }
      
      if (length(temp) < 2) {
        
        temp <- c(1,1,1)
        
      }
      
      # create a vector of initial values for growth rate calculation. correct for ID column at end of data.
      start.vals <- temp[1:(length(temp) - 1)]
      
      # create a vector of final values for growth rate calculation. correct for ID column at end of data.
      final.vals <- temp[2:(length(temp))]
      
      # calculate the growth rates
      g.rate <- log(final.vals / start.vals)
      
      g.rate[g.rate >= 2.302585] <- 2.302585
      
      g.rate[g.rate <= -2.302585] <- -2.302585
      
      # convert to a vector
      g.rate.vec <- as.vector(as.matrix(g.rate))
      
      # find the maximum growth rate from the given row.
      maxr[j] <- max(g.rate.vec)
      
      # find the mean growth rate from the given row (geometric mean).
      geomeanr[j] <- exp(1) ^ mean(log(g.rate.vec))
      armeanr[j] <- mean(g.rate.vec)
      
      # find the mimimum growth rate from the given row.
      minr[j] <- min(g.rate.vec)
      
      # find the standard deviation of growth rates for a given row
      sd.mr[j] <- sd(g.rate.vec)
      
      gratelist[[j]] <- g.rate.vec
      
    }
    
    maxr2[counter] <- max(maxr)
    
    # arithmetic mean of population mean growth rates for the species
    geomeanr2[counter] <- mean(geomeanr, na.rm=TRUE)
    armeanr2[counter] <- mean(armeanr, na.rm=TRUE)
    
    minr2[counter] <- min(minr)
    
    # standard deviation of population mean growth rates for the species
    sd.mr2[counter] <- sd(geomeanr, na.rm=TRUE)
    
    gratelist.temp <- unlist(gratelist)
    
    specgr[counter] <- mean(gratelist.temp)
    
    specsd[counter] <- sd(gratelist.temp)
    
    gratelist2[[counter]] <- gratelist.temp
    
    # list of population mean growth rates by species
    popgeogrlist2[[counter]] <- geomeanr[!is.na(geomeanr) & is.finite(geomeanr)]
    popargrlist2[[counter]] <- armeanr[!is.na(armeanr) & is.finite(armeanr)]
    
    # list of population growth rate standard deviations by species
    popsdlist2[[counter]] <- sd.mr[!is.na(sd.mr) & is.finite(sd.mr)]
    
    counter <- counter + 1
    
  }
  
  allgrates <- unlist(gratelist2)
  
  sd.all <- sd(allgrates)
  
  popargrates <- unlist(popargrlist2)
  popargrates <- popargrates[!is.na(popargrates) & is.finite(popargrates)]
  
  popgeogrates <- unlist(popgeogrlist2)
  popgeogrates <- popgeogrates[!is.na(popgeogrates) & is.finite(popgeogrates)]
  
  # arithmetic mean of population mean growth rates
  armeanpopargrates <- mean(popargrates)
  armeanpopgeogrates <- mean(popgeogrates)
  
  # geometric mean of population mean growth rates
  geomeanpopargrates <- exp(1) ^ mean(log(popargrates))
  geomeanpopgeogrates <- exp(1) ^ mean(log(popgeogrates))
  
  # population standard deviations of growth rates
  popsds <- unlist(popsdlist2)
  
  # arithmetic mean of population standard deviations of growth rates
  meanpopsds <- mean(popsds, na.rm=TRUE)
  # standard deviation of population standard deviations of growth rates
  sdpopsds <- sd(popsds, na.rm=TRUE)
  
  # standard deviation of population mean growth rates
  sd.geopop <- sd(popgeogrates)
  sd.arpop <- sd(popargrates)
  
  maxr2 <- maxr2[!is.na(maxr2) & is.finite(maxr2)]
  
  # arithmetic means of population geometric mean growth rates for each species
  geomeanr2 <- geomeanr2[!is.na(geomeanr2) & is.finite(geomeanr2)]
  
  minr2 <- minr2[!is.na(minr2) & is.finite(minr2)]
  
  # standard deviations of population geometric mean growth rates for each species
  #sd.mr2 <- sd.mr2[!is.na(sd.mr2) & is.finite(sd.mr2)]
  
  mean.sdmr2 <- mean(sd.mr2, na.rm=TRUE)
  
  maxr.mean <- mean(maxr2, na.rm=TRUE)
  
  maxr.sd <- sd(maxr2, na.rm=TRUE)
  
  maxr.var <- maxr.sd / sqrt(length(maxr2))
  
  maxr.max <- max(maxr2)
  
  maxr.min <- min(maxr2)
  
  maxr.geomean <- exp(1) ^ mean(log(maxr2))
  
  minr.mean <- mean(minr2, na.rm=TRUE)
  
  minr.sd <- sd(minr2, na.rm=TRUE)
  
  minr.var <- minr.sd / sqrt(length(minr2))
  
  minr.max <- max(minr2)
  
  minr.min <- min(minr2)
  
  minr.geomean <- exp(1) ^ mean(log(minr2))
  
  meanr.mean <- mean(geomeanr2)
  
  meanr.sd <- sd(geomeanr2)
  
  meanr.var <- meanr.sd / sqrt(length(geomeanr2))
  
  meanr.max <- max(geomeanr2)
  
  meanr.min <- min(geomeanr2)
  
  meanr.geomean <- exp(1) ^ mean(log(geomeanr2))
  
  return.vals <- list(armeanpopargrates, sd.arpop, meanpopsds, sdpopsds, popsds, popargrlist2)#, geomeanr2, sd.mr2, popgeogrlist2, popsdlist2)
  
  names(return.vals) <- c("Mean of Population Mean Log Growth Rates",
                          "Standard Deviation of Population Mean Log Growth Rates",
                          "Mean of Population Growth Rate Standard Deviations",
                          "Standard Deviation of Population Growth Rate Standard Deviations",
                          "Population Standard Deviations of Growth Rates",
                          "Population Mean Growth Rates by Species")#
  
  return(return.vals)
  
}