
# plot the number of data points per year in a subset of the LPD
LPIyeardata <- function(groupdata, first_year = 1950, max_year = 2019, plot_title = "") {
  groupdata_years <- groupdata[,((65+first_year-1950)):(134-(2019-max_year))] # extract only the columns with yearly counts
  years_sum_all <- vector(length = 0)
  for (i in 1:length(groupdata_years)) {
    years_sum <- sum(groupdata_years[i]!="NULL" & groupdata_years[i] > 0, na.rm = TRUE)
    years_sum_all[i] <- years_sum
  }
  LPI_years <- c(first_year:max_year)
  years_matrix<-matrix(c(years_sum_all,LPI_years), ncol=2, byrow=FALSE)
  plot(years_matrix[,2],years_matrix[,1], main=plot_title, xlab="Year", ylab="Number of Observations")
  return(years_matrix)
}

LPItslength <- function(groupdata, first_year = 1970, max_year = 2019, plot_title = "") {
  groupdata_years <- groupdata[,((65+first_year-1950)):(134-(2019-max_year))] # extract only the columns with yearly counts
  tslength_sum_all <- vector(length = 0)
  for (i in 1:nrow(groupdata_years)) {
    tslength_sum <- sum(groupdata_years[i,]!="NULL" & groupdata_years[i,] > 0, na.rm = TRUE)
    tslength_sum_all[i] <- tslength_sum
  }
  hist(tslength_sum_all, freq=TRUE, main=plot_title)
  median_tslength <- median(tslength_sum_all)
  mean_tslength <- mean(tslength_sum_all)
  print(paste("The median time series length is", median_tslength, "data points."))
  print(paste("The mean time series length is", mean_tslength, "data points."))
}

LPIclean <- function(groupdata, minm = 3) {
  cleaned <- groupdata[rowSums(groupdata[,65:134] > 0, na.rm = TRUE) >= minm,]
  return(cleaned)
}

LPIclean2 <- function(groupdata, minm = 10) {
  groupdata_years <- groupdata[,65:134]
  names(groupdata_years) <- as.integer(sub(pattern='X', replacement='', x = names(groupdata_years)))
  groupdata_years[is.na(groupdata_years)]<-0
  cleaned <- groupdata[FALSE,]
  # tswidth_all <- groupdata_years[FALSE,]
  # tswidth_all <- vector(length = 0)
  for (i in 1:nrow(groupdata_years)) {
    tswidth <- max(as.integer(colnames(groupdata_years)[groupdata_years[i,]>0], na.rm = TRUE)) - min(as.integer(colnames(groupdata_years)[groupdata_years[i,]>0], na.rm = TRUE))
    if (tswidth < minm) {
      cleaned <- rbind(cleaned, groupdata[i,])
    }
  }
  return(cleaned)
}

LPItswidth <- function(groupdata, first_year = 1970, max_year = 2019, plot_title = "") {
  groupdata_years <- groupdata[,((65+first_year-1950)):(134-(2019-max_year))] # extract only the columns with yearly counts
  names(groupdata_years) <- as.integer(sub(pattern='X', replacement='', x = names(groupdata_years)))
  groupdata_years[groupdata_years=="NULL"]<-0
  tswidth_all <- vector(length = 0)
  tslength_sum_all <- vector(length = 0)
  for (i in 1:nrow(groupdata_years)) {
    tswidth <- max(as.integer(colnames(groupdata_years)[groupdata_years[i,]>0], na.rm = TRUE)) - min(as.integer(colnames(groupdata_years)[groupdata_years[i,]>0], na.rm = TRUE))
    tswidth_all[i] <- tswidth
    tslength_sum <- sum(groupdata_years[i,] > 0, na.rm = TRUE)
    tslength_sum_all[i] <- tslength_sum
  }
  hist(tswidth_all, freq=TRUE, main=plot_title)
  median_tswidth <- median(tswidth_all)
  print(paste("The median time series width is", median_tswidth, "years."))
  mean_tswidth <- mean(tswidth_all)
  print(paste("The mean time series width is", median_tswidth, "years."))
  hist(tslength_sum_all, freq=TRUE, main=plot_title)
  median_tslength <- median(tslength_sum_all)
  print(paste("The median time series length is", median_tslength, "data points."))
  mean_tslength <- mean(tslength_sum_all)
  print(paste("The mean time series length is", mean_tslength, "data points."))
  plot(tswidth_all, tslength_sum_all)
  plot(tslength_sum_all, tswidth_all)
}

# calculate growth rates
growth_rate_byrow <- function(merged.matrix, model=FALSE) {
  
  # create empty vectors to store max and mean growth rates
  maxr <- vector()
  meanr <- vector()
  minr <- vector()
  
  # loop to get max and mean growth rates
  for (i in 1:nrow(merged.matrix)) {
    
    temp <- merged.matrix[i, which(!is.na(merged.matrix[i,]))] # get non-na values
    
    if (model==TRUE) {
      
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
    g.rate <- final.vals / start.vals
    
    g.rate[g.rate >= 10] <- 10
    
    g.rate[g.rate <= -10] <- -10
    
    # convert to a vector
    g.rate.vec <- as.vector(as.matrix(g.rate))
    
    # find the maximum growth rate from the given row.
    maxr[i] <- max(g.rate.vec)
    
    # find the mean log growth rate from the given row (geometric mean).
    #meanr[i] <- exp(1) ^ mean(log(g.rate.vec))
    meanr[i] <- mean(log(g.rate.vec))
    
    # find the mimimum growth rate from the given row.
    minr[i] <- min(g.rate.vec)
  }
  #maxr <- maxr[!is.na(maxr) & maxr > 0 & is.finite(maxr)]
  maxr <- maxr[!is.na(maxr) & is.finite(maxr)]
  
  #meanr <- meanr[!is.na(meanr) & meanr > 0 & is.finite(meanr)]
  meanr <- meanr[!is.na(meanr) & is.finite(meanr)]
  
  minr <- minr[!is.na(minr) & is.finite(minr)]
  
  maxr.mean <- mean(maxr, na.rm=TRUE)
  
  maxr.sd <- sd(maxr, na.rm=TRUE)
  
  maxr.var <- maxr.sd / sqrt(length(maxr))
  
  maxr.max <- max(maxr)
  
  maxr.min <- min(maxr)
  
  maxr.geomean <- exp(1) ^ mean(log(maxr))
  
  minr.mean <- mean(minr, na.rm=TRUE)
  
  minr.sd <- sd(minr, na.rm=TRUE)
  
  minr.var <- minr.sd / sqrt(length(minr))
  
  minr.max <- max(minr)
  
  minr.min <- min(minr)
  
  minr.geomean <- exp(1) ^ mean(log(minr))
  
  meanr.mean <- mean(meanr)
  
  meanr.sd <- sd(meanr)
  
  meanr.var <- meanr.sd / sqrt(length(meanr))
  
  meanr.max <- max(meanr)
  
  meanr.min <- min(meanr)
  
  #meanr.geomean <- exp(1) ^ mean(log(meanr))
  
  #return(list(c(maxr.mean, maxr.sd, maxr.var, maxr.min, maxr.max, maxr.geomean), 
  #            c(meanr.mean, meanr.sd, meanr.var, meanr.min, meanr.max, meanr.geomean)))
  
  #return.vals <- c(meanr.mean, meanr.sd, meanr.max, meanr.min, meanr.var, maxr.max, minr.min)
  
  #names(return.vals) <- c("mean.gr", "sd.gr", "max.mean.gr", "min.mean.gr", "var.mean.gr", "max.max.gr", "min.min.gr")
  
  return.vals <- list(maxr, minr, meanr)
  
  names(return.vals) <- c("maxr", "minr", "meanr")
  
  return(return.vals)
  
}

tsnumobs_byrow <- function(groupdata) {
  tsnumobs_sum_all <- vector(length = 0)
  for (i in 1:nrow(groupdata)) {
    tsnumobs_sum <- sum(groupdata[i,] > 0, na.rm = TRUE)
    tsnumobs_sum_all[i] <- tsnumobs_sum
  }
  return(tsnumobs_sum_all)
}

tslength_byrow <- function(groupdata) {
  tslength_all <- vector(length = 0)
  for (i in 1:nrow(groupdata)) {
    tslength <- max(as.integer(colnames(groupdata)[groupdata[i,]>0]), na.rm = TRUE) - min(as.integer(colnames(groupdata)[groupdata[i,]>0]), na.rm = TRUE)
    tslength_all[i] <- tslength
  }
  return(tslength_all)
}
