# function to shorten and degrade time series
degrade_ts_fn <- function(all_pops_index, c, mlength=10, numobs=5) {
  
  # create a Poisson distribution around the mean length of time series
  # restricted between 2 and the starting length of the t.s.
  length_dist <- rpois(1000, mlength)
  
  # create a Poisson distribution around the mean number of observations of t.s.
  # restricted between 2 and the new mean length of t.s.
  numobs_dist <- rpois(1000, numobs)
  
  # loop through time series
  for (i in 1:nrow(all_pops_index)) {
    
    # randomly set length of time series from Poisson distribution
    length_ts <- sample(length_dist[(length_dist>=2) & (length_dist<=c)], 1)
    
    # randomly set number of observations of t.s. from Poisson distribution,
    numobs_ts <- sample(numobs_dist[(numobs_dist>=2) & (numobs_dist<=length_ts)], 1)
    
    #randomly set start year for the t.s. based on length
    start_ts_year <- sample(1:(c-length_ts+1), 1)
    
    #calculate end year for the t.s. based on start year and length
    end_ts_year <- start_ts_year + length_ts - 1
    
    # cut the time series to the correct length
    cut_ts <- all_pops_index[i,start_ts_year:end_ts_year]
    
    # calculate how many observations to replace with NA
    remove_obs <- length_ts - numobs_ts
    
    # degrade time series to the randomly selected number of observations
    cut_ts[which(cut_ts %in% sample(cut_ts, remove_obs))] <- NA
    
    # fill time series with NAs
    all_pops_index[i,1:c] <- NA
    
    # put cut ts back in
    all_pops_index[i,start_ts_year:end_ts_year] <- cut_ts
    
  }
  
  return(all_pops_index)
  
}

# function to shorten and degrade time series, with all time series in two groups,
# either starting at the beginning or ending at the end
degrade_ts_fn_endpoints <- function(all_pops_index, c, mlength=10, numobs=5) {
  
  # create a Poisson distribution around the mean length of time series
  # restricted between 2 and the starting length of the t.s.
  length_dist <- rpois(1000, mlength)
  
  # create a Poisson distribution around the mean number of observations of t.s.
  # restricted between 2 and the new mean length of t.s.
  numobs_dist <- rpois(1000, numobs)
  
  
  
  # loop through time series
  for (i in 1:nrow(all_pops_index)) {
    
    
    # randomly set length of time series from Poisson distribution
    length_ts <- sample(length_dist[(length_dist>=2) & (length_dist<=c)], 1)
    
    # randomly set number of observations of t.s. from Poisson distribution,
    numobs_ts <- sample(numobs_dist[(numobs_dist>=2) & (numobs_dist<=length_ts)], 1)
    
    # randomly choose whether to start the t.s. at year 1 or end it at the final year
    version <- sample(c(1,2), 1)
    
    # a 1 means start at year 1
    if (version==1) {
      
      # set start year to year 1
      start_ts_year <- 1
      
      #calculate end year for the t.s. based on start year and length
      end_ts_year <- start_ts_year + length_ts - 1
      
      # a 2 means end at final year
    } else if (version==2) {
      
      # set end year to final year
      end_ts_year <- c
      
      # calculate start year for the t.s. based on end year and length
      start_ts_year <- end_ts_year - length_ts + 1
      
    }

    # cut the time series to the correct length
    cut_ts <- all_pops_index[i,start_ts_year:end_ts_year]
    
    # calculate how many observations to replace with NA
    remove_obs <- length_ts - numobs_ts
    
    # degrade time series to the randomly selected number of observations
    cut_ts[which(cut_ts %in% sample(cut_ts, remove_obs))] <- NA
    
    # fill time series with NAs
    all_pops_index[i,1:c] <- NA
    
    # put cut ts back in
    all_pops_index[i,start_ts_year:end_ts_year] <- cut_ts
    
  }
  
  return(all_pops_index)
  
}

# function to shorten and degrade time series, with time series divided into clusters
# each cluster has the same length and randomized start and end years
degrade_ts_fn_clust <- function(all_pops_index, c, mlength=10, numobs=5, numclust=10) {
  
  # create a Poisson distribution around the mean length of time series
  # restricted between 2 and the starting length of the t.s.
  length_dist <- rpois(1000, mlength)
  
  # create a Poisson distribution around the mean number of observations of t.s.
  # restricted between 2 and the new mean length of t.s.
  numobs_dist <- rpois(1000, numobs)
  
  # create vectors
  length_ts_clust <- vector()
  numobs_ts_clust <- vector()
  start_ts_year_clust <- vector()
  end_ts_year_clust <- vector()
  
  # loop through clusters
  for (h in 1:numclust) {
    
    # randomly set length of cluster from Poisson distribution
    length_ts_clust[h] <- sample(length_dist[(length_dist>=2) & (length_dist<=c)], 1)
    
    # randomly set number of observations of cluster from Poisson distribution,
    numobs_ts_clust[h] <- sample(numobs_dist[(numobs_dist>=2) & (numobs_dist<=length_ts_clust[h])], 1)
    
    #randomly set start year for the cluster based on length
    start_ts_year_clust[h] <- sample(1:(c-length_ts_clust[h]+1), 1)
    
    #calculate end year for the cluster based on start year and length
    end_ts_year_clust[h] <- start_ts_year_clust[h] + length_ts_clust[h] - 1
    
  }
  
  # loop through time series
  for (i in 1:nrow(all_pops_index)) {
    
    # randomly choose whether to assign to the t.s. to a cluster (probably is 2 out of 3)
    version <- sample(c(1,1,2), 1)
    
    # a 1 means assign to a cluster
    if (version==1) {
      
      # assign row to a cluster
      clust_num <- sample(1:numclust, 1)
      
      # cut the time series to the correct length
      cut_ts <- all_pops_index[i,start_ts_year_clust[clust_num]:end_ts_year_clust[clust_num]]
      
      # calculate how many observations to replace with NA
      remove_obs <- length_ts_clust[clust_num] - numobs_ts_clust[clust_num]
      
      # degrade time series to the randomly selected number of observations
      cut_ts[which(cut_ts %in% sample(cut_ts, remove_obs))] <- NA
      
      # fill time series with NAs
      all_pops_index[i,1:c] <- NA
      
      # put cut ts back in
      all_pops_index[i,start_ts_year_clust[clust_num]:end_ts_year_clust[clust_num]] <- cut_ts
      
      # a 2 means do not assign to a cluster
    } else if (version==2) {
      
      # randomly set length of time series from Poisson distribution
      length_ts <- sample(length_dist[(length_dist>=2) & (length_dist<=c)], 1)
      
      # randomly set number of observations of t.s. from Poisson distribution,
      numobs_ts <- sample(numobs_dist[(numobs_dist>=2) & (numobs_dist<=length_ts)], 1)
      
      #randomly set start year for the t.s. based on length
      start_ts_year <- sample(1:(c-length_ts+1), 1)
      
      #calculate end year for the t.s. based on start year and length
      end_ts_year <- start_ts_year + length_ts - 1
      
      # cut the time series to the correct length
      cut_ts <- all_pops_index[i,start_ts_year:end_ts_year]
      
      # calculate how many observations to replace with NA
      remove_obs <- length_ts - numobs_ts
      
      # degrade time series to the randomly selected number of observations
      cut_ts[which(cut_ts %in% sample(cut_ts, remove_obs))] <- NA
      
      # fill time series with NAs
      all_pops_index[i,1:c] <- NA
      
      # put cut ts back in
      all_pops_index[i,start_ts_year:end_ts_year] <- cut_ts
      
    }
    

    
  }
  
  return(all_pops_index)
  
}

# function to shorten and degrade time series
degrade_ts_fn_clustend <- function(all_pops_index, c, mlength=10, numobs=5, clustlength=3) {
  
  # create a Poisson distribution around the mean length of time series
  # restricted between 2 and the starting length of the t.s.
  length_dist <- rpois(1000, mlength)
  
  # create a Poisson distribution around the mean number of observations of t.s.
  # restricted between 2 and the new mean length of t.s.
  numobs_dist <- rpois(1000, numobs)
  

  # set length of endpoint cluster
  length_ts_clust <- clustlength
    
  # set number of observations of endpoint cluster
  numobs_ts_clust <- clustlength
    
  #set end year for the cluster
  end_ts_year_clust <- c
    
  #calculate start year for the cluster based on length
  start_ts_year_clust <- c - length_ts_clust + 1

  # loop through time series
  for (i in 1:nrow(all_pops_index)) {
    
    # randomly choose whether to assign to the t.s. to a cluster 
    # (probably is the number of 1s divided by number of 2s)
    version <- sample(c(1,2,2,2), 1)
    
    # a 1 means assign to a cluster
    if (version==1) {

      # cut the time series to the correct length
      cut_ts <- all_pops_index[i,start_ts_year_clust:end_ts_year_clust]
      
      # calculate how many observations to replace with NA
      remove_obs <- length_ts_clust - numobs_ts_clust
      
      # degrade time series to the randomly selected number of observations
      cut_ts[which(cut_ts %in% sample(cut_ts, remove_obs))] <- NA
      
      # fill time series with NAs
      all_pops_index[i,1:c] <- NA
      
      # put cut ts back in
      all_pops_index[i,start_ts_year_clust:end_ts_year_clust] <- cut_ts
      
      # a 2 means do not assign to a cluster
    } else if (version==2) {
      
      # randomly set length of time series from Poisson distribution
      length_ts <- sample(length_dist[(length_dist>=2) & (length_dist<=c)], 1)
      
      # randomly set number of observations of t.s. from Poisson distribution,
      numobs_ts <- sample(numobs_dist[(numobs_dist>=2) & (numobs_dist<=length_ts)], 1)
      
      #randomly set start year for the t.s. based on length
      start_ts_year <- sample(1:(c-length_ts+1), 1)
      
      #calculate end year for the t.s. based on start year and length
      end_ts_year <- start_ts_year + length_ts - 1
      
      # cut the time series to the correct length
      cut_ts <- all_pops_index[i,start_ts_year:end_ts_year]
      
      # calculate how many observations to replace with NA
      remove_obs <- length_ts - numobs_ts
      
      # degrade time series to the randomly selected number of observations
      cut_ts[which(cut_ts %in% sample(cut_ts, remove_obs))] <- NA
      
      # fill time series with NAs
      all_pops_index[i,1:c] <- NA
      
      # put cut ts back in
      all_pops_index[i,start_ts_year:end_ts_year] <- cut_ts
      
    }
    
    
    
  }
  
  return(all_pops_index)
  
}

# function to shorten and degrade time series
degrade_ts_fn_endreveal <- function(all_pops_index, c, mlength=10, numobs=5, endlength=5) {
  
  # sample rows for random endreveal
  rev_rows <- sample(1:nrow(all_pops_index), ceiling(0.5*nrow(all_pops_index)))
  # save x final values, where x is endlength
  reveal <- all_pops_index[,(c-endlength+1):c]
  
  # create a Poisson distribution around the mean length of time series
  # restricted between 2 and the starting length of the t.s.
  length_dist <- rpois(1000, mlength)
  
  # create a Poisson distribution around the mean number of observations of t.s.
  # restricted between 2 and the new mean length of t.s.
  numobs_dist <- rpois(1000, numobs)
  
  # loop through time series
  for (i in 1:nrow(all_pops_index)) {
    
    # randomly set length of time series from Poisson distribution
    length_ts <- sample(length_dist[(length_dist>=2) & (length_dist<=c)], 1)
    
    # randomly set number of observations of t.s. from Poisson distribution,
    numobs_ts <- sample(numobs_dist[(numobs_dist>=2) & (numobs_dist<=length_ts)], 1)
    
    #randomly set start year for the t.s. based on length
    start_ts_year <- sample(1:(c-length_ts+1), 1)
    
    #calculate end year for the t.s. based on start year and length
    end_ts_year <- start_ts_year + length_ts - 1
    
    # cut the time series to the correct length
    cut_ts <- all_pops_index[i,start_ts_year:end_ts_year]
    
    # calculate how many observations to replace with NA
    remove_obs <- length_ts - numobs_ts
    
    # degrade time series to the randomly selected number of observations
    cut_ts[which(cut_ts %in% sample(cut_ts, remove_obs))] <- NA
    
    # fill time series with NAs
    all_pops_index[i,1:c] <- NA
    
    # put cut ts back in
    all_pops_index[i,start_ts_year:end_ts_year] <- cut_ts
    
    # if current row is chosen to have its end value(s) revealed
    if (i %in% rev_rows) {
      
      # restore x final values, where x is endlength
      all_pops_index[i,(c-endlength+1):c] <- reveal[i,]
    }

  }
  
  return(all_pops_index)
  
}

# function to shorten and degrade time series. this one removes all except the final x data points,
# where x is clustlength
degrade_ts_fn_clustend2 <- function(all_pops_index, c, clustlength=3) {
  

  # set length of endpoint cluster
  length_ts_clust <- clustlength
  
  # set number of observations of endpoint cluster
  numobs_ts_clust <- clustlength
  
  #set end year for the cluster
  end_ts_year_clust <- c
  
  #calculate start year for the cluster based on length
  start_ts_year_clust <- c - length_ts_clust + 1
  
  # loop through time series
  for (i in 1:nrow(all_pops_index)) {

      # cut the time series to the correct length
      cut_ts <- all_pops_index[i,start_ts_year_clust:end_ts_year_clust]
      
      # calculate how many observations to replace with NA
      remove_obs <- length_ts_clust - numobs_ts_clust
      
      # degrade time series to the randomly selected number of observations
      cut_ts[which(cut_ts %in% sample(cut_ts, remove_obs))] <- NA
      
      # fill time series with NAs
      all_pops_index[i,1:c] <- NA
      
      # put cut ts back in
      all_pops_index[i,start_ts_year_clust:end_ts_year_clust] <- cut_ts

  }
  
  return(all_pops_index)
  
}
