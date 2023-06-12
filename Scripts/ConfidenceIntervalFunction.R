######################
### Author: Shawn Dove
######################

# Function to calculate confidence intervals

ci_fn <- function(index, c, m_colnames, boots=100) {
  
  # check whether data is in list format
  if (inherits(index, "list")) {
    
    # if it is, convert the list to a data frame
    index.matrix <- do.call(rbind, index)
    
  } else {
    
    # if it is not a list, copy the data frame directly
    index.matrix <- index
    
  }
  
  # group ID set to 1
  grp_ids <- 1
  
  # loop through groups
  for (grp in grp_ids) {
    
    # select the whole dataset
    grp_specdata <- index.matrix
    
    # check if these are already lambda values
    # if they are, there will be one less column than in m_colnames
    # so we test if the final element of m_colnames is missing from the column names
    if (!(tail(m_colnames, 1) %in% colnames(grp_specdata))) {
      
      # if they are already lambda values, remove ID tags
      grp_specdata1.4 <- grp_specdata[,1:(c-1)]
      
      # if they are not already lambda values
    } else {
      
      # remove ID tags
      grp_specdata1.1 <- grp_specdata[,1:c]
      
      # update variable name
      grp_specdata1.3 <- grp_specdata1.1
      
      # convert to log10
      grp_specdata1.4 <- log10(grp_specdata1.3)
      
    }
    
    # create matrix to hold sampled species lambdas/indices
    grp_samp_list <- list()
    
    # loop to bootstrap sampling process
    for (i in 1:boots) {
      
      # create matrix to hold bootstrapped sample
      grp_samp_temp <- matrix(NA, nrow=nrow(grp_specdata1.4), ncol=ncol(grp_specdata1.4))
      
      # loop to sample from each interval
      for (j in 1:ncol(grp_specdata1.4)) {
        
        # get observed values from interval i
        temp <- grp_specdata1.4[!is.na(grp_specdata1.4[,j]),j]
        
        #if there are no observed values in the column
        if(length(temp)==0) {
          
          # put NA into temp2 to avoid an empty vector
          temp2 <- NA
          
          # if there are observed values
        } else {
          
          # sample n observed values with replacement
          temp2 <- sample(temp, replace=TRUE)
          
        }
        
        # put the sample into the sample matrix
        grp_samp_temp[(1:length(temp2)),j] <- temp2
        
      }
      
      # take the mean of the natural log of the species lambdas
      grp_samp_temp1.4 <- colMeans(grp_samp_temp, na.rm=TRUE)
      
      
      # if there are weights...
      # convert NaN values to NA
      grp_samp_temp1.4[is.nan(grp_samp_temp1.4)] <- NA
      
      
      # back convert from log10
      grp_samp_temp1.5 <- 10^(grp_samp_temp1.4)
      
      # if these were lambda values
      if (!(tail(m_colnames, 1) %in% colnames(grp_specdata))) {
        
        # convert to index values
        grp_samp_temp1.7 <- append(grp_samp_temp1.5, 100, after = 0)
        grp_samp_temp1.8 <- cumprod(grp_samp_temp1.7)
        
        # if they were not lambda values
      } else {
        
        # set each index to start at 100
        grp_samp_temp1.8 <- (grp_samp_temp1.5 / grp_samp_temp1.5[1] * 100)
        
      }
      
      # add bootstrapped sample to list
      grp_samp_list[[i]] <- grp_samp_temp1.8
      
    }
    
    # convert list to a data frame
    grp_samp <- do.call(rbind, grp_samp_list)
    
    
    # create matrix to hold confidence intervals
    grp_samp_ci <- matrix(NA, nrow=2, ncol=ncol(grp_samp))
    
    # loop to get confidence intervals for each year
    for (i in 1:ncol(grp_samp)) {
      
      # order the index values for year i
      temp <- sort(as.vector(grp_samp[,i]))
      
      # get the lower bound for year i
      lower_bound <- quantile(temp, 0.025, names=FALSE)
      
      # get the upper bound for year i
      upper_bound <- quantile(temp, 0.975, names=FALSE)
      
      # put the lower bounding confidence interval into the matrix
      grp_samp_ci[1,i] <- lower_bound
      
      # put the upper bounding confidence interval into the matrix
      grp_samp_ci[2,i] <- upper_bound
      
    }
    
    # add years as column names
    colnames(grp_samp_ci) <- m_colnames
    
    # restructure as a data frame
    grp_samp_ci <- as.data.frame(grp_samp_ci)
    
    # rename as confidence interval list
    grp_ci.list <- grp_samp_ci
    
    # print information to show that the msi has been completed
    print(paste("completed confidence intervals for msi"))
    
  }
  
  # exit function
  return(grp_ci.list)
  
}