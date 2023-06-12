######################
### Author: Shawn Dove
######################

# Function to introduce observation error to time series

obs_error_fn <- function(all_pops_index, m_colnames, mean_cv, cv_sd) {
  
  # create matrix to hold new index values
  all_pops_error <- matrix(data = NA, nrow = nrow(all_pops_index), ncol= length(m_colnames))
  
  # generate distribution of cv values
  cv_dist <- rnorm(1000, mean = mean_cv, sd = cv_sd)
  
  # remove negative values
  cv_dist <- cv_dist[cv_dist >= 0]
  
  # sample cv value for each row of the dataset, with replacement
  row_cv <- sample(cv_dist, size = nrow(all_pops_index), replace = TRUE)
  
  # loop through rows
  for (i in 1:nrow(all_pops_index)) {
    
    # loop through observations
    for (j in 1:(length(m_colnames))) {
      
      # create a normal distribution of index values with actual index value as mean and std_val as st. dev.
      # note that this uses the absolute value of the normal distribution to avoid negative values
      distribution <- rnorm(1000, mean = all_pops_index[i,j], sd = all_pops_index[i,j] * row_cv[i])
      
      distribution_nz <- distribution[distribution > 0] # remove zeros
      
      new_ival <- sample(distribution_nz, 1) # sample from distribution to get new index value with error
      
      all_pops_error[i,j] <- new_ival # put new index value into matrix
      
    }
    
  }
  
  # adjust all indices to begin at a base value of 100
  all_pops_error <- all_pops_error / all_pops_error[,1] * 100
  
  all_pops_error <- as.data.frame(all_pops_error) # convert to data frame
  
  names(all_pops_error) <- m_colnames # get column names from original data frame
  
  all_pops_error$SpecID <- all_pops_index$SpecID # get species ID values from original data frame
  
  all_pops_error$PopID <- all_pops_index$PopID # get population ID values from original data frame
  
  all_pops_error$GrpID <- all_pops_index$GrpID # get group ID values from original data frame
  
  return(all_pops_error)
}