######################
### Author: Shawn Dove
######################

# function to remove time series which have fewer than the number of counts specified
# in 'count_thres' and which are shorter than the length specified in 'min_ts_length'

cull_fn <- function(grp_data, count_thres, min_ts_length, c) {
  
  # create a vector to hold row names of time series that meet the thresholds
  notna <- vector()
  
  # set counter to record number of time series that meet the thresholds
  counter <- 1
  
  # loop over each row
  for (i in 1:nrow(grp_data)) {
    
    # check whether the number of observations meets the threshold
    if (length(which(!is.na(grp_data[i,1:c]))) >= count_thres) {
      
      # check whether the length of the time series (from first to last observation) meets the threshold
      if ((max(as.numeric(colnames(grp_data[i,1:c][which(!is.na(grp_data[i,1:c]))]))) - 
           min(as.numeric(colnames(grp_data[i,1:c][which(!is.na(grp_data[i,1:c]))])))) >= (min_ts_length - 1)) {
        
        # if the time series passes both tests, add its row number to the vector notna
        notna[counter] <- rownames(grp_data[i,1:c])
        
        # increase the counter
        counter <- counter + 1
        
      }
      
    }
    
  }
  
  # create a new data frame, including only time series that meet the thresholds
  grp_data_culled <- grp_data[rownames(grp_data) %in% notna,]
  
  return(grp_data_culled)
  
}