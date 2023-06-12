######################
### Author: Shawn Dove
######################

## function to interpolate missing years in time series using log-linear interpolation

interpolate_fn <- function(grp_data_culled, c, m_colnames, sample_pop_ids=NA, lessthansix=TRUE, calcsd=FALSE) {
  
  # set forecast to TRUE if you want to project the time series backwards and forewards
  
  if (!is.na(sample_pop_ids)) {
    
    # select all copies of the populations listed in the sample
    grp_data_culled <- grp_data_culled[grp_data_culled$PopID %in% sample_pop_ids,]
    
  }
  
  # check for GAM quality status column
  if ("gqfail" %in% colnames(grp_data_culled)) {
    
    POSTGAM <- TRUE
    
  } else {
    
    POSTGAM <- FALSE
    
  }
  
  if(nrow(grp_data_culled)==0) {
    
    return(grp_data_culled)
    
  }
  
  ## INTERPOLATION
  
  # create a matrix to hold the interpolated time series
  new.grp_data <- as.data.frame(matrix(NA, nrow = nrow(grp_data_culled[,1:c]), ncol = ncol(grp_data_culled[,1:c])))
  
  # create a vector to hold the rows
  rowsmat <- vector()
  
  # create a vector to record whether time series are copied or interpolated
  copied <- vector()
  
  # create a counter to track how many populations have been interpolated
  counter1 <- 1
  
  # begin the interpolation loop
  for (i in 1:nrow(grp_data_culled)) {
    
    # put a single time series (row) into a vector
    new_pop_counts <- as.matrix(grp_data_culled[i,m_colnames])
    
    # put the row name of the time series into a vector
    rownum <- rownames(grp_data_culled[i,])
    
    # check if GAMs have already been performed
    if (POSTGAM==TRUE) {
      
      # if so, check if this population passed quality check
      if (grp_data_culled$gqfail[i]==0) {
        
        # copy the time series into the new matrix
        new.grp_data[i,] <- new_pop_counts
        
        # put the row number of the time series into the row numbers vector
        rowsmat[i] <- rownum
        
        print(paste("copied population ", counter1, sep=""))
        
        # increase the counter
        counter1 <- counter1 + 1
        
        next
        
      }
      
    }
    
    # check if there are any zeros
    if (length(which(new_pop_counts==0))>0 & calcsd==FALSE) {
      
      # check if all non-NA observations are zero
      if (mean(new_pop_counts[which(!is.na(new_pop_counts))])==0) {
        
        # if so, set zero adjust to a very small value
        zero_adjust <- 1e-17
        
      } else {
        
        # otherwise, calculate 1% of the mean of the observed values (excluding zeros)
        # if there are any zeros, this will be added to every observation to avoid issues with log of zero
        zero_adjust <- 0.01 * mean(new_pop_counts[which(new_pop_counts>0)], na.rm=TRUE)
        
      }
      
      # add the zero_adjust value
      new_pop_counts = new_pop_counts + zero_adjust
      
    }
    
    # get the missing counts from the time series
    pop_count_blanks <- new_pop_counts[,which(is.na(new_pop_counts))]
    
    # the lines below are a fix for a problem that occurs when there is only a single NA value
    new_pop_count_names <- as.integer(colnames(new_pop_counts))
    names(pop_count_blanks) <- new_pop_count_names[which(new_pop_counts %in% pop_count_blanks)]
    
    # if there are no missing counts, copy the time series directly and move on
    if (length(pop_count_blanks) == 0) {
      
      # copy the time series into the new matrix
      new.grp_data[i,] <- new_pop_counts
      
      # put the row number of the time series into the row numbers vector
      rowsmat[i] <- rownum
      
      # put 1 into the copied vector to record that this time series was not interpolated
      copied[i] <- 1
      
      print(paste("copied population ", counter1, sep=""))
      
      # increase the counter
      counter1 <- counter1 + 1
      
      next
      
    }
    
    # if this flag is set, only interpolate populations with less than 6 data points
    if (lessthansix==TRUE) {
      
      # if there are 6 or more non-missing counts...
      if (length(new_pop_counts[,which(!is.na(new_pop_counts))]) >= 6) {
        
        # check if GAMs have already been performed
        if (POSTGAM==FALSE) {
          
          # if not, copy the time series into the new matrix
          new.grp_data[i,] <- new_pop_counts
          
          # put the row number of the time series into the row numbers vector
          rowsmat[i] <- rownum
          
          # put 1 into the copied vector to record that this time series was not interpolated
          copied[i] <- 1
          
          print(paste("copied population ", counter1, sep=""))
          
          # increase the counter
          counter1 <- counter1 + 1
          
          next
          
        }
        
      }
      
    } 
    
    
    # get the non-missing counts
    pop_count_filled <- new_pop_counts[,which(!is.na(new_pop_counts))]
    
    # get the column numbers of the missing counts
    pop_count_blank_cols <- which(new_pop_counts %in% pop_count_blanks)
    
    # get the column numbers of the non-missing counts
    pop_count_filled_cols <- which(new_pop_counts %in% pop_count_filled)
    
    # get the years of the missing counts
    pop_count_blank_years <- as.integer(names(pop_count_blanks))
    
    # get the years of the non-missing counts
    pop_count_filled_years <- as.integer(names(pop_count_filled))
    
    # get the year of the first non-missing count
    first_filled_year <- min(pop_count_filled_years)
    
    # get the year of the last non-missing count
    last_filled_year <- max(pop_count_filled_years)
    
    # get the column number of the first non-missing count
    first_filled_col <- min(pop_count_filled_cols)
    
    # get the column number of the last non-missing count
    last_filled_col <- max(pop_count_filled_cols)
    
    # make a vector of all columns between and including the first to the last non-missing counts of the time series
    # this will be used to determine years which must be interpolated. Projection will be done later.
    actual_ts_cols <- first_filled_col:last_filled_col
    
    # make a vector of years that match the columns
    actual_ts_years <- first_filled_year:last_filled_year
    
    # make a vector of columns to be interpolated by matching the above columns vector with the non-missing counts
    actual_ts_blank_cols <- pop_count_blank_cols[pop_count_blank_cols %in% actual_ts_cols]
    
    # get the years of the missing counts within the time series period
    actual_ts_blank_years <- actual_ts_years[actual_ts_years %in% pop_count_blank_years]
    
    # get column numbers of existing counts within the time series period
    actual_ts_filled_cols <- actual_ts_cols[actual_ts_cols %in% pop_count_filled_cols]
    
    # get existing counts for use in interpolation
    actual_ts_filled_vals <- new_pop_counts[,actual_ts_filled_cols]
    
    # get the years of the non-missing counts within the period of the time series
    actual_ts_filled_years <- actual_ts_years[actual_ts_years %in% pop_count_filled_years]
    
    #interpolate to fill in missing values within the time series period (not projected), using log-linear interpolation
    pop_interp_temp <- approx(actual_ts_filled_years, log(actual_ts_filled_vals), actual_ts_blank_years)
    
    # back convert from log
    pop_interp <- exp(pop_interp_temp$y)
    
    # create a counter for use when adding the interpolated values into the time series
    counter2 <- 1
    
    # loop for adding the interpolated values into the time series
    for (j in actual_ts_blank_cols) {
      
      # add interpolated values. counter is used for pop_interp because it has fewer values than new_pop_counts
      new_pop_counts[,j] <- pop_interp[counter2]
      
      # increase the counter each time a value is added to new_pop_counts
      counter2 <- counter2 + 1
      
    }
    
    # put the interpolated time series into the new matrix
    new.grp_data[i,] <- as.matrix(new_pop_counts)
    
    # put the row number of the time series into the row numbers vector
    rowsmat[i] <- rownum
    
    # put 0 into the copied vector to record that this time series was interpolated
    copied[i] <- 0
    
    print(paste("completed interpolation of population ", counter1, sep=""))
    
    counter1 <- counter1 + 1
    
  }
  
  if(length(grp_data_culled)>c) {
    
    # add id tags back
    new.grp_data[,(c+1):length(grp_data_culled)] <- grp_data_culled[,(c+1):length(grp_data_culled)]
    
  }
  
  # put column names into the new matrix
  colnames(new.grp_data) <- colnames(grp_data_culled)
  
  # put row names into the new matrix
  rownames(new.grp_data) <- rowsmat
  
  # check if GAMs have already been performed
  if (POSTGAM==FALSE) {
    
    # if not, add "copied" vector as a column to show the gam function which time series to ignore
    new.grp_data$copied <- copied
    
  }
  
  return(new.grp_data)
    
}