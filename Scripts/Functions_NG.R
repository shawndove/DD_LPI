### functions ###

# remove time series which have fewer than the number of counts specified in 'count_thres' and which are shorter than
# the length specified in 'min_ts_length'
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


## function to project time series forward and backward and interpolate missing years so all are complete
## this version only uses data from within a species to project other populations of that species
complete_time_series <- function(grp_data_culled, c, m_colnames, sample_pop_ids=NA, lambda=FALSE, gam_all=FALSE, calcsd=FALSE) {
  
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
    
    # if the lambda flag is set, only interpolate populations with less than 6 data points
    if (lambda==TRUE) {

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

  # stop here if using lambda method or if calcsd is set to true
  if (lambda==TRUE | calcsd==TRUE) {
    
    return(new.grp_data)
    
  }
  
  ## PROJECTION
  
  # create a matrix to hold the projected time series
  new.grp_data2 <- as.data.frame(matrix(NA, nrow = nrow(new.grp_data[,1:c]), ncol = ncol(new.grp_data[,1:c])))
  
  # create a vector to hold the rows
  rowsmat <- vector()
  
  # create a counter to track how many populations have been projected
  counter2 <- 1
  
  # begin the projection loop
  for (i in 1:nrow(new.grp_data)) {
    
    grp <- new.grp_data[i,"GrpID"]
    
    grp_all <- new.grp_data[new.grp_data$GrpID==grp,]
    
    rownum <- rownames(new.grp_data[i,])
    
    new_pop_counts <- new.grp_data[i,1:c]
    
    # get the non-missing values
    pop_counts <- new_pop_counts[which(!is.na(new_pop_counts))]
    
    # get the missing values
    pop_counts_missing <- new_pop_counts[which(is.na(new_pop_counts))]
    
    # if there are no missing values, copy the time series directly and move on
    if (length(pop_counts_missing) == 0) {
      
      # copy time series into new data frame
      new.grp_data2[i,] <- new.grp_data[i,1:c]
      
      # add the row number to rows vector
      rowsmat[i] <- rownum
      
      print(paste("completed forecasting of population ", counter2, sep=""))
      
      # increase the counter
      counter2 <- counter2 + 1
      
      next
    }
    
    # get all the years from the data
    pop_all_years <- as.numeric(colnames(new_pop_counts))
    
    # get all the column numbers
    pop_all_cols <- which(new_pop_counts %in% new_pop_counts)
    
    # get the first non-missing value
    pop_count_first <- pop_counts[1]
    
    # get the last non-missing value
    pop_count_last <- pop_counts[length(pop_counts)]
    
    # get the years containing actual counts
    pop_count_years <- as.numeric(colnames(pop_counts))
    
    # get the first year containing an actual count
    pop_count_first_year <- min(pop_count_years)
    
    # get the last year containing an actual count
    pop_count_last_year <- max(pop_count_years)
    
    # get the columns containing actual counts
    pop_count_cols <- which(new_pop_counts %in% pop_counts)
    
    # get the first column containing an actual count
    pop_count_first_col <- min(pop_count_cols)
    
    # get the last column containing an actual count
    pop_count_last_col <- max(pop_count_cols)
    
    # get the years which need to be projected backwards
    pop_backproject_years <- min(pop_all_years):pop_count_first_year
    
    # get the columns corresponding to these years
    pop_backproject_cols <- min(pop_all_cols):pop_count_first_col
    
    # get the columns which need to be projected forwards
    pop_forproject_years <- pop_count_last_year:max(pop_all_years)
    
    # get the columns corresponding to these years
    pop_forproject_cols <- pop_count_last_col:max(pop_all_cols)
    
    # put the first count into a temporary vector for use in the backward projection loop
    pop_count_temp_back <- pop_count_first
    
    # put the last count into a temporary vector for use in the forward projection loop
    pop_count_temp_for <- pop_count_last
    
    # only proceed with backward projection if there are actually years to project
    if (length(pop_backproject_years) > 1) {
      
      # begin backward projection loop
      for (j in (length(pop_backproject_years) - 1):1) {
        
        if (!is.na(new_pop_counts[,j])) {
          
          next
          
        }
        
        # check which populations have non-NA values for the years which need to be projected
        b.int1_year_rows <- which(!is.na(grp_all[,pop_backproject_cols[j + 1]]))
        b.int2_year_rows <- which(!is.na(grp_all[,pop_backproject_cols[j]]))
        
        if (length(b.int2_year_rows)==0) {
          
            if (j > 1) { # check if there is more than 1 year left to be projected

              for (k in (j-1):1) {
                
                # check which populations have non-NA values for the years which need to be projected
                b.int1_year_rows <- which(!is.na(grp_all[,pop_backproject_cols[j + 1]]))
                b.int2_year_rows <- which(!is.na(grp_all[,pop_backproject_cols[k]]))
                
                if (length(b.int2_year_rows)==0) {
                  
                  next
                  
                } else {
                  
                  # take only the populations which have a count for BOTH years
                  b.rows_intersect <- intersect(b.int1_year_rows, b.int2_year_rows)
                  b.int1_year_counts <- grp_all[b.rows_intersect, pop_backproject_cols[j + 1]]
                  b.int2_year_counts <- grp_all[b.rows_intersect, pop_backproject_cols[k]]
                  
                  if (length(b.rows_intersect)==0) {
                    
                    # take populations which have a count for EITHER year
                    b.int1_year_counts <- grp_all[b.int1_year_rows, pop_backproject_cols[j + 1]]
                    b.int2_year_counts <- grp_all[b.int2_year_rows, pop_backproject_cols[k]]
                    
                  }
                  
                  # take the geometric means from the two years and divide
                  b.int1_2_multiplier <- (exp(1) ^ mean(log(b.int2_year_counts + 10e-17))) / 
                    (exp(1) ^ mean(log(b.int1_year_counts + 10e-17)))
                  
                  # get first non-NA value for the time series being projected
                  b.int1_year_value <- pop_count_temp_back
                  
                  # get year of first non-NA value
                  b.int1_year <- pop_backproject_years[j + 1]
                  
                  # calculate previous year value
                  b.int2_year_value <- b.int1_year_value * b.int1_2_multiplier
                  
                  # get previous year
                  b.int2_year <- pop_backproject_years[k]
                  
                  # put projected values into the time series
                  new_pop_counts[,pop_backproject_cols[j + 1]] <- b.int1_year_value
                  new_pop_counts[,pop_backproject_cols[k]] <- b.int2_year_value
                  names(new_pop_counts[,pop_backproject_cols[j + 1]]) <- b.int1_year
                  names(new_pop_counts[,pop_backproject_cols[k]]) <- b.int2_year
                  
                  # fill in any missing years using log-linear interpolation
                  k_mis_years <- (k+1):j
                  k_present_vals <- c(new_pop_counts[,k],new_pop_counts[,j+1])
                  
                  k_interp_temp <- approx(c(k,(j+1)), log(k_present_vals + 10e-17), k_mis_years)
                  # back convert from log
                  new_pop_counts[,k_mis_years] <- exp(k_interp_temp$y)
                  
                  # replace temporary count vector with the next one
                  pop_count_temp_back <- b.int2_year_value
                  
                  
                }
                
              }
                
            } else { # if there is only one year left to be projected
              
              next
              
            }
            
        } else {
          
          # take only the populations which have a count for BOTH years
          b.rows_intersect <- intersect(b.int1_year_rows, b.int2_year_rows)
          b.int1_year_counts <- grp_all[b.rows_intersect, pop_backproject_cols[j + 1]]
          b.int2_year_counts <- grp_all[b.rows_intersect, pop_backproject_cols[j]]
          
          
          if (length(b.rows_intersect)==0) {
            
            # take populations which have a count for EITHER year
            b.int1_year_counts <- grp_all[b.int1_year_rows, pop_backproject_cols[j + 1]]
            b.int2_year_counts <- grp_all[b.int2_year_rows, pop_backproject_cols[j]]
            
          }
          
          # take the geometric means from the two years and divide
          b.int1_2_multiplier <- (exp(1) ^ mean(log(b.int2_year_counts + 10e-17))) / 
            (exp(1) ^ mean(log(b.int1_year_counts + 10e-17)))
          
          # get first non-NA value for the time series being projected
          b.int1_year_value <- pop_count_temp_back
          
          # get year of first non-NA value
          b.int1_year <- pop_backproject_years[j + 1]
          
          # calculate previous year value
          b.int2_year_value <- b.int1_year_value * b.int1_2_multiplier
          
          # get previous year
          b.int2_year <- pop_backproject_years[j]
          
          # put projected values into the time series
          new_pop_counts[,pop_backproject_cols[j + 1]] <- b.int1_year_value
          new_pop_counts[,pop_backproject_cols[j]] <- b.int2_year_value
          names(new_pop_counts[,pop_backproject_cols[j + 1]]) <- b.int1_year
          names(new_pop_counts[,pop_backproject_cols[j]]) <- b.int2_year
          
          # replace temporary count vector with the next one
          pop_count_temp_back <- b.int2_year_value
          
        }

      }
      
    }
    
    # only proceed if there are actually years to project forward
    if (length(pop_forproject_years) > 1) {
      
      # begin forward projection loop
      for(j in 1:(length(pop_forproject_years) - 1)) {

        if (!is.na(new_pop_counts[,pop_forproject_cols[j + 1]])) {
          
          next
          
        }
        
        # check which populations have non-NA values for the years which need to be projected
        f.int1_year_rows <- which(!is.na(grp_all[,pop_forproject_cols[j]]))
        f.int2_year_rows <- which(!is.na(grp_all[,pop_forproject_cols[j + 1]]))
        
        if (length(f.int2_year_rows)==0) {
          
          if (j < (length(pop_forproject_years) - 1)) { # check if there is more than 1 year left to be projected
            
            for (k in (j+2):length(pop_forproject_years)) {
              
              # check which populations have non-NA values for the years which need to be projected
              f.int1_year_rows <- which(!is.na(grp_all[,pop_forproject_cols[j]]))
              f.int2_year_rows <- which(!is.na(grp_all[,pop_forproject_cols[k]]))
              
              if (length(f.int2_year_rows)==0) {
                
                next
                
              } else {
                
                # take only the populations which have a count for BOTH years
                f.rows_intersect <- intersect(f.int1_year_rows, f.int2_year_rows)
                f.int1_year_counts <- grp_all[f.rows_intersect, pop_forproject_cols[j]]
                f.int2_year_counts <- grp_all[f.rows_intersect, pop_forproject_cols[k]]
                
                if (length(f.rows_intersect)==0) {
                  
                  # take populations which have a count for EITHER year
                  f.int1_year_counts <- grp_all[f.int1_year_rows, pop_forproject_cols[j]]
                  f.int2_year_counts <- grp_all[f.int2_year_rows, pop_forproject_cols[k]]
                  
                }

                # take the geometric means from the two years and divide
                f.int1_2_multiplier <- (exp(1) ^ mean(log(f.int2_year_counts + 10e-17))) / 
                  (exp(1) ^ mean(log(f.int1_year_counts + 10e-17)))
                
                # get first non-NA value for the time series being projected
                f.int1_year_value <- pop_count_temp_for
                
                # get year of first non-NA value
                b.int1_year <- pop_forproject_years[j]
                
                # calculate next year value
                f.int2_year_value <- f.int1_year_value * f.int1_2_multiplier
                
                # get next year
                f.int2_year <- pop_forproject_years[k]
                
                # put projected values into the time series
                new_pop_counts[,pop_forproject_cols[j]] <- f.int1_year_value
                new_pop_counts[,pop_forproject_cols[k]] <- f.int2_year_value
                names(new_pop_counts[,pop_forproject_cols[j]]) <- f.int1_year
                names(new_pop_counts[,pop_forproject_cols[k]]) <- f.int2_year
                
                # fill in any missing years using log-linear interpolation
                kf_mis_years <- pop_forproject_cols[(j+1):(k-1)]
                kf_present_vals <- c(new_pop_counts[,pop_forproject_cols[j]], 
                                     new_pop_counts[,pop_forproject_cols[k]])
                
                kf_interp_temp <- approx(c(pop_forproject_cols[j],pop_forproject_cols[k]), 
                                                                             log(kf_present_vals + 10e-17), 
                                                                             kf_mis_years)
                
                # back convert from log
                new_pop_counts[,kf_mis_years] <- exp(kf_interp_temp$y)
                
                # replace temporary count vector with the next one
                pop_count_temp_for <- f.int2_year_value
                
                
              }
              
            }
            
          } else { # if there is only one year left to be projected
            
            next
            
          }
          
        } else {
        
          # take only the populations which have a count for BOTH years
          f.rows_intersect <- intersect(f.int1_year_rows, f.int2_year_rows)
          f.int1_year_counts <- grp_all[f.rows_intersect, pop_forproject_cols[j]]
          f.int2_year_counts <- grp_all[f.rows_intersect, pop_forproject_cols[j + 1]]
          
          if (length(f.rows_intersect)==0) {
            
            # take populations which have a count for EITHER year
            f.int1_year_counts <- grp_all[f.int1_year_rows, pop_forproject_cols[j]]
            f.int2_year_counts <- grp_all[f.int2_year_rows, pop_forproject_cols[j + 1]]
            
          }
          
          # take the geometric means from the two years and divide
          f.int1_2_multiplier <- (exp(1) ^ mean(log(f.int2_year_counts + 10e-17))) / 
            (exp(1) ^ mean(log(f.int1_year_counts + 10e-17)))
          
          # get first non-NA value for the time series being projected
          f.int1_year_value <- pop_count_temp_for
          
          # get year of first non-NA value
          f.int1_year <- pop_forproject_years[j]
          
          # calculate next year value
          f.int2_year_value <- f.int1_year_value * f.int1_2_multiplier
          
          # get next year
          f.int2_year <- pop_forproject_years[j + 1]
          
          # put projected values into the time series
          new_pop_counts[,pop_forproject_cols[j]] <- f.int1_year_value
          new_pop_counts[,pop_forproject_cols[j + 1]] <- f.int2_year_value
          names(new_pop_counts[,pop_forproject_cols[j]]) <- f.int1_year
          names(new_pop_counts[,pop_forproject_cols[j + 1]]) <- f.int2_year
          pop_count_temp_for <- f.int2_year_value
          
        }
        
      }
        
    }
    
    # put the projected counts into the new matrix
    new.grp_data2[i,] <- as.matrix(new_pop_counts)
    
    # put the row number into the row numbers vector
    rowsmat[i] <- rownum

    print(paste("completed forecasting of population ", counter2, sep=""))
    
    counter2 <- counter2 + 1
    
  }
  
  # add id tags back
  new.grp_data2[,(c+1):length(grp_data_culled)] <- grp_data_culled[,(c+1):length(grp_data_culled)]
    
  # put column names into the new matrix
  colnames(new.grp_data2) <- colnames(grp_data_culled)

  # put row names into the new matrix
  rownames(new.grp_data2) <- rowsmat
  
  # add "copied" vector as a column to show the gam function which time series to ignore
  new.grp_data2$copied <- copied

  return(new.grp_data2)
  
}


pop_gam_fn <- function(new.grp_data, c, m_colnames, n=NA, lambda=FALSE, resample=FALSE, forecast=FALSE, quality=FALSE) {
  
  # create a list to put resampled populations into
  gam_poplist <- list()
  
  # create a vector of population IDs
  pop_ids <- new.grp_data$PopID
  
  # create a vector of row numbers
  rows <- 1:nrow(new.grp_data)
  
  # create vector to record populations which fail GAM quality check
  gqfail <- vector()
  
  # reorganize the data into a long format that works with the population resampling code
  trim.mat <- matrix(NA, nrow = nrow(new.grp_data) * c, ncol = 3)
  
  # name columns
  colnames(trim.mat) <- c("population", "year", "count")
  
  # convert to data frame
  trim.mat <- as.data.frame(trim.mat)
  
  # fill columns with data
  trim.mat[,1] <- as.factor(rep(pop_ids, each = c)) # population ID
  
  trim.mat[,2] <- as.numeric(rep(m_colnames, times = nrow(new.grp_data))) # year
  
  trim.mat[,3] <- as.numeric(as.vector(t(new.grp_data[m_colnames])))  # count
  
  # if resampling flag is turned off...
  if (resample==FALSE) {
    
    # create copied vector to check which populations have already been interpolated
    copied <- new.grp_data$copied
    
    # create vector of rows which have been interpolated
    copied_rows <- which(copied[rows]==0)
    
    # create vector of rows which have not been interpolated
    gam_rows <- which(copied[rows]==1)
    
    for (row in copied_rows) {
      
      # copy population into completed list, without the "copied" column
      gam_poplist[[row]] <- new.grp_data[row,1:c]
      
      # flag the GAM quality status as pass (because pop was not GAM'd)
      gqfail[row] <- 0
      
      print(paste("copied population ", row, sep=""))
      
    } 
    
  } else {
    
    # otherwise, GAM everything
    gam_rows <- rows
    
  }
  

  for (row in gam_rows) {
    # Get pop data
    pop_data = subset(trim.mat, population == pop_ids[row])
    
    # check if there are any zeros
    if (length(which(pop_data$count==0))>0) {
      
      # check if all non-NA observations are zero
      if (mean(pop_data$count[which(!is.na(pop_data$count))])==0) {

        # if so, set zero adjust to a very small value
        zero_adjust <- 1e-17
        
      } else {
        
        # otherwise, calculate 1% of the mean of the observed values (excluding zeros)
        # if there are any zeros, this will be added to every observation to avoid issues with log of zero
        zero_adjust <- 0.01 * mean(pop_data$count[which(pop_data$count>0)], na.rm=TRUE)
        
      }
      
      # add the zero_adjust value
      pop_data$log_popvalue = log(pop_data$count + zero_adjust)
      
    } else {
      
      # add log column
      pop_data$log_popvalue = log(pop_data$count)
      
    }
    
    # K is half of the number of non-NA values in pop_data$count (if it is an odd number, K will be rounded down)
    K = round(length(which(!is.na(pop_data$count)))/2)

    # Make GAM
    #b <- gam(log_popvalue~s(year,k = K, bs = "gp"),data=pop_data)
    b <- gam(log_popvalue ~ s(year, k = K), data = pop_data)
    
    # if we are reproducing the LPI method...
    if (quality==TRUE) {
      
      # check the model fit
      # first, get the residuals of the GAM model
      resid <- residuals(b)
      
      # change K to the full number of non-NA values in pop_data$count 
      K = length(which(!is.na(pop_data$count)))

      # set years for residuals
      resid.years <- pop_data$year[which(!is.na(pop_data$count))]
      
      # then GAM the residuals (using same GAM settings as LPI)
      resid.gam <- gam(resid ~ s(resid.years, k = K, bs = "cs"), gamma = 1.4)
      
      # finally, check whether the sum of the estimated degrees of freedom is close to 1
      if ((abs(sum(resid.gam$edf) - 1)) < 0.01) {
        
        # flag the GAM quality status as a pass
        gqfail[row] <- 0
        
      } else {
        
        # copy the original data for using the chain method later
        gam_poplist[[row]] <- new.grp_data[row,1:c]
        
        # flag the GAM quality status as a fail
        gqfail[row] <- 1
        
        print(paste0("GAM of population ", row, " failed quality check.", sep=""))
        
        next
        
      }
      
      
    }
    
    # create matrix to hold GAM'd population
    pred.a <- matrix(NA, nrow=1, ncol=c)
      
    # add column names
    colnames(pred.a) = paste(m_colnames)
    
    if (resample==TRUE) {
      
      if (lambda==TRUE) {
        
        # predict all values between the first and last missing values using GAM
        # in this case we are interpolating using the GAM, as in the LPI
        startGAM <- min(which(!is.na(pop_data$count)))
        
        endGAM <- max(which(!is.na(pop_data$count)))
        
        # map coefs to fitted curves
        Xp <- predict(b, pop_data[startGAM:endGAM,], type="lpmatrix")
        
      } else if (lambda==FALSE) {
        
        # map coefs to fitted curves
        Xp <- predict(b, pop_data, type="lpmatrix")
        
      }
      
      # posterior mean and cov of coefs - (capture error of the model)
      beta <- coef(b);Vb <- vcov(b)
      
      # Samples from a multivariate normal distribution, where beta is the means of the variables (coefficients)
      # and Vb is the variance-covariance matrix of the coefficients
      br <- MASS::mvrnorm(n,beta,Vb) ## simulate n rep coef vectors from post.
      
      # convert to data frame
      pred.a <- as.data.frame(pred.a)
      
      # loop to get trough to peak diff for each sim
      for (i in 1:n) {
      
        if (lambda==TRUE) {
        
          # curve for this replicate
          pred.a[i,startGAM:endGAM] <- as.data.frame(t(Xp%*%br[i,]))
          
        } else if (lambda==FALSE) {
        
          # curve for this replicate
          pred.a[i,] <- as.data.frame(t(Xp%*%br[i,]))
          
        }
        
      }
      
    } else if (resample==FALSE) {
      
      if (lambda==TRUE) {
        
        # predict non-missing values using GAM
        # this is only to be used after log-linear interpolation has been done
        #pred.a[,which(!is.na(pop_data$count))] <- t(predict(b, pop_data[which(!is.na(pop_data$count)),]))
        
        # predict all values between the first and last missing values using GAM
        # in this case we are interpolating using the GAM, as in the LPI
        startGAM <- min(which(!is.na(pop_data$count)))
        
        endGAM <- max(which(!is.na(pop_data$count)))
        
        # predict the missing values from the GAM
        pred.a[,startGAM:endGAM] <- t(predict(b, pop_data[startGAM:endGAM,]))
        
      } else if (lambda==FALSE & forecast==TRUE) {
        
        # predict all values between the first and last missing values using GAM
        # in this case we are interpolating using the GAM, as in the LPI
        startGAM <- min(which(is.na(pop_data$count)))-1
        
        endGAM <- max(which(is.na(pop_data$count)))+1
        
        pred.a <- log(new.grp_data[row,1:c])
        
        # predict the missing values from the GAM
        pred.a[,startGAM:endGAM] <- t(predict(b, pop_data[startGAM:endGAM,]))
        
      } else if (lambda==FALSE & forecast==FALSE) {
        
        # if not using lambda method, GAM the whole time series
        pred.a <- t(predict(b, pop_data))
      
      }
      
    }
    
    # convert to matrix
    pred.a <- as.matrix(pred.a)
    
    # convert back to index values
    pred.a <- exp(pred.a)
    
    # convert any negative values to 0s
    pred.a <- ifelse(pred.a < 0, 0, pred.a)
    
    # convert to data frame
    pred.a <- as.data.frame(pred.a)
    
    # add column names
    colnames(pred.a) <- colnames(new.grp_data[,m_colnames])
    
    if (resample==TRUE) {
      
      # add resample ID
      #pred.a$ResID <- 1:nrow(pred.a)
      
      # add other IDs back in
      pred.a[,(c+1):(length(new.grp_data))] <- new.grp_data[row,(c+1):(length(new.grp_data))]
      
    }
    
    # add GAM'd population to list
    gam_poplist[[row]] <- pred.a
    
    print(paste("completed GAM of population ", row, sep=""))
    
  }
  
  # convert from list to data frame
  gam_popmat <- do.call(rbind, gam_poplist)
  
  if (is.null(gam_popmat)) {
    
    gam_popmat <- as.data.frame(matrix(NA, nrow=1, ncol=(ncol(new.grp_data)-1)))
    
    gqfail <- 0
    
  }
  
  if (resample==FALSE) {
    
    # add extra columns back in from the original data frame
    gam_popmat[,(c+1):(length(new.grp_data)-1)] <- new.grp_data[,(c+1):(length(new.grp_data)-1)]
    
    # add column names back in
    colnames(gam_popmat) <- colnames(new.grp_data[,1:(length(new.grp_data)-1)])
    
  } else {
    
    # add column names back in
    colnames(gam_popmat) <- colnames(new.grp_data[,1:(length(new.grp_data))])
    
  }
  

  
  # if we are reproducing the LPI method...
  if (quality==TRUE) {
    
    # add the GAM quality fail status to the data frame
    gam_popmat$gqfail <- gqfail
    
  }

  return(gam_popmat)
  
}



species_index_fn <- function(resamp_popmat, c, n=NA, n_boot=NA, lambda=FALSE, random_pops=TRUE, resample=FALSE) {
  
  # make a vector of all unique species IDs in the sample
  spec_ids <- unique(resamp_popmat$SpecID)
  
  # create list to hold the species indices
  spec_index.list <- list()
  
  # start a counter to iterate species
  counter1 <- 1
  
  if(!any(!is.na(resamp_popmat[,1:c])) & nrow(resamp_popmat)==1) {

    # get species data
    sample_mat5 <- as.data.frame(resamp_popmat[,1:c])

    # add a species ID column
    sample_mat5$SpecID <- spec_ids[1]

    # add a group ID column
    sample_mat5$GrpID <- resamp_popmat$GrpID[1]

    # add all copies of the species index to the list of species indices
    spec_index.list[[counter1]] <- sample_mat5

    # print information to show that a species has been completed
    print(paste("completed index for species ", counter1, sep=""))

    return(spec_index.list)

  }
  
  # loop to create species indices
  for (spec in spec_ids) {
    
    # select all populations that belong to a particular species
    spec_popdata <- subset(resamp_popmat, SpecID==spec)
    
    if (resample==FALSE) {
      
      # get species data
      sample_mat1.2 <- as.matrix(spec_popdata[,1:c])
      
      if (nrow(spec_popdata)==1) {
        
        if (lambda==FALSE) {
          
          # adjust all indices to begin at a base value of 100
          sample_mat1.5 <- (sample_mat1.2 / sample_mat1.2[,1] * 100)
          
          # restructure as data frame
          sample_mat5 <- as.data.frame(sample_mat1.5)
          
        } else if (lambda==TRUE) {
          
          # if using lambda method, convert to lambda values
          sample_mat1.3_pt1 <- sample_mat1.2[,1:(ncol(sample_mat1.2)-1)]
          sample_mat1.3_pt2 <- sample_mat1.2[,2:ncol(sample_mat1.2)]
          sample_mat1.5 <- sample_mat1.3_pt2 / (sample_mat1.3_pt1)
          
          # restructure as data frame
          sample_mat1.7 <- as.data.frame(t(sample_mat1.5))
          
          # convert time series values to natural log
          sample_mat1.8 <- log10(sample_mat1.7)
          
          # remove all change values with natural log values outside of -2.3:2.3 (log10 -1:1)
          #sample_mat5 <- as.data.frame(t(apply(sample_mat1.8, 2, function(i) {ifelse(i > -1, ifelse(i < 1, i, NA), NA)})))
          sample_mat5 <- as.data.frame(t(apply(sample_mat1.8, 2, function(i) {ifelse(i > -1, ifelse(i < 1, i, 1), -1)})))
          
        }
        
        # add years as column names
        colnames(sample_mat5) <- colnames(spec_popdata[,1:ncol(sample_mat5)])

        # add a species ID column
        sample_mat5$SpecID <- spec
        
        # add a group ID column
        sample_mat5$GrpID <- spec_popdata$GrpID[1]
        
      } else {
        
        # adjust all time series to begin at a base value of 100
        #sample_mat1.2 <- sample_mat / sample_mat[,1] * 100
        
        if (lambda==TRUE) {
          
          # if using lambda method, convert to lambda values
          sample_mat1.3_pt1 <- sample_mat1.2[,1:(ncol(sample_mat1.2)-1)]
          sample_mat1.3_pt2 <- sample_mat1.2[,2:ncol(sample_mat1.2)]
          sample_mat1.5 <- sample_mat1.3_pt2 / (sample_mat1.3_pt1)
          
        } 
        else if (lambda==FALSE) {
          
          # if not using lambda method, update the variable name
          sample_mat1.5 <- sample_mat1.2
          
        }
        
        # convert time series values to log10
        sample_mat1.7 <- log10(sample_mat1.5)
        
        if (lambda==TRUE) {
          
          # remove all change values with natural log values outside of -2.3:2.3 (log10 -1:1)
          #sample_mat1.8 <- apply(sample_mat1.7, 2, function(i) {ifelse(i > -1, ifelse(i < 1, i, NA), NA)})
          sample_mat1.8 <- apply(sample_mat1.7, 2, function(i) {ifelse(i > -1, ifelse(i < 1, i, 1), -1)})
          
        } else if (lambda==FALSE) {
          
          # update variable name
          sample_mat1.8 <- sample_mat1.7
          
        }

        # take mean of population index values at each time point
        sample_mat2 <- colMeans(sample_mat1.8, na.rm=TRUE)
        
        if (lambda==FALSE) {
          
          # convert back from natural log
          sample_mat3 <- 10^(sample_mat2)
          
          # if not using lambda method, set species index to start at 100
          sample_mat4 <- sample_mat3 / sample_mat3[1] * 100
          
        } else if (lambda==TRUE) {
          
          # convert NaN values to NA
          sample_mat2[is.nan(sample_mat2)] <- NA
          
          # update variable name
          sample_mat4 <- sample_mat2
          
        }
        
        # transpose species index as a matrix
        sample_mat5 <- t(as.matrix(sample_mat4))
        
        # add column names
        colnames(sample_mat5) <- colnames(spec_popdata[,1:ncol(sample_mat5)])
        
        # add rownames
        rownames(sample_mat5) <- rownames(sample_mat2)
        
        # restructure as data frame
        sample_mat5 <- as.data.frame(sample_mat5)
        
        # add a species ID column
        sample_mat5$SpecID <- spec
        
        # add a group ID column
        sample_mat5$GrpID <- spec_popdata$GrpID[1]
        
      }
      
    } 
    else if (resample==TRUE) {
      
      # if there is only one population in the species
      if (nrow(spec_popdata) <= n) {
        
        # randomly sample row numbers from all resampled time series x times, where x is specified by n_boot
        sample_nums <- sample(1:nrow(spec_popdata), n_boot, replace = TRUE)
        
        sample_mat <- spec_popdata[match(sample_nums, row(spec_popdata)),1:c]
        
        if (lambda==FALSE) {
          
          # adjust all indices to begin at a base value of 100
          sample_mat1.2 <- (sample_mat / sample_mat[,1] * 100)
          
          # restructure as data frame
          sample_mat5 <- as.data.frame(sample_mat1.2)
          
        } else if (lambda==TRUE) {
          
          # if using lambda method, convert to lambda values
          sample_mat1.1_pt1 <- sample_mat[,1:(ncol(sample_mat)-1)]
          sample_mat1.1_pt2 <- sample_mat[,2:ncol(sample_mat)]
          sample_mat1.2 <- sample_mat1.1_pt2 / (sample_mat1.1_pt1)
          
          # restructure as data frame
          sample_mat1.5 <- as.data.frame(sample_mat1.2)
          
          # convert time series values to log10
          sample_mat1.7 <- log10(sample_mat1.5)
          
          # remove all change values with log10 values outside of -1:1
          #sample_mat5 <- apply(sample_mat1.7, 2, function(i) {ifelse(i > -1, ifelse(i < 1, i, NA), NA)})
          sample_mat5 <- apply(sample_mat1.7, 2, function(i) {ifelse(i > -1, ifelse(i < 1, i, 1), -1)})
          
        }
        
        # add years as column names
        colnames(sample_mat5) <- colnames(spec_popdata[,1:ncol(sample_mat5)])
        
        # restructure as data frame
        sample_mat5 <- as.data.frame(sample_mat5)
        
        # add a species ID column
        sample_mat5$SpecID <- spec
        
        # add a group ID column
        sample_mat5$GrpID <- spec_popdata$GrpID[1]
        
      } else {
        
        if (random_pops==TRUE) {
          
          ### bootstrapping the species index, with replacement ###
          
          # randomly sample from all resampled populations x times, where x is specified by n_boot times the number
          # of populations that belong to the species
          sample_nums <- sample(1:nrow(spec_popdata), n_boot*nrow(spec_popdata) / n, replace = TRUE)
          
          # get data from spec_popdata by matching to row numbers in sample_nums
          sample_mat <- spec_popdata[match(sample_nums, row(spec_popdata)),]
          
        } else if (random_pops==FALSE) {
          
          ### bootstrapping the species index, with replacement ###
          
          # make a vector of all unique population IDs in the species
          pop_ids <- unique(spec_popdata$PopID)
          
          # create a counter to keep track of the number of populations sampled
          p_counter <- 1
          
          # create a vector to hold time series row numbers
          sample_mat <- data.frame()
          
          for (pop in pop_ids) {
            
            # select all resampled time series that belong to a particular population
            popdata <- subset(spec_popdata, PopID==pop)
            
            # randomly sample row numbers from all resampled time series x times, where x is specified by n_boot
            sample_nums <- sample(1:nrow(popdata), n_boot, replace = TRUE)
            
            # get data from popdata by matching to row numbers in sample_nums
            sample_mat[(1+(p_counter*n_boot)-n_boot):(p_counter*n_boot),1:ncol(popdata)] <- popdata[match(sample_nums, row(popdata)),]
            
            # advance the population counter
            p_counter <- p_counter + 1
          }
          
        }
        
        
        # remove the id, popid and specid columns from sample_mat
        sample_mat1.2 <- sample_mat[,1:c]
        
        # adjust all indices to begin at a base value of 100
        #sample_mat1.2 <- (sample_mat1.1 / sample_mat1.1[,1] * 100)
        
        if (lambda==TRUE) {
          
          # if using lambda method, convert to lambda values
          sample_mat1.3_pt1 <- sample_mat1.2[,1:(ncol(sample_mat1.2)-1)]
          sample_mat1.3_pt2 <- sample_mat1.2[,2:ncol(sample_mat1.2)]
          sample_mat1.5 <- sample_mat1.3_pt2 / (sample_mat1.3_pt1)
          
        } else if (lambda==FALSE) {
          
          # if not using lambda method, update the variable name
          sample_mat1.5 <- sample_mat1.2
          
        }
        
        # convert to natural log and convert to matrix format for faster processing
        sample_mat1.7 <- log10(as.matrix(sample_mat1.5))
        
        if (lambda==TRUE) {
          
          # remove all change values with log10 values outide of -1:1
          #sample_mat1.75 <- apply(sample_mat1.7, 2, function(i) {ifelse(i > -1, ifelse(i < 1, i, NA), NA)})
          sample_mat1.75 <- apply(sample_mat1.7, 2, function(i) {ifelse(i > -1, ifelse(i < 1, i, 1), -1)})
          
        } else if (lambda==FALSE) {
          
          sample_mat1.75 <- sample_mat1.7
          
        }

        # transpose data and convert to vector
        sample_mat1.8 <- as.vector(t(sample_mat1.75))
        
        # create a new matrix in wide format, so column means can be used to create species indices
        sample_mat1.9 <- matrix(sample_mat1.8, nrow=length(unique(spec_popdata$PopID)), 
                                ncol=(n_boot*ncol(sample_mat1.7)), byrow=TRUE)
        
        # take the column means. the .colMeans function requires extra information but is very fast
        sample_mat2 <- .colMeans(sample_mat1.9, nrow(sample_mat1.9), ncol(sample_mat1.9), na.rm=TRUE)
        
        # convert back to long format, with a species index in each row
        sample_mat3 <- matrix(sample_mat2, nrow = n_boot, ncol = ncol(sample_mat1.7), byrow=TRUE)
        
        if (lambda==FALSE) {
          
          # convert back from natural log to index values
          sample_mat3 <- 10^(sample_mat3)
          
          # set each species index to start at 100
          sample_mat4 <- (sample_mat3 / sample_mat3[,1] * 100)
          
        } else if (lambda==TRUE) {
          
          # convert NaN values to NA
          sample_mat3[is.nan(sample_mat3)] <- NA
          
          # update variable name
          sample_mat4 <- sample_mat3
          
        } 
        
        # convert to a data frame
        sample_mat5 <- as.data.frame(sample_mat4)
        
        # add years as column names
        colnames(sample_mat5) <- colnames(spec_popdata[,1:ncol(sample_mat5)])
        
        # add a species id column
        sample_mat5$SpecID <- rep(spec, times=nrow(sample_mat5))
        
        # add a group id column
        sample_mat5$GrpID <- rep(spec_popdata$GrpID[1], times=nrow(sample_mat5))
        
      }
      
    }
    
    # add all resampled copies of the species index to the list of species indices
    spec_index.list[[counter1]] <- sample_mat5
    
    # print information to show that a species has been completed
    print(paste("completed index for species ", counter1, sep=""))
    
    # increase counter value to keep track of the number of species completed
    counter1 <- counter1 + 1
    
  }
  
  return(spec_index.list)
  
}


group_index_fn <- function(grp_index.list, c, m_colnames, n=NA, n_boot=NA, weights=NA, stay_lambda=FALSE) {
  
  # check if a value is provided for n
  if (!is.na(n)) {
    
    # check if the number of rows for each group index is equal to n_boot
    if (nrow(grp_index.list[[1]]) == n_boot) {
      
      # if so, set the resample flag to TRUE
      resample<-TRUE
      
    } else {
      
      # otherwise something is wrong, so stop with an error message
      stop("function parameters are not set correctly")
      
    }
    
  } else {
    
    # if no value is provided for n, set the resample flag to FALSE
    resample<-FALSE
    
  }

  # check whether data is in list format
  if (inherits(grp_index.list, "list")) {
    
    # if so, convert group index list to data frame
    grp_index <- do.call(rbind, grp_index.list)
    
  } else {
    
    # otherwise, copy group index directly
    grp_index <- grp_index.list
    
  }
    
  # check whether the final column from the original dataset is missing
  # if so, turn the lambda FLAG on. As lambdas are year-to-year change values,
  # there should one fewer lambda value than the number of observation years
  if (!(tail(m_colnames, 1) %in% colnames(grp_index))) {

    lambda<-TRUE

  # check whether there are weightings provided
  # if so, turn the lambda flag on
  } else if (!is.na(weights)) {
    
    lambda<-TRUE 
    
  } else {
    
    lambda<-FALSE
    
  }
  
  # check for species ID column. if it exists, set the group flag to TRUE
  # because the data contain species indices intended to form a group index
  # otherwise set group flag to FALSE, as the intention is to make an msi
  if ("SpecID" %in% colnames(grp_index) & "GrpID" %in% colnames(grp_index)) {
    
    group<-TRUE
    
  } else {
    
    group<-FALSE
    
  }
  
  if (nrow(grp_index)==1 & !any(!is.na(grp_index))) {
    
    # if in lambda format...
    if (!(tail(m_colnames, 1) %in% colnames(grp_index))) {
      
      # remove id tags
      grp_index1.1 <- grp_index[,1:(c-1)]
      
      # if NOT tagged to stay in lambda format
      if (stay_lambda==FALSE) {
        
        # convert to index values
        grp_index1.5 <- append(grp_index1.1, 100, after = 0)
        grp_index2 <- cumprod(grp_index1.5)
        
        # if tagged to stay in lambda format
      } else if (stay_lambda==TRUE) {
        
        # update variable name
        grp_index2 <- grp_index1.1
        
      }

      
    # if tagged to stay in lambda format
    } else if (stay_lambda==TRUE) {
      
      # remove id tags
      grp_index2 <- grp_index[,1:(c-1)]
      
    } else {
      
      # remove id tags
      grp_index2 <- grp_index[,1:c]
      
    }
    
    if (group==TRUE) {
      
      # add a group id column
      grp_index2$GrpID <- grp_index$GrpID[1]
      
      # create list to hold completed group indices
      grp.list <- list()
      
      # add completed group index to group index list
      grp.list[[1]] <- grp_index2
      
      # print information to show that a species has been completed
      print(paste("completed index for group 1"))

    } else {
      
      print(paste("completed msi"))
      
      grp.list <- grp_index2
      
    }
    
    return(grp.list)
    
  }
  
  if (lambda==TRUE & resample==TRUE) {
    
    # if the lambda flag and the resample flag are both set
    # create a single-column matrix of 100s for later conversion to index values
    first_col <- matrix(100, nrow = n_boot, ncol = 1)
    
  }
  
  if (group==TRUE) {
    
    # create vector of group ids
    grp_ids <- unique(grp_index$GrpID)
    
    # set counter to track completed group indices
    counter1 <- 1
    
    # create list to hold completed group indices
    grp.list <- list()
    
  } else {
    
    grp_ids <- 1
    
  }
  
  for (grp in grp_ids) {
    
    if (group==TRUE) {
      
      # select all species indices that belong to a particular group
      ms_grpdata <- subset(grp_index, GrpID==grp)
      
    } else {
      
      # otherwise, copy all data
      ms_grpdata <- grp_index
      
    }
    
    if (lambda==TRUE) {
        
      # remove ID tags
      ms_grpdata1.4 <- ms_grpdata[,1:(c-1)]
        
    } else {
      
      # remove ID tags
      ms_grpdata1.1 <- ms_grpdata[,1:c]
      
      # convert to log10
      ms_grpdata1.4 <- log10(ms_grpdata1.1)
      
    }
  

    if (resample==FALSE) {
      
      # if weights are provided
      if (!is.na(weights)) {
        
        # create a list to hold weighted values
        w_grpdata <- list()
        
        # apply weighting to the group lambdas
        for (i in 1:length(weights)) {
          
          temp <- ms_grpdata1.4[i,] * weights[i]
          
          w_grpdata[[i]] <- temp
          
        }
        
        # convert to data frame
        ms_grpdata1.4 <- as.data.frame(do.call(rbind, w_grpdata))
        
        # take the sum of the group lambdas
        ms_grpdata1.5 <- colSums(ms_grpdata1.4, na.rm=TRUE)

      } else {
        
        if (lambda==TRUE) {
 
          # take the mean of the group lambdas
          ms_grpdata1.5 <- colMeans(ms_grpdata1.4, na.rm=TRUE)

        } else {

          # take the mean of the group lambdas
          ms_grpdata1.5 <- colMeans(ms_grpdata1.4, na.rm=TRUE)
          
        }

      }

      if (lambda==TRUE & stay_lambda==FALSE) {
        
        # back convert from natural log
        ms_grpdata1.6 <- 10^(ms_grpdata1.5)
        
        # convert to index values
        ms_grpdata1.7 <- append(ms_grpdata1.6, 100, after = 0)
        ms_grpdata1.8 <- cumprod(ms_grpdata1.7)
        
      } else if (lambda==FALSE & stay_lambda==FALSE) {

        # back convert from natural log
        ms_grpdata1.6 <- 10^(ms_grpdata1.5)
        
        # ensure index starts at 100
        ms_grpdata1.8 <- (ms_grpdata1.6 / ms_grpdata1.6[1] * 100)
      
      } else if (lambda==TRUE & stay_lambda==TRUE) {
        
        # update variable name
        ms_grpdata1.8 <- ms_grpdata1.5
        
      } else if (lambda==FALSE & stay_lambda==TRUE) {
        
        # update variable name
        ms_grpdata1.8 <- ms_grpdata1.5[1:(length(ms_grpdata1.5)-1)]
        
      }
    
      # transpose group index as a matrix
      ms_grpdata2.1 <- t(as.matrix(ms_grpdata1.8))
      
      # convert NaN values to NA
      ms_grpdata2.1[is.nan(ms_grpdata2.1)] <- NA
    
      if (stay_lambda==FALSE) {
        
        # add column names
        colnames(ms_grpdata2.1) <- m_colnames
        
      } else {
        
        if (any(!is.na(ms_grpdata2.1))) {
          
          # add column names
          colnames(ms_grpdata2.1) <- m_colnames[1:(length(m_colnames)-1)]
          
        }
        
      }

      # restructure as data frame
      ms_grpdata2.1 <- as.data.frame(ms_grpdata2.1)
      
      if (group==TRUE) {
        
        # add a group id column
        ms_grpdata2.1$GrpID <- rep(ms_grpdata$GrpID[1], times=nrow(ms_grpdata2.1))
        
        # add completed group index to group index list
        grp.list[[counter1]] <- ms_grpdata2.1
        
        # print information to show that a species has been completed
        print(paste("completed index for group ", counter1, sep=""))
        
        # advance counter
        counter1 <- counter1 + 1
        
      } else {
        
        print(paste("completed msi"))
        
        grp.list <- ms_grpdata2.1
        
      }

    } else if (resample==TRUE) {
    
      # transpose data and convert to vector
      ms_grpdata1.5 <- as.vector(t(ms_grpdata1.4))
      
      if (group==TRUE) {
        
        # create a new matrix in wide format, so column means can be used to create group indices
        ms_grpdata1.6 <- matrix(ms_grpdata1.5, nrow=length(unique(ms_grpdata$SpecID)), 
                                ncol=(n_boot*ncol(ms_grpdata1.4)), byrow=TRUE)
        
      } else {
    
        # create a new matrix in wide format, so column means can be used to create group indices
        ms_grpdata1.6 <- matrix(ms_grpdata1.5, nrow=length(unique(grp_index$GrpID)), 
                                ncol=(n_boot*ncol(ms_grpdata1.4)), byrow=TRUE)
      
      }
    
      # take the column means. the .colMeans function requires extra information but is very fast
      ms_grpdata1.7 <- .colMeans(ms_grpdata1.6, nrow(ms_grpdata1.6), ncol(ms_grpdata1.6), na.rm=TRUE)
    
      # convert back to long format, with a species index in each row
      ms_grpdata1.8 <- matrix(ms_grpdata1.7, nrow = n_boot, ncol = ncol(ms_grpdata1.4), byrow=TRUE)
      
      # convert NaN values to NA
      ms_grpdata1.8[is.nan(ms_grpdata1.8)] <- NA
    
      if (lambda==TRUE & stay_lambda==FALSE) {
      
        # back convert from natural log
        ms_grpdata1.9 <- 10^(ms_grpdata1.8)
        
        # convert to index values
        ms_grpdata2 <- cbind(first_col, ms_grpdata1.9)
        ms_grpdata2.1 <- rowCumprods(ms_grpdata2)
      
      } else if (lambda==FALSE & stay_lambda==FALSE) {
      
        # back convert from natural log
        ms_grpdata1.9 <- 10^(ms_grpdata1.8)
        
        # set each index to start at 100
        ms_grpdata2.1 <- (ms_grpdata1.9 / ms_grpdata1.9[,1] * 100)
      
      } else if (lambda==TRUE & stay_lambda==TRUE) {
        
        # update variable name
        ms_grpdata2.1 <- ms_grpdata1.8
        
      } else if (lambda==FALSE & stay_lambda==TRUE) {
        
        # update variable name
        ms_grpdata2.1 <- ms_grpdata1.8[,(1:length(ms_grpdata1.8)-1)]
        
      }
      
      # convert NaN values to NA
      ms_grpdata2.1[is.nan(ms_grpdata2.1)] <- NA

      # convert to a data frame
      ms_grpdata2.1 <- as.data.frame(ms_grpdata2.1)
      
      if (stay_lambda==FALSE) {
        
        
        # add years as column names
        colnames(ms_grpdata2.1) <- m_colnames
        
      } else {
      
        if (any(!is.na(ms_grpdata2.1))) {
          
          # add years as column names
          colnames(ms_grpdata2.1) <- m_colnames[1:(length(m_colnames)-1)]
          
        }

      }
    
      # take mean of all msi bootstraps
      #ms_index <- colMeans(ms_grpdata2.1, na.rm=TRUE)
      
      if (group==TRUE) {
        
        # add a group id column
        ms_grpdata2.1$GrpID <- rep(ms_grpdata$GrpID[1], times=nrow(ms_grpdata2.1))
        
        # add completed group index to group index list
        grp.list[[counter1]] <- ms_grpdata2.1
        
        # print information to show that a species has been completed
        print(paste("completed index for group ", counter1, sep=""))
        
        # advance counter
        counter1 <- counter1 + 1
        
      } else {
        
        print(paste("completed msi"))
        
        grp.list <- ms_grpdata2.1
        
      }
      
    }
    
  }
  
  return(grp.list)
  
}

# note: grouplevel = 1 to 4 for LPI, 1 to 3 for sampled LPI
ci_fn <- function(index, c, m_colnames, boots=100, n=NA, weights=NA, lambda=FALSE, savedata=FALSE, grouplevel) {
  
  if (grouplevel==1) {
    
    # check whether data is in list format
    if (inherits(index, "list")) {
      
      # if it is, convert the list to a data frame
      index.matrix <- do.call(rbind, index)
      
      # check whether there is more than 1 row per list element
      if (nrow(index[[1]]) > 1) {
        
        # if there is, set resample flag to TRUE
        resample<-TRUE
        
      } else {
        
        # if not, set resample flag to FALSE
        resample<-FALSE
        
      }
      
    } else {
      
      # if it is not a list, copy the data frame directly
      index.matrix <- index
      
      # check whether there is more than 1 row in the data frame
      if (!is.na(n)) {
        
        # if there is, set resample flag to TRUE
        resample<-TRUE
        
      } else {
        
        # if not, set resample flag to FALSE
        resample<-FALSE
        
      }
      
    }
    
    if (resample==FALSE) {
      
      # check whether there is a species id column AND a group id column
      if ("SpecID" %in% colnames(index.matrix) & "GrpID" %in% colnames(index.matrix)) {
        
        # if so, set the group flag to TRUE
        group<-TRUE
        
      } else {
        
        # if not, set the group flag to FALSE
        group<-FALSE
        
      }
      
      if (group==TRUE) {
        
        # check whether the data frame contains only a row of NAs (no data)
        if (nrow(index.matrix)==1 & !any(!is.na(index.matrix))) {
          
          # create matrix to hold confidence intervals
          grp_ci.list <- list()
          
          # upper and lower confidence intervals will be the same as the index
          grp_samp_ci1 <- index.matrix[,1:c]
          grp_samp_ci2 <- index.matrix[,1:c]
          
          # bind the confidence intervals together into a data frame
          grp_samp_ci <- as.data.frame(rbind(grp_samp_ci1, grp_samp_ci2))
          
          # add year column names
          colnames(grp_samp_ci) <- m_colnames
          
          # add a group id column
          grp_samp_ci$GrpID <- rep(index.matrix$GrpID[1], times=2)
          
          # add completed group index to group index list
          grp_ci.list[[1]] <- grp_samp_ci
          
          # print information to show that a species has been completed
          print(paste("completed confidence intervals for group 1"))
          
          return(grp_ci.list)
          
        }
        
        # if there is data...
        
        # create vector of group ids
        grp_ids <- unique(index.matrix$GrpID)
        
        # set counter to track completed group indices
        counter1 <- 1
        
        # create list to hold completed group indices
        grp_ci.list <- list()
        
        # if group flag is FALSE
      } else {
        
        # group ID set to 1
        grp_ids <- 1
        
      }
      
      # loop through groups
      for (grp in grp_ids) {
        
        if (group==TRUE) {
          
          # select all species indices that belong to a particular group
          grp_specdata <- subset(index.matrix, GrpID==grp)
          
          # if group flag is FALSE
        } else {
          
          # select the whole dataset
          grp_specdata <- index.matrix
          
        }
        
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
          
          # check if the lambda flag is set to TRUE (LPI mode)
          if (lambda==TRUE) {
            
            # if so, convert to lambda values
            grp_specdata1.2_pt1 <- grp_specdata1.1[,1:(ncol(grp_specdata1.1)-1)]
            grp_specdata1.2_pt2 <- grp_specdata1.1[,2:ncol(grp_specdata1.1)]
            grp_specdata1.3 <- grp_specdata1.2_pt2 / (grp_specdata1.2_pt1)
            
            # if lambda flag is set to FALSE
          } else {
            
            # update variable name
            grp_specdata1.3 <- grp_specdata1.1
            
          }
          
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
          
          # if there are no weights
          if (is.na(weights)) {
            
            # take the mean of the natural log of the species lambdas
            grp_samp_temp1.4 <- colMeans(grp_samp_temp, na.rm=TRUE)
            
          }
          
          # if there are weights...
          # convert NaN values to NA
          grp_samp_temp1.4[is.nan(grp_samp_temp1.4)] <- NA
          
          # if savedata flag is set to TRUE
          if (savedata==TRUE) {
            
            # add bootstrapped sample to list
            grp_samp_list[[i]] <- grp_samp_temp1.4
            
            # if savedata flag is set to FALSE
          } else {
            
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
          
        }
        
        # convert list to a data frame
        grp_samp <- do.call(rbind, grp_samp_list)
        
        # if savedata flag is set to TRUE
        if (savedata==TRUE) {
          
          # add years as column names
          colnames(grp_samp) <- m_colnames
          
          # restructure as a data frame
          grp_samp <- as.data.frame(grp_samp)
          
          # if group flag is set to TRUE
          if (group==TRUE) {
            
            # add a group id column
            grp_samp$GrpID <- rep(grp_specdata$GrpID[1], times=nrow(grp_samp))
            
            # add data frame to confidence interval list
            grp_ci.list[[counter1]] <- grp_samp
            
            # print information to show that a group index has been completed
            print(paste("saved confidence interval data for group ", counter1, sep=""))
            
            # update counter
            counter1 <- counter1 + 1
            
            # if this is the final group
            if (grp==tail(grp_ids, 1)) {
              
              # exit the function
              return(grp_ci.list)
              
              # if it is not the final group
            } else {
              
              # proceed to the next group
              next
              
            }
            
            # if group flag is set to FALSE
          } else{
            
            # rename as confidence interval list
            grp_ci.list <- grp_samp
            
            # exit function
            return(grp_ci.list)
            
          }
          
        }
        
        # if savedata is FALSE...
        
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
        
        # if group flag set to TRUE
        if (group==TRUE) {
          
          # add a group id column
          grp_samp_ci$GrpID <- rep(grp_specdata$GrpID[1], times=2)
          
          # add completed group index to group index list
          grp_ci.list[[counter1]] <- grp_samp_ci
          
          # print information to show that a species has been completed
          print(paste("completed confidence intervals for group ", counter1, sep=""))
          
          # advance counter
          counter1 <- counter1 + 1
          
          # if group flag is set to FALSE  
        } else {
          
          # rename as confidence interval list
          grp_ci.list <- grp_samp_ci
          
          # print information to show that the msi has been completed
          print(paste("completed confidence intervals for msi"))
          
        }
        
      }
      
      # exit function
      return(grp_ci.list)
      
      # if resample is set to TRUE
    } else if (resample==TRUE) {
      
      # load GET package
      library(GET)
      
      # check for group ID column
      if ("GrpID" %in% colnames(index.matrix)) {
        
        # if there is one, set group flag to TRUE
        group<-TRUE
        
        # if there is not one
      } else {
        
        # set group flag to FALSE
        group<-FALSE
        
      }
      
      # if group flag set to TRUE
      if (group==TRUE) {
        
        # create matrix to hold confidence intervals
        grp_ci <- matrix(NA, nrow=2, ncol=(ncol(index.matrix)-1))
        
        # create vector of group ids
        grp_ids <- unique(index.matrix$GrpID)
        
        # set counter to track completed group indices
        counter1 <- 1
        
        # create list to hold completed group indices
        grp_ci.list <- list()
        
        # if group flag set to FALSE
      } else {
        
        # create matrix to hold confidence intervals
        grp_ci <- matrix(NA, nrow=2, ncol=ncol(index.matrix))
        
        # set group id as 1
        grp_ids <- 1
        
      }
      
      # loop through groups
      for (grp in grp_ids) {
        
        # if group flag is set to TRUE
        if (group==TRUE) {
          
          # select all species indices that belong to a particular group, omitting NAs
          grp.data <- na.omit(subset(index.matrix, GrpID==grp)[,(1:ncol(index.matrix)-1)])
          
          # if group flag is set to FALSE
        } else {
          
          # select the whole dataset, omitting NAs
          grp.data <- na.omit(index.matrix)
          
        }
        
        # calculate final msi from bootstraps
        grp.final <- colMeans(grp.data, na.rm = TRUE)
        
        ## rank envelope method
        
        # convert msi into a vector
        grp.final.vec <- as.vector(grp.final)
        
        # transpose matrix of msi bootstraps
        grp.data.t <- t(grp.data)
        
        # create a list for the create_curve_set function
        c1 <- list(m_colnames, grp.final.vec, grp.data.t)
        
        # name the list elements appropriately for the function  
        names(c1) <- c("r", "obs", "sim_m")
        
        # create curve set for the rank envelope function
        c2 <- create_curve_set(c1)
        
        # create confidence intervals from the curve set
        res <- rank_envelope(c2)
        
        # extract the lower bounding confidence interval
        grp_ci[1,] <- res$lo
        
        # extract the upper bounding confidence interval
        grp_ci[2,] <- res$hi
        
        # restructure as data frame for export
        grp_ci <- as.data.frame(grp_ci)
        
        # restore column names
        colnames(grp_ci) <- m_colnames
        
        # if group flag is set to TRUE
        if (group==TRUE) {
          
          # add a group id column
          grp_ci$GrpID <- rep(grp.data$GrpID[1], times=2)
          
          # add group confidence intervals to list
          grp_ci.list[[counter1]] <- grp_ci
          
          # print information to show that a group index has been completed
          print(paste("completed confidence intervals for group ", counter1, sep=""))
          
          # update counter
          counter1 <- counter1 + 1
          
          # if group flag is set to FALSE
        } else {
          
          # rename as confidence interval list
          grp_ci.list <- grp_ci
          
          # print information to show that the msi has been completed
          print(paste("completed confidence intervals for msi"))
          
        }
        
      }
      
      # export the confidence intervals
      return(grp_ci.list)
      
    }
    
    # if grouplevel is set to 2
  } else if (grouplevel==2) {
    
    # convert the list to a data frame
    index.matrix <- do.call(rbind, index)
    
    # get group ids
    grp_ids <- unique(index.matrix$GrpID)
    
    # create matrix to hold sampled species lambdas/indices
    grp_samp_list <- list()
    
    # check whether n is set to NA
    if (!is.na(n)) {
      
      # if it isn't, set resample flag to TRUE
      resample<-TRUE
      
    } else {
      
      # if it is, set resample flag to FALSE
      resample<-FALSE
      
    }
    
    # if weights are provided
    if (!is.na(weights)) {
      
      # check if the correct number of weights are provided
      if (length(weights) != length(grp_ids)) {
        
        # if not, stop with an error message
        stop('Length of weightings vector does not match number of indices to weight.', call.=FALSE)
        
      }
      
      # create a list to hold weighted values
      w_grpdata <- list()
      
      # apply weighting to the group lambdas
      for (i in 1:length(grp_ids)) {
        
        # copy all rows from the group
        temp <- index.matrix[index.matrix$GrpID %in% grp_ids[i],]
        
        # weight the rows
        temp2 <- temp[,1:c] * weights[i]
        
        # add weighted rows to list
        w_grpdata[[i]] <- temp2
        
      }
      
      
      # bind list into a dataframe
      w_grpdata.df <- do.call(rbind, w_grpdata)
      
      if (resample==TRUE) {
        
        # create a new matrix in wide format, so column sums can be used
        w_grpdata.df2 <- matrix(w_grpdata.df, nrow=length(grp_ids), 
                                ncol=(ncol(w_grpdata.df) * nrow(w_grpdata.df) / length(grp_ids)), byrow=TRUE)
        
        # take the column sums. the .colSums function requires extra information but is very fast
        w_grpdata.df3 <- .colSums(w_grpdata.df2, nrow(w_grpdata.df2), ncol(w_grpdata.df2), na.rm=TRUE)
        
        # convert back to long format, with a species index in each row
        w_grpdata.df4 <- matrix(w_grpdata.df3, nrow = (nrow(w_grpdata.df) / length(grp_ids)), ncol = ncol(w_grpdata.df), byrow=TRUE)
        
      } else {
        
        # update name and convert to matrix
        w_grpdata.df4 <- matrix(w_grpdata.df)
        
      }
      
      # convert NaN values to NA
      w_grpdata.df4[is.nan(w_grpdata.df4)] <- NA
      
      # if no weights are provided
    } else {
      
      # check if these are already lambda values
      # if they are, there will be one less column than in m_colnames
      # so we test if the final element of m_colnames is missing from the column names
      if (!(tail(m_colnames, 1) %in% colnames(index.matrix))) {
        
        # copy unweighted dataframe as a matrix
        w_grpdata.df <- as.matrix(index.matrix[,1:(c-1)])
        
      } else {
        
        # copy unweighted dataframe as a matrix
        w_grpdata1.1 <- as.matrix(index.matrix[,1:c])
        
        # check if the lambda flag is set to TRUE (LPI mode)
        if (lambda==TRUE) {
          
          # if so, convert to lambda values
          w_grpdata1.2_pt1 <- w_grpdata1.1[,1:(ncol(w_grpdata1.1)-1)]
          w_grpdata1.2_pt2 <- w_grpdata1.1[,2:ncol(w_grpdata1.1)]
          w_grpdata1.3 <- w_grpdata1.2_pt2 / (w_grpdata1.2_pt1)
          
          # if lambda flag is set to FALSE
        } else {
          
          # update variable name
          w_grpdata1.3 <- w_grpdata1.1
          
        }
        
        # convert to log10
        w_grpdata.df <- log10(w_grpdata1.3)
        
      }
      
    }
    
    
    if (resample==TRUE) {
      
      # create a new matrix in wide format, so column means can be used
      w_grpdata.df2 <- matrix(w_grpdata.df, nrow=length(grp_ids), 
                              ncol=(ncol(w_grpdata.df) * nrow(w_grpdata.df) / length(grp_ids)), byrow=TRUE)
      
      # take the column sums. the .colMeans function requires extra information but is very fast
      w_grpdata.df3 <- .colMeans(w_grpdata.df2, nrow(w_grpdata.df2), ncol(w_grpdata.df2), na.rm=TRUE)
      
      # convert back to long format, with a species index in each row
      w_grpdata.df4 <- matrix(w_grpdata.df3, nrow = (nrow(w_grpdata.df) / length(grp_ids)), ncol = ncol(w_grpdata.df), byrow=TRUE)
      
    } else {
      
      # update name
      w_grpdata.df4 <- w_grpdata.df
      
    }
    
    # convert NaN values to NA
    w_grpdata.df4[is.nan(w_grpdata.df4)] <- NA
    
    # }
    
    # if savedata flag set to TRUE
    if (savedata==TRUE) {
      
      # add bootstrapped sample to list
      grp_samp_list[[i]] <- w_grpdata.df4
      
      # if savedata flag set to FALSE
    } else {
      
      # back convert from log10
      grp_samp_temp1.5 <- 10^(w_grpdata.df4)
      
      # if lambda flag set to TRUE (lpi method)
      if (lambda==TRUE) {
        
        # create a single-column matrix of 100s for later conversion to index values
        first_col <- matrix(100, nrow = nrow(grp_samp_temp1.5), ncol = 1)
        
        # convert to index values
        grp_samp_temp1.7 <- cbind(first_col, grp_samp_temp1.5)
        grp_samp_temp1.8 <- rowCumprods(grp_samp_temp1.7)
        
        # convert to index values
        #grp_samp_temp1.7 <- append(grp_samp_temp1.5, 100, after = 0)
        #grp_samp_temp1.8 <- cumprod(grp_samp_temp1.7)
        
        # if lambda flag set to FALSE
      } else {
        
        # set each index to start at 100
        grp_samp_temp1.8 <- (grp_samp_temp1.5 / grp_samp_temp1.5[1] * 100)
        
      }
      
      # add bootstrapped sample to list
      grp_samp_list[[i]] <- grp_samp_temp1.8
      
    }
    
    # convert list to a data frame
    grp_samp <- do.call(rbind, grp_samp_list)
    
    # if savedata flag is set to TRUE
    if (savedata==TRUE) {
      
      # add years as column names
      colnames(grp_samp) <- m_colnames
      
      # restructure as a data frame
      grp_samp <- as.data.frame(grp_samp)
      
      # rename as confidence interval list
      grp_ci.list <- grp_samp
      
      # exit function
      return(grp_ci.list)
      
    }
    
    # if savedata is FALSE...
    
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
    
    return(grp_ci.list)
    
  }
  
}

# function, with x and y as time series, and alpha as moving avg window length
devdist <- function(x, y, alpha=10) {
  
  # take difference of two time series
  z = x - y
  
  # set limits for moving average
  dist <- alpha - 1
  
  # create vector to hold standard deviation values
  t <- vector()
  
  # loop to get moving average of the standard deviation
  for (i in 1:(length(z)-dist)) {
    
    t[i] <- sd(z[i:(i+dist)])
    
  }
  
  # take the mean of the standard deviation values
  dev_val <- mean(t)
  
  # return the mean
  return(dev_val)
  
}


# find the average width of the confidence intervals
ci_width_fn <- function(msi_ci, c) {
  
  # subtract lower c.i. from upper c.i. for all time points and take the mean
  width <- mean(as.matrix(msi_ci[2,2:c] - msi_ci[1,2:c]), na.rm=TRUE)
  
  if (is.nan(width)) {
    
    width <- NA
    
  }
  
  return(width)
  
}


# test how much of the real trend is within the confidence intervals of the sampled trend
real_trend_within_ci_fn <- function(msi_ci, true) {
  
  # create a vector to count how many times the true trend falls outside of the confidence intervals
  false.count <- 0
  
  # loop over each time point
  for (i in 1:length(msi_ci)) {
    
    # put upper confidence interval into a vector
    test.upper <- msi_ci[2,][i]
    
    # put lower confidence interval into a vector
    test.lower <- msi_ci[1,][i]
    
    # put true trend into a vector
    test <- true[i]
    
    # check for NAs in confidence intervals
    if (!is.na(test.upper) & !is.na(test.lower)) {
      
      # check if the true trend is above the upper c.i.
      if (test > test.upper) {
        
        # if so, increase the count
        false.count <- false.count + 1
        
        # in this case we do not need to test the lower c.i., so move to next time point
        next
        
      # check if the true trend is below the lower c.i.
      } else if (test.lower > test) {

        # if so, increase the count
        false.count <- false.count + 1

      }
      
    # if there is an NA value in the upper or lower confidence interval...  
    } else {
      
      # move to the next time point
      next
      
    }

  }
  
  # calculate the percentage of non-NA time points where the true trend was within the c.i.'s
  # we remove 1 from the length because the first time point has no confidence intervals
  p.cent <- 100 - (false.count / (length(!is.na(msi_ci[1,]))-1) * 100)
  
  return(p.cent)
  
}


# create a database of artificial populations
# randomized change points ensure trends are not linear
# this version runs the creation process separately for each species
pgrowth4.2 <- function(tpops, tmax, gr_mean, gr_sd_vec, popspec) {
  
  all_pops <- matrix(NA, nrow = tpops, ncol = tmax+1) # create matrix to put all the population counts
  
  spec_id <- sample(1:(tpops/popspec), tpops, replace = TRUE) # assign each population to a species
  
  all_pops[,tmax+1] <- t(spec_id) # fill the final column of the matrix with the species IDs
  
  for (i in unique(spec_id)) { # loop over each species
  
    num_cp <- sample(1:(tmax/5), 1) # number of times to change the growth rate
    
    popmean <- seq(0.00, 2*sample(gr_mean, 1), length.out = 1000) # sequence of mean growth rates to sample from

    cp <- vector("logical", length = tmax) # create change points vector
    
    cp[c(sample(1:tmax, num_cp))] <- "TRUE" # assign change points to randomly chosen time points
  
    for (j in 1:tmax) { # loop for all time points
    
      if (cp[j] == "TRUE") { # check if each time point is a change point
    
        popmean <- seq(0.00, 2 * sample(gr_mean, 1), length.out = 1000) # set new mean growth rate distribution

      } 
      
      popmean1 <- sample(popmean, 1) # choose a mean growth rate for this species
      
      if (is.vector(all_pops[all_pops[,tmax+1]==i,])) { # check if there is only one pop. assigned to this species
      
        temp_pops <- 10^popmean1 # if so, assign it the mean growth rate, as no distribution is needed

      } else { # but if there is more than one population assigned to this species..
      
        temp_pops <- rlnorm(nrow(all_pops[all_pops[,tmax+1]==i,]), 
                            meanlog = popmean1, sdlog = gr_sd_vec) # create a growth rate distr. for this species
      
      }
      
      all_pops[all_pops[,tmax+1]==i,j] <- temp_pops # assign growth rates to each population in this species
    
    }  
  
  }
  
  for (i in 1:tpops) {
    
    ts <- cumprod(all_pops[i,1:tmax])
    
    ts <- (ts/ts[1])*100
    
    all_pops[i,1:tmax] <- ts
  
  }
  
  colnames(all_pops) <- c(1:(ncol(all_pops)-1), "SpecID")
  
  all_pops <- as.data.frame(all_pops)
  
  all_pops <- all_pops[complete.cases(all_pops),]
  
  return(all_pops)

}

# create a database of artificial populations
# randomized change points ensure trends are not linear
# this version runs the creation process separately for each species
pgrowth4.3 <- function(tpops, tmax, gr_mean, gr_sd_vec, popspec) {
  
  all_pops <- matrix(NA, nrow = tpops, ncol = tmax+1) # create matrix to put all the population counts
  
  spec_id <- sample(1:(tpops/popspec), tpops, replace = TRUE) # assign each population to a species
  
  all_pops[,tmax+1] <- t(spec_id) # fill the final column of the matrix with the species IDs
  
  num_cp <- sample(1:(tmax/10), 1) # number of times to change the growth rate
  
  cp <- vector("logical", length = tmax) # create change points vector
  
  cp[c(sample(1:tmax, num_cp))] <- "TRUE" # assign change points to randomly chosen time points
  
  for (i in unique(spec_id)) { # loop over each species
    
    popmean <- sample(rnorm(1000, mean=gr_mean, sd=gr_sd_vec), 1) # randomly select mean growth rate from inputs
    
    for (j in 1:tmax) { # loop for all time points
      
      if (cp[j] == "TRUE") { # check if each time point is a change point
        
        popmean <- sample(rnorm(1000, mean=gr_mean, sd=gr_sd_vec), 1) # randomly select new mean growth rate from inputs
        
      } 
      
      #popmean1 <- sample(popmean, 1) # choose a mean growth rate for this species
      
      if (is.vector(all_pops[all_pops[,tmax+1]==i,])) { # check if there is only one pop. assigned to this species
        
        temp_pops <- 10^popmean # if so, assign it the mean growth rate, as no distribution is needed
        
      } else { # but if there is more than one population assigned to this species..
        
        temp_pops <- rlnorm(nrow(all_pops[all_pops[,tmax+1]==i,]), 
                            meanlog = popmean, sdlog = gr_sd_vec) # create a growth rate distr. for this species
        
        #temp_pops <- rnorm(nrow(all_pops[all_pops[,tmax+1]==i,]), 
        #                                       mean = exp(popmean), sd = gr_sd_vec)
        
      }
      
      all_pops[all_pops[,tmax+1]==i,j] <- temp_pops # assign growth rates to each population in this species
      
    }  
    
  }
  
  for (i in 1:tpops) {
    
    ts <- cumprod(all_pops[i,1:tmax])
    
    ts <- (ts/ts[1])*100
    
    all_pops[i,1:tmax] <- ts
    
  }
  
  colnames(all_pops) <- c(1:(ncol(all_pops)-1), "SpecID")
  
  all_pops <- as.data.frame(all_pops)
  
  all_pops <- all_pops[complete.cases(all_pops),]
  
  return(all_pops)
  
}


# function to randomly remove observations from synthetic dataset
remove_vals_fn <- function(all_pops_index, c, rmin=0.65, rmax=0.95) {
  
  # minimum number of observations to remove
  a <- round(c * rmin)
  
  # maximum number of observations to remove
  b <- round(c * rmax)
  
  # sequence with all possible numbers of observations that can be removed
  num <- seq(a,b)
  
  # loop for each population
  for (i in 1:nrow(all_pops_index)) {
    
    # sample from the sequence to choose how many observations to remove
    d <- sample(num, 1)
    
    # replace d randomly chosen observations with NA
    all_pops_index[i,][sample(ncol(all_pops_index[1:(ncol(all_pops_index)-3)]), d)] <- NA 
    
  }
  
  return(all_pops_index)
  
}

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


# function to introduce observation error to time series
error_intr_fn <- function(all_pops_index, m_colnames) {
  
  # create matrix to hold new index values
  all_pops_error <- matrix(data = NA, nrow = nrow(all_pops_index), ncol= length(m_colnames))
  
  mult_seq <- seq(0.03, 0.3, length.out=1000) # create values for setting standard deviation at 3% to 30% of mean
  
  for (i in 1:nrow(all_pops_index)) {
    
    for (j in 1:(length(m_colnames))) {
      
      multiplier <- sample(mult_seq, 1) # randomly select a multiplier value to set the standard deviation
      
      std_val <- all_pops_index[i,j] * multiplier # set standard deviation as mean * multiplier
      
      # create a normal distribution of index values with actual index value as mean and std_val as st. dev.
      # note that this uses the absolute value of the normal distribution to avoid negative values
      distribution <- abs(rnorm(1000, mean = all_pops_index[i,j], sd = std_val))
      
      new_ival <- sample(distribution, 1) # sample from distribution to get new index value with error
      
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

# function to introduce observation error to time series
error_intr_fn2 <- function(all_pops_index, m_colnames, mean_cv, cv_sd) {
  
  # remove zeros and convert to log
  #all_pops_index2 <- log(all_pops_index[,1:length(m_colnames)] + 0.0000001)
  
  # create matrix to hold new index values
  all_pops_error <- matrix(data = NA, nrow = nrow(all_pops_index), ncol= length(m_colnames))
  
  # find the difference between the mean and minimum
  #mean_min_diff <- mean_cv - min_cv
  
  # find the difference between the maximum and the mean
  #mean_max_diff <- max_cv - mean_cv
  
  # the standard deviation for the cv distribution is calculated as one third of
  # the difference between the max and mean or mean and min, whichever is smaller
  # because a normal distribution of 1000 values should give max and min values 
  # of roughly 3x the standard deviation (this is subject to chance)
  #cv_sd <- min(c(mean_min_diff, mean_max_diff)) / 3
  
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
      
     # new_ival <- sample(distribution, 1) # sample from distribution to get new index value with error
      
      all_pops_error[i,j] <- new_ival # put new index value into matrix
      
    }
    
  }
  
  # back convert from log
  #all_pops_error <- exp(all_pops_error)
  
  # adjust all indices to begin at a base value of 100
  all_pops_error <- all_pops_error / all_pops_error[,1] * 100
  
  all_pops_error <- as.data.frame(all_pops_error) # convert to data frame
  
  names(all_pops_error) <- m_colnames # get column names from original data frame
  
  all_pops_error$SpecID <- all_pops_index$SpecID # get species ID values from original data frame
  
  all_pops_error$PopID <- all_pops_index$PopID # get population ID values from original data frame
  
  all_pops_error$GrpID <- all_pops_index$GrpID # get group ID values from original data frame
  
  return(all_pops_error)
}

# calculate growth rates
growth_rate_calc_fn <- function(merged.matrix, model=FALSE) {
  
  # create empty vectors to store max, min, and mean growth rates
  maxr <- vector()
  meanr <- vector()
  minr <- vector()
  sd.mr <- vector()
  gratelist <- list()
  
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
    
    g.rate[g.rate <= 0.1] <- 0.1
    
    # convert to a vector
    g.rate.vec <- as.vector(as.matrix(g.rate))
    
    # find the maximum growth rate from the given row.
    maxr[i] <- max(g.rate.vec)
    
    # find the mean growth rate from the given row (geometric mean).
    meanr[i] <- exp(1) ^ mean(log(g.rate.vec))
    #meanr[i] <- mean(g.rate.vec)
    
    # find the mimimum growth rate from the given row.
    minr[i] <- min(g.rate.vec)
    
    sd.mr[i] <- sd(g.rate.vec)
    
    gratelist[[i]] <- g.rate.vec
  }
  
  allgrates <- unlist(gratelist)
  
  sd.all <- sd(allgrates)
  
  sd.mr <- sd.mr[!is.na(sd.mr) & is.finite(sd.mr)]
  
  maxr <- maxr[!is.na(maxr) & is.finite(maxr)]
  
  meanr <- meanr[!is.na(meanr) & is.finite(meanr)]
  
  minr <- minr[!is.na(minr) & is.finite(minr)]
  
  mean.sdmr <- mean(sd.mr, na.rm=TRUE)
  
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
  
  meanr.geomean <- exp(1) ^ mean(log(meanr))
  
  return.vals <- c(meanr.geomean, sd.all, meanr.sd, meanr.max, meanr.min, meanr.var, maxr.max, minr.min, meanr.mean)
  
  names(return.vals) <- c("mean.geomean.gr", "sd.all.gr", "sd.mean.gr", "max.mean.gr", "min.mean.gr", "var.mean.gr", "max.max.gr", "min.min.gr", "mean.mean.gr")
  
  return(return.vals)
  
}


# calculate species mean growth rates
growth_rate_calc_fn2 <- function(merged.matrix, c2=50, model=FALSE) {
  
  # create empty vectors to store max, min, and mean growth rates
  maxr2 <- vector()
  meanr2 <- vector()
  minr2 <- vector()
  sd.mr2 <- vector()
  gratelist2 <- list()
  popgrlist2 <- list()
  
  # loop to get max and mean growth rates
  for (i in unique(merged.matrix$SpecID)) {
    
    # create empty vectors to store max, min, and mean growth rates
    maxr <- vector()
    meanr <- vector()
    minr <- vector()
    sd.mr <- vector()
    gratelist <- list()

    tempa <- merged.matrix[merged.matrix$SpecID==i,1:c2] # get non-na values
    
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
      g.rate <- final.vals / start.vals
      
      g.rate[g.rate >= 10] <- 10
      
      g.rate[g.rate <= 0.1] <- 0.1
      
      # convert to a vector
      g.rate.vec <- as.vector(as.matrix(g.rate))
      
      # find the maximum growth rate from the given row.
      maxr[j] <- max(g.rate.vec)
      
      # find the mean growth rate from the given row (geometric mean).
      meanr[j] <- exp(1) ^ mean(log(g.rate.vec))
      #meanr[j] <- mean(g.rate.vec)
      
      # find the mimimum growth rate from the given row.
      minr[j] <- min(g.rate.vec)
      
      sd.mr[j] <- sd(g.rate.vec)
      
      gratelist[[j]] <- g.rate.vec
      
    }
    
    maxr2[i] <- max(maxr)
    
    meanr2[i] <- mean(meanr)
    
    minr2[i] <- min(minr)
    
    sd.mr2[i] <- sd(meanr)
    
    gratelist.temp <- unlist(gratelist)
    
    gratelist2[[i]] <- gratelist.temp
    
    popgrlist2[[i]] <- meanr
   
  }
  
  allgrates <- unlist(gratelist2)
  
  sd.all <- sd(allgrates)
  
  popgrates <- unlist(popgrlist2)
  
  sd.pop <- sd(popgrates)

  maxr2 <- maxr2[!is.na(maxr2) & is.finite(maxr2)]
  
  meanr2 <- meanr2[!is.na(meanr2) & is.finite(meanr2)]
  
  minr2 <- minr2[!is.na(minr2) & is.finite(minr2)]
  
  sd.mr2 <- sd.mr2[!is.na(sd.mr2) & is.finite(sd.mr2)]
  
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
  
  meanr.mean <- mean(meanr2)
  
  meanr.sd <- sd(meanr2)
  
  meanr.var <- meanr.sd / sqrt(length(meanr2))
  
  meanr.max <- max(meanr2)
  
  meanr.min <- min(meanr2)
  
  meanr.geomean <- exp(1) ^ mean(log(meanr2))
  
  return.vals <- c(meanr.geomean, sd.all, sd.pop, meanr.sd, meanr.max, meanr.min, meanr.var, maxr.max, minr.min, meanr.mean)
  
  names(return.vals) <- c("mean.geomean.gr", "sd.all.gr", "sd.pop.gr", "sd.mean.gr", "max.mean.gr", "min.mean.gr", "var.mean.gr", "max.max.gr", "min.min.gr", "mean.mean.gr")
  
  return(return.vals)
  
}

# calculate species mean growth rates
growth_rate_calc_fn3 <- function(merged.matrix, c2=50, model=FALSE) {
  
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

    tempa <- merged.matrix[merged.matrix$SpecID==i,1:c2] # get non-na values
    
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
  
#  return.vals <- c(meanr.geomean, sd.all, sd.pop, meanr.sd, meanr.max, meanr.min, meanr.var, maxr.max, minr.min, meanr.mean)
  
#  names(return.vals) <- c("mean.geomean.gr", "sd.all.gr", "sd.pop.gr", "sd.mean.gr", "max.mean.gr", "min.mean.gr", "var.mean.gr", "max.max.gr", "min.min.gr", "mean.mean.gr")
 
  return.vals <- list(armeanpopargrates, sd.arpop, meanpopsds, sdpopsds, popsds, popargrlist2)#, geomeanr2, sd.mr2, popgeogrlist2, popsdlist2)
  
  names(return.vals) <- c("Mean of Population Mean Log Growth Rates",
                          "Standard Deviation of Population Mean Log Growth Rates",
                          "Mean of Population Growth Rate Standard Deviations",
                          "Standard Deviation of Population Growth Rate Standard Deviations",
                          "Population Standard Deviations of Growth Rates",
                          "Population Mean Growth Rates by Species")#,
                          #"Species Arithmetic Means of Population Geometric Mean Growth Rates",
                          #"Species Standard Deviations of Population Geometric Mean Growth Rates",
                          #"Pop Geo Mean Growth Rates by Species",
                          #"Pop Growth Rate St Dev by Species")
   
  return(return.vals)
  
}

# create a database of artificial populations
# randomized change points ensure trends are not linear
# this version runs the creation process separately for each species
pgrowth4.4 <- function(tpops, tmax, gr_mean, gr_sd_vec, popspec) {
  
  all_pops <- matrix(NA, nrow = tpops, ncol = tmax+1) # create matrix to put all the population counts
  
  all_pops[,tmax+1] <- 1:tpops # fill the final column of the matrix with the species IDs
  
  for (i in 1:nrow(all_pops)) { # loop over each species
    
    popmean <- sample(rnorm(1000, mean=gr_mean, sd=gr_sd_vec), 1)
    
      temp_pops <- rlnorm(tmax, 
                            meanlog = popmean, sdlog = gr_sd_vec) # create a growth rate distr. for this species
        
      all_pops[i,1:tmax] <- temp_pops # assign growth rates to each population in this species
      
  }
  
  for (i in 1:tpops) {
    
    ts <- cumprod(all_pops[i,1:tmax])
    
    ts <- (ts/ts[1])*100
    
    all_pops[i,1:tmax] <- ts
    
  }
  
  colnames(all_pops) <- c(1:(ncol(all_pops)-1), "SpecID")
  
  all_pops <- as.data.frame(all_pops)
  
  all_pops <- all_pops[complete.cases(all_pops),]
  
  return(all_pops)
  
}

# create a database of artificial populations
pgrowth4.5 <- function(tpops, tmax, gr_mean, gr_sd_vec, popspec, sd_mean) {
  
  all_pops <- matrix(NA, nrow = tpops, ncol = tmax+1) # create matrix to put all the population counts
  
  spec_id <- sample(1:(tpops/popspec), tpops, replace = TRUE) # assign each population to a species
  
  all_pops[,tmax+1] <- t(spec_id) # fill the final column of the matrix with the species IDs
  
  # create normal distribution of species mean growth rates
  specmeangrdist <- rnorm(length(unique(spec_id)), mean=gr_mean, sd=gr_sd_vec)
  
  # create exponential distribution of population standard deviations of growth rates
  popsddist <- rexp(nrow(all_pops), rate=1/sd_mean)

  for (i in unique(spec_id)) { # loop over each species
    
    if (is.null(nrow(all_pops[all_pops[,tmax+1]==i,]))) { # check if there are no populations assigned to this species
      
      next
      
    }
    
    # create temporary matrix to fill with species growth rates
    spec_data <- all_pops[all_pops[,tmax+1]==i,1:tmax]
    
    # randomly select species mean growth rate from normal distribution
    specmean <- sample(specmeangrdist, 1)
    
    # randomly select standard deviation of mean growth rate from a sequence between 0 and dataset st.dev
    specsd <- sample(seq(0, gr_sd_vec, by=0.01), 1)
    
    # create normal distribution of population mean growth rates
    popmeangrdist <- rnorm(nrow(spec_data), mean=specmean, sd=specsd)
    
    for (j in 1:nrow(spec_data)) { # loop over all populations in species i

      # randomly select a mean growth rate for this population from normal distribution
      popmean <- sample(popmeangrdist, 1) 
      
      # randomly select a standard deviation of growth rates for this population from exponential dist
      popsd <- sample(popsddist, 1)

      # assign growth rates to each population in this species
      spec_data[j,] <- rlnorm(tmax, meanlog=popmean, sdlog=popsd)
      
    }
    
    # put the species growth rates into the population matrix
    all_pops[all_pops[,tmax+1]==i,1:tmax] <- spec_data
    
  }
  
  # loop over each population
  for (i in 1:tpops) {
    
    # calculate cumulative product of growth rates
    ts <- cumprod(all_pops[i,1:tmax])
    
    # convert to index, starting at value of 100. year 1 is the base year
    ts <- (ts/ts[1])*100
    
    # put the index values into the matrix
    all_pops[i,1:tmax] <- ts
    
  }
  
  # name the columns
  colnames(all_pops) <- c(1:(ncol(all_pops)-1), "SpecID")
  
  # convert to data frame
  all_pops <- as.data.frame(all_pops)
  
  # remove any NA values
  all_pops <- all_pops[complete.cases(all_pops),]
  
  return(all_pops)
  
}

