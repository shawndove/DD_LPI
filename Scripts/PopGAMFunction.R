######################
### Author: Shawn Dove
######################

# Function to model each population time series with a GAM

pop_gam_fn <- function(new.grp_data, c, m_colnames, quality=FALSE) {
  
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
    
    # predict all values between the first and last missing values using GAM
    # in this case we are interpolating using the GAM, as in the LPI
    startGAM <- min(which(!is.na(pop_data$count)))
    
    endGAM <- max(which(!is.na(pop_data$count)))
    
    # predict the missing values from the GAM
    pred.a[,startGAM:endGAM] <- t(predict(b, pop_data[startGAM:endGAM,]))

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
  
  # add extra columns back in from the original data frame
  gam_popmat[,(c+1):(length(new.grp_data)-1)] <- new.grp_data[,(c+1):(length(new.grp_data)-1)]
  
  # add column names back in
  colnames(gam_popmat) <- colnames(new.grp_data[,1:(length(new.grp_data)-1)])
    

  # if we are reproducing the LPI method...
  if (quality==TRUE) {
    
    # add the GAM quality fail status to the data frame
    gam_popmat$gqfail <- gqfail
    
  }
  
  return(gam_popmat)
  
}
