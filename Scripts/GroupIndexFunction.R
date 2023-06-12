######################
### Author: Shawn Dove
######################

# Function to calculate group indices

group_index_fn <- function(grp_index.list, c, m_colnames) {
  
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
    
  } else {
    
    lambda<-FALSE
    
  }
  

  if (nrow(grp_index)==1 & !any(!is.na(grp_index))) {
    
    # if in lambda format...
    if (!(tail(m_colnames, 1) %in% colnames(grp_index))) {
      
      # remove id tags
      grp_index1.1 <- grp_index[,1:(c-1)]
      
        # convert to index values
        grp_index1.5 <- append(grp_index1.1, 100, after = 0)
        grp_index2 <- cumprod(grp_index1.5)
        
      # if tagged to stay in lambda format
    } else {
      
      # remove id tags
      grp_index2 <- grp_index[,1:c]
      
    }
    
    print(paste("completed msi"))
    
    grp.list <- grp_index2
      
    return(grp.list)
    
  }
  
  grp_ids <- 1
  
  for (grp in grp_ids) {
    
    # otherwise, copy all data
    ms_grpdata <- grp_index
      
    if (lambda==TRUE) {
      
      # remove ID tags
      ms_grpdata1.4 <- ms_grpdata[,1:(c-1)]
      
    } else {
      
      # remove ID tags
      ms_grpdata1.1 <- ms_grpdata[,1:c]
      
      # convert to log10
      ms_grpdata1.4 <- log10(ms_grpdata1.1)
      
    }
    
    if (lambda==TRUE) {
      
      # take the mean of the group lambdas
      ms_grpdata1.5 <- colMeans(ms_grpdata1.4, na.rm=TRUE)
      
    } else {
      
      # take the mean of the group lambdas
      ms_grpdata1.5 <- colMeans(ms_grpdata1.4, na.rm=TRUE)
      
    }
    
    if (lambda==TRUE) {
      
      # back convert from natural log
      ms_grpdata1.6 <- 10^(ms_grpdata1.5)
      
      # convert to index values
      ms_grpdata1.7 <- append(ms_grpdata1.6, 100, after = 0)
      ms_grpdata1.8 <- cumprod(ms_grpdata1.7)
      
    } else {
      
      # back convert from natural log
      ms_grpdata1.6 <- 10^(ms_grpdata1.5)
      
      # ensure index starts at 100
      ms_grpdata1.8 <- (ms_grpdata1.6 / ms_grpdata1.6[1] * 100)
      
    } 
    
    # transpose group index as a matrix
    ms_grpdata2.1 <- t(as.matrix(ms_grpdata1.8))
    
    # convert NaN values to NA
    ms_grpdata2.1[is.nan(ms_grpdata2.1)] <- NA
    
    # add column names
    colnames(ms_grpdata2.1) <- m_colnames
    
    # restructure as data frame
    ms_grpdata2.1 <- as.data.frame(ms_grpdata2.1)
    
    print(paste("completed msi"))
    
    grp.list <- ms_grpdata2.1

  }
  
  return(grp.list)
  
}