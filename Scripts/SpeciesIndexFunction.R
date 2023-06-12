######################
### Author: Shawn Dove
######################

# Function to calculate species indices

species_index_fn <- function(resamp_popmat, c) {
  
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
    
    
    # get species data
    sample_mat1.2 <- as.matrix(spec_popdata[,1:c])
    
    if (nrow(spec_popdata)==1) {
      
      # if using lambda method, convert to lambda values
      sample_mat1.3_pt1 <- sample_mat1.2[,1:(ncol(sample_mat1.2)-1)]
      sample_mat1.3_pt2 <- sample_mat1.2[,2:ncol(sample_mat1.2)]
      sample_mat1.5 <- sample_mat1.3_pt2 / (sample_mat1.3_pt1)
      
      # restructure as data frame
      sample_mat1.7 <- as.data.frame(t(sample_mat1.5))
      
      # convert time series values to natural log
      sample_mat1.8 <- log10(sample_mat1.7)
      
      # remove all change values with natural log values outside of -2.3:2.3 (log10 -1:1)
      sample_mat5 <- as.data.frame(t(apply(sample_mat1.8, 2, function(i) {ifelse(i > -1, ifelse(i < 1, i, 1), -1)})))
      
      # add years as column names
      colnames(sample_mat5) <- colnames(spec_popdata[,1:ncol(sample_mat5)])
      
      # add a species ID column
      sample_mat5$SpecID <- spec
      
      # add a group ID column
      sample_mat5$GrpID <- spec_popdata$GrpID[1]
      
    } else {
      
      # adjust all time series to begin at a base value of 100
      #sample_mat1.2 <- sample_mat / sample_mat[,1] * 100
      
      # if using lambda method, convert to lambda values
      sample_mat1.3_pt1 <- sample_mat1.2[,1:(ncol(sample_mat1.2)-1)]
      sample_mat1.3_pt2 <- sample_mat1.2[,2:ncol(sample_mat1.2)]
      sample_mat1.5 <- sample_mat1.3_pt2 / (sample_mat1.3_pt1)
      
      # convert time series values to log10
      sample_mat1.7 <- log10(sample_mat1.5)
      
      
      # remove all change values with natural log values outside of -2.3:2.3 (log10 -1:1)
      sample_mat1.8 <- apply(sample_mat1.7, 2, function(i) {ifelse(i > -1, ifelse(i < 1, i, 1), -1)})
      
      # take mean of population index values at each time point
      sample_mat2 <- colMeans(sample_mat1.8, na.rm=TRUE)
      
      # convert NaN values to NA
      sample_mat2[is.nan(sample_mat2)] <- NA
      
      # update variable name
      sample_mat4 <- sample_mat2
      
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
    
    # add all resampled copies of the species index to the list of species indices
    spec_index.list[[counter1]] <- sample_mat5
    
    # print information to show that a species has been completed
    print(paste("completed index for species ", counter1, sep=""))
    
    # increase counter value to keep track of the number of species completed
    counter1 <- counter1 + 1
    
  }
  
  return(spec_index.list)
  
}
