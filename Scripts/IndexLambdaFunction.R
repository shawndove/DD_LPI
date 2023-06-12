######################
### Author: Shawn Dove
######################

# Function to generate and plot a trend from a dataset using the LPI's lambda method
#
# This function has 8 arguments, explained below. It is intended to be called by 
# MainFunction.R, therefore arguments should be supplied by MainFunction.R.
#
# @dataset - The dataset from which to calculate and plot a trend.
# @sample_pop_id_list - A list of time series to be included in each sample (should
#                       have been produced by MainFunction.R).
# @true_trend - The "true" trend of the dataset (should have been produced by
#               MainFunction.R).
# @iter_num - This is an identification number which will be added to saved files
#             which will be placed within a subdirectory titled with this number.
#             This should be carried over from MainFunction.R.
# @dir_name - The directory to store the data. The directory should have been
#             created by MainFunction.R. Files will be saved in a numbered subdirectory 
#             (see @iter_num) within this directory.
# @samp_size - Number of time series to sample from the dataset.
# @tmax - The number of years for which to generate data.
# @resamp_size - Number of times to sample from the dataset.


indexlambda_fn <- function(dataset,
                           sample_pop_id_list, 
                           true_trend,
                           iter_num,
                           dir_name,
                           samp_size, 
                           tmax,
                           resamp_size) {
  
  library(philentropy) # load the philentropy package containing the jaccard distance measure function
  
  m_colnames <- 1:tmax # generate a vector of years of the dataset
                       # this is used when ignoring columns without observations

  cat(paste0("Interpolating time series.\n"))
    
  # call a function to log-linear interpolate populations shorter than 6 years
  # populations 6 years or longer will be interpolated with a GAM in the next step
  new.grp_data <- interpolate_fn(dataset, tmax, m_colnames)
  #saveRDS(new.grp_data, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_completed_time_series.RData", sep=""))
    
  cat(paste0("GAMing populations.\n"))
    
  # call a function to GAM the population indices
  gam_popmat <- pop_gam_fn(new.grp_data, tmax, m_colnames)
  saveRDS(gam_popmat, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_completed_time_series.RData", sep=""))
  
    
  print(paste("All population GAMs complete."))
    
  cat(paste0("Creating species indices.\n"))
  
  # call a function to calculate species indices from the population indices
  full_spec_index <- species_index_fn(gam_popmat, tmax)
  saveRDS(full_spec_index, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_species_indices.RData", sep=""))
  
  print(paste("All species indices complete."))
  
  # free up memory, as gam_popmat is a very large object
  #remove(gam_popmat)
  
  cat(paste0("Creating final index.\n"))
  
  # call a function to calculate multi species indices from the group indices
  msi_final <- group_index_fn(full_spec_index, tmax, m_colnames)
  saveRDS(msi_final, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_msi.RData", sep=""))
  
  print(paste("Final index complete."))
  
  # call a function to create final confidence intervals
  msi_ci_final <- ci_fn(full_spec_index, tmax, m_colnames)
  saveRDS(msi_ci_final, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_msi_ci.RData", sep=""))
    
  print(paste("Confidence intervals complete."))
    
  # create sampled multi species indices from the population indices (species indices created but not saved)
  samp_list <- lapply(sample_pop_id_list, function(i) {
    
    # log-linear interpolate populations with less than 6 data points
    pop_samp <- interpolate_fn(dataset, tmax, m_colnames, sample_pop_ids=sort(as.numeric(unlist(i))))
      
    # use GAM to interpolate populations with 6 or more data points
    gam_samp <- pop_gam_fn(pop_samp, tmax, m_colnames)
      
    # calculate species indices
    spec_samp <- species_index_fn(gam_samp, tmax)
    
    # calculate final index
    msi_final_samp <- group_index_fn(spec_samp, tmax, m_colnames)
    
    # create msi confidence intervals
    msi_ci <- ci_fn(spec_samp, tmax, m_colnames)
      
    return(list(msi_ci, msi_final_samp))
    
  })
  
  # create lists to store sampled data
  msi_ci_list <- list()
  msi_final_samp_list <- list()
  
  # loop over samples
  for (i in 1:length(samp_list)) {
    
    # add data to lists
    msi_ci_list[[i]] <- samp_list[[i]][[1]]
    msi_final_samp_list[[i]] <- samp_list[[i]][[2]]
    
  }
  
  # save all the sampled data
  saveRDS(msi_ci_list, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_msi_sampled_list_ci_", samp_size, ".RData", sep=""))
  saveRDS(msi_final_samp_list, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_msi_final_sampled_list_", samp_size, ".RData", sep=""))
  
  # calculate how much the trend being tested deviates from the real trend using the jaccard distance measure
  trend.dev <- jaccard(as.numeric(msi_final), as.numeric(true_trend), testNA=FALSE)
  write.csv(trend.dev, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_trend_dev_full.csv", sep=""))
  
  # calculate how much the sampled trends deviate from the real trend
  trend.dev.list <- lapply(msi_final_samp_list, function(i) {if(anyNA(i)) {NA} else {jaccard(as.numeric(i), as.numeric(true_trend), testNA=FALSE)}})
  trend.dev.matrix <- do.call(rbind, trend.dev.list)
  colnames(trend.dev.matrix) <- "MSI"
  write.csv(trend.dev.matrix, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_trend_dev_sampled_list.csv", sep=""))
  
  mean.trend.dev <- mean(trend.dev.matrix, na.rm=TRUE)
  write.csv(mean.trend.dev, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_trend_dev_sampled_mean.csv", sep=""))
  
  # plot trend
  plot_msi <- msi_final
  plot_ci <- msi_ci_final
  
  # open a png file to save the plot to
  png(file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_msi.png", sep=""),
      width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
  
  # build the plot
  par(mar=c(6,6,6,2) + 0.1)
  plot(m_colnames, plot_msi, xlab="", ylab="", type="l", lty=1,  ylim=c(max(c(0, min(plot_ci, na.rm=TRUE)-10)), max(plot_ci, na.rm=TRUE)+10), lwd=3, frame.plot=TRUE,
       cex.lab=1.5, cex.axis=1.5,
       panel.first = polygon(c(min(m_colnames):max(m_colnames), max(m_colnames):min(m_colnames)), c(plot_ci[2,], rev(plot_ci[1,])), col="grey90", border = NA))
  lines(m_colnames, plot_ci[2,], lty=2, lwd=2)
  lines(m_colnames, plot_ci[1,], lty=2, lwd=2)
  grid(col="grey85", lty=2)
  mtext(text = "Year", side = 1, line = 4, cex = 1.5)
  mtext(text = "Index", side = 2, line = 3.5, cex = 1.5)
  
  dev.off() # close the png file
  
  # plot trend with true trend
  plot_msi <- msi_final
  plot_ci <- msi_ci_final
  
  # open png file
  png(file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_msi_with_true.png", sep=""),
      width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
  
  # build plot
  par(mar=c(6,6,6,2) + 0.1)
  plot(m_colnames, plot_msi, xlab="", ylab="", type="l", lty=1,  ylim=c(max(c(0, min(plot_ci, na.rm=TRUE)-10)), max(plot_ci, na.rm=TRUE)+10), lwd=3, frame.plot=TRUE,
       cex.lab=1.5, cex.axis=1.5,
       panel.first = polygon(c(min(m_colnames):max(m_colnames), max(m_colnames):min(m_colnames)), c(plot_ci[2,], rev(plot_ci[1,])), col="grey90", border = NA))
  lines(m_colnames, plot_ci[2,], lty=2, lwd=2)
  lines(m_colnames, plot_ci[1,], lty=2, lwd=2)
  lines(m_colnames, true_trend, lty=1, lwd=3, col="red")
  grid(col="grey85", lty=2)
  mtext(text = "Year", side = 1, line = 4, cex = 1.5)
  mtext(text = "Index", side = 2, line = 3.5, cex = 1.5)
  
  dev.off() # close png file
  
  #plot sampled trends
  for (i in 1:resamp_size) {
    
    plot_msi <- msi_final_samp_list[[i]]
    plot_ci <- msi_ci_list[[i]]
    
    # create a directory to store the file (if it does not exist)
    dir.create(paste(dir_name, "/", iter_num, "/Samp_Graphs_", iter_num, sep=""))
    
    # open a png file to save the plot to
    png(file=paste(dir_name, "/", iter_num, "/Samp_Graphs_", iter_num, "/saved_synth_", iter_num, "_msi_samp_", samp_size, "_", i, ".png", sep=""),
        width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
    
    # build the plot
    par(mar=c(6,6,6,2) + 0.1)
    plot(m_colnames, plot_msi, xlab="", ylab="", type="l", lty=1,  ylim=c(max(c(0, min(plot_ci, na.rm=TRUE)-10)), max(plot_ci, na.rm=TRUE)+10), lwd=3, frame.plot=TRUE,
         cex.lab=1.5, cex.axis=1.5,
         panel.first = polygon(c(min(m_colnames):max(m_colnames), max(m_colnames):min(m_colnames)), c(plot_ci[2,], rev(plot_ci[1,])), col="grey90", border = NA))
    lines(m_colnames, plot_ci[2,], lty=2, lwd=2)
    lines(m_colnames, plot_ci[1,], lty=2, lwd=2)
    grid(col="grey85", lty=2)
    mtext(text = "Year", side = 1, line = 4, cex = 1.5)
    mtext(text = "Index", side = 2, line = 3.5, cex = 1.5)
    
    dev.off() # close the png file
    
    
    # plot sample trends with true trend
    
    # open png file
    png(file=paste(dir_name, "/", iter_num, "/Samp_Graphs_", iter_num, "/saved_synth_", iter_num, "_msi_samp_vs_full_", samp_size, "_", i, ".png", sep=""),
        width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
    
    # build plot
    par(mar=c(6,6,6,2) + 0.1)
    plot(m_colnames, plot_msi, xlab="", ylab="", type="l", lty=1,  ylim=c(max(c(0, min(plot_ci, na.rm=TRUE)-10)), max(plot_ci, na.rm=TRUE)+10), lwd=3, frame.plot=TRUE,
         cex.lab=1.5, cex.axis=1.5,
         panel.first = polygon(c(min(m_colnames):max(m_colnames), max(m_colnames):min(m_colnames)), c(plot_ci[2,], rev(plot_ci[1,])), col="grey90", border = NA))
    lines(m_colnames, plot_ci[2,], lty=2, lwd=2)
    lines(m_colnames, plot_ci[1,], lty=2, lwd=2)
    lines(m_colnames, true_trend, lty=1, lwd=3, col="red")
    grid(col="grey85", lty=2)
    mtext(text = "Year", side = 1, line = 4, cex = 1.5)
    mtext(text = "Index", side = 2, line = 3.5, cex = 1.5)
    
    dev.off() # close png file
    
  }
  
  # garbage cleaning
  gc()
  
  return(1)
  
}
