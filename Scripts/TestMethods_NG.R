# function to perform the methods
method_fn <- function(grp_data_culled, sample_pop_id_list, msi_real, c, m_colnames, n, n_boot, iter_num, samp_size, bootstrap_size, method) {
  
  if (method=="lambda" | method=="lambda2") {
    
    lambda<-TRUE
    
    resample<-FALSE
    
  } else if (method=="nores") {
    
    lambda<-FALSE
    
    resample<-FALSE
    
  } else if (method=="lr") {
    
    lambda<-TRUE
    
    resample<-TRUE
    
  } else if (method=="my") {
    
    lambda<-FALSE
    
    resample<-TRUE
    
  }
  
  if (method=="lambda2") {
    
    quality<-TRUE
    
  }
  
  if (method=="nores") {
    
    cat(paste0("Interpolating time series for the ", method, " method.\n"))
    
    # log-linear interpolate populations shorter than 6 years
    # populations 6 years or longer will be interpolated with a GAM in the next step
    # if method is set to nores, populations will be forecast as well
    new.grp_data <- complete_time_series(grp_data_culled, c, m_colnames, lambda=FALSE)
    saveRDS(new.grp_data, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_completed_time_series_", method, ".RData", sep=""))
    
    cat(paste0("GAMing populations for the ", method, " method.\n"))
    
    # GAM the population indices
    gam_popmat <- pop_gam_fn(new.grp_data, c, m_colnames, n=n, lambda=FALSE, resample=resample)
    
    print(paste("All population GAMs for", method, "method complete."))
    
  } else if (method=="lambda") {
    
    cat(paste0("Interpolating time series for the ", method, " method.\n"))
    
    # log-linear interpolate populations shorter than 6 years
    # populations 6 years or longer will be interpolated with a GAM in the next step
    new.grp_data <- complete_time_series(grp_data_culled, c, m_colnames, lambda=TRUE)
    saveRDS(new.grp_data, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_completed_time_series_", method, ".RData", sep=""))
    
    cat(paste0("GAMing populations for the ", method, " method.\n"))
    
    # GAM the population indices
    gam_popmat <- pop_gam_fn(new.grp_data, c, m_colnames, n=n, lambda=TRUE, resample=resample)
    
    print(paste("All population GAMs for", method, "method complete."))
    
    #cat(paste0("GAMing populations that failed quality check for the ", method, " method.\n"))
    
    # log-linear interpolate populations that failed the GAM check
    #gam_popmat <- complete_time_series(gam_popmat, c, m_colnames, lambda=TRUE)
    
  } else if (method=="lambda2") {
    
    cat(paste0("Interpolating time series for the ", method, " method.\n"))
    
    # log-linear interpolate populations shorter than 6 years
    # populations 6 years or longer will be interpolated with a GAM in the next step
    new.grp_data <- complete_time_series(grp_data_culled, c, m_colnames, lambda=TRUE)
    saveRDS(new.grp_data, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_completed_time_series_", method, ".RData", sep=""))
    
    cat(paste0("GAMing populations for the ", method, " method.\n"))
    
    # GAM the population indices
    gam_popmat <- pop_gam_fn(new.grp_data, c, m_colnames, n=n, lambda=TRUE, resample=resample, quality=TRUE)
    
    print(paste("All population GAMs for", method, "method complete."))
    
    #cat(paste0("GAMing populations that failed quality check for the ", method, " method.\n"))
    
    # log-linear interpolate populations that failed the GAM check
    gam_popmat <- complete_time_series(gam_popmat, c, m_colnames, lambda=TRUE)
    
  } else {
    
    # update variable name
    new.grp_data <- grp_data_culled
    
    cat(paste0("GAMing populations for the ", method, " method.\n"))
    
    # GAM the population indices
    # model fit checks in LPI have not been implemented as it is not automated
    # brief testing with synthetic data suggested poorly fitting models are rare
    gam_popmat <- pop_gam_fn(new.grp_data, c, m_colnames, n=n, lambda=lambda, resample=resample)
    
    print(paste("All population GAMs for", method, "method complete."))
    
  }

  cat(paste0("Creating species indices for the ", method, " method.\n"))
  
  # create species indices from the population indices
  full_spec_index <- species_index_fn(gam_popmat, c, n=n, n_boot=n_boot, lambda=lambda, resample=resample)
  saveRDS(full_spec_index, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_species_indices_", method, ".RData", sep=""))
  
  print(paste("All species indices for", method, "method complete."))
  
  remove(gam_popmat)
  
  cat(paste0("Creating msi for the ", method, " method.\n"))
  
  # create multi species indices from the group indices
  msi_full <- group_index_fn(full_spec_index, c, m_colnames, n=n, n_boot=n_boot)
  saveRDS(msi_full, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_", method, ".RData", sep=""))
  
  print(paste("Final MSI for", method, "method complete."))
  
  if (method=="lr" | method=="my") {
    
    # create final msi as mean of all msi bootstraps
    msi_final <- colMeans(msi_full, na.rm=TRUE)
    saveRDS(msi_final, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_final_", method, ".RData", sep=""))
    
  } else {
    
    msi_final <- msi_full
    
  }
  
  if (method=="lr" | method=="my") {
    
    # create final confidence intervals
    msi_ci_full <- ci_fn(full_spec_index, c, m_colnames, n=n, grouplevel=1)
    saveRDS(msi_ci_full, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_ci_", method, ".RData", sep=""))
    
    print(paste("MSI confidence intervals for", method, "method complete."))
    
  } else if (method=="lambda" | method=="lambda2") {
    
    # create final confidence intervals
    msi_ci_full <- ci_fn(full_spec_index, c, m_colnames, grouplevel=1)
    saveRDS(msi_ci_full, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_ci_", method, ".RData", sep=""))
    
    print(paste("MSI confidence intervals for", method, "method complete."))
    
  } else if (method=="nores") {
    
    # create confidence intervals for the multi species indices
    msi_ci_full <- ci_fn(full_spec_index, c, m_colnames, grouplevel=1)
    saveRDS(msi_ci_full, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_ci_", method, ".RData", sep=""))
    
    print(paste("MSI confidence intervals for", method, "method complete."))
    
  }
  
  
  # create sampled multi species indices from the population indices (species indices created but not saved)
  samp_list <- lapply(sample_pop_id_list, function(i) {
    
    if (method=="nores") {
      # interpolate and forecast time series
      pop_samp <- complete_time_series(grp_data_culled, c, m_colnames, sort(as.numeric(unlist(i))), lambda=FALSE)
      
      # GAM all populations
      gam_samp <- pop_gam_fn(pop_samp, c, m_colnames, n=n, lambda=FALSE, resample=resample)
      
    } else if (method=="lambda") {
      
      # log-linear interpolate populations with less than 6 data points
      pop_samp <- complete_time_series(grp_data_culled, c, m_colnames, sort(as.numeric(unlist(i))), lambda=TRUE)
      
      # use GAM to interpolate populations with 6 or more data points
      gam_samp <- pop_gam_fn(pop_samp, c, m_colnames, n=n, lambda=TRUE, resample=resample)
      
      # log-linear interpolate populations that failed GAM model quality check
      #gam_samp <- complete_time_series(gam_samp, c, m_colnames, lambda=TRUE)
      
    } else if (method=="lambda2") {
      
      # log-linear interpolate populations with less than 6 data points
      pop_samp <- complete_time_series(grp_data_culled, c, m_colnames, sort(as.numeric(unlist(i))), lambda=TRUE)
      
      # use GAM to interpolate populations with 6 or more data points
      gam_samp <- pop_gam_fn(pop_samp, c, m_colnames, n=n, lambda=TRUE, resample=resample, quality=TRUE)
      
      # log-linear interpolate populations that failed GAM model quality check
      gam_samp <- complete_time_series(gam_samp, c, m_colnames, lambda=TRUE)
      
    } else {
      
      # copy sample data into pop_samp but do not interpolate or forecast
      pop_samp <- grp_data_culled[grp_data_culled$PopID %in% sort(as.numeric(unlist(i))),]
      
      # use GAM to interpolate all populations, then resample from the GAM. Do not check model quality.
      gam_samp <- pop_gam_fn(pop_samp, c, m_colnames, n=n, lambda=lambda, resample=resample)
      
    }
    
    # create species indices
    spec_samp <- species_index_fn(gam_samp, c, n=n, n_boot=n_boot, lambda=lambda, resample=resample)
    
    # create final index
    msi_samp <- group_index_fn(spec_samp, c, m_colnames, n=n, n_boot=n_boot)
    
    if (method=="lr" | method=="my") {
      
      # create final msi as mean of all msi bootstraps
      msi_final_samp <- colMeans(msi_samp, na.rm=TRUE)
      
    } else {
      
      # update variable names
      msi_final_samp <- msi_samp
      
    }
    
    if (method=="lr" | method=="my") {
      
      # create msi confidence intervals
      msi_ci <- ci_fn(spec_samp, c, m_colnames, n=n, grouplevel=1)
      
    } else if (method=="lambda" | method=="lambda2") {
      
      # create msi confidence intervals
      msi_ci <- ci_fn(spec_samp, c, m_colnames, savedata=FALSE, grouplevel=1)
      
    } else if (method=="nores") {
      
      # create msi confidence intervals
      msi_ci <- ci_fn(spec_samp, c, m_colnames, grouplevel=1)
      
    }
    
    return(list(msi_samp, msi_ci, msi_final_samp))
    
    
  })
  
  # create lists to store sampled data
  msi_samp_list <- list()
  msi_ci_list <- list()
  msi_final_samp_list <- list()
  
  # loop over samples
  for (i in 1:length(samp_list)) {
    
    # add data to lists
    msi_samp_list[[i]] <- samp_list[[i]][[1]]
    msi_ci_list[[i]] <- samp_list[[i]][[2]]
    msi_final_samp_list[[i]] <- samp_list[[i]][[3]]
    
  }
  
  # save all the sampled data
  saveRDS(msi_samp_list, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_sampled_list_", samp_size, "_", method, ".RData", sep=""))
  saveRDS(msi_ci_list, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_sampled_list_ci_", samp_size, "_", method, ".RData", sep=""))
  saveRDS(msi_final_samp_list, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_final_sampled_list_", samp_size, "_", method, ".RData", sep=""))
  
  
  # test how much the trend being tested deviates from the real trend
  trend.dev <- EuclideanDistance(as.numeric(msi_final), as.numeric(msi_real))
  write.csv(trend.dev, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_trend_dev_full_", method, ".csv", sep=""))
  
  # test how much of the real trend is within the confidence intervals of the trend being tested
  within.ci <- real_trend_within_ci_fn(msi_ci_full, msi_real)
  write.csv(within.ci, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_within_ci_full_", method, ".csv", sep=""))
  
  # check the width of the confidence intervals
  ci.width <- ci_width_fn(msi_ci_full, c)
  write.csv(ci.width, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_ci_width_", method, ".csv", sep=""))
  
  # test how much the sampled trends deviate from the real trend
  trend.dev.list <- lapply(msi_final_samp_list, function(i) {if(anyNA(i)) {NA} else {EuclideanDistance(as.numeric(i), as.numeric(msi_real))}})
  trend.dev.matrix <- do.call(rbind, trend.dev.list)
  colnames(trend.dev.matrix) <- "MSI"
  write.csv(trend.dev.matrix, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_trend_dev_sampled_list_", method, ".csv", sep=""))
  
  mean.trend.dev <- mean(trend.dev.matrix, na.rm=TRUE)
  write.csv(mean.trend.dev, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_trend_dev_sampled_mean_", method, ".csv", sep=""))
  
  # test how much of the real trend is within the confidence intervals of the sampled trends
  within.ci.list <- lapply(msi_ci_list, real_trend_within_ci_fn, msi_real)
  within.ci.matrix <- do.call(rbind, within.ci.list)
  colnames(within.ci.matrix) <- "MSI"
  write.csv(within.ci.matrix, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_within_ci_sampled_list_", method, ".csv", sep=""))
  
  mean.within.ci <- mean(within.ci.matrix, na.rm=TRUE)
  write.csv(mean.within.ci, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_within_ci_sampled_mean_", method, ".csv", sep=""))
  
  # check the width of the confidence intervals of the sampled trends
  ci.width.list <- lapply(msi_ci_list, ci_width_fn, c)
  ci.width.matrix <- do.call(rbind, ci.width.list)
  colnames(ci.width.matrix) <- "MSI"
  write.csv(ci.width.matrix, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_ci_width_sampled_list_", method, ".csv", sep=""))
  
  mean.ci.width <- mean(ci.width.matrix, na.rm=TRUE)
  write.csv(mean.ci.width, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_ci_width_sampled_mean_", method, ".csv", sep=""))
  
  # plot trend
  plot_msi <- msi_final
  plot_ci <- msi_ci_full
  
  png(file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_", method, ".png", sep=""),
      width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
  
  par(mar=c(6,6,6,2) + 0.1)
  plot(m_colnames, plot_msi, xlab="", ylab="", type="l", lty=1,  ylim=c(max(c(0, min(plot_ci, na.rm=TRUE)-10)), max(plot_ci, na.rm=TRUE)+10), lwd=3, frame.plot=TRUE,
       cex.lab=1.5, cex.axis=1.5,
       panel.first = polygon(c(min(m_colnames):max(m_colnames), max(m_colnames):min(m_colnames)), c(plot_ci[2,], rev(plot_ci[1,])), col="grey90", border = NA))
  lines(m_colnames, plot_ci[2,], lty=2, lwd=2)
  lines(m_colnames, plot_ci[1,], lty=2, lwd=2)
  grid(col="grey85", lty=2)
  mtext(text = "Year", side = 1, line = 4, cex = 1.5)
  mtext(text = "Index", side = 2, line = 3.5, cex = 1.5)
  
  dev.off()
  
  # plot trend with true trend
  plot_msi <- msi_final
  plot_ci <- msi_ci_full
  
  png(file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_", method, "_with_true.png", sep=""),
      width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
  
  par(mar=c(6,6,6,2) + 0.1)
  plot(m_colnames, plot_msi, xlab="", ylab="", type="l", lty=1,  ylim=c(max(c(0, min(plot_ci, na.rm=TRUE)-10)), max(plot_ci, na.rm=TRUE)+10), lwd=3, frame.plot=TRUE,
       cex.lab=1.5, cex.axis=1.5,
       panel.first = polygon(c(min(m_colnames):max(m_colnames), max(m_colnames):min(m_colnames)), c(plot_ci[2,], rev(plot_ci[1,])), col="grey90", border = NA))
  lines(m_colnames, plot_ci[2,], lty=2, lwd=2)
  lines(m_colnames, plot_ci[1,], lty=2, lwd=2)
  lines(m_colnames, msi_real, lty=1, lwd=3, col="red")
  grid(col="grey85", lty=2)
  mtext(text = "Year", side = 1, line = 4, cex = 1.5)
  mtext(text = "Index", side = 2, line = 3.5, cex = 1.5)
  
  dev.off()
  
  #plot sampled trends
  for (i in 1:bootstrap_size) {
    
    plot_msi <- msi_final_samp_list[[i]]
    plot_ci <- msi_ci_list[[i]]
    
    dir.create(paste("TestData/", iter_num, "/Samp_Graphs_", method, "_", iter_num, sep=""))
    png(file=paste("TestData/", iter_num, "/Samp_Graphs_", method, "_", iter_num, "/saved_synth_", iter_num, "_msi_samp_", samp_size, "_", i, "_", method, ".png", sep=""),
        width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
    
    par(mar=c(6,6,6,2) + 0.1)
    plot(m_colnames, plot_msi, xlab="", ylab="", type="l", lty=1,  ylim=c(max(c(0, min(plot_ci, na.rm=TRUE)-10)), max(plot_ci, na.rm=TRUE)+10), lwd=3, frame.plot=TRUE,
         cex.lab=1.5, cex.axis=1.5,
         panel.first = polygon(c(min(m_colnames):max(m_colnames), max(m_colnames):min(m_colnames)), c(plot_ci[2,], rev(plot_ci[1,])), col="grey90", border = NA))
    lines(m_colnames, plot_ci[2,], lty=2, lwd=2)
    lines(m_colnames, plot_ci[1,], lty=2, lwd=2)
    grid(col="grey85", lty=2)
    mtext(text = "Year", side = 1, line = 4, cex = 1.5)
    mtext(text = "Index", side = 2, line = 3.5, cex = 1.5)
    
    dev.off()
    
    #plot sample trends with true trend
    png(file=paste("TestData/", iter_num, "/Samp_Graphs_", method, "_", iter_num, "/saved_synth_", iter_num, "_msi_samp_vs_full_", samp_size, "_", i, "_", method, ".png", sep=""),
        width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
    
    par(mar=c(6,6,6,2) + 0.1)
    plot(m_colnames, plot_msi, xlab="", ylab="", type="l", lty=1,  ylim=c(max(c(0, min(plot_ci, na.rm=TRUE)-10)), max(plot_ci, na.rm=TRUE)+10), lwd=3, frame.plot=TRUE,
         cex.lab=1.5, cex.axis=1.5,
         panel.first = polygon(c(min(m_colnames):max(m_colnames), max(m_colnames):min(m_colnames)), c(plot_ci[2,], rev(plot_ci[1,])), col="grey90", border = NA))
    lines(m_colnames, plot_ci[2,], lty=2, lwd=2)
    lines(m_colnames, plot_ci[1,], lty=2, lwd=2)
    lines(m_colnames, msi_real, lty=1, lwd=3, col="red")
    grid(col="grey85", lty=2)
    mtext(text = "Year", side = 1, line = 4, cex = 1.5)
    mtext(text = "Index", side = 2, line = 3.5, cex = 1.5)
    
    dev.off()
    
  }
  
  print(paste(method, "method complete"))
  
  # garbage cleaning
  gc()
  
  return(1)
  
}
