######################
### Author: Shawn Dove
######################

# Main function to create simulated datasets for modelling reliability of the LPI
#
# This function has 23 arguments, explained below. The first 2 have no defaults.
#
# @iter_num - This is an identification number which will be added to saved files
#             which will be placed within a subdirectory titled with this number.
#             There is no default.
# @dir_name - The directory to store the data. The directory will be created if
#             it does not exist. Files will be saved in a numbered subdirectory 
#             (see @iter_num) within this directory. There is no default.
# @popvar - The standard deviation of the population mean growth rates.
#           (Default is 0.2)
# @popmean - The mean of the population mean growth rates (default is 0).
# @sdmean - The mean of the population standard deviations in growth rates.
#           (Default is 0.2)
# @pgrowthx - This selects between different time series generators.
#             7 uses a normal distribution of the numbers of populations per species,
#             while 8 uses a discretized exponential distribution and 9 uses a
#             negative binomial distribution (default is 7).
# @tmax - The number of years for which to generate data (default is 50).
# @tpops - The total number of time series to generate (default is 1,000).
# @popspec - The mean number of populations per species (default is 10).
# @count_thres - Remove any time series with fewer observations than this (default is 2).
# @min_ts_length - Remove any time series shorter than this (default is 2).
# @samp_size - Number of time series to sample from the dataset (default is 200).
# @mlength - Mean time series length (default is 20).
# @numobs - Mean number of observations per time series (default is 10).
# @resamp_size - Number of times to sample from the dataset (default is 20).
# @error - Whether or not to add observation error (default is TRUE).
# @degrade - Method of degrading time series. This changes the way observations
#            are distributed across time series. There are 5 options.
#            "normal" randomly removes observations randomly
#            "endpoints" forces all time series to either start at the first year 
#                        or end at the final year of the dataset
#            "cluster" divides time series into a specified number of clusters,
#                      with time series in each cluster having the same start and 
#                      end years (but randomized for each cluster)
#            "clustend" removes all but the final x observations (where x is set 
#                       by the user) for a certain portion of the sampled time 
#                       series (the portion is also set by the user), while the
#                       rest are degraded using the "normal" method.
#            "endreveal" avoids removing the final x observations (where x is set
#                        by the user) of a certain portion of the sampled time
#                        series (the portion is also set by the user), while the
#                        rest are degraded using the "normal" method.
#            (default is "normal")
# @mean_cv - Mean coefficient of variation for observation error. This setting is
#            ignored if error is set to FALSE (default is 0.15).
# @cv_sd - Standard deviation of coefficient of variation for observation error.
#          This setting is ignored if error set to FALSE (default is 0.1).
# @clustlength - The number of final observations to leave intact for "clustend"
#                degradation method. This setting is ignored for all other methods.
#                (Default is 10).
# @endlength - The number of final observations to leave intact for "endreveal"
#             degradation method. This setting is ignored for all other methods.
#             (Default is 1).
# @endpops_ratio - The portion of the sampled time series to degrade using the
#                  "clustend" method. The rest are degraded using the "normal" 
#                  method (default is 0.2). This setting is ignored for all except
#                  the "clustend" method.
# @endreveal_ratio - The portion of the sampled time series for which the final x
#                    observations will be left intact during degradation. The rest
#                    are degraded using the "normal" method (default is 0.2). This
#                    setting is ignored for all except the "endreveal" method.



main_fn <- function(iter_num, 
                    dir_name,
                    popvar=0.2, 
                    popmean=0, 
                    sdmean=0.2,
                    pgrowthx=7, 
                    tmax=50, 
                    tpops=1000, 
                    popspec=10, 
                    count_thres=2, 
                    min_ts_length=2, 
                    samp_size=200, 
                    mlength=20, 
                    numobs=10, 
                    resamp_size=20, 
                    error=TRUE,
                    degrade="normal",
                    mean_cv=0.15,
                    cv_sd=0.1,
                    clustlength=10,
                    endlength=1,
                    endpops_ratio=0.2,
                    endreveal_ratio=0.2) {
  
  # create directory to store files
  if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))}
  if(!dir.exists(paste(dir_name, "/", iter_num, sep=""))) {dir.create(paste(dir_name, "/", iter_num, sep=""))}
  
  m_colnames <- 1:tmax # generate a vector of years of the dataset
                       # this is used when ignoring columns without observations
  
  # create synthetic populations, assigned to different species
  if (pgrowthx == 7) {
    
    # create synthetic populations, assigned to different species
    all_pops_index <- pgrowth4.5_fn(tpops, tmax, popmean, popvar, popspec, sdmean)
    
    
  } else if (pgrowthx == 8) {
    
    # create synthetic populations, assigned to different species
    all_pops_index <- pgrowth4.7_fn(tpops, tmax, popmean, popvar, popspec, sdmean)
    
  } else if (pgrowthx == 9) {
    
    # create synthetic populations, assigned to different species
    all_pops_index <- pgrowth4.8_fn(tpops, tmax, popmean, popvar, popspec, sdmean)
    
  } else("Check your synthetic data generator function input setting.")
  
  # save raw synthetic dataset
  saveRDS(all_pops_index, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_raw.RData", sep=""))
  
  cat(paste("Constructing dataset:", "\nPopulations:", tpops, "\nSpecies:", length(unique(all_pops_index$SpecID)), "\n\n"))
  
  cat(paste0("Plotting geometric mean of the dataset.\n"))
  
  # create geometric mean and plot it
  index_geomean <- exp(colMeans(log(all_pops_index[,1:tmax])))
  
  png(file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_geometric_mean_raw.png", sep=""),
      width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
  
  par(mar=c(6,6,6,2) + 0.1)
  plot(m_colnames, index_geomean[1:tmax], type="l", lty=1, lwd=2, xlab="", ylab="", cex.axis=1.5)
  grid(col="grey85", lty=2)
  mtext(text = "Year", side = 1, line = 4, cex = 2)
  mtext(text = "Index", side = 2, line = 3.5, cex = 2)
  
  dev.off()
  
  cat(paste0("Plotting a GAM of the geometric mean.\n"))
  
  # GAM the geometric mean and plot it
  g <- gam(index_geomean[1:tmax]~s(m_colnames,k = round(tmax/2)))$fitted.values
  g.gam <- g / g[1] * 100
  y.low <- min(g.gam) - 5
  y.high <- max(g.gam) + 5
  
  png(file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_raw_gam_geometric_mean.png", sep=""),
      width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
  
  plot(m_colnames, g.gam, xlab="", ylab="", type="l", lty=1, ylim=c(y.low,y.high), lwd=2, frame.plot=TRUE,
       cex.axis=1.5)
  grid(col="grey85", lty=2)
  mtext(text = "Year", side = 1, line = 4, cex = 2)
  mtext(text = "Index", side = 2, line = 3.5, cex = 2)
  
  dev.off()
  
  #all_pops_index <- readRDS(file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_raw.RData", sep=""))
  
  # add a population ID column
  all_pops_index$PopID <- rownames(all_pops_index)
  
  if (error==TRUE) {
    
    cat(paste("Adding random sampling error to dataset.\n"))
    
    # add sampling error to dataset
    all_pops_error <- obs_error_fn(all_pops_index, m_colnames, mean_cv, cv_sd)
    saveRDS(all_pops_error, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_error.RData", sep=""))
    
    cat(paste0("Degrading dataset.\n"))
    
    if (degrade=="normal") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn(all_pops_error, tmax, mlength, numobs)
      
    } else if (degrade=="endpoints") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_endpoints(all_pops_error, tmax, mlength, numobs)
      
    } else if (degrade=="cluster") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_clust(all_pops_error, tmax, mlength, numobs)
      
    } else if (degrade=="clustend") {
      
      if (samp_size >= ceiling((1-endpops_ratio)*nrow(all_pops_error))) {
        
        stop(paste("Sample size must be less than", 100 * (1 - endpops_ratio), "% of total populations to use clustend method.", sep=""))
        
      }
      
      # select populations to make endpoints
      endpops_ids <- sample(all_pops_error$PopID, ceiling(endpops_ratio*nrow(all_pops_index)))
      endpops <- all_pops_error[all_pops_error$PopID %in% endpops_ids,]
      
      # remove endpops from pops index
      all_pops_error2 <- all_pops_error[!(all_pops_error$PopID %in% endpops_ids),]
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn(all_pops_error2, tmax, mlength, numobs)
      
      # degraded endpops
      endpops_degraded <- degrade_ts_fn_clustend2(endpops, tmax, clustlength)
      saveRDS(endpops_degraded, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_endpops_degraded.RData", sep=""))
      
      
    } else if (degrade=="endreveal") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_endreveal(all_pops_error, tmax, mlength, numobs, endlength, endreveal_ratio)
      
    }
    
    saveRDS(all_pops_degraded, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_degraded.RData", sep=""))
    
  } else {
    
    cat(paste0("Degrading dataset.\n"))
    
    if (degrade=="normal") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn(all_pops_index, tmax, mlength, numobs)
      
    } else if (degrade=="endpoints") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_endpoints(all_pops_index, tmax, mlength, numobs)
      
    } else if (degrade=="cluster") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_clust(all_pops_index, tmax, mlength, numobs)
      
    } else if (degrade=="clustend") {
      
      if (samp_size >= ceiling((1-endpops_ratio)*nrow(all_pops_index))) {
        
        stop(paste("Sample size must be less than", 100 * (1 - endpops_ratio), "% of total populations to use clustend method.", sep=""))
        
      }
      
      # select populations to make endpoints
      endpops_ids <- sample(all_pops_index$PopID, ceiling(endpops_ratio*nrow(all_pops_index)))
      endpops <- all_pops_index[all_pops_index$PopID %in% endpops_ids,]
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn(all_pops_index, tmax, mlength, numobs)
      #all_pops_degraded <- readRDS(file=paste(dir_name, "/", iter_num-40, "/saved_synth_", iter_num-40, "_degraded.RData", sep=""))
      
      # remove endpops from pops index
      all_pops_degraded <- all_pops_degraded[!(all_pops_degraded$PopID %in% endpops_ids),]
      
      # degraded endpops
      endpops_degraded <- degrade_ts_fn_clustend2(endpops, tmax, clustlength)
      saveRDS(endpops_degraded, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_endpops_degraded.RData", sep=""))
      
    } else if (degrade=="endreveal") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_endreveal(all_pops_index, tmax, mlength, numobs, endlength, endreveal_ratio)
      
    }
    
    saveRDS(all_pops_degraded, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_degraded.RData", sep=""))
    
  }
  
  cat(paste0("Removing time series with less than ", count_thres, " total observations and/or with a length of less than ", min_ts_length, ".\n"))
  
  # remove time series which are too short or have too few counts
  grp_data_culled <- cull_fn(all_pops_degraded, count_thres, min_ts_length, tmax)
  saveRDS(grp_data_culled, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_culled.RData", sep=""))
  
  cat(paste0("Creating ", resamp_size, " randomly sampled subsets of ", samp_size, " populations each.\n"))
  
  # if sample size is too large, set it to the total number of populations - 1
  if (nrow(grp_data_culled) < samp_size) {
    samp_size <- (nrow(grp_data_culled)-1)
  }
  
  # create list of x population ID matrices at a given sample size, where x is the number of times to sample
  sample_pop_id_list <- list()
  for (i in 1:resamp_size) {
    sample_pop_id_list[[i]] <- sample(grp_data_culled$PopID, samp_size)
  }
  #sample_pop_id_list <- readRDS(file=paste(dir_name, "/", iter_num-40, "/saved_synth_", iter_num-40, "_sample_pop_id_list.RData", sep=""))
  saveRDS(sample_pop_id_list, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_sample_pop_id_list.RData", sep=""))
  
  # if using clustend method, add selected endpops to all the samples and merge them back to the time series dataframe
  if (degrade=="clustend") {
    for (i in 1:resamp_size) {
      sample_pop_id_list[[i]] <- c(sample_pop_id_list[[i]], endpops_ids)
    }
    saveRDS(sample_pop_id_list, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_sample_pop_id_list.RData", sep=""))
    grp_data_culled <- rbind(grp_data_culled, endpops_degraded)
    saveRDS(grp_data_culled, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_culled.RData", sep=""))
  }
  
  
  ## TRUE TREND ##
  
  cat(paste0("Creating ", length(unique(all_pops_index$SpecID)), " species indices for true trend.\n"))
  
  # create species indices from the population indices
  full_spec_true <- species_index_fn(all_pops_index, tmax)
  saveRDS(full_spec_true, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_species_indices_TrueTrend.RData", sep=""))
  
  cat(paste0("Creating multi-species index for true trend.\n"))
  
  # create multi species indices from the group indices
  msi_true <- group_index_fn(full_spec_true, tmax, m_colnames)
  saveRDS(msi_true, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_msi_TrueTrend.RData", sep=""))
  
  cat(paste0("Creating msi confidence intervals for true trend.\n"))
  
  # create confidence intervals for the multi species index
  msi_ci_true <- ci_fn(full_spec_true, tmax, m_colnames)
  saveRDS(msi_ci_true, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_msi_ci_TrueTrend.RData", sep=""))
  
  temp <- indexlambda_fn(dataset=grp_data_culled, 
                         sample_pop_id_list=sample_pop_id_list, 
                         true_trend=msi_true,
                         iter_num=iter_num,
                         dir_name=dir_name,
                         samp_size=samp_size, 
                         tmax=tmax,
                         resamp_size=resamp_size)
  
  ## DATA GATHERING ##
  
  # get mean and standard deviation of the mean growth rate
  gr.stats.raw <- growth_rate_calc_fn(all_pops_index, tmax)
  
  all_pops_completed <- readRDS(file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_completed_time_series.RData", sep=""))

  gr.stats.degraded <- growth_rate_calc_fn(all_pops_completed, tmax, model=TRUE)
  
  # create csv file with important info about the dataset
  info.dat <- data.frame(row.names = "1")
  info.dat$ID <- iter_num
  info.dat$synth_version <- pgrowthx
  info.dat$num_pops <- tpops
  info.dat$num_years <- tmax
  info.dat$min_ts_length <- min_ts_length + 1
  info.dat$mean_gr_raw <- gr.stats.raw[[1]]
  info.dat$gr_sd_raw <- gr.stats.raw[[2]]
  info.dat$mean_sd_raw <- gr.stats.raw[[3]]
  info.dat$mean_gr_degraded <- gr.stats.degraded[[1]]
  info.dat$gr_sd_degraded <- gr.stats.degraded[[2]]
  info.dat$mean_sd_degraded <- gr.stats.degraded[[3]]
  info.dat$samp_size <- samp_size
  info.dat$pops_per_species <- popspec
  info.dat$mean_ts_length <- sum(!is.na(as.vector(all_pops_completed[,1:tmax]))) / nrow(all_pops_completed)
  info.dat$mean_num_obs <- sum(!is.na(as.vector(all_pops_degraded[,1:tmax]))) / nrow(all_pops_degraded)
  info.dat$degrade_type <- degrade
  info.dat$error <- error
  info.dat$mean_cv <- mean_cv
  info.dat$cv_sd <- cv_sd
  
  write.csv(info.dat, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_info.csv", sep=""))
  
  gr.stats.samples <- list()
  tslength.samples <- vector()
  numobs.samples <- vector()
  for (i in 1:length(sample_pop_id_list)) {
    
    temp <- all_pops_degraded[all_pops_degraded$PopID %in% sample_pop_id_list[[i]],]
    temp2 <- interpolate_fn(temp, tmax, m_colnames, calcsd=TRUE)
    gr.stats.samples[[i]] <- growth_rate_calc_fn(temp2, tmax, model=TRUE)
    tslength.samples[i] <- sum(!is.na(as.vector(temp2[,1:tmax]))) / nrow(temp2)
    numobs.samples[i] <- sum(!is.na(as.vector(temp[,1:tmax]))) / nrow(temp)
    
  }
  
  info.samples.meangr <- list()
  info.samples.grsd <- list()
  info.samples.meansd <- list()
  info.samples.meantslength <- vector()
  info.samples.meannumobs <- vector()
  for (i in 1:length(sample_pop_id_list)) {
    
    info.samples.meangr[[i]] <- gr.stats.samples[[i]][[1]]
    info.samples.grsd[[i]] <- gr.stats.samples[[i]][[2]]
    info.samples.meansd[[i]] <- gr.stats.samples[[i]][[3]]
    info.samples.meantslength[i] <- tslength.samples[i]
    info.samples.meannumobs[i] <- numobs.samples[i]
    
  }
  
  info.samples <- data.frame(mean_gr_sample = do.call(rbind, info.samples.meangr),
                             gr_sd_sample = do.call(rbind, info.samples.grsd),
                             mean_sd_sample = do.call(rbind, info.samples.meansd),
                             mean_ts_length = info.samples.meantslength,
                             mean_num_obs = info.samples.meannumobs)
  
  write.csv(info.samples, file=paste(dir_name, "/", iter_num, "/saved_synth_", iter_num, "_gr_samples.csv", sep=""))
  
}
