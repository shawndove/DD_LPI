### admin ###
#
# Shawn Dove
# May 17, 2021
#
#

## load packages ----

library(plyr)
library(ggplot2)
library(dplyr)
library(mgcv)
library(GET)
library(MASS)
library(reshape2)
library(matrixStats)
library(foreach)
library(doSNOW)

## load external functions ----

source("Functions.R")
source("Methods.R")

## Main Function ----

all_fn <- function(popvar, popmean, pgrowthx, iter_num, tmax, tpops, popspec, n, n_boot, ngrps,
                   count_thres, min_ts_length, c, samp_size, m_colnames, rmin, rmax, bootstrap_size, error=FALSE) {
  
  # create synthetic populations, assigned to different species
  if (pgrowthx == 5) {
    
    # create data frame to hold the populations
    all_pops_index <- as.data.frame(matrix(data=NA, nrow=tpops, ncol=c+2))
    
    # loop for each group
    for (i in 1:ngrps) {
      
      # create group of populations
      temp <- pgrowth4.2(tpops/ngrps, tmax, popmean, popvar, popspec)
      
      # put them into the data frame
      all_pops_index[(((tpops/ngrps)*i)-(tpops/ngrps)+1):((tpops/ngrps)*i),1:(c+1)] <- temp
      
      # add a group ID
      all_pops_index[(((tpops/ngrps)*i)-(tpops/ngrps)+1):((tpops/ngrps)*i),(c+2)] <- i
      
    }
    
    # force back to data frame
    all_pops_index <- as.data.frame(all_pops_index)
    
    # add column names
    colnames(all_pops_index) <- c(1:c,"SpecID","GrpID")
    
    # get approx number of species per group
    num_spec <- (tpops / (ngrps * popspec))
    
    # make species IDs unique
    all_pops_index$SpecID <- all_pops_index$SpecID + (all_pops_index$GrpID*num_spec) - num_spec
    
  } else("Check your synthetic data generator function input setting.")
  
  # save raw synthetic dataset
  saveRDS(all_pops_index, file=paste("saved_synth_", iter_num, "_raw.RData", sep=""))
  
  cat(paste("Constructing dataset:", "\nPopulations:", tpops, "\nSpecies:", length(unique(all_pops_index$SpecID)), "\nGroups:", ngrps, "\n\n"))
  
  cat(paste0("Plotting geometric mean of the dataset.\n"))
  
  # create geometric mean and plot it
  index_geomean <- exp(colMeans(log(all_pops_index[,1:c])))
  
  png(file=paste("saved_synth_", iter_num, "_geometric_mean_raw.png", sep=""),
      width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
  
  par(mar=c(6,6,6,2) + 0.1)
  plot(m_colnames, index_geomean[1:c], type="l", lty=1, lwd=2, xlab="", ylab="", cex.axis=1.5)
  grid(col="grey85", lty=2)
  mtext(text = "Year", side = 1, line = 4, cex = 2)
  mtext(text = "Index", side = 2, line = 3.5, cex = 2)
  
  dev.off()
  
  cat(paste0("Plotting a GAM of the geometric mean.\n"))
  
  # GAM the geometric mean and plot it
  g <- gam(index_geomean[1:c]~s(m_colnames,k = round(c/2)))$fitted.values
  g.gam <- g / g[1] * 100
  y.low <- min(g.gam) - 5
  y.high <- max(g.gam) + 5 
  
  png(file=paste("saved_synth_", iter_num, "_raw_gam_geometric_mean.png", sep=""),
      width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
  
  plot(m_colnames, g.gam, xlab="", ylab="", type="l", lty=1, ylim=c(y.low,y.high), lwd=2, frame.plot=TRUE,
       cex.axis=1.5)
  grid(col="grey85", lty=2)
  mtext(text = "Year", side = 1, line = 4, cex = 2)
  mtext(text = "Index", side = 2, line = 3.5, cex = 2)
  
  dev.off()
  
  # add a population ID column
  all_pops_index$PopID <- rownames(all_pops_index)
  
  if (error==TRUE) {
    
    cat(paste("Adding random sampling error to dataset.\n"))
    
    # add sampling error to dataset
    all_pops_error <- error_intr_fn(all_pops_index, m_colnames)
    saveRDS(all_pops_error, file=paste("saved_synth_", iter_num, "_error.RData", sep=""))
    
    cat(paste0("Degrading dataset by ", mean(c(rmin, rmax))*100, "%.\n"))
    
    # degrade dataset
    all_pops_degraded <- remove_vals_fn(all_pops_error, c, rmin, rmax)
    saveRDS(all_pops_degraded, file=paste("saved_synth_", iter_num, "_degraded.RData", sep=""))
    
  } else {
    
    cat(paste0("Degrading dataset by ", mean(c(rmin, rmax))*100, "%.\n"))
    
    # degrade dataset
    all_pops_degraded <- remove_vals_fn(all_pops_index, c, rmin, rmax)
    saveRDS(all_pops_degraded, file=paste("saved_synth_", iter_num, "_degraded.RData", sep=""))
    
  }
  
  cat(paste0("Removing time series with less than ", count_thres, " total observations and/or with a length of less than ", min_ts_length, ".\n"))
  
  # remove time series which are too short or have too few counts
  grp_data_culled <- cull_fn(all_pops_degraded, count_thres, min_ts_length, c)
  saveRDS(grp_data_culled, file=paste("saved_synth_", iter_num, "_culled.RData", sep=""))
  
  cat(paste0("Creating ", bootstrap_size, " randomly sampled subsets of ", samp_size, " populations each.\n"))
  
  # create list of x population ID matrices at a given sample size, where x is the number of sample bootstraps
  sample_pop_id_list <- list()
  for (i in 1:bootstrap_size) {
    sample_pop_id_list[[i]] <- sample(grp_data_culled$PopID, samp_size)
  }
  saveRDS(sample_pop_id_list, file=paste("saved_synth_", iter_num, "_sample_pop_id_list.RData", sep=""))
  
  
  ## TRUE TREND ##
  
  cat(paste0("Creating ", length(unique(all_pops_index$SpecID)), " species indices for true trend.\n"))
  
  # create species indices from the population indices
  full_spec_real <- species_index_fn(all_pops_index, c)
  saveRDS(full_spec_real, file=paste("saved_synth_", iter_num, "_species_indices_TrueTrend.RData", sep=""))
  
  cat(paste0("Creating ", ngrps, " group indices for true trend.\n"))
  
  # create group indices from the species indices
  grp_real <- group_index_fn(full_spec_real, c, m_colnames)
  saveRDS(grp_real, file=paste("saved_synth_", iter_num, "_group_indices_TrueTrend.RData", sep=""))
  
  cat(paste0("Creating ", ngrps, " group confidence intervals for true trend.\n"))
  
  # create confidence intervals for the group indices
  grp_ci_real <- ci_fn(full_spec_real, c, m_colnames)
  saveRDS(grp_ci_real, file=paste("saved_synth_", iter_num, "_grp_ci_TrueTrend.RData", sep=""))
  
  cat(paste0("Creating msi for true trend.\n"))
  
  # create multi species indices from the group indices
  msi_real <- group_index_fn(grp_real, c, m_colnames)
  saveRDS(msi_real, file=paste("saved_synth_", iter_num, "_msi_TrueTrend.RData", sep=""))
  
  cat(paste0("Creating msi confidence intervals for true trend.\n"))
  
  # create confidence intervals for the multi species indices
  msi_ci_real <- ci_fn(grp_real, c, m_colnames)
  saveRDS(msi_ci_real, file=paste("saved_synth_", iter_num, "_msi_ci_TrueTrend.RData", sep=""))
  
  cat(paste0("Plotting true trend.\n"))
  
  # plot trend
  plot_msi <- msi_real
  plot_ci <- msi_ci_real
  
  png(file=paste("saved_synth_", iter_num, "_msi_TrueTrend.png", sep=""),
      width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
  
  par(mar=c(6,6,6,2) + 0.1)
  plot(m_colnames, plot_msi, xlab="", ylab="", type="l", lty=1, ylim=c(0,max(c(200, max(plot_ci)+10))), lwd=3, frame.plot=TRUE,
       cex.lab=1.5, cex.axis=1.5,
       panel.first = polygon(c(min(m_colnames):max(m_colnames), max(m_colnames):min(m_colnames)), c(plot_ci[2,], rev(plot_ci[1,])), col="grey90", border = NA))
  lines(m_colnames, plot_ci[2,], lty=2, lwd=2)
  lines(m_colnames, plot_ci[1,], lty=2, lwd=2)
  lines(m_colnames, g.gam, lty=1, lwd=3, col="red")
  grid(col="grey85", lty=2)
  mtext(text = "Year", side = 1, line = 4, cex = 1.5)
  mtext(text = "Index", side = 2, line = 3.5, cex = 1.5)
  
  dev.off()
  
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, grp_real, c, m_colnames, n=NA, n_boot=NA, iter_num, samp_size, bootstrap_size, method="nores")
  
  temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, grp_real, c, m_colnames, n=NA, n_boot=NA, iter_num, samp_size, bootstrap_size, method="lambda")
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, grp_real, c, m_colnames, n=n, n_boot=n_boot, iter_num, samp_size, bootstrap_size, method="lr")
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, grp_real, c, m_colnames, n=n, n_boot=n_boot, iter_num, samp_size, bootstrap_size, method="my")
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, grp_real, c, m_colnames, n=NA, n_boot=NA, iter_num, samp_size, bootstrap_size, method="lambda2")
  
  
  ## DATA GATHERING ##
  
  # get mean and standard deviation of the mean growth rate
  gr.stats.raw <- growth_rate_calc_fn(all_pops_index[,1:c])
  
  #gr.stats.rebuilt <- growth_rate_calc_fn(new.grp_data[,1:c])
  
  # create csv file with important info about the dataset
  info.dat <- data.frame(row.names = "1")
  info.dat$ID <- iter_num
  info.dat$synth_version <- pgrowthx
  info.dat$num_pops <- tpops
  info.dat$num_years <- tmax
  info.dat$num_resamps <- n
  info.dat$num_bootstraps <- n_boot
  info.dat$min_ts_length <- min_ts_length + 1
  info.dat$mean_gr_raw <- gr.stats.raw[1]
  info.dat$gr_sd_raw <- gr.stats.raw[2]
  
  #info.dat$mean_gr_rebuilt <- gr.stats.rebuilt[1]
  #info.dat$gr_sd_rebuilt <- gr.stats.rebuilt[2]
  #info.dat$pops_per_species <- nrow(new.grp_data) / length(spec_ids) # average populations per species in culled data
  #info.dat$pops_per_species_raw <- nrow(all_pops_index) / length(unique(SpecID_vec)) # avg pops per species in unculled data
  
  write.csv(info.dat, file=paste("saved_synth_", iter_num, "_info.csv", sep=""))
  
}




## set parameters ----

n <- 100 # number of GAM resamples per population
n_boot <- 3000 # number of index bootstraps per species
bootstrap_size <- 20 # number of samples
samp_size <- 500 # number of populations to select for each sample
count_thres <- 3 # threshold at which number of counts is too low and we should not include this population
min_ts_length <- 3 # minimum length of time series to be included
startyear <- 1 # column name for first data point
endyear <- 50 # column name for final data point
tmax <- length(startyear:endyear) # number of years (data points per population)
tpops <- 10000 # total number of populations
popspec <- 20 # mean number of populations per species
ngrps <- 10 # number of groups to divide species into
m_colnames <- startyear:endyear # column names
c <- length(m_colnames) # number of years or columns (same as tmax)

## Setup Testing ----

iter_num <- 1200
gr_mean_a <- c(-0.2, -0.1, -0.05, 0.05, 0.1, 0.2)
gr_sd_vec_a <- c(0.1, 0.15, 0.2, 0.25, 0.3)
popspec <- 10
rmin <- 0.8
rmax <- 0.95
samp_size <- c(100, 200, 500, 1000, 2000)
tmax <- 50
c <- tmax
tpops <- 10000
m_colnames <- 1:c
j_choice_a <- rep(gr_sd_vec_a, 5)
k_choice_a <- rep(samp_size, each=5)
setwd("TestData/") # set the working directory

no_cores <- 5 # the number of cores to be used for parallel processing
cl <- makeCluster(no_cores, outfile="output.txt") # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk), library(dplyr), library(MASS), library(GET), library(mgcv), library(reshape2), 
                   library(matrixStats))) # send necessary functions to the cluster

  foreach(i = 1:25) %dopar% {  # loop for parallel processing
    all_fn(j_choice_a[[i]], gr_mean_a, pgrowthx=5, iter_num=(iter_num+i), tmax,
           tpops, popspec, n, n_boot, ngrps, count_thres, min_ts_length, c, k_choice_a[[i]], m_colnames, rmin, rmax, bootstrap_size, error=FALSE)
  }
  
stopCluster(cl) # stop the cluster

