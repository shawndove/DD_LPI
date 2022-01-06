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
library(TSdist)

## load external functions ----

source("Scripts/Functions_NG.R")
source("Scripts/TestMethods_NG.R")

## Main Function ----

all_fn <- function(popvar, 
                   popmean, 
                   pgrowthx, 
                   iter_num, 
                   tmax, 
                   tpops, 
                   popspec, 
                   n, 
                   n_boot, 
                   ngrps,
                   count_thres, 
                   min_ts_length, 
                   c, 
                   samp_size, 
                   m_colnames, 
                   mlength, 
                   numobs, 
                   bootstrap_size, 
                   error=FALSE) {
  
  # create directory to store files
  if(!dir.exists("TestData/")) {dir.create("TestData/")}
  if(!dir.exists(paste("TestData/", iter_num, sep=""))) {dir.create(paste("TestData/", iter_num, sep=""))}

  # create synthetic populations, assigned to different species
  if (pgrowthx == 5) {
    
    # create synthetic populations, assigned to different species
    all_pops_index <- pgrowth4.2(tpops, tmax, popmean, popvar, popspec)

  } else("Check your synthetic data generator function input setting.")
  
  # save raw synthetic dataset
  saveRDS(all_pops_index, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_raw.RData", sep=""))
  
  cat(paste("Constructing dataset:", "\nPopulations:", tpops, "\nSpecies:", length(unique(all_pops_index$SpecID)), "\n\n"))
  
  cat(paste0("Plotting geometric mean of the dataset.\n"))
  
  # create geometric mean and plot it
  index_geomean <- exp(colMeans(log(all_pops_index[,1:c])))

  png(file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_geometric_mean_raw.png", sep=""),
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
  
  png(file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_raw_gam_geometric_mean.png", sep=""),
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
    saveRDS(all_pops_error, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_error.RData", sep=""))
    
    cat(paste0("Degrading dataset.\n"))
    
    # degrade dataset
    all_pops_degraded <- degrade_ts_fn(all_pops_error, c, mlength, numobs)
    saveRDS(all_pops_degraded, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_degraded.RData", sep=""))
    
  } else {
    
    cat(paste0("Degrading dataset.\n"))
    
    # degrade dataset
    all_pops_degraded <- degrade_ts_fn(all_pops_index, c, mlength, numobs)
    saveRDS(all_pops_degraded, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_degraded.RData", sep=""))
    
  }
  
  cat(paste0("Removing time series with less than ", count_thres, " total observations and/or with a length of less than ", min_ts_length, ".\n"))
  
  # remove time series which are too short or have too few counts
  grp_data_culled <- cull_fn(all_pops_degraded, count_thres, min_ts_length, c)
  saveRDS(grp_data_culled, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_culled.RData", sep=""))
  
  cat(paste0("Creating ", bootstrap_size, " randomly sampled subsets of ", samp_size, " populations each.\n"))
  
  # if sample size is too large, set it to 90% of the total number of populations
  if (nrow(grp_data_culled) < samp_size) {
    samp_size <- ceiling(0.9*nrow(grp_data_culled))
  }
  
  # create list of x population ID matrices at a given sample size, where x is the number of sample bootstraps
  sample_pop_id_list <- list()
  for (i in 1:bootstrap_size) {
    sample_pop_id_list[[i]] <- sample(grp_data_culled$PopID, samp_size)
  }
  saveRDS(sample_pop_id_list, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_sample_pop_id_list.RData", sep=""))
  
  
  ## TRUE TREND ##
  
  cat(paste0("Creating ", length(unique(all_pops_index$SpecID)), " species indices for true trend.\n"))
  
  # create species indices from the population indices
  full_spec_real <- species_index_fn(all_pops_index, c)
  saveRDS(full_spec_real, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_species_indices_TrueTrend.RData", sep=""))
  
  cat(paste0("Creating msi for true trend.\n"))
  
  # create multi species indices from the group indices
  msi_real <- group_index_fn(full_spec_real, c, m_colnames)
  saveRDS(msi_real, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_TrueTrend.RData", sep=""))
  
  cat(paste0("Creating msi confidence intervals for true trend.\n"))
  
  # create confidence intervals for the multi species index
  msi_ci_real <- ci_fn(full_spec_real, c, m_colnames, grouplevel=1)
  saveRDS(msi_ci_real, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_ci_TrueTrend.RData", sep=""))
  
  cat(paste0("Plotting true trend.\n"))
  
  # plot trend
  plot_msi <- msi_real
  plot_ci <- msi_ci_real
  
  png(file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_msi_TrueTrend.png", sep=""),
      width = 1280, height = 720, units = "px", pointsize = 12, bg = "white")
  
  par(mar=c(6,6,6,2) + 0.1)
  plot(m_colnames, plot_msi, xlab="", ylab="", type="l", lty=1, ylim=c(max(c(0, min(plot_ci)-10)),max(plot_ci)+10), lwd=3, frame.plot=TRUE,
       cex.lab=1.5, cex.axis=1.5,
       panel.first = polygon(c(min(m_colnames):max(m_colnames), max(m_colnames):min(m_colnames)), c(plot_ci[2,], rev(plot_ci[1,])), col="grey90", border = NA))
  lines(m_colnames, plot_ci[2,], lty=2, lwd=2)
  lines(m_colnames, plot_ci[1,], lty=2, lwd=2)
  lines(m_colnames, g.gam, lty=1, lwd=3, col="red")
  grid(col="grey85", lty=2)
  mtext(text = "Year", side = 1, line = 4, cex = 1.5)
  mtext(text = "Index", side = 2, line = 3.5, cex = 1.5)
  
  dev.off()
  
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, c, m_colnames, n=NA, n_boot=NA, iter_num, samp_size, bootstrap_size, method="nores")
  
  temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, c, m_colnames, n=NA, n_boot=NA, iter_num, samp_size, bootstrap_size, method="lambda")
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, c, m_colnames, n=n, n_boot=n_boot, iter_num, samp_size, bootstrap_size, method="lr")
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, c, m_colnames, n=n, n_boot=n_boot, iter_num, samp_size, bootstrap_size, method="my")
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, c, m_colnames, n=NA, n_boot=NA, iter_num, samp_size, bootstrap_size, method="lambda2")
  
  
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
  info.dat$samp_size <- samp_size
  info.dat$pops_per_species <- popspec
  info.dat$mean_ts_length <- mlength
  info.dat$mean_num_obs <- numobs
  
  #info.dat$mean_gr_rebuilt <- gr.stats.rebuilt[1]
  #info.dat$gr_sd_rebuilt <- gr.stats.rebuilt[2]
  #info.dat$pops_per_species <- nrow(new.grp_data) / length(spec_ids) # average populations per species in culled data
  #info.dat$pops_per_species_raw <- nrow(all_pops_index) / length(unique(SpecID_vec)) # avg pops per species in unculled data
  
  write.csv(info.dat, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_info.csv", sep=""))
  
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

iter_num <- 11700
gr_mean_a <- 0
gr_sd_vec_a <- 0.2
popspec <- 10
mlength <- 15
numobs <- 7
samp_size <- rep(c(20, 40, 60, 80, 100, 125, 150, 175, 200, 250, 300, 350), 8)
tmax <- 50
c <- tmax
tpops <- rep(rep(c(100, 200, 400), each=4), 8)
m_colnames <- 1:c
#j_choice_a <- rep(numobs, each=20)
#k_choice_a <- rep(mlength, 20)

## load external functions ----

source("Scripts/Functions_NG.R")
source("Scripts/TestMethods_NG.R")


no_cores <- 8 # the number of cores to be used for parallel processing
cl <- makeCluster(no_cores, outfile="TestData/output.txt") # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

sink(file="console_output.txt", split=TRUE)
# call the main function
foreach(i = 1:96) %dopar% {  # loop for parallel processing
  all_fn(popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         pgrowthx = 5, # which time series generator to use
         iter_num = (iter_num+i), 
         tmax = tmax, # number of years
         tpops = tpops[i], # total number of time series
         popspec = popspec, # mean number of populations per species
         n = n, # number of GAM resamples 
         n_boot = n_boot, # number of index bootstraps for each species
         ngrps = ngrps, # number of groups to divide time series into
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         c = c, # number of columns (years: same as tmax)
         samp_size = samp_size[i], # number of time series in each sample
         m_colnames = m_colnames, # column names
         mlength = mlength, # mean length of time series 
         numobs = numobs, # mean number of observations in each time series
         bootstrap_size = bootstrap_size, # number of samples
         error = FALSE) # add sampling error
}
sink()
stopCluster(cl) # stop the cluster
