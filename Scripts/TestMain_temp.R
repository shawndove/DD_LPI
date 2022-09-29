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
                   sdmean,
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
                   error=FALSE,
                   degrade,
                   mean_cv,
                   cv_sd,
                   clustlength,
                   endlength,
                   endpops_ratio,
                   endreveal_ratio) {
  
  # create directory to store files
  if(!dir.exists("TestData/")) {dir.create("TestData/")}
  if(!dir.exists(paste("TestData/", iter_num, sep=""))) {dir.create(paste("TestData/", iter_num, sep=""))}

  # create synthetic populations, assigned to different species
  if (pgrowthx == 5) {

    # create synthetic populations, assigned to different species
    all_pops_index <- pgrowth4.2(tpops, tmax, popmean, popvar, popspec)
    #all_pops_index <- readRDS(file=paste("TestData/", iter_num-40, "/saved_synth_", iter_num-40, "_raw.RData", sep=""))

  } else if (pgrowthx == 6) {

    # create synthetic populations, assigned to different species
    all_pops_index <- pgrowth4.3(tpops, tmax, popmean, popvar, popspec)
    
  } else if (pgrowthx == 7) {
    
    # create synthetic populations, assigned to different species
    all_pops_index <- pgrowth4.5(tpops, tmax, popmean, popvar, popspec, sdmean)
    
    
  } else if (pgrowthx == 8) {
    
    # create synthetic populations, assigned to different species
    all_pops_index <- pgrowth4.7(tpops, tmax, popmean, popvar, popspec, sdmean)
    
  } else if (pgrowthx == 9) {
    
    # create synthetic populations, assigned to different species
    all_pops_index <- pgrowth4.8(tpops, tmax, popmean, popvar, popspec, sdmean)
    
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
  
  #all_pops_index <- readRDS(file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_raw.RData", sep=""))
  
  # add a population ID column
  all_pops_index$PopID <- rownames(all_pops_index)
  
  if (error==TRUE) {
    
    cat(paste("Adding random sampling error to dataset.\n"))
    
    # add sampling error to dataset
    all_pops_error <- error_intr_fn2(all_pops_index, m_colnames, mean_cv, cv_sd)
    saveRDS(all_pops_error, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_error.RData", sep=""))
    
    cat(paste0("Degrading dataset.\n"))
    
    if (degrade=="normal") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn(all_pops_error, c, mlength, numobs)

    } else if (degrade=="endpoints") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_endpoints(all_pops_error, c, mlength, numobs)
      
    } else if (degrade=="cluster") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_clust(all_pops_error, c, mlength, numobs)
      
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
      all_pops_degraded <- degrade_ts_fn(all_pops_error2, c, mlength, numobs)
      
      # degraded endpops
      endpops_degraded <- degrade_ts_fn_clustend2(endpops, c, clustlength)
      saveRDS(endpops_degraded, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_endpops_degraded.RData", sep=""))
      
      
    } else if (degrade=="endreveal") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_endreveal(all_pops_error, c, mlength, numobs, endlength, endreveal_ratio)
      
    }
    
    saveRDS(all_pops_degraded, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_degraded.RData", sep=""))
    
  } else {
    
    cat(paste0("Degrading dataset.\n"))
    
    if (degrade=="normal") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn(all_pops_index, c, mlength, numobs)
      
    } else if (degrade=="endpoints") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_endpoints(all_pops_index, c, mlength, numobs)
      
    } else if (degrade=="cluster") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_clust(all_pops_index, c, mlength, numobs)
      
    } else if (degrade=="clustend") {
      
      if (samp_size >= ceiling((1-endpops_ratio)*nrow(all_pops_index))) {
        
        stop(paste("Sample size must be less than", 100 * (1 - endpops_ratio), "% of total populations to use clustend method.", sep=""))
        
      }
      
      # select populations to make endpoints
      endpops_ids <- sample(all_pops_index$PopID, ceiling(endpops_ratio*nrow(all_pops_index)))
      endpops <- all_pops_index[all_pops_index$PopID %in% endpops_ids,]
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn(all_pops_index, c, mlength, numobs)
      #all_pops_degraded <- readRDS(file=paste("TestData/", iter_num-40, "/saved_synth_", iter_num-40, "_degraded.RData", sep=""))
      
      # remove endpops from pops index
      all_pops_degraded <- all_pops_degraded[!(all_pops_degraded$PopID %in% endpops_ids),]
      
      # degraded endpops
      endpops_degraded <- degrade_ts_fn_clustend2(endpops, c, clustlength)
      saveRDS(endpops_degraded, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_endpops_degraded.RData", sep=""))
      
    } else if (degrade=="endreveal") {
      
      # degrade dataset
      all_pops_degraded <- degrade_ts_fn_endreveal(all_pops_index, c, mlength, numobs, endlength, endreveal_ratio)
      
    }
    
    saveRDS(all_pops_degraded, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_degraded.RData", sep=""))
    
  }
  
  cat(paste0("Removing time series with less than ", count_thres, " total observations and/or with a length of less than ", min_ts_length, ".\n"))
  
  # remove time series which are too short or have too few counts
  grp_data_culled <- cull_fn(all_pops_degraded, count_thres, min_ts_length, c)
  saveRDS(grp_data_culled, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_culled.RData", sep=""))
  
  cat(paste0("Creating ", bootstrap_size, " randomly sampled subsets of ", samp_size, " populations each.\n"))
  
  # if sample size is too large, set it to the total number of populations - 1
  if (nrow(grp_data_culled) < samp_size) {
    samp_size <- (nrow(grp_data_culled)-1)
  }
  
  # create list of x population ID matrices at a given sample size, where x is the number of sample bootstraps
  sample_pop_id_list <- list()
  for (i in 1:bootstrap_size) {
    sample_pop_id_list[[i]] <- sample(grp_data_culled$PopID, samp_size)
  }
  #sample_pop_id_list <- readRDS(file=paste("TestData/", iter_num-40, "/saved_synth_", iter_num-40, "_sample_pop_id_list.RData", sep=""))
  saveRDS(sample_pop_id_list, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_sample_pop_id_list.RData", sep=""))
  
  # if using clustend method, add selected endpops to all the samples and merge them back to the time series dataframe
  if (degrade=="clustend") {
    for (i in 1:bootstrap_size) {
      sample_pop_id_list[[i]] <- c(sample_pop_id_list[[i]], endpops_ids)
    }
    saveRDS(sample_pop_id_list, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_sample_pop_id_list.RData", sep=""))
   grp_data_culled <- rbind(grp_data_culled, endpops_degraded)
   saveRDS(grp_data_culled, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_culled.RData", sep=""))
  }
  
  
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
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, c, m_colnames, n=NA, n_boot=NA, iter_num, samp_size, bootstrap_size, method="nores")
  
  temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, c, m_colnames, n=NA, n_boot=NA, iter_num, samp_size, bootstrap_size, method="lambda")
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, c, m_colnames, n=n, n_boot=n_boot, iter_num, samp_size, bootstrap_size, method="lr")
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, c, m_colnames, n=n, n_boot=n_boot, iter_num, samp_size, bootstrap_size, method="my")
  
  #temp <- method_fn(grp_data_culled, sample_pop_id_list, msi_real, c, m_colnames, n=NA, n_boot=NA, iter_num, samp_size, bootstrap_size, method="lambda2")
  
  
  ## DATA GATHERING ##
  
  # get mean and standard deviation of the mean growth rate
  gr.stats.raw <- growth_rate_calc_fn3(all_pops_index, c)
  
  all_pops_completed <- complete_time_series(all_pops_degraded, c, m_colnames, calcsd=TRUE)
  
  gr.stats.degraded <- growth_rate_calc_fn3(all_pops_completed, c, model=TRUE)
  
  # create csv file with important info about the dataset
  info.dat <- data.frame(row.names = "1")
  info.dat$ID <- iter_num
  info.dat$synth_version <- pgrowthx
  info.dat$num_pops <- tpops
  info.dat$num_years <- tmax
  info.dat$num_resamps <- n
  info.dat$num_bootstraps <- n_boot
  info.dat$min_ts_length <- min_ts_length + 1
  info.dat$mean_gr_raw <- gr.stats.raw[[1]]
  info.dat$gr_sd_raw <- gr.stats.raw[[2]]
  info.dat$mean_sd_raw <- gr.stats.raw[[3]]
  info.dat$mean_gr_degraded <- gr.stats.degraded[[1]]
  info.dat$gr_sd_degraded <- gr.stats.degraded[[2]]
  info.dat$mean_sd_degraded <- gr.stats.degraded[[3]]
  info.dat$samp_size <- samp_size
  info.dat$pops_per_species <- popspec
  info.dat$mean_ts_length <- sum(!is.na(as.vector(all_pops_completed[,1:c]))) / nrow(all_pops_completed)
  info.dat$mean_num_obs <- sum(!is.na(as.vector(all_pops_degraded[,1:c]))) / nrow(all_pops_degraded)
  info.dat$degrade_type <- degrade
  info.dat$error <- error
  info.dat$mean_cv <- mean_cv
  info.dat$cv_sd <- cv_sd
  
  write.csv(info.dat, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_info.csv", sep=""))
  
  gr.stats.samples <- list()
  tslength.samples <- vector()
  numobs.samples <- vector()
  for (i in 1:length(sample_pop_id_list)) {
    
    temp <- all_pops_degraded[all_pops_degraded$PopID %in% sample_pop_id_list[[i]],]
    temp2 <- complete_time_series(temp, c, m_colnames, calcsd=TRUE)
    gr.stats.samples[[i]] <- growth_rate_calc_fn3(temp2, c, model=TRUE)
    tslength.samples[i] <- sum(!is.na(as.vector(temp2[,1:c]))) / nrow(temp2)
    numobs.samples[i] <- sum(!is.na(as.vector(temp[,1:c]))) / nrow(temp)
    
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
  write.csv(info.samples, file=paste("TestData/", iter_num, "/saved_synth_", iter_num, "_gr_samples.csv", sep=""))

}




## set parameters ----

n <- 100 # number of GAM resamples per population
n_boot <- 3000 # number of index bootstraps per species
bootstrap_size <- 20 # number of samples
samp_size <- 500 # number of populations to select for each sample
count_thres <- 2 # threshold at which number of counts is too low and we should not include this population
min_ts_length <- 2 # minimum length of time series to be included
startyear <- 1 # column name for first data point
endyear <- 50 # column name for final data point
tmax <- length(startyear:endyear) # number of years (data points per population)
tpops <- 10000 # total number of populations
popspec <- 20 # mean number of populations per species
ngrps <- 10 # number of groups to divide species into
m_colnames <- startyear:endyear # column names
c <- length(m_colnames) # number of years or columns (same as tmax)

## Setup Testing ----

iter_num <- 70000
gr_mean_a <- runif(3000, min = -0.08, max = 0.08)
gr_sd_vec_a <- runif(3000, min = 0.05, max = 0.5)
sd_mean <- runif(3000, min = 0.05, max = 2)
popspec <- 10
mlength_ <- round(runif(3000, min = 6, max = 40))
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- round(exp(runif(3000, min = log(50), max = log(10000))))
tmax <- 50
c <- tmax
tpops <- 10000
m_colnames <- 1:c
mean_cv <- 0.15
cv_sd <- 0.1
pgrowthx <- 7

#d_method <- rep(c("normal", "clustend", "endreveal"), each=50)
#cl_length <- rep(c(1, 3, 6, 10, 3, 
 #                  6, 10, 3, 6, 10, 
  #                 1, 1, 1), each=100)
#e_ratio <-rep(c(1, 0.2, 1), each=50)
#er_ratio <- rep(c(1, 1, 1, 1, 1, 
 #                 1, 1, 1, 1, 1, 
  #                0.25, 0.5, 1), each=100)

#iter_num <- 45040
#gr_mean_a <- rep(c(-0.1, -0.05, -0, 0.05, 0.1), 448)
#gr_sd_vec_a <- rep(rep(c(0.1, 0.3, 0.5, 0.7), each=5), 112)
#sd_mean <- rep(rep(c(0.15, 0.35, 0.55, 0.75), each=20), 28)
#popspec <- 10
#mlength <- rep(rep(c(10, 20, 30, 40), each=80), 7)
#numobs <- ceiling(0.5*mlength)
#samp_size <- rep(c(50, 100, 200, 500, 1000, 2000, 5000), each=320)
#tmax <- 50
#c <- tmax
#tpops <- 10000
#m_colnames <- 1:c
#mean_cv <- 0.15
#cv_sd <- 0.1
#pgrowthx <- 7
#j_choice_a <- rep(numobs, each=20)
#k_choice_a <- rep(mlength, 20)

## load external functions ----

source("Scripts/Functions_NG.R")
source("Scripts/TestMethods_NG.R")
source("Scripts/DegradeTSVersions_NG.R")


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
foreach(i = 1:3000) %dopar% {  # loop for parallel processing
  all_fn(popvar = gr_sd_vec_a[i], # variance in mean growth rate
         popmean = gr_mean_a[i], # mean growth rate
         sdmean =  sd_mean[i], # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         iter_num = (iter_num+i), 
         tmax = tmax, # number of years
         tpops = tpops, # total number of time series
         popspec = popspec, # mean number of populations per species
         n = n, # number of GAM resamples 
         n_boot = n_boot, # number of index bootstraps for each species
         ngrps = ngrps, # number of groups to divide time series into
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         c = c, # number of columns (years: same as tmax)
         samp_size = samp_size_[i], # number of time series in each sample
         m_colnames = m_colnames, # column names
         mlength = mlength_[i], # mean length of time series 
         numobs = numobs_[i], # mean number of observations in each time series
         bootstrap_size = bootstrap_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series to be added for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # fraction of populations to add for clustend degradation method
         endreveal_ratio=1)
}
sink()
stopCluster(cl) # stop the cluster


#########
# test different total trend lengths
#########

## Setup Testing ----

iter_num <- 83000
gr_mean_a <- 0
gr_sd_vec_a <- 0.4
sd_mean <- 0.6
popspec <- 20
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 200
tmax <- rep(c(30, 40, 50, 60, 70, 80, 90, 100), each=20)
c <- tmax
tpops <- 1000
m_colnames <- 1:c
mean_cv <- 0
cv_sd <- 0
pgrowthx <- 7

## load external functions ----

source("Scripts/Functions_NG.R")
source("Scripts/TestMethods_NG.R")
source("Scripts/DegradeTSVersions_NG.R")


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
foreach(i = 1:160) %dopar% {  # loop for parallel processing
  c <- tmax[i]
  m_colnames <- 1:c
  all_fn(popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         iter_num = (iter_num+i), 
         tmax = tmax[i], # number of years
         tpops = tpops, # total number of time series
         popspec = popspec, # mean number of populations per species
         n = n, # number of GAM resamples 
         n_boot = n_boot, # number of index bootstraps for each species
         ngrps = ngrps, # number of groups to divide time series into
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         c = c, # number of columns (years: same as tmax)
         samp_size = samp_size_, # number of time series in each sample
         m_colnames = m_colnames, # column names
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         bootstrap_size = bootstrap_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series to be added for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # fraction of populations to add for clustend degradation method
         endreveal_ratio=1)
}
sink()
stopCluster(cl) # stop the cluster


#########
# test different numbers of pops/species
#########

## Setup Testing ----

iter_num <- 84200
gr_mean_a <- 0
gr_sd_vec_a <- 0.2
sd_mean <- 0.2
popspec <- rep(c(5, 10, 15, 20, 30, 50), each=20)
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 200
tmax <- 50
c <- tmax
tpops <- 1000
m_colnames <- 1:c
mean_cv <- 0.15
cv_sd <- 0.1
pgrowthx <- 7

## load external functions ----

source("Scripts/Functions_NG.R")
source("Scripts/TestMethods_NG.R")
source("Scripts/DegradeTSVersions_NG.R")


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
foreach(i = 1:120) %dopar% {  # loop for parallel processing
  all_fn(popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = pgrowthx, # which time series generator to use
         iter_num = (iter_num+i), 
         tmax = tmax, # number of years
         tpops = tpops, # total number of time series
         popspec = popspec[i], # mean number of populations per species
         n = n, # number of GAM resamples 
         n_boot = n_boot, # number of index bootstraps for each species
         ngrps = ngrps, # number of groups to divide time series into
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         c = c, # number of columns (years: same as tmax)
         samp_size = samp_size_, # number of time series in each sample
         m_colnames = m_colnames, # column names
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         bootstrap_size = bootstrap_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series to be added for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # fraction of populations to add for clustend degradation method
         endreveal_ratio=1)
}
sink()
stopCluster(cl) # stop the cluster


#########
# test different timepoint distributions
#########

## Setup Testing ----

iter_num <- 84400
gr_mean_a <- 0
gr_sd_vec_a <- 0.2
sd_mean <- 0.2
popspec <- 20
mlength_ <- 25
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 200
tmax <- 50
c <- tmax
tpops <- 1000
m_colnames <- 1:c
mean_cv <- 0.15
cv_sd <- 0.1
pgrowthx <- 7

d_method <- rep(c("normal", "endpoints"), each=40)

## load external functions ----

source("Scripts/Functions_NG.R")
source("Scripts/TestMethods_NG.R")
source("Scripts/DegradeTSVersions_NG.R")


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
foreach(i = 1:80) %dopar% {  # loop for parallel processing
  all_fn(popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         iter_num = (iter_num+i), 
         tmax = tmax, # number of years
         tpops = tpops, # total number of time series
         popspec = popspec, # mean number of populations per species
         n = n, # number of GAM resamples 
         n_boot = n_boot, # number of index bootstraps for each species
         ngrps = ngrps, # number of groups to divide time series into
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         c = c, # number of columns (years: same as tmax)
         samp_size = samp_size_, # number of time series in each sample
         m_colnames = m_colnames, # column names
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         bootstrap_size = bootstrap_size, # number of samples
         error = TRUE,  # add observation error
         degrade = d_method[i], # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series to be added for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # fraction of populations to add for clustend degradation method
         endreveal_ratio=1)
}
sink()
stopCluster(cl) # stop the cluster


#########
# test different dataset sizes
#########

## Setup Testing ----

iter_num <- 82800
gr_mean_a <- 0
gr_sd_vec_a <- 0.4
sd_mean <- 0.6
popspec <- 10
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 50
tmax <- 50
c <- tmax
tpops <- rep(c(50, 100, 200, 500, 1000, 2000, 5000, 10000), each=10)
m_colnames <- 1:c
mean_cv <- 0
cv_sd <- 0
pgrowthx <- 7

## load external functions ----

source("Scripts/Functions_NG.R")
source("Scripts/TestMethods_NG.R")
source("Scripts/DegradeTSVersions_NG.R")


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
foreach(i = 1:80) %dopar% {  # loop for parallel processing
  all_fn(popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
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
         samp_size = samp_size_, # number of time series in each sample
         m_colnames = m_colnames, # column names
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         bootstrap_size = bootstrap_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series to be added for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # fraction of populations to add for clustend degradation method
         endreveal_ratio=1)
}
sink()
stopCluster(cl) # stop the cluster



#########
# test different amounts of observation error
#########

## Setup Testing ----

iter_num <- 80800
gr_mean_a <- 0
gr_sd_vec_a <- 0.2
sd_mean <- 0.2
popspec <- 20
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 200
tmax <- 50
c <- tmax
tpops <- 1000
m_colnames <- 1:c
mean_cv <- rep(c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2), each=20)
cv_sd <- rep(c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2), each=20)
pgrowthx <- 7

#d_method <- rep(c("normal", "clustend", "endreveal"), each=50)
#cl_length <- rep(c(1, 3, 6, 10, 3, 
#                  6, 10, 3, 6, 10, 
#                 1, 1, 1), each=100)
#e_ratio <-rep(c(1, 0.2, 1), each=50)
#er_ratio <- rep(c(1, 1, 1, 1, 1, 
#                 1, 1, 1, 1, 1, 
#                0.25, 0.5, 1), each=100)

## load external functions ----

source("Scripts/Functions_NG.R")
source("Scripts/TestMethods_NG.R")
source("Scripts/DegradeTSVersions_NG.R")


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
foreach(i = 1:160) %dopar% {  # loop for parallel processing
  all_fn(popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         iter_num = (iter_num+i), 
         tmax = tmax, # number of years
         tpops = tpops, # total number of time series
         popspec = popspec, # mean number of populations per species
         n = n, # number of GAM resamples 
         n_boot = n_boot, # number of index bootstraps for each species
         ngrps = ngrps, # number of groups to divide time series into
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         c = c, # number of columns (years: same as tmax)
         samp_size = samp_size_, # number of time series in each sample
         m_colnames = m_colnames, # column names
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         bootstrap_size = bootstrap_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv[i], # mean coefficient of variation for observation error
         cv_sd = cv_sd[i],# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series to be added for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # fraction of populations to add for clustend degradation method
         endreveal_ratio=1)
}
sink()
stopCluster(cl) # stop the cluster
