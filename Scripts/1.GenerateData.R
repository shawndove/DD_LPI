######################
### Author: Shawn Dove
######################

# This script contains code to generate datasets datasets for building a model of 
# trend reliability in the LPI and for testing various parameters of the data for 
# effects on trend reliability.
#
# The script calls MainFunction.R to generate sets of data for different types of testing.
# The process can be quite slow, especially if each dataset contains thousands of time series
# and there are many datasets. Parallel processing* (see note below) is used to speed things 
# up, but generating the 3000 datasets used to build the accuracy model will likely take
# several days, even with a fast computer.
#
# *Note: If you do want to use parallel processing, skip the "parallel processing"
# section before each function call and edit the "foreach" line before the function call 
# to change %dopar% to %do%.
# Example: "foreach(i = 1:80) %dopar% {" becomes "foreach(i = 1:80) %do% {"
# R will then execute the loop in series instead of in parallel. But be aware that
# this will massively increase the time it takes to complete the loop.
#
# The first step is to load all packages and external functions, as these are needed
# by MainFunction.R and are not loaded by the function itself.
#
# Next, check that the parameters and directory to store files are appropriate and run
# these lines. The directory will be created as a subdirectory within your working
# directory if it does not already exist. If you are using parallel processing, be sure 
# to alter the "no_cores" setting to reflect the number of CPU cores you have or want to use.
#
# Find the section heading that corresponds to the type of datasets you want to generate.
# Check the parameter settings. "iter_num" can be any number you want to identify you datsets. 
# The first dataset generated will have an ID number equal to iter_num + 1, and the number
# will increase by 1 with each subsequent dataset. "dir_name" should be an appropriately
# titled subdirectory under the name you chose for "test_dir". Parameter settings are briefly
# commented in the trend accuracy model section and on each function call. However, they are
# explained more extensively in the MainFunction.R script.
# 
# When the process is finished, the output within R will be "NULL". All the data is saved to
# files, and the necessary information for analysis will be retrieved in the next script.


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

source("Scripts/MainFunction.R")
source("Scripts/GenerateTSFunction.R")
source("Scripts/ObservationErrorFunction.R")
source("Scripts/DegradeTSFunction.R")
source("Scripts/CullFunction.R")
source("Scripts/IndexLambdaFunction.R")
source("Scripts/InterpolateFunction.R")
source("Scripts/PopGAMFunction.R")
source("Scripts/SpeciesIndexFunction.R")
source("Scripts/GroupIndexFunction.R")
source("Scripts/ConfidenceIntervalFunction.R")
source("Scripts/GrowthRateCalcFunction.R")


## set parameters ----

resamp_size <- 20 # number of times to sample each dataset
count_thres <- 2 # threshold at which number of counts is too low and we should not include this population
min_ts_length <- 2 # minimum length of time series to be included
no_cores <- 8 # the number of cores to be used for parallel processing

## set directory to store datasets ----

test_dir <- "TestData"
if(!dir.exists(paste(test_dir, "/", sep=""))) {dir.create(paste(test_dir, "/", sep=""))} # create directory

#########
## Generate datasets for trend accuracy model ----
#########

## set parameters ----

iter_num <- 70000 # ID number to include in filenames, iterated for each dataset
gr_mean_a <- runif(3000, min = -0.08, max = 0.08) # mean growth rate
gr_sd_vec_a <- runif(3000, min = 0.05, max = 0.5) # growth rate standard deviation
sd_mean <- runif(3000, min = 0.05, max = 2) # st. dev. of the mean growth rate
popspec <- 10 # mean number of populations per species
mlength_ <- round(runif(3000, min = 6, max = 40)) # mean time series length
numobs_ <- ceiling(0.5*mlength_) # mean number of observations per time series
samp_size_ <- round(exp(runif(3000, min = log(50), max = log(10000)))) # sample size
tmax <- 50 # total number of years in dataset (overall trend length)
tpops <- 10000 # total number of populations in dataset
mean_cv <- 0.15 # mean observation error
cv_sd <- 0.1 # st. dev. of observation error
pgrowthx <- 7 # ID for time series generation function
dir_name <- "TestData/ModelData" # name of directory to store the datasets

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

## parallel processing ----

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:3000) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
         dir_name = dir_name, # name of directory to save the dataset(s)
         popvar = gr_sd_vec_a[i], # variance in mean growth rate
         popmean = gr_mean_a[i], # mean growth rate
         sdmean =  sd_mean[i], # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         tmax = tmax, # number of years
         tpops = tpops, # total number of time series
         popspec = popspec, # mean number of populations per species
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         samp_size = samp_size_[i], # number of time series in each sample
         mlength = mlength_[i], # mean length of time series 
         numobs = numobs_[i], # mean number of observations in each time series
         resamp_size = resamp_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # portion of populations to apply clustend degradation method
         endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster



#########
# test total trend length ----
#########

## set parameters ----

iter_num <- 80000
gr_mean_a <- 0
gr_sd_vec_a <- 0.2
sd_mean <- 0.2
popspec <- 20
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 200
tmax <- rep(c(30, 40, 50, 60, 70, 80, 90, 100), each=20)
tpops <- 1000
mean_cv <- 0.15
cv_sd <- 0.1
dir_name <- "TestData/TrendLength"

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

## parallel processing ----

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:160) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
         dir_name = dir_name, # name of directory to save the dataset(s)
         popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         tmax = tmax[i], # number of years
         tpops = tpops, # total number of time series
         popspec = popspec, # mean number of populations per species
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         samp_size = samp_size_, # number of time series in each sample
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         resamp_size = resamp_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # portion of populations to apply clustend degradation method
         endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster



#########
# test observation error
#########

## set parameters ----

iter_num <- 81000
gr_mean_a <- 0
gr_sd_vec_a <- 0.2
sd_mean <- 0.2
popspec <- 20
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 200
tmax <- 50
tpops <- 1000
mean_cv <- rep(c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2), each=20)
cv_sd <- rep(c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2), each=20)
dir_name <- "TestData/ObsError"

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

## parallel processing ----

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:160) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
         dir_name = dir_name, # name of directory to save the dataset(s)
         popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         tmax = tmax, # number of years
         tpops = tpops, # total number of time series
         popspec = popspec, # mean number of populations per species
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         samp_size = samp_size_, # number of time series in each sample
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         resamp_size = resamp_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv[i], # mean coefficient of variation for observation error
         cv_sd = cv_sd[i],# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # portion of populations to apply clustend degradation method
         endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster


#########
# test dataset size and sample size ----
#########


### part 1 of 4: sample size 50

## set parameters ----

iter_num <- 82000
gr_mean_a <- 0
gr_sd_vec_a <- 0.3
sd_mean <- 0.4
popspec <- 20
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 50
tmax <- 50
tpops <- rep(c(50, 100, 200, 500, 1000, 2000, 5000, 10000), each=10)
mean_cv <- 0.15
cv_sd <- 0.1
dir_name <- "TestData/DSSize"

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

## parallel processing ----

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:80) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
          dir_name = dir_name, # name of directory to save the dataset(s)
          popvar = gr_sd_vec_a, # variance in mean growth rate
          popmean = gr_mean_a, # mean growth rate
          sdmean =  sd_mean, # mean of standard deviations in growth rates
          pgrowthx = 7, # which time series generator to use
          tmax = tmax, # number of years
          tpops = tpops[i], # total number of time series
          popspec = popspec, # mean number of populations per species
          count_thres = count_thres, # minimum number of population counts
          min_ts_length = min_ts_length, # minimum time series length
          samp_size = samp_size_, # number of time series in each sample
          mlength = mlength_, # mean length of time series 
          numobs = numobs_, # mean number of observations in each time series
          resamp_size = resamp_size, # number of samples
          error = TRUE,  # add observation error
          degrade = "normal", # data degradation method
          mean_cv = mean_cv, # mean coefficient of variation for observation error
          cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
          clustlength=10, # length of time series for clustend degradation method
          endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
          endpops_ratio=1,  # portion of populations to apply clustend degradation method
          endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster


### part 2 of 4: sample size 100 ----

## set parameters ----

iter_num <- 82100
gr_mean_a <- 0
gr_sd_vec_a <- 0.3
sd_mean <- 0.4
popspec <- 20
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 100
tmax <- 50
tpops <- rep(c(100, 200, 500, 1000, 2000, 5000, 10000), each=10)
mean_cv <- 0.15
cv_sd <- 0.1
dir_name <- "TestData/DSSize"

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

## parallel processing ----

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:70) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
         dir_name = dir_name, # name of directory to save the dataset(s)
         popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         tmax = tmax, # number of years
         tpops = tpops[i], # total number of time series
         popspec = popspec, # mean number of populations per species
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         samp_size = samp_size_, # number of time series in each sample
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         resamp_size = resamp_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # portion of populations to apply clustend degradation method
         endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster


### part 3 of 4: sample size 200 ----

## set parameters ----

iter_num <- 82200
gr_mean_a <- 0
gr_sd_vec_a <- 0.3
sd_mean <- 0.4
popspec <- 20
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 200
tmax <- 50
tpops <- rep(c(200, 500, 1000, 2000, 5000, 10000), each=10)
mean_cv <- 0.15
cv_sd <- 0.1
dir_name <- "TestData/DSSize"

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

## parallel processing ----

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:60) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
         dir_name = dir_name, # name of directory to save the dataset(s)
         popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         tmax = tmax, # number of years
         tpops = tpops[i], # total number of time series
         popspec = popspec, # mean number of populations per species
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         samp_size = samp_size_, # number of time series in each sample
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         resamp_size = resamp_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # portion of populations to apply clustend degradation method
         endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster


### part 4 of 4: sample size 500 ----

## set parameters ----

iter_num <- 82300
gr_mean_a <- 0
gr_sd_vec_a <- 0.3
sd_mean <- 0.4
popspec <- 20
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 500
tmax <- 50
tpops <- rep(c(500, 1000, 2000, 5000, 10000), each=10)
mean_cv <- 0.15
cv_sd <- 0.1
dir_name <- "TestData/DSSize"

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:50) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
         dir_name = dir_name, # name of directory to save the dataset(s)
         popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         tmax = tmax, # number of years
         tpops = tpops[i], # total number of time series
         popspec = popspec, # mean number of populations per species
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         samp_size = samp_size_, # number of time series in each sample
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         resamp_size = resamp_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # portion of populations to apply clustend degradation method
         endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster


#########
# test mean and distribution of the number of populations per species ----
#########

### part 1 of 3: negative bionomial distribution ----

## set parameters ----

iter_num <- 83000
gr_mean_a <- 0
gr_sd_vec_a <- 0.2
sd_mean <- 0.2
popspec <- rep(c(5, 10, 15, 20, 30, 50), each=20)
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 200
tmax <- 50
tpops <- 1000
mean_cv <- 0.15
cv_sd <- 0.1
pgrowthx <- 9 # negative binomial distribution
dir_name <- "TestData/PPSDist"

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

## parallel processing ----

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:120) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
         dir_name = dir_name, # name of directory to save the dataset(s)
         popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = pgrowthx, # which time series generator to use
         tmax = tmax, # number of years
         tpops = tpops, # total number of time series
         popspec = popspec[i], # mean number of populations per species
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         samp_size = samp_size_, # number of time series in each sample
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         resamp_size = resamp_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # portion of populations to apply clustend degradation method
         endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster


### part 2 of 3: discretized exponential distribution ----

## set parameters ----

iter_num <- 83200
gr_mean_a <- 0
gr_sd_vec_a <- 0.2
sd_mean <- 0.2
popspec <- rep(c(5, 10, 15, 20, 30, 50), each=20)
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 200
tmax <- 50
tpops <- 1000
mean_cv <- 0.15
cv_sd <- 0.1
pgrowthx <- 8 # discretized exponential distribution
dir_name <- "TestData/PPSDist"

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

## parallel processing ----

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:120) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
         dir_name = dir_name, # name of directory to save the dataset(s)
         popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = pgrowthx, # which time series generator to use
         tmax = tmax, # number of years
         tpops = tpops, # total number of time series
         popspec = popspec[i], # mean number of populations per species
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         samp_size = samp_size_, # number of time series in each sample
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         resamp_size = resamp_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # portion of populations to apply clustend degradation method
         endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster


### part 3 of 3: normal distribution ----

## set parameters ----

iter_num <- 83400
gr_mean_a <- 0
gr_sd_vec_a <- 0.2
sd_mean <- 0.2
popspec <- rep(c(5, 10, 15, 20, 30, 50), each=20)
mlength_ <- 20
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 200
tmax <- 50
tpops <- 1000
mean_cv <- 0.15
cv_sd <- 0.1
pgrowthx <- 7 # normal distribution
dir_name <- "TestData/PPSDist"

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

## parallel processing ----

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:120) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
         dir_name = dir_name, # name of directory to save the dataset(s)
         popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = pgrowthx, # which time series generator to use
         tmax = tmax, # number of years
         tpops = tpops, # total number of time series
         popspec = popspec[i], # mean number of populations per species
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         samp_size = samp_size_, # number of time series in each sample
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         resamp_size = resamp_size, # number of samples
         error = TRUE,  # add observation error
         degrade = "normal", # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # portion of populations to apply clustend degradation method
         endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster


#########
# test different time point distributions ----
#########

## set parameters ----

iter_num <- 84000
gr_mean_a <- 0
gr_sd_vec_a <- 0.2
sd_mean <- 0.2
popspec <- 20
mlength_ <- 25
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- 200
tmax <- 50
tpops <- 1000
mean_cv <- 0.15
cv_sd <- 0.1
d_method <- rep(c("normal", "endpoints"), each=40)
dir_name <- "TestData/TPDist"

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

## parallel processing ----

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:80) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
         dir_name = dir_name, # name of directory to save the dataset(s)
         popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         tmax = tmax, # number of years
         tpops = tpops, # total number of time series
         popspec = popspec, # mean number of populations per species
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         samp_size = samp_size_, # number of time series in each sample
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         resamp_size = resamp_size, # number of samples
         error = TRUE,  # add observation error
         degrade = d_method[i], # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=1,  # portion of populations to apply clustend degradation method
         endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster


#########
# solutions testing ----
#########

## set parameters ----

iter_num <- 85000
gr_mean_a <- 0
gr_sd_vec_a <- 0.3
sd_mean <- 0.4
popspec <- 10
mlength_ <- 14
numobs_ <- ceiling(0.5*mlength_)
samp_size_ <- rep(c(200, 200, 200, 400), each=20)
tmax <- 50
tpops <- 1000
mean_cv <- 0.15
cv_sd <- 0.1
d_method <- rep(c("normal", "clustend", "endreveal", "normal"), each=20) # degradation method
ep_ratio <- rep(c(1, 0.2, 1, 1), each=20) # endpops ratio (portion to apply "clustend" to)
dir_name <- "TestData/Solutions"

if(!dir.exists(paste(dir_name, "/", sep=""))) {dir.create(paste(dir_name, "/", sep=""))} # create directory

## parallel processing ----

cl <- makeCluster(no_cores, outfile=paste(dir_name, "/output.txt", sep="")) # create cluster for parallel processing
registerDoSNOW(cl) # register the cluster
clusterEvalQ(cl, c(library(tcltk),  # send necessary functions to the cluster
                   library(dplyr), 
                   library(MASS), 
                   library(GET), 
                   library(mgcv), 
                   library(reshape2), 
                   library(matrixStats),
                   library(TSdist)))

## call main function ----

foreach(i = 1:80) %dopar% {  # loop for parallel processing
  main_fn(iter_num = (iter_num+i), # ID number for dataset
         dir_name = dir_name, # name of directory to save the dataset(s)
         popvar = gr_sd_vec_a, # variance in mean growth rate
         popmean = gr_mean_a, # mean growth rate
         sdmean =  sd_mean, # mean of standard deviations in growth rates
         pgrowthx = 7, # which time series generator to use
         tmax = tmax, # number of years
         tpops = tpops, # total number of time series
         popspec = popspec, # mean number of populations per species
         count_thres = count_thres, # minimum number of population counts
         min_ts_length = min_ts_length, # minimum time series length
         samp_size = samp_size_[i], # number of time series in each sample
         mlength = mlength_, # mean length of time series 
         numobs = numobs_, # mean number of observations in each time series
         resamp_size = resamp_size, # number of samples
         error = TRUE,  # add observation error
         degrade = d_method[i], # data degradation method
         mean_cv = mean_cv, # mean coefficient of variation for observation error
         cv_sd = cv_sd,# standard deviation of coefficient of variation for observation error
         clustlength=10, # length of time series for clustend degradation method
         endlength=1, # how many observations to reveal at the final year(s) for endreveal degradation method
         endpops_ratio=ep_ratio[i],  # portion of populations to apply clustend degradation method
         endreveal_ratio=1) # portion of populations to apply endreveal degradation method
}

stopCluster(cl) # stop the cluster

