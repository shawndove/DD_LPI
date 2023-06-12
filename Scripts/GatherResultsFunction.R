######################
### Author: Shawn Dove
######################

# This functions gathers and compiles data from test datasets in the directory the user provides.
# It returns a list with two data frames. The first contains information from all individual
# dataset samples and is only needed for building the accuracy model. The second contains
# mean data and is used for all analyses. To get only the second data frame, add [[2]]
# at the end of the function call, e.g. testdata <- datacompile_fn(dir_name)[[2]].

gather_results_fn <- function(dir_name) {
  
  dir_names <- list.dirs(path=dir_name, 
                         full.names = TRUE, 
                         recursive = FALSE)
  dir_names <- paste(dir_names, "/", sep="")
  
  
  # create vectors to hold file lists
  info_list <- vector()
  gr_samples_list <- vector()
  td_list <- vector()
  tdm_list <- vector()
  tsd_list <- vector()
  tscd_list <- vector()
  samp_list <- vector()
  
  
  # get lists of csv files in the test directory
  for (i in 1:length(dir_names)) {
    
    info_list[i] <- list.files(paste(dir_names[i], 
                                     sep=""), 
                               pattern = "info",
                               full.names = TRUE,
                               recursive = FALSE)
    
    gr_samples_list[i] <- list.files(paste(dir_names[i], 
                                           sep=""), 
                                     pattern = "gr_samples",
                                     full.names = TRUE,
                                     recursive = FALSE)
    
    tdm_list[i] <- list.files(paste(dir_names[i], 
                                    sep=""), 
                              pattern = "trend_dev_sampled_mean",
                              full.names = TRUE,
                              recursive = FALSE)
    
    td_list[i] <- list.files(paste(dir_names[i], 
                                   sep=""), 
                             pattern = "trend_dev_sampled_list",
                             full.names = TRUE,
                             recursive = FALSE)
    
    tsd_list[i] <- list.files(paste(dir_names[i],
                                    sep=""),
                              pattern = "culled",
                              full.names = TRUE,
                              recursive = FALSE)
    
    tscd_list[i] <- list.files(paste(dir_names[i],
                                     sep=""),
                               pattern = "completed_time_series",
                               full.names = TRUE,
                               recursive = FALSE)
    
    samp_list[i] <- list.files(paste(dir_names[i],
                                     sep=""),
                               pattern = "sample_pop_id_list",
                               full.names = TRUE,
                               recursive = FALSE)
    
  }
  
  
  
  # create a set of temporary lists and vectors to hold various data from saved files
  trenddev_tl <- list()
  trenddev_m_tl <- list()
  samp_sdgr_tl <- list()
  samp_meangr_tl <- list()
  samp_meansd_tl <- list()
  samp_meants_tl <- list()
  samp_sdgr_m_tl <- list()
  samp_meangr_m_tl <- list()
  samp_meansd_m_tl <- list()
  samp_meants_m_tl <- list()
  sampspecsize_tl <- list()
  sampspecsize_m_tl <- vector()
  tyears_tl <- vector()
  tpops_tl <- vector()
  boots_tl <- vector()
  sampsize_tl <- vector()
  samppercent_tl <- vector()
  popspec_tl <- vector()
  sdgr_tl <- vector()
  meangr_tl <- vector()
  meansd_tl <- vector()
  iternum_tl <- vector()
  tsgenver_tl <- vector()
  tspec_tl <- vector()
  meanobs_tl <- vector()
  meanlength_tl <- vector()
  degrade_type_tl <- vector()
  
  # get data from files
  for (i in 1:length(info_list)) {
    
    # import data from csv files
    info <- read.csv(info_list[i])
    gr_samples <- read.csv(gr_samples_list[i])
    trenddev <- read.csv(td_list[i])
    trenddev_m <- read.csv(tdm_list[i])
    tsdata <- readRDS(tsd_list[i])
    tscdata <- readRDS(tscd_list[i])
    sampdata <- readRDS(samp_list[i])
    
    # read data into lists
    iternum_tl[i] <- info$ID[1]
    tsgenver_tl[i] <- info$synth_version[1]
    tyears_tl[i] <- info$num_years[1]
    meangr_tl[i] <- info$mean_gr_raw[1]
    sdgr_tl[i] <- info$gr_sd_raw[1]
    meansd_tl[i] <- info$mean_sd_raw[1]
    degrade_type_tl[i] <- info$degrade_type[1]
    trenddev_tl[[i]] <- trenddev$MSI
    trenddev_m_tl[[i]] <- trenddev_m[[2]]
    samp_sdgr_tl[[i]] <- gr_samples$gr_sd_sample
    samp_meangr_tl[[i]] <- gr_samples$mean_gr_sample
    samp_meansd_tl[[i]] <- gr_samples$mean_sd_sample
    samp_meants_tl[[i]] <- gr_samples$mean_ts_length
    samp_sdgr_m_tl[[i]] <- mean(gr_samples$gr_sd_sample, na.rm=TRUE)
    samp_meangr_m_tl[[i]] <- mean(gr_samples$mean_gr_sample, na.rm=TRUE)
    samp_meansd_m_tl[[i]] <- mean(gr_samples$mean_sd_sample, na.rm=TRUE)
    samp_meants_m_tl[[i]] <- mean(gr_samples$mean_ts_length, na.rm=TRUE)
    boots_tl[i] <- length(trenddev_tl[[i]])
    tpops_tl[i] <- nrow(tsdata)
    tspec_tl[i] <- length(unique(tsdata$SpecID))
    popspec_tl[i] <- round((tpops_tl[i] / tspec_tl[i]), digits = 1)
    sampsize_tl[i] <- length(sampdata[[1]])
    samppercent_tl[i] <- (sampsize_tl[i] / tpops_tl[i]) * 100
    meanobs_tl[i] <- sum(!is.na(as.vector(tsdata[,1:tyears_tl[i]]))) / tpops_tl[i]
    # meanlength_tl[i] <- sum(!is.na(as.vector(tscdata[,1:tyears_tl[i]]))) / tpops_tl[i]
    meanlength_tl[i] <- info$mean_ts_length[1]
    
    sampspecsize_temp <- vector()
    for (j in 1:boots_tl[i]) {
      
      sampspecsize_temp[j] <- length(unique(tsdata$SpecID[tsdata$PopID %in% sampdata[[j]]]))
      
    }
    sampspecsize_tl[[i]] <- sampspecsize_temp
    sampspecsize_m_tl[i] <- mean(sampspecsize_temp, na.rm=TRUE)
    
    
  }
  
  ## create data frame to hold results
  test_results <- data.frame(matrix(NA, ncol = 20, nrow = length(info_list)*20))
  test_results_m <- data.frame(matrix(NA, ncol = 20, nrow = length(info_list)))
  
  # name columns
  colnames(test_results) <- c("ID",
                              "MeanGR",
                              "SDGR",
                              "MeanSD",
                              "MeanGRSamp",
                              "SDGRSamp",
                              "MeanSDSamp",
                              "MeanTSSamp",
                              "SampSize",
                              "SampPercent",
                              "SampSpecSize",
                              "MeanTSLength",
                              "MeanNumObs",
                              "TotalPops",
                              "TotalSpec",
                              "PopSpec",
                              "TotalYears",
                              "TSGenVersion",
                              "TrendDev",
                              "DegradeType")
  
  colnames(test_results_m) <- c("ID",
                                "MeanGR",
                                "SDGR",
                                "MeanSD",
                                "MeanGRSamp",
                                "SDGRSamp",
                                "MeanSDSamp",
                                "MeanTSSamp",
                                "SampSize",
                                "SampPercent",
                                "SampSpecSize",
                                "MeanTSLength",
                                "MeanNumObs",
                                "TotalPops",
                                "TotalSpec",
                                "PopSpec",
                                "TotalYears",
                                "TSGenVersion",
                                "TrendDev",
                                "DegradeType")
  
  counter <- 1
  counter_m <- 1
  
  for (i in 1:length(info_list)) {
    
    for (j in 1:boots_tl[i]) {
      
      test_results$ID[counter] <- paste(iternum_tl[i], j, sep="_")
      test_results$MeanGR[counter] <- meangr_tl[i]
      test_results$SDGR[counter] <- sdgr_tl[i]
      test_results$MeanSD[counter] <- meansd_tl[i]
      test_results$MeanGRSamp[counter] <- samp_meangr_tl[[i]][j]
      test_results$SDGRSamp[counter] <- samp_sdgr_tl[[i]][j]
      test_results$MeanSDSamp[counter] <- samp_meansd_tl[[i]][j]
      test_results$MeanTSSamp[counter] <- samp_meants_tl[[i]][j]
      test_results$SampSize[counter] <- sampsize_tl[i]
      test_results$SampPercent[counter] <- samppercent_tl[i]
      test_results$SampSpecSize[counter] <- sampspecsize_tl[[i]][j]
      test_results$MeanTSLength[counter] <- meanlength_tl[i]
      test_results$MeanNumObs[counter] <- meanobs_tl[i]
      test_results$TotalPops[counter] <- tpops_tl[i]
      test_results$TotalSpec[counter] <- tspec_tl[i]
      test_results$PopSpec[counter] <- popspec_tl[i]
      test_results$TotalYears[counter] <- tyears_tl[i]
      test_results$TSGenVersion[counter] <- tsgenver_tl[i]
      test_results$TrendDev[counter] <- trenddev_tl[[i]][j]
      test_results$DegradeType[counter] <- degrade_type_tl[i]
      
      counter <- counter + 1
      
    }
    
    test_results_m$ID[counter_m] <- iternum_tl[i]
    test_results_m$MeanGR[counter_m] <- meangr_tl[i]
    test_results_m$SDGR[counter_m] <- sdgr_tl[i]
    test_results_m$MeanSD[counter_m] <- meansd_tl[i]
    test_results_m$MeanGRSamp[counter_m] <- samp_meangr_m_tl[[i]]
    test_results_m$SDGRSamp[counter_m] <- samp_sdgr_m_tl[[i]]
    test_results_m$MeanSDSamp[counter_m] <- samp_meansd_m_tl[[i]]
    test_results_m$MeanTSSamp[counter_m] <- samp_meants_m_tl[[i]]
    test_results_m$SampSize[counter_m] <- sampsize_tl[i]
    test_results_m$SampPercent[counter_m] <- samppercent_tl[i]
    test_results_m$SampSpecSize[counter_m] <- sampspecsize_m_tl[i]
    test_results_m$MeanTSLength[counter_m] <- meanlength_tl[i]
    test_results_m$MeanNumObs[counter_m] <- meanobs_tl[i]
    test_results_m$TotalPops[counter_m] <- tpops_tl[i]
    test_results_m$TotalSpec[counter_m] <- tspec_tl[i]
    test_results_m$PopSpec[counter_m] <- popspec_tl[i]
    test_results_m$TotalYears[counter_m] <- tyears_tl[i]
    test_results_m$TSGenVersion[counter_m] <- tsgenver_tl[i]
    test_results_m$TrendDev[counter_m] <- trenddev_m_tl[[i]]
    test_results_m$DegradeType[counter_m] <- degrade_type_tl[i]
    
    counter_m <- counter_m + 1
    
  }
  
  return(list(test_results, test_results_m))
  
}