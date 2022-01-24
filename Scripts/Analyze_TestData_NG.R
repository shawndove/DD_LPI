
# save the test directory name to a variable
dir_name <- "TestData/Testing2/"
dir_name2 <- "TestData/Testing2/fixed3/"

dir_names <- list.dirs(path="TestData", 
                       full.names = TRUE, 
                       recursive = FALSE)

# remove unwanted directories
dir_names <- dir_names[4283:4498]
#dir_names <- dir_names[1773:2668]
dir_names <- paste(dir_names, "/", sep="")

# get lists of csv files in the test directory
info_list <- list.files(paste(dir_name, 
                              sep=""), 
                        pattern = "info",
                        full.names = TRUE,
                        recursive = FALSE)

wci_list <- list.files(paste(dir_name, 
                             sep=""), 
                       pattern = "within_ci_sampled_list_lambda",
                       full.names = TRUE,
                       recursive = FALSE)

ciw_list <- list.files(paste(dir_name, 
                             sep=""), 
                       pattern = "ci_width_sampled_list_lambda",
                       full.names = TRUE,
                       recursive = FALSE)

tsd_list <- list.files(paste(dir_name,
                             sep=""),
                       pattern = "culled",
                       full.names = TRUE,
                       recursive = FALSE)

tscd_list <- list.files(paste(dir_name,
                              sep=""),
                        pattern = "completed_time_series",
                        full.names = TRUE,
                        recursive = FALSE)

samp_list <- list.files(paste(dir_name,
                              sep=""),
                        pattern = "sample_pop_id_list",
                        full.names = TRUE,
                        recursive = FALSE)

# copy files to new directory
file.copy(c(info_list,
            wci_list,
            ciw_list,
            tsd_list,
            tscd_list,
            samp_list), 
          dir_name2,
          overwrite = TRUE)

# reset previously used vectors to empty
info_list <- vector()
wci_list <- vector()
wcim_list <- vector()
ciw_list <- vector()
ciwm_list <- vector()
tsd_list <- vector()
tscd_list <- vector()
samp_list <- vector()

for (i in 1:length(dir_names)) {
  
  info_list[i] <- list.files(paste(dir_names[i], 
                                sep=""), 
                          pattern = "info",
                          full.names = TRUE,
                          recursive = FALSE)
  
  wcim_list[i] <- list.files(paste(dir_names[i], 
                               sep=""), 
                         pattern = "within_ci_sampled_mean_lambda",
                         full.names = TRUE,
                         recursive = FALSE)
  
  wci_list[i] <- list.files(paste(dir_names[i], 
                                  sep=""), 
                            pattern = "within_ci_sampled_list_lambda",
                            full.names = TRUE,
                            recursive = FALSE)
  
  ciwm_list[i] <- list.files(paste(dir_names[i], 
                               sep=""), 
                         pattern = "ci_width_sampled_mean_lambda",
                         full.names = TRUE,
                         recursive = FALSE)
  
  ciw_list[i] <- list.files(paste(dir_names[i], 
                                  sep=""), 
                            pattern = "ci_width_sampled_list_lambda",
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

# copy files to new directory
file.copy(c(info_list,
            wci_list,
            wcim_list,
            ciw_list,
            ciwm_list,
            tsd_list,
            tscd_list,
            samp_list), 
          dir_name2,
          overwrite = TRUE)

# reset previously used vectors to empty again
info_list <- vector()
td_list <- vector()
tdm_list <- vector()
wci_list <- vector()
wcim_list <- vector()
ciw_list <- vector()
ciwm_list <- vector()
tsd_list <- vector()
tscd_list <- vector()
samp_list <- vector()

# get lists of csv files in the test directory
info_list <- list.files(paste(dir_name2, 
                              sep=""), 
                        pattern = "info",
                        full.names = TRUE,
                        recursive = FALSE)

tdm_list <- list.files(paste(dir_name2, 
                            sep=""), 
                      pattern = "trend_dev_sampled_mean_lambda",
                      full.names = TRUE,
                      recursive = FALSE)

td_list <- list.files(paste(dir_name2, 
                            sep=""), 
                      pattern = "trend_dev_sampled_list_lambda",
                      full.names = TRUE,
                      recursive = FALSE)

wcim_list <- list.files(paste(dir_name2, 
                             sep=""), 
                       pattern = "within_ci_sampled_mean_lambda",
                       full.names = TRUE,
                       recursive = FALSE)

wci_list <- list.files(paste(dir_name2, 
                             sep=""), 
                       pattern = "within_ci_sampled_list_lambda",
                       full.names = TRUE,
                       recursive = FALSE)

ciwm_list <- list.files(paste(dir_name2, 
                             sep=""), 
                       pattern = "ci_width_sampled_mean_lambda",
                       full.names = TRUE,
                       recursive = FALSE)

ciw_list <- list.files(paste(dir_name2, 
                             sep=""), 
                       pattern = "ci_width_sampled_list_lambda",
                       full.names = TRUE,
                       recursive = FALSE)

tsd_list <- list.files(paste(dir_name2,
                             sep=""),
                       pattern = "culled",
                       full.names = TRUE,
                       recursive = FALSE)

tscd_list <- list.files(paste(dir_name2,
                             sep=""),
                       pattern = "completed_time_series",
                       full.names = TRUE,
                       recursive = FALSE)

samp_list <- list.files(paste(dir_name2,
                              sep=""),
                        pattern = "sample_pop_id_list",
                        full.names = TRUE,
                        recursive = FALSE)

# create a set of temporary lists and vectors to hold various data from saved files
trenddev_tl <- list()
trenddev_m_tl <- list()
withinci_tl <- list()
withinci_m_tl <- list()
ciwidth_tl <- list()
ciwidth_m_tl <- list()
tyears_tl <- vector()
tpops_tl <- vector()
boots_tl <- vector()
sampsize_tl <- vector()
popspec_tl <- vector()
sdgr_tl <- vector()
meangr_tl <- vector()
iternum_tl <- vector()
tsgenver_tl <- vector()
tspec_tl <- vector()
meanobs_tl <- vector()
meanlength_tl <- vector()

# get data from files
for (i in 1:length(info_list)) {
  
  # import data from csv files
  info <- read.csv(info_list[i])
  trenddev <- read.csv(td_list[i])
  trenddev_m <- read.csv(tdm_list[i])
  withinci <- read.csv(wci_list[i])
  withinci_m <- read.csv(wcim_list[i])
  ciwidth <- read.csv(ciw_list[i])
  ciwidth_m <- read.csv(ciwm_list[i])
  tsdata <- readRDS(tsd_list[i])
  tscdata <- readRDS(tscd_list[i])
  sampdata <- readRDS(samp_list[i])
  
  # read data into lists
  iternum_tl[i] <- info$ID[1]
  tsgenver_tl[i] <- info$synth_version[1]
  tyears_tl[i] <- info$num_years[1]
  meangr_tl[i] <- info$mean_gr_raw[1]
  sdgr_tl[i] <- info$gr_sd_raw[1]
  trenddev_tl[[i]] <- trenddev$MSI
  trenddev_m_tl[[i]] <- trenddev_m[[2]]
  withinci_tl[[i]] <- withinci$MSI
  withinci_m_tl[[i]] <- withinci_m[[2]]
  ciwidth_tl[[i]] <- ciwidth$MSI
  ciwidth_m_tl[[i]] <- ciwidth_m[[2]]
  boots_tl[i] <- length(trenddev_tl[[i]])
  tpops_tl[i] <- nrow(tsdata)
  tspec_tl[i] <- length(unique(tsdata$SpecID))
  popspec_tl[i] <- round((tpops_tl[i] / tspec_tl[i]), digits = 1)
  sampsize_tl[i] <- length(sampdata[[1]])
  meanobs_tl[i] <- sum(!is.na(as.vector(tsdata[,1:tyears_tl[i]]))) / tpops_tl[i]
 # meanlength_tl[i] <- sum(!is.na(as.vector(tscdata[,1:tyears_tl[i]]))) / tpops_tl[i]
  meanlength_tl[i] <- info$mean_ts_length[1]

}

# create data frame to hold results
test_results <- data.frame(matrix(NA, ncol = 14, nrow = length(info_list)*20))
test_results_m <- data.frame(matrix(NA, ncol = 14, nrow = length(info_list)))

# name columns
colnames(test_results) <- c("ID",
                            "MeanGR",
                            "SDGR",
                            "SampSize",
                            "MeanTSLength",
                            "MeanNumObs",
                            "TotalPops",
                            "TotalSpec",
                            "PopSpec",
                            "TotalYears",
                            "TSGenVersion",
                            "TrendDev",
                            "WithinCI",
                            "CIWidth")

colnames(test_results_m) <- c("ID",
                            "MeanGR",
                            "SDGR",
                            "SampSize",
                            "MeanTSLength",
                            "MeanNumObs",
                            "TotalPops",
                            "TotalSpec",
                            "PopSpec",
                            "TotalYears",
                            "TSGenVersion",
                            "TrendDev",
                            "WithinCI",
                            "CIWidth")

counter <- 1
counter_m <- 1

for (i in 1:length(info_list)) {
  
  for (j in 1:boots_tl[i]) {
    
    test_results$ID[counter] <- paste(iternum_tl[i], j, sep="_")
    test_results$MeanGR[counter] <- meangr_tl[i]
    test_results$SDGR[counter] <- sdgr_tl[i]
    test_results$SampSize[counter] <- sampsize_tl[i]
    test_results$MeanTSLength[counter] <- meanlength_tl[i]
    test_results$MeanNumObs[counter] <- meanobs_tl[i]
    test_results$TotalPops[counter] <- tpops_tl[i]
    test_results$TotalSpec[counter] <- tspec_tl[i]
    test_results$PopSpec[counter] <- popspec_tl[i]
    test_results$TotalYears[counter] <- tyears_tl[i]
    test_results$TSGenVersion[counter] <- tsgenver_tl[i]
    test_results$TrendDev[counter] <- trenddev_tl[[i]][j]
    test_results$WithinCI[counter] <- withinci_tl[[i]][j]
    test_results$CIWidth[counter] <- ciwidth_tl[[i]][j]
    
    counter <- counter + 1
    
  }
  
  test_results_m$ID[counter_m] <- iternum_tl[i]
  test_results_m$MeanGR[counter_m] <- meangr_tl[i]
  test_results_m$SDGR[counter_m] <- sdgr_tl[i]
  test_results_m$SampSize[counter_m] <- sampsize_tl[i]
  test_results_m$MeanTSLength[counter_m] <- meanlength_tl[i]
  test_results_m$MeanNumObs[counter_m] <- meanobs_tl[i]
  test_results_m$TotalPops[counter_m] <- tpops_tl[i]
  test_results_m$TotalSpec[counter_m] <- tspec_tl[i]
  test_results_m$PopSpec[counter_m] <- popspec_tl[i]
  test_results_m$TotalYears[counter_m] <- tyears_tl[i]
  test_results_m$TSGenVersion[counter_m] <- tsgenver_tl[i]
  test_results_m$TrendDev[counter_m] <- trenddev_m_tl[[i]]
  test_results_m$WithinCI[counter_m] <- withinci_m_tl[[i]]
  test_results_m$CIWidth[counter_m] <- ciwidth_m_tl[[i]]
  
  counter_m <- counter_m + 1
  
}
testm <- lm(log(TrendDev)~(MeanGR^2+MeanGR), data=test_results)
summary(lm(log(TrendDev)~log(SampSize)+log(SDGR)+MeanGR+MeanTSLength, data=test_results))
summary(lm(log(TrendDev)~log(SampSize)*log(SDGR)*MeanGR+MeanTSLength, data=test_results_m))
summary(lm(log(TrendDev)~log(SampSize)+log(SDGR)+MeanGR+MeanTSLength+SampSize*SDGR+SampSize*MeanGR+SDGR*MeanGR+SampSize*SDGR*MeanGR+MeanTSLength, data=test_results))
summary(nls(log(TrendDev)~(a*log(SampSize))+(b*log(SDGR))+(c*MeanGR)+(d*MeanTSLength)+(e*SampSize*SDGR)+(f*SampSize*MeanGR)+(g*SDGR*MeanGR)+k, data=test_results_m, start=list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, k=0)))

test_results_m2 <- test_results_m2[test_results_m2$SampSize > 20,]

trmean <- test_results_m2 %>%
  group_by(SampSize) %>%
  summarise(TrendDev = mean(TrendDev, na.rm=TRUE))

trmean_lm <- lm(log(TrendDev)~log(SampSize), trmean)
summary(trmean_lm)
nn <- data.frame(SampSize = seq(50, 5000, 50))

trmean_fit <- data.frame(TrendDev=exp(predict(trmean_lm, newdata=nn, SampSize=nn$SampSize)), SampSize=nn$SampSize)

ggplot(test_results, aes(x=SampSize, y=TrendDev))+
  geom_point(size=2)+
  geom_line(data=trmean_fit, aes(x=SampSize, y=TrendDev))+
  scale_y_reverse()

