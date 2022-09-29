library(tidyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)

# save the test directory name to a variable
dir_name <- "TestData/Testing2/"
dir_name2 <- "TestData/Testing2/fixed36/"

dir_names <- list.dirs(path="TestData", 
                       full.names = TRUE, 
                       recursive = FALSE)

# remove unwanted directories
dir_names <- dir_names[26890:26969]
#dir_names <- dir_names[21390:24389]
#dir_names <- dir_names[11226:15705]
#dir_names <- dir_names[11226:13465]
#dir_names <- dir_names[5696:7695]
#dir_names <- dir_names[5616:5655]
#dir_names <- dir_names[5676:5695]
#dir_names <- dir_names[797:4498]
#dir_names <- dir_names[1773:2668]
dir_names <- paste(dir_names, "/", sep="")

# get lists of csv files in the test directory
# info_list <- list.files(paste(dir_name, 
#                               sep=""), 
#                         pattern = "info",
#                         full.names = TRUE,
#                         recursive = FALSE)
# 
# gr_samples_list <- list.files(paste(dir_name, 
#                               sep=""), 
#                              pattern = "gr_samples",
#                              full.names = TRUE,
#                              recursive = FALSE)
# 
# wci_list <- list.files(paste(dir_name, 
#                              sep=""), 
#                        pattern = "within_ci_sampled_list_lambda",
#                        full.names = TRUE,
#                        recursive = FALSE)
# 
# ciw_list <- list.files(paste(dir_name, 
#                              sep=""), 
#                        pattern = "ci_width_sampled_list_lambda",
#                        full.names = TRUE,
#                        recursive = FALSE)
# 
# tsd_list <- list.files(paste(dir_name,
#                              sep=""),
#                        pattern = "culled.",
#                        full.names = TRUE,
#                        recursive = FALSE)
# 
# tscd_list <- list.files(paste(dir_name,
#                               sep=""),
#                         pattern = "completed_time_series",
#                         full.names = TRUE,
#                         recursive = FALSE)
# 
# samp_list <- list.files(paste(dir_name,
#                               sep=""),
#                         pattern = "sample_pop_id_list",
#                         full.names = TRUE,
#                         recursive = FALSE)
# 
# # copy files to new directory
# file.copy(c(info_list,
#             gr_samples_list,
#             wci_list,
#             ciw_list,
#             tsd_list,
#             tscd_list,
#             samp_list), 
#           dir_name2,
#           overwrite = TRUE)

# reset previously used vectors to empty
info_list <- vector()
gr_samples_list <- vector()
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
  
  
  gr_samples_list[i] <- list.files(paste(dir_names[i], 
                                   sep=""), 
                                   pattern = "gr_samples",
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
                            pattern = "culled.",
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
            gr_samples_list,
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
gr_samples_list <- vector()
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

gr_samples_list <- list.files(paste(dir_name2, 
                              sep=""), 
                              pattern = "gr_samples",
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
sdgrdeg_tl <- vector()
meangrdeg_tl <- vector()
meansddeg_tl <- vector()
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
  meansd_tl[i] <- info$mean_sd_raw[1]
  meangrdeg_tl[i] <- info$mean_gr_degraded[1]
  sdgrdeg_tl[i] <- info$gr_sd_degraded[1]
  meansddeg_tl[i] <- info$mean_sd_degraded[1]
  degrade_type_tl[i] <- info$degrade_type[1]
  trenddev_tl[[i]] <- trenddev$MSI
  trenddev_m_tl[[i]] <- trenddev_m[[2]]
  withinci_tl[[i]] <- withinci$MSI
  withinci_m_tl[[i]] <- withinci_m[[2]]
  ciwidth_tl[[i]] <- ciwidth$MSI
  ciwidth_m_tl[[i]] <- ciwidth_m[[2]]
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
test_results <- data.frame(matrix(NA, ncol = 25, nrow = length(info_list)*20))
test_results_m <- data.frame(matrix(NA, ncol = 25, nrow = length(info_list)))

# name columns
colnames(test_results) <- c("ID",
                            "MeanGR",
                            "SDGR",
                            "MeanSD",
                            "MeanGRDeg",
                            "SDGRDeg",
                            "MeanSDDeg",
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
                            "DegradeType",
                            "WithinCI",
                            "CIWidth")

colnames(test_results_m) <- c("ID",
                            "MeanGR",
                            "SDGR",
                            "MeanSD",
                            "MeanGRDeg",
                            "SDGRDeg",
                            "MeanSDDeg",
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
                            "DegradeType",
                            "WithinCI",
                            "CIWidth")

counter <- 1
counter_m <- 1

for (i in 1:length(info_list)) {
  
  for (j in 1:boots_tl[i]) {
    
    test_results$ID[counter] <- paste(iternum_tl[i], j, sep="_")
    test_results$MeanGR[counter] <- meangr_tl[i]
    test_results$SDGR[counter] <- sdgr_tl[i]
    test_results$MeanSD[counter] <- meansd_tl[i]
    test_results$MeanGRDeg[counter] <- meangrdeg_tl[i]
    test_results$SDGRDeg[counter] <- sdgrdeg_tl[i]
    test_results$MeanSDDeg[counter] <- meansddeg_tl[i]
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
    test_results$WithinCI[counter] <- withinci_tl[[i]][j]
    test_results$CIWidth[counter] <- ciwidth_tl[[i]][j]
    test_results$DegradeType[counter] <- degrade_type_tl[i]
    
    counter <- counter + 1
    
  }
  
  test_results_m$ID[counter_m] <- iternum_tl[i]
  test_results_m$MeanGR[counter_m] <- meangr_tl[i]
  test_results_m$SDGR[counter_m] <- sdgr_tl[i]
  test_results_m$MeanSD[counter_m] <- meansd_tl[i]
  test_results_m$MeanGRDeg[counter_m] <- meangrdeg_tl[i]
  test_results_m$SDGRDeg[counter_m] <- sdgrdeg_tl[i]
  test_results_m$MeanSDDeg[counter_m] <- meansddeg_tl[i]
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
  test_results_m$WithinCI[counter_m] <- withinci_m_tl[[i]]
  test_results_m$CIWidth[counter_m] <- ciwidth_m_tl[[i]]
  test_results_m$DegradeType[counter_m] <- degrade_type_tl[i]
  
  counter_m <- counter_m + 1
  
}

#saveRDS(test_results, file="test_results.RData")
#saveRDS(test_results_m, file="test_results_m.RData")

# remove datasets outside LPD parameter range
test_results_mc <- test_results_m[test_results_m$SDGRSamp < 0.63 
                             & test_results_m$SDGRSamp > 0.12 
                             & test_results_m$MeanSDSamp > 0.16 
                             & test_results_m$MeanSDSamp < 0.89
                             & test_results_m$MeanGRSamp > -0.19
                             & test_results_m$MeanGRSamp < 0.16
                             & test_results_m$MeanTSLength > 6.0
                             & test_results_m$MeanTSLength < 39,]

saveRDS(test_results_mc, file="test_results_mc.RData")

# randomly select data for training the model
model_train_IDs <- sample(test_results_mc$ID, size=0.67*nrow(test_results_mc), replace=FALSE)
model_train_df <- test_results_mc[test_results_mc$ID %in% model_train_IDs,]
saveRDS(model_train_df, file="model_train_df.RData")
model_train_samp_IDs_list <- list()
for (i in 1:length(model_train_IDs)) {
  model_train_samp_IDs_list[[i]] <- grep(model_train_IDs[[i]], test_results$ID)
}
model_train_samp_IDs_vec <- unlist(model_train_samp_IDs_list)
model_train_samp_df <- test_results[model_train_samp_IDs_vec,]

# use the remaining data for testing the model
model_test_IDs <- test_results_mc$ID[!test_results_mc$ID %in% model_train_IDs]
model_test_df_mean <- test_results_mc[test_results_mc$ID %in% model_test_IDs,]
model_test_IDs_list <- list()
for (i in 1:length(model_test_IDs)) {
  model_test_IDs_list[[i]] <- grep(model_test_IDs[[i]], test_results$ID)
}
model_test_IDs_vec <- unlist(model_test_IDs_list)
model_test_df <- test_results[model_test_IDs_vec,]
saveRDS(model_test_df, file="model_test_df.RData")

# check parameter ranges
range(model_test_df$MeanSDSamp, na.rm=TRUE)
range(model_train_samp_df$MeanSDSamp, na.rm=TRUE)
range(model_test_df$SDGRSamp, na.rm=TRUE)
range(model_train_samp_df$SDGRSamp, na.rm=TRUE)
range(model_test_df$MeanGRSamp, na.rm=TRUE)
range(model_train_samp_df$MeanGRSamp, na.rm=TRUE)
range(model_test_df$MeanTSSamp, na.rm=TRUE)
range(model_train_samp_df$MeanTSSamp, na.rm=TRUE)
range(model_test_df$SampSize, na.rm=TRUE)
range(model_train_samp_df$SampSize, na.rm=TRUE)

range(model_test_df$MeanSD, na.rm=TRUE)
range(model_train_samp_df$MeanSD, na.rm=TRUE)
range(model_test_df$SDGR, na.rm=TRUE)
range(model_train_samp_df$SDGR, na.rm=TRUE)
range(model_test_df$MeanGR, na.rm=TRUE)
range(model_train_samp_df$MeanGR, na.rm=TRUE)
range(model_test_df$MeanTSLength, na.rm=TRUE)
range(model_train_samp_df$MeanTSLength, na.rm=TRUE)
range(model_test_df$SampSize, na.rm=TRUE)
range(model_train_samp_df$SampSize, na.rm=TRUE)

# build the model
model_pops <- lm(log(TrendDev) 
                 ~ log(SampSize) 
                 + log(SDGR) 
                 + MeanGR 
                 + MeanSD 
                 + MeanTSLength, 
                 data = model_train_samp_df)
summary(model_pops)
saveRDS(model_pops, file="model_pops.RData")

# get beta coefficients
library(QuantPsyc)
model_pops_beta <- lm.beta(model_pops)
model_pops_beta
saveRDS(model_pops_beta, file="model_pops_beta.RData")

#model_species <- lm(I(TrendDev^(1/6))~I(SampSpecSize^(0.9))+I(SDGR^(1/4))+MeanGR+MeanTSLength+MeanSD, data=test_results)
#model_species <- lm(log(TrendDev) ~ SampSpecSize + log(SDGR) + MeanGR + MeanSD + MeanTSLength, data = test_results)
#summary(model_species)
#saveRDS(model_species, file="model_species.RData")

##############################

test_results_m$ObsErr <- NA
test_results$ObsErr <- NA
test_results_m$ObsErr <- rep(c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 
                               0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8), each=20)
test_results$ObsErr <- rep(c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 
                             0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8), each=400)
test_results_1 <- test_results[test_results$MeanGR < 1.02,]
test_results_2 <- test_results[test_results$MeanGR > 1.02,]
boxplot(log(test_results$TrendDev)~test_results$GType)
boxplot(log(test_results_m$TrendDev)~test_results_m$GType)
summary(aov(test_results_m$TrendDev~test_results_m$ObsErr))
boxplot(test_results$TrendDev~test_results$ObsErr)
summary(aov(log(test_results$TrendDev)~test_results$ObsErr))

testm <- lm(log(TrendDev)~(MeanGR^2+MeanGR), data=test_results)
rawmodel <- lm(log(TrendDev)~log(SampSize)+log(SDGR)+MeanGR+MeanTSLength+MeanSD, data=test_results)
summary(rawmodel)
hist(residuals(rawmodel), prob=TRUE, main="Histogram of residuals", xlab="Residuals", ylab="Density")
lines(density(residuals(rawmodel)), col="red", lwd=2)

rawmodel2 <- lm(log(TrendDev)~log(SampSize)+log(SDGR)+MeanGR+MeanTSLength+log(TotalPops)+SampPercent, data=test_results_m)
summary(rawmodel2)

summary(lm(log(TrendDev)~log(SampSize)*log(SDGR)*MeanGR+MeanTSLength, data=test_results_m2))
summary(lm(log(TrendDev)~log(SampSize)+log(SDGR)+MeanGR+MeanTSLength+SampSize*SDGR+SampSize*MeanGR+SDGR*MeanGR+SampSize*SDGR*MeanGR+MeanTSLength, data=test_results))
summary(nls(log(TrendDev)~(a*log(SampSize))+(b*log(SDGR))+(c*MeanGR)+(d*MeanTSLength)+(e*SampSize*SDGR)+(f*SampSize*MeanGR)+(g*SDGR*MeanGR)+k, data=test_results_m, start=list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, k=0)))

test_results_m2 <- test_results_m[test_results_m$SampSize > 20,]

trmean <- test_results_m2 %>%
  group_by(SampSize) %>%
  summarise(TrendDev = mean(TrendDev, na.rm=TRUE))

trmean_lm <- lm(log(TrendDev)~log(SampSize), trmean)
summary(trmean_lm)
nn <- data.frame(SampSize = seq(50, 5000, 50))

trmean_fit <- data.frame(TrendDev=exp(predict(trmean_lm, newdata=nn, SampSize=nn$SampSize)), SampSize=nn$SampSize)

ggplot(test_results, aes(x=SampSize, y=TrendDev))+
  geom_boxplot(aes(group=SampSize))+
  geom_line(data=trmean_fit, aes(x=SampSize, y=TrendDev))+
  scale_y_reverse()



#######

test_results_m$cl_length <- rep(c(1, 3, 6, 10, 3, 6, 10, 3, 6, 10, 1, 1, 1), each=100)
test_results$cl_length <- rep(c(1, 3, 6, 10, 3, 6, 10, 3, 6, 10, 1, 1, 1), each=2000)

test_results_m$e_ratio <- rep(c(1, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 1, 1, 1), each=100)
test_results$e_ratio <- rep(c(1, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 1, 1, 1), each=2000)

test_results_m$er_ratio <- rep(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.25, 0.5, 1), each=100)
test_results$er_ratio <- rep(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.25, 0.5, 1), each=2000)

test_results_m$group <- rep(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"), each=100)
test_results$group <- rep(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"), each=2000)

library(ggplot2)
#library(ggsignif)


ggplot(test_results_m, aes(x=group, y=TrendDev, fill=factor(group)))+
  geom_boxplot(aes())

ggsave("solutions_testing.jpg",
       plot=last_plot(),
       device=jpeg)

#######

test_results_m$group <- rep(c(200, 400, 600, 800, 1000), each=20)

ggplot(test_results_m, aes(x=group, y=TrendDev, fill=factor(group)))+
         geom_boxplot()

ggsave("samp_size_testing.jpg",
       plot=last_plot(),
       device=jpeg)

############

test_results_m$e_ratio <- rep(c(1, 0.2, 0.4, 0.6, 0.79), each=20)
test_results_m$group <- rep(c("normal", "add 200", "add 400", "add 600", "add 790"), each=20)                              

ggplot(test_results_m, aes(x=group, y=TrendDev, fill=factor(group)))+
  geom_boxplot()

ggsave("new_data_testing_10yrs.jpg",
       plot=last_plot(),
       device=jpeg)

##########################################

# build solutions plot -- data in fixed31 directory
# to load data, comment out lines of code involving mean ts length from samples

test_results_m$e_ratio <- rep(c(1, 0.2, 1, 1), each=50)
test_results_m$group <- rep(c("Control", "Group A", "Group B", "Group C"), each=50)

test_results_m$group <- factor(test_results_m$group, 
                               levels=c("Control", "Group A", "Group B", "Group C"),
                               labels=c("Control", "Group A", "Group B", "Group C"))
saveRDS(test_results_m, "comparison_data.RData")
ggplot(test_results_m, aes(x=group, y=TrendDev, fill=factor(group)))+
  geom_boxplot(show.legend=FALSE)+
  ylab("TDV")+
  xlab("")+
  theme_bw()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("solutions.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw")


##########################################

## build supplementary plots for trend length, number of pops per species,
## clustered vs random distribution of time points, size of data set, and
## amount of observation error -- data in fixed36 directory

## prepare data for trend length graph
trendlength <- test_results_m[test_results_m$ID<80200,]
#trendlength <- test_results_m[test_results_m$ID>83000,]
trendlength$TotalYears <- factor(trendlength$TotalYears)

## prepare data for populations per species graph
# popsperspec <- test_results_m[test_results_m$ID>80200 & 
#                                 test_results_m$ID<80400 & 
#                                 test_results_m$PopSpec>6 &
#                                 test_results_m$PopSpec<150,]

# popsperspec <- test_results_m[test_results_m$ID>83200 & 
#                                 test_results_m$ID<83400 &
#                                 test_results_m$PopSpec<80,]
# 
# popsperspec <- test_results_m[test_results_m$ID>83400 & 
#                                 test_results_m$ID<83600 &
#                                 test_results_m$PopSpec<80,]
# 
# popsperspec <- test_results_m[test_results_m$ID>83600 & 
#                                 test_results_m$ID<83800 &
#                                 test_results_m$PopSpec<80,]

pps_norm <- test_results_m[test_results_m$ID>83800 & 
                                test_results_m$ID<84000,]

pps_norm$SampPopSpec <- 200/pps_norm$SampSpecSize

pps_norm <- pps_norm %>%
  mutate(PopSpec = cut(PopSpec, breaks=c(3,7,12,18,25,40,75),
                       labels=c(5,10,15,20,30,50)))

pps_exp <- test_results_m[test_results_m$ID>84000 & 
                                test_results_m$ID<84200,]

pps_exp$SampPopSpec <- 200/pps_exp$SampSpecSize

pps_exp <- pps_exp %>%
  mutate(PopSpec = cut(PopSpec, breaks=c(3,7,12,18,25,40,75),
                       labels=c(5,10,15,20,30,50)))

pps_nb <- test_results_m[test_results_m$ID>84200 & 
                                test_results_m$ID<84400,]

pps_nb$SampPopSpec <- 200/pps_nb$SampSpecSize
 
pps_nb <- pps_nb %>%
  mutate(PopSpec = cut(PopSpec, breaks=c(3,7,12,18,25,40,75),
                       labels=c(5,10,15,20,30,50)))
 
#popsperspec2$SampPopSpec <- factor(popsperspec2$SampPopSpec)

#popsperspec2 <- popsperspec2 %>%
#  mutate(SampPopSpec = cut(SampPopSpec, breaks=c(1,2,2.8,3.5,5,8,15),
#                       labels=c(1.6,2.3,3.2,4.1,6.1,10)))

## prepare data for timepoint distribution graph
#tpdist <- test_results_m[test_results_m$ID>80400 & test_results_m$ID<80600,]

tpdist <- test_results_m[test_results_m$ID>84400 & test_results_m$ID<84500,]

tpdist$DegradeType <- factor(tpdist$DegradeType)

## prepare data for dataset size graph
#dssize <- test_results_m[test_results_m$ID>80600 & test_results_m$ID<80800,]
#dssize <- test_results_m[test_results_m$ID>81000 & test_results_m$ID<82000,]
dssize <- test_results_m[test_results_m$ID>82100 & test_results_m$ID<82200,]
dssize <- test_results_m[test_results_m$ID>82200 & test_results_m$ID<82300,]
dssize <- test_results_m[test_results_m$ID>82300 & test_results_m$ID<82400,]
dssize <- test_results_m[test_results_m$ID>82400 & test_results_m$ID<82500,]
dssize <- test_results_m[test_results_m$ID>82100 & test_results_m$SampSize==50,]
#dssize <- test_results_m[(test_results_m$ID>80600 & test_results_m$ID<80800) | (test_results_m$ID>81000),]
dssize$TotalPops <- rep(rep(c(50, 100, 200, 500, 1000, 2000, 5000, 10000), each=10),2)
dssize$TotalPops <- factor(dssize$TotalPops)
dssize <- dssize[dssize$TotalPops!=200,]
dssize50 <- dssize

dssize <- test_results_m[test_results_m$ID>82100 & test_results_m$SampSize==100,]
dssize$TotalPops <- factor(dssize$TotalPops)
dssize <- dssize[dssize$TotalPops!=150,]
dssize100 <- dssize

dssize <- test_results_m[test_results_m$ID>82100 & test_results_m$SampSize==500,]
dssize$TotalPops <- c(rep(c(500, 1000, 2000, 5000, 10000), each=10), 
                      rep(c(500, 1000, 2000, 5000, 10000), each=20),
                      rep(c(500, 1000, 2000, 5000, 10000), each=20))
dssize$TotalPops <- factor(dssize$TotalPops)
dssize500 <- dssize

dssize <- test_results_m[test_results_m$ID>82100 & test_results_m$SampSize==200,]
dssize$TotalPops <- c(rep(c(200, 500, 1000, 2000, 5000, 10000), each=10), 
                      rep(c(200, 500, 1000, 2000, 5000, 10000), each=15))
dssize$TotalPops <- factor(dssize$TotalPops)
dssize200 <- dssize


## prepare data for observation error graph
obserror <- test_results_m[test_results_m$ID>80800 & test_results_m$ID<81000,]

obserror$ObsError <- rep(c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2), each=20)

obserror$ObsError <- factor(obserror$ObsError)

## trend length

ggplot(trendlength, aes(x=TotalYears, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
  ylim(c(0,1))+
  ylab("Trend Deviation Value")+
  xlab("Number of Years Modelled")+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("yearsmodelled.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw")

## pops per species

# uniform distribution
pps_norm_plot <- ggplot(pps_norm, aes(x=PopSpec, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
  #scale_x_discrete(labels=c("5", "10", "15", "20", "30", "50", "100", "200"))+
  ylim(c(0,1))+
  ylab("Trend Deviation Value")+
  xlab("Mean Populations Per Species")+
  ggtitle("Normal")+
  theme_bw()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("popsperspec_normal.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw")

# discretized exponential distribution
pps_exp_plot <- ggplot(pps_exp, aes(x=PopSpec, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
  #scale_x_discrete(labels=c("5", "10", "15", "20", "30", "50", "100", "200"))+
  ylim(c(0,1))+
  ylab("Trend Deviation Value")+
  xlab("Mean Populations Per Species")+
  ggtitle("Discretized Exponential")+
  theme_bw()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("popsperspec_exponential.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw")

# zero-truncated negative binomial distribution
pps_nb_plot <- ggplot(pps_nb, aes(x=PopSpec, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
  #scale_x_discrete(labels=c("5", "10", "15", "20", "30", "50", "100", "200"))+
  ylim(c(0,1))+
  ylab("Trend Deviation Value")+
  xlab("Mean Populations Per Species")+
  ggtitle("Negative Binomial")+
  theme_bw()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("popsperspec_negbin.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw")

ppsplot <- ggarrange(pps_norm_plot + rremove("ylab") + rremove("xlab"),
                    pps_exp_plot + rremove("ylab") + rremove("xlab"),
                    pps_nb_plot + rremove("ylab") + rremove("xlab"),
                    ncol=3)

annotate_figure(ppsplot,
                left = textGrob("Trend Deviation Value",
                                rot = 90, vjust = 0.5,
                                gp = gpar(cex = 1.3)),
                bottom = textGrob("Mean Populations Per Species",
                                  hjust = 0.5,
                                  gp = gpar(cex = 1.3)))

ggsave("popsperspec_all.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw",
       width = 9000,
       height = 3000,
       units = "px")

## time point distribution

ggplot(tpdist, aes(x=DegradeType, y=TrendDev, fill=factor(DegradeType)))+
  geom_boxplot(show.legend=FALSE)+
  scale_x_discrete(labels=c("Endpoint Clustering", "Random"))+
  scale_fill_manual(values=c("skyblue", "orange"))+
  ylim(c(0,0.4))+
  ylab("Trend Deviation Value")+
  xlab("Time Series Distribution")+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("timepointdistribution.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw",
       width = 9000,
       height = 4000,
       units = "px")

## dataset size
ds50 <- ggplot(dssize50, aes(x=TotalPops, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
  scale_y_continuous(limits=c(0,1))+
  ylab("TDV")+
  xlab("Total Populations")+
  ggtitle("Sample Size: 50")+
  theme_bw()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("totalpops50.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw")

ds100 <- ggplot(dssize100, aes(x=TotalPops, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
  scale_y_continuous(limits=c(0,1))+
  ylab("TDV")+
  xlab("Total Populations")+
  ggtitle("Sample Size: 100")+
  theme_bw()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("totalpops100.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw")

ds200 <- ggplot(dssize200, aes(x=TotalPops, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
  scale_y_continuous(limits=c(0,1))+
  ylab("TDV")+
  xlab("Total Populations")+
  ggtitle("Sample Size: 200")+
  theme_bw()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("totalpops200.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw")

ds500 <- ggplot(dssize500, aes(x=TotalPops, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
  scale_y_continuous(limits=c(0,1))+
  ylab("TDV")+
  xlab("Total Populations")+
  ggtitle("Sample Size: 500")+
  theme_bw()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("totalpops500.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw")

dsplot <- ggarrange(ds50 + rremove("ylab") + rremove("xlab"),
                    ds100 + rremove("ylab") + rremove("xlab"),
                    ds200 + rremove("ylab") + rremove("xlab"),
                    ds500 + rremove("ylab") + rremove("xlab"),
                    ncol=2, nrow=2)

annotate_figure(dsplot,
                left = textGrob("Trend Deviation Value",
                                rot = 90, vjust = 0.5,
                                gp = gpar(cex = 1.3)),
                bottom = textGrob("Total Populations",
                                  hjust = 0.5,
                                  gp = gpar(cex = 1.3)))

ggsave("totalpopsall.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw")

## observation error

ggplot(obserror, aes(x=ObsError, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
  scale_x_discrete(labels=c("0%", "10%", "20%", "40%", "60%", "80%", "100%", "200%"))+
  ylim(c(0,0.6))+
  ylab("Trend Deviation Value")+
  xlab("Percent Observation Error")+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("observationerror.tiff",
       device = tiff,
       dpi = 1000,
       compression = "lzw",
       width = 10000,
       height = 5000,
       units = "px")
