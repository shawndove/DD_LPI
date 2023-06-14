######################
### Author: Shawn Dove
######################

# This script finds the minimum trend deviation value for LPI trends and plots it
# on a sqrt-log model of the upper 90% confidence interval of trend accuracy of the
# test data.


# load packages

library(plyr)
library(ggplot2)
library(dplyr)
library(boot)
library(simpleboot)


# set directory paths for model data, .RData, and plots
md_name <- "TestData/ModelData" # directory where data used for the accuracy model is stored
rd_name <- "RDataFiles" # directory to where .RData files are stored and should be saved
pd_name <- "Plots2" # directory where plots should be saved

# Check if directory to save plots exists. If not, create it.
if(!dir.exists(paste(pd_name, "/", sep=""))) {dir.create(paste(pd_name, "/", sep=""))} # create directory


# get IDs of testing datasets to use for testing the model
# (in the analysis script, datasets were divided into training and testing data)
model_test_df_mean <- readRDS(file=paste(rd_name, "/", "model_test_df_mean.RData", sep=""))
testIDs <- model_test_df_mean$ID

# create lists and vectors
TD_actual <- list()
TD_actual_90 <- vector()
TD_actual_mean <- vector()
TD_expected <- vector()
Percent_diff_TD <- list()
Percent_diff_TD_90 <- vector()
Percent_diff_SS <- list()
Percent_diff_SS_90 <- vector()
SS_expected <- list()
SS_expected_90 <- vector()
M_sampsize <- vector()
M_meangr <- vector()
M_sdgr <- vector()
M_tslength <- vector()
M_psd <- vector()


counter <- 1
# get TEST data and calculate values
for (i in testIDs) {
  
  M_info <- read.csv(file=paste(md_name, "/", i, "/saved_synth_", i, "_info.csv", sep="")) # copy the info csv to a variable
  
  M_sampsize[counter] <- M_info$samp_size # get sample size
  
  M_meangr[counter] <- M_info$mean_gr_raw # get mean growth rate
  
  M_sdgr[counter] <- M_info$gr_sd_raw # get standard deviation in the growth rate
  
  M_tslength[counter] <- M_info$mean_ts_length # get mean time series length
  
  M_psd[counter] <- M_info$mean_sd_raw # get mean of the standard deviation
  
  # get the actual trend deviation values for each dataset
  TD_actual[[counter]] <- read.csv(file=paste(md_name, "/", i, "/saved_synth_", i, "_trend_dev_sampled_list_lambda.csv", sep=""))$MSI
   
  # calculate expected trend deviation value for each dataset using the model built in the analysis script
  TD_expected[counter] <- exp( 
      (log(M_sampsize[counter]) * model_pops$coefficients[2])
      + (log(M_sdgr[counter]) * model_pops$coefficients[3])
      + (M_meangr[counter] * model_pops$coefficients[4])
      + (M_tslength[counter] * model_pops$coefficients[6])
      + (M_psd[counter] * model_pops$coefficients[5])
      + model_pops$coefficients[1])
 
  SSe_tmp <- vector()
  PdTD_tmp <- vector()
  PdSS_tmp <- vector()
  for (j in 1:length(TD_actual[[counter]])) {
    
    SSe_tmp[j] <- exp(
      (
        - (log(TD_actual[[counter]][j]))
           + (log(M_sdgr[counter]) * model_pops$coefficients[3])
           + (M_meangr[counter] * model_pops$coefficients[4])
           + (M_tslength[counter] * model_pops$coefficients[6])
           + (M_psd[counter] * model_pops$coefficients[5])
           + model_pops$coefficients[1])
        / (-model_pops$coefficients[2]))
    
    # squared error
    PdTD_tmp[j] <- (TD_actual[[counter]][j] - TD_expected[counter])^2
       
    # percentage error
    PdSS_tmp[j] <- (abs(SSe_tmp[j] - M_sampsize[counter]) / abs(SSe_tmp[j])) * 100

  }
  
  SS_expected[[counter]] <- SSe_tmp
  
  Percent_diff_TD[[counter]] <- PdTD_tmp
  
  Percent_diff_SS[[counter]] <- PdSS_tmp
  
  counter <- counter + 1
  
}


# save TDV data
saveRDS(TD_actual, file = paste(rd_name, "/", "TD_actual.RData", sep=""))

mean_PdTD <- mean(unlist(Percent_diff_TD), na.rm=TRUE)
mean_PdTD
mean_abs_PdTD <- mean(abs(unlist(Percent_diff_TD)), na.rm=TRUE)
mean_abs_PdTD

# Root Mean Squared Error
sqrt(mean_PdTD)

# Relative Root Mean Squared Error
RRMSEP <- mean_PdTD / sd(unlist(TD_actual), na.rm=TRUE)
RRMSEP
# save RRMSEP
saveRDS(RRMSEP, file = paste(rd_name, "/", "RRMSEP.RData", sep=""))
1 / RRMSEP

mean_PdSS <- mean(unlist(Percent_diff_SS), na.rm=TRUE)
mean_PdSS
mean_abs_PdSS <- mean(abs(unlist(Percent_diff_SS)), na.rm=TRUE)
mean_abs_PdSS


# get upper 90% bootstrap confidence interval using the BCa method to correct for skewness
# this doesn't make assumptions of a normal distribution and should work fine with beta

TD_train_actual <- list()

# get TRAINING data
counter <- 1
for (i in testIDs) {
  
  M_info <- read.csv(file=paste(md_name, "/", i, "/saved_synth_", i, "_info.csv", sep=""))
  
  M_sampsize[counter] <- M_info$samp_size
  
  TD_train_actual[[counter]] <- read.csv(file=paste(md_name, "/", i, "/saved_synth_", i, "_trend_dev_sampled_list_lambda.csv", sep=""))$MSI

  counter <- counter + 1
  
}


# calculate 90% confidence intervals
TD_bootstrap_CI <- list()

for (i in 1:length(TD_train_actual)) {

  if (all(is.na(TD_train_actual[[i]]))) {

    TD_bootstrap_CI[[i]] <- NA

  } else if (length(unique(TD_train_actual[[i]][!is.na(TD_train_actual[[i]])])) <= 1) {

    TD_bootstrap_CI[[i]] <- NA

  } else {

    nona_temp <- TD_train_actual[[i]][which(!is.na(TD_train_actual[[i]]))]
    
    bs_temp <- one.boot(nona_temp, mean, R=10000)
    
    ci_temp <- boot.ci(bs_temp, conf = 0.90, type="bca")
    
    TD_bootstrap_CI[[i]] <- ci_temp$bca[5]

  }
  
}


# save confidence intervals
saveRDS(TD_bootstrap_CI, file = paste(rd_name, "/", "TD_bootstrap_CI.RData", sep=""))


# build log-log model of upper 90% confidence interval of TDV vs sample size
TD_actual_90 <- TD_bootstrap_CI # get confidence intervals

td90_tempdf <- data.frame(td90 = unlist(TD_actual_90), ss = M_sampsize) # create data frame

saveRDS(td90_tempdf, file = paste(rd_name, "/", "td90_tempdf.RData", sep="")) # save data frame

td90_templm <- lm(sqrt(td90_tempdf$td90)~log(td90_tempdf$ss)) # make log-log model from the data

summary(td90_templm)

saveRDS(td90_templm, file = paste(rd_name, "/", "td90_templm.RData", sep="")) # save log-log model


# create data frame of grouped means
td90 <- td90_tempdf %>%
  group_by(ss) %>%
  dplyr::summarise(td90 = mean(td90, na.rm=TRUE))

td90_lm <- lm(sqrt(td90)~log(ss), td90) # make new sqrt-log model from grouped means

summary(td90_lm) # summarize model

saveRDS(td90_lm, file = paste(rd_name, "/", "td90_lm.RData", sep="")) # save model

ss90 <- data.frame(ss = seq(50, 10000, 10)) # create sequence of sample sizes


# use the sqrt-log model to predict tdv values for the sample sizes in ss90
td90_fit <- data.frame(td90=(predict(td90_lm, newdata=ss90, ss=ss90$ss))^2, ss=ss90$ss)


# calculate cut point (maximum TDV value) from the predicted values using the concordance probability method
cutpoint_df <- data.frame(sens = 1 - td90_fit$td90, spec = 1 - (td90_fit$ss/10000)) # create data frame

cutpoint_df$cz <- cutpoint_df$sens * cutpoint_df$spec # calculate sens * spec for cutpoint

cz_cutpoint <- max(cutpoint_df$cz) # calculate optimal cut point

max_tdv <- 1 - cutpoint_df$sens[cutpoint_df$cz==cz_cutpoint] # find tdv at optimal cut point

max_tdv

saveRDS(max_tdv, file = paste(rd_name, "/", "max_tdv.RData", sep="")) #save optimal cut point


# find sample size at optimal cut point
ss_at_cutpoint <- (1 - cutpoint_df$spec[cutpoint_df$cz==cz_cutpoint])*10000

ss_at_cutpoint

# plot 90% confidence interval of TDV vs sample size, with a fitted sqrt-log model
ggplot(td90_tempdf, aes(x=ss, y=td90))+
  geom_point(aes(group=ss), size=0.3)+
  geom_line(data=td90_fit, aes(x=ss, y=td90), size=0.5, color="blue")+
  scale_y_reverse()+
  geom_vline(xintercept = 1340, size=0.5, color="red")+
  ggtitle("TDV (trend deviation value) curve with optimal cut point")+
  ylab("TDV")+
  xlab("Sample Size")+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))


# save the plot
ggsave(filename=paste(pd_name, "/", "TDV_curve.tiff", sep=""),
       plot=last_plot(),
       dpi=1000,
       width=5000,
       height=3000,
       units="px",
       compression="lzw")
