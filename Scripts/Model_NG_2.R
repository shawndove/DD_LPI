library(plyr)
library(ggplot2)
library(dplyr)
library(boot)
library(simpleboot)

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
for (i in 70001:73000) {
  
  if (!file.exists(paste("TestData/", i, "/saved_synth_", i, "_info.csv", sep=""))) {
    
    next
    
  }
  
  if (!i %in% model_test_df_mean$ID) {
    
    next
    
  }
  
  M_info <- read.csv(file=paste("TestData/", i, "/saved_synth_", i, "_info.csv", sep=""))
  
  M_sampsize[counter] <- M_info$samp_size
  
  M_meangr[counter] <- M_info$mean_gr_raw
  
  M_sdgr[counter] <- M_info$gr_sd_raw
  
  M_tslength[counter] <- M_info$mean_ts_length
  
  M_psd[counter] <- M_info$mean_sd_raw
  
#  TD_actual_mean[counter] <- read.csv(file=paste("TestData/Testing2/fixed4/saved_synth_", i, "_trend_dev_sampled_mean_lambda.csv", sep=""))[,2]

  TD_actual[[counter]] <- read.csv(file=paste("TestData/Testing2/fixed35/saved_synth_", i, "_trend_dev_sampled_list_lambda.csv", sep=""))$MSI
   
 # TD_actual_90[counter] <- max(TD_actual[[counter]][TD_actual[[counter]]!=max(TD_actual[[counter]], na.rm=TRUE)], na.rm=TRUE)
   
  #TD_expected[counter] <- exp((log(M_sampsize[counter]) * -0.951) + (log(M_sdgr[counter]) * 1.857) + (M_meangr[counter] * 4.609) + (M_tslength[counter] * -0.0437) + 3.518)
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
    
    # percent difference
   # PdTD_tmp[j] <- ((TD_actual[[counter]][j] - TD_expected[counter]) / ((TD_actual[[counter]][j] + TD_expected[counter]) / 2)) * 100
    
    # percentage error
   # PdTD_tmp[j] <- (abs(TD_actual[[counter]][j] - TD_expected[counter]) / abs(TD_actual[[counter]][j])) * 100
    
    # squared error
    PdTD_tmp[j] <- (TD_actual[[counter]][j] - TD_expected[counter])^2
       
    # percent difference
   # PdSS_tmp[j] <- ((M_sampsize[counter] - SSe_tmp[j]) / ((M_sampsize[counter] + SSe_tmp[j]) / 2)) * 100
   
    # percentage error
    PdSS_tmp[j] <- (abs(SSe_tmp[j] - M_sampsize[counter]) / abs(SSe_tmp[j])) * 100
    
    # squared error
    #PdSS_tmp[j] <- (SSe_tmp[j] - M_sampsize[counter])^2
    
  }
  
  SS_expected[[counter]] <- SSe_tmp
  
  Percent_diff_TD[[counter]] <- PdTD_tmp
  
  Percent_diff_SS[[counter]] <- PdSS_tmp
  
  # SS_expected_90[counter] <- exp(
  #     (
  #       - (log(TD_actual_90[counter]))
  #       + (log(M_sdgr[counter]) * model_pops$coefficients[3])
  #       + (M_meangr[counter] * model_pops$coefficients[4])
  #       + (M_tslength[counter] * model_pops$coefficients[6])
  #       + (M_psd[counter] * model_pops$coefficients[5])
  #       + model_pops$coefficients[1])
  #     / (-model_pops$coefficients[2]))
  
 # Percent_diff_TD_90[counter] <- ((TD_actual_90[counter] - TD_expected[counter]) / ((TD_actual_90[counter] + TD_expected[counter]) / 2)) * 100
  
#  Percent_diff_SS_90[counter] <- ((M_sampsize[counter] - SS_expected_90[counter]) / ((M_sampsize[counter] + SS_expected_90[counter]) / 2)) * 100 
  
  counter <- counter + 1
  
}
# save TDV data
saveRDS(TD_actual, file="TD_actual.RData")

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
saveRDS(RRMSEP, file = "RRMSEP.RData")
1 / RRMSEP

mean_PdSS <- mean(unlist(Percent_diff_SS), na.rm=TRUE)
mean_PdSS
mean_abs_PdSS <- mean(abs(unlist(Percent_diff_SS)), na.rm=TRUE)
mean_abs_PdSS

# get upper 90% bootstrap confidence interval using the BCa method to correct for skewness
# this doesn't make assumptions of a normal distribution and should work fine with beta
# first get TRAINING data and use that to get confidence intervals
counter <- 1
TD_train_actual <- list()
for (i in 70001:73000) {
  if (!file.exists(paste("TestData/", i, "/saved_synth_", i, "_info.csv", sep=""))) {
    next
  }
  if (!i %in% model_train_df$ID) {
    next
  }
  M_info <- read.csv(file=paste("TestData/", i, "/saved_synth_", i, "_info.csv", sep=""))
  M_sampsize[counter] <- M_info$samp_size
  TD_train_actual[[counter]] <- read.csv(file=paste("TestData/Testing2/fixed35/saved_synth_", i, "_trend_dev_sampled_list_lambda.csv", sep=""))$MSI
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
saveRDS(TD_bootstrap_CI, file = "TD_bootstrap_CI.RData")

# build log-log model of upper 90% confidence interval of TDV vs sample size
TD_actual_90 <- TD_bootstrap_CI # get confidence intervals
td90_tempdf <- data.frame(td90 = unlist(TD_actual_90), ss = M_sampsize) # create data frame
# save data frame
save(td90_tempdf, file="td90_tempdf.RData")
td90_templm <- lm(sqrt(td90_tempdf$td90)~log(td90_tempdf$ss)) # make log-log model from the data
summary(td90_templm)
# save log-log model
save(td90_templm, file="td90_templm.RData")

# create data frame of grouped means
td90 <- td90_tempdf %>%
  group_by(ss) %>%
  dplyr::summarise(td90 = mean(td90, na.rm=TRUE))
td90_lm <- lm(sqrt(td90)~log(ss), td90) # make new log-log model from grouped means
summary(td90_lm) # summarize model
# save model
saveRDS(td90_lm, file="td90_lm.RData")
ss90 <- data.frame(ss = seq(50, 10000, 10)) # create sequence of sample sizes

# use the log-log grouped means model to predict tdv values for the sample sizes in ss90
td90_fit <- data.frame(td90=(predict(td90_lm, newdata=ss90, ss=ss90$ss))^2, ss=ss90$ss)

# calculate cut point (maximum TDV value) from the predicted values using the concordance probability method
cutpoint_df <- data.frame(sens = 1 - td90_fit$td90, spec = 1 - (td90_fit$ss/10000)) # create data frame
cutpoint_df$cz <- cutpoint_df$sens * cutpoint_df$spec # calculate sens * spec for cutpoint
cz_cutpoint <- max(cutpoint_df$cz) # calculate optimal cut point
max_tdv <- 1 - cutpoint_df$sens[cutpoint_df$cz==cz_cutpoint] # find tdv at optimal cut point
max_tdv
#save optimal cut point
saveRDS(max_tdv, file="max_tdv.RData")
# find sample size at optimal cut point
ss_at_cutpoint <- (1 - cutpoint_df$spec[cutpoint_df$cz==cz_cutpoint])*10000
ss_at_cutpoint

# plot 90% confidence interval of TDV vs sample size, with a fitted linear log-log model
ggplot(td90_tempdf, aes(x=ss, y=td90))+
  #geom_point(size=2)+
  geom_point(aes(group=ss), size=0.3)+
  geom_line(data=td90_fit, aes(x=ss, y=td90), size=0.5, color="blue")+
  scale_y_reverse()+
  geom_vline(xintercept = 1340, size=0.5, color="red")+
  #geom_hline(yintercept = 0.0906, size=0.5, color="black")+
  #geom_point(aes(x=1500, y=0.045), size=2, shape=1, stroke=1, color="black")+
  #annotate("segment", x = 0, xend = 3000, y = 0.045, yend = 0.002, color="black", size = 0.5)+
  #annotate("point", x = 520, y = 0.058, color="red", size=4)+
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
ggsave(filename="TDV_curve3.tiff",
       plot=last_plot(),
       dpi=1000,
       height=5000,
       width=8000,
       units="px",
       compression="lzw")
