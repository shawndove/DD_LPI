library(plyr)
library(ggplot2)
library(dplyr)

TD_actual <- list()
TD_actual_95 <- vector()
TD_actual_mean <- vector()
TD_expected <- vector()
Percent_diff_TD <- list()
Percent_diff_TD_95 <- vector()
Percent_diff_SS <- list()
Percent_diff_SS_95 <- vector()
SS_expected <- list()
SS_expected_95 <- vector()
M_sampsize <- vector()
M_meangr <- vector()
M_sdgr <- vector()
M_tslength <- vector()
counter <- 1
for (i in 21001:21549) {
  
  M_info <- read.csv(file=paste("TestData/", i, "/saved_synth_", i, "_info.csv", sep=""))
  
  M_sampsize[counter] <- M_info$samp_size
  
  M_meangr[counter] <- M_info$mean_gr_raw
  
  M_sdgr[counter] <- M_info$gr_sd_raw
  
  M_tslength[counter] <- M_info$mean_ts_length
  
#  TD_actual_mean[counter] <- read.csv(file=paste("TestData/Testing2/fixed4/saved_synth_", i, "_trend_dev_sampled_mean_lambda.csv", sep=""))[,2]

  TD_actual[[counter]] <- read.csv(file=paste("TestData/Testing2/fixed4/saved_synth_", i, "_trend_dev_sampled_list_lambda.csv", sep=""))$MSI
   
  TD_actual_95[counter] <- max(TD_actual[[counter]][TD_actual[[counter]]!=max(TD_actual[[counter]], na.rm=TRUE)], na.rm=TRUE)
   
  TD_expected[counter] <- exp((log(M_sampsize[counter]) * -0.951) + (log(M_sdgr[counter]) * 1.857) + (M_meangr[counter] * 4.609) + (M_tslength[counter] * -0.0437) + 3.518)
 
  SSe_tmp <- vector()
  PdTD_tmp <- vector()
  PdSS_tmp <- vector()
  for (j in 1:length(TD_actual[[counter]])) {
    
    SSe_tmp[j] <- exp(((log(M_sdgr[counter]) * 1.857) - log(TD_actual[[counter]][j]) + (M_meangr[counter] * 4.609) + (M_tslength[counter] * -0.0437) + 3.518) / 0.951)
    
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
    PdSS_tmp[j] <- (SSe_tmp[j] - M_sampsize[counter])^2
    
  }
  
  SS_expected[[counter]] <- SSe_tmp
  
  Percent_diff_TD[[counter]] <- PdTD_tmp
  
  Percent_diff_SS[[counter]] <- PdSS_tmp
  
  SS_expected_95[counter] <- exp(((log(M_sdgr[counter]) * 1.857) - log(TD_actual_95[counter]) + (M_meangr[counter] * 4.609) + (M_tslength[counter] * -0.0437) + 3.518) / 0.951)
  
  Percent_diff_TD_95[counter] <- ((TD_actual_95[counter] - TD_expected[counter]) / ((TD_actual_95[counter] + TD_expected[counter]) / 2)) * 100
  
  Percent_diff_SS_95[counter] <- ((M_sampsize[counter] - SS_expected_95[counter]) / ((M_sampsize[counter] + SS_expected_95[counter]) / 2)) * 100 
  
  counter <- counter + 1
  
}

mean_PdTD <- mean(unlist(Percent_diff_TD), na.rm=TRUE)
mean_PdTD
mean_abs_PdTD <- mean(abs(unlist(Percent_diff_TD)), na.rm=TRUE)
mean_abs_PdTD

# Root Mean Squared Error
sqrt(mean_PdTD)

# Relative Root Mean Squared Error
RRMSEP <- mean_PdTD / sd(unlist(TD_actual), na.rm=TRUE)
RRMSEP
# 1 / RRMSEP
1 / RRMSEP

mean_PdSS <- mean(unlist(Percent_diff_SS), na.rm=TRUE)
mean_PdSS
mean_abs_PdSS <- mean(abs(unlist(Percent_diff_SS)), na.rm=TRUE)
mean_abs_PdSS

# get upper 90% bootstrap confidence interval using the BCa method to correct for skewness
# this doesn't make assumptions of a normal distribution and should work fine with beta
TD_bootstrap_CI <- list()
for (i in 1:length(TD_actual)) {
  
  nona_temp <- TD_actual[[i]][which(!is.na(TD_actual[[i]]))]
  bs_temp <- one.boot(nona_temp, mean, R=10000)
  ci_temp <- boot.ci(bs_temp, conf = 0.9, type="bca")
  TD_bootstrap_CI[[i]] <- ci_temp$bca[5]
  
}

TD_actual_95 <- lapply(TD_actual, function(i) {
  
  temp <- i[i!=max(i, na.rm=TRUE)]
  
 # temp2 <- temp[temp!=max(temp, na.rm=TRUE)]
  
  actual_95 <- max(temp, na.rm=TRUE)
  
  return(actual_95)
  
})

TD_actual_95 <- TD_bootstrap_CI

td95_tempdf <- data.frame(td95 = unlist(TD_actual_95), ss = M_sampsize)

td95_templm <- lm(log(td95_tempdf$td95)~log(td95_tempdf$ss))

td95 <- td95_tempdf %>%
  group_by(ss) %>%
  summarise(td95 = mean(td95, na.rm=TRUE))

td95_lm <- lm(log(td95)~log(ss), td95)
summary(td95_lm)
ss95 <- data.frame(ss = seq(50, 5000, 10))

td95_fit <- data.frame(td95=exp(predict(td95_lm, newdata=ss95, ss=ss95$ss)), ss=ss95$ss)

ggplot(td95_tempdf, aes(x=ss, y=td95))+
  #geom_point(size=2)+
  geom_boxplot(aes(group=ss), size=0.3, outlier.size=0.7)+
  geom_line(data=td95_fit, aes(x=ss, y=td95), size=0.5, color="blue")+
  scale_y_reverse()+
  geom_vline(xintercept = 520, size=0.5, color="red")+
  #geom_hline(yintercept = 0.058, size=0.5, color="black")+
  #geom_point(aes(x=1500, y=0.045), size=2, shape=1, stroke=1, color="black")+
  #annotate("segment", x = 0, xend = 3000, y = 0.045, yend = 0.002, color="black", size = 0.5)+
  #annotate("point", x = 520, y = 0.058, color="red", size=4)+
  ggtitle("TDV (trend deviation value) curve")+
  ylab("TDV")+
  xlab("Sample Size")+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))



ggsave(filename="test4.tiff",
       plot=last_plot(),
       dpi=1000,
       height=5000,
       width=8000,
       units="px",
       compression="lzw")
