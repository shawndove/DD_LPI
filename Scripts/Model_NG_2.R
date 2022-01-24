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
for (i in 19001:19550) {
  
  M_info <- read.csv(file=paste("TestData/", i, "/saved_synth_", i, "_info.csv", sep=""))
  
  M_sampsize[counter] <- M_info$samp_size
  
  M_meangr[counter] <- M_info$mean_gr_raw
  
  M_sdgr[counter] <- M_info$gr_sd_raw
  
  M_tslength[counter] <- M_info$mean_ts_length
  
#  TD_actual_mean[counter] <- read.csv(file=paste("TestData/Testing2/fixed4/saved_synth_", i, "_trend_dev_sampled_mean_lambda.csv", sep=""))[,2]

  TD_actual[[counter]] <- read.csv(file=paste("TestData/Testing2/fixed4/saved_synth_", i, "_trend_dev_sampled_list_lambda.csv", sep=""))$MSI
   
  TD_actual_95[counter] <- max(TD_actual[[counter]][TD_actual[[counter]]!=max(TD_actual[[counter]], na.rm=TRUE)], na.rm=TRUE)
   
  TD_expected[counter] <- exp((log(M_sampsize[counter]) * -0.951) + (log(M_sdgr[counter]) * 1.857) + (M_meangr[counter] * 4.609) + (M_tslength[counter] * -0.0437) + 3.518)
 
#  SSe_tmp <- vector()
#  PdTD_tmp <- vector()
#  PdSS_tmp <- vector()
#  for (j in 1:length(TD_actual[[counter]])) {
#    
#    SSe_tmp[j] <- exp(((log(M_sdgr[counter]) * 1.857) - log(TD_actual[[counter]][j]) + (M_meangr[counter] * 4.609) + (M_tslength[counter] * -0.0437) + 3.518) / 0.951)
#    
#    PdTD_tmp[j] <- ((TD_actual[[counter]][j] - TD_expected[counter]) / ((TD_actual[[counter]][j] + TD_expected[counter]) / 2)) * 100
#    
#    PdSS_tmp[j] <- ((M_sampsize[counter] - SSe_tmp[j]) / ((M_sampsize[counter] + SSe_tmp[j]) / 2)) * 100
#    
#  }
#  
#  SS_expected[[counter]] <- SSe_tmp
#  
#  Percent_diff_TD[[counter]] <- PdTD_tmp
#  
#  Percent_diff_SS[[counter]] <- PdSS_tmp
  
  SS_expected_95[counter] <- exp(((log(M_sdgr[counter]) * 1.857) - log(TD_actual_95[counter]) + (M_meangr[counter] * 4.609) + (M_tslength[counter] * -0.0437) + 3.518) / 0.951)
  
  Percent_diff_TD_95[counter] <- ((TD_actual_95[counter] - TD_expected[counter]) / ((TD_actual_95[counter] + TD_expected[counter]) / 2)) * 100
  
  Percent_diff_SS_95[counter] <- ((M_sampsize[counter] - SS_expected_95[counter]) / ((M_sampsize[counter] + SS_expected_95[counter]) / 2)) * 100 
  
  counter <- counter + 1
  
}

mean_PdTD <- mean(unlist(Percent_diff_TD), na.rm=TRUE)
mean_PdTD
mean_abs_PdTD <- mean(abs(unlist(Percent_diff_TD)), na.rm=TRUE)
mean_abs_PdTD

mean_PdSS <- mean(unlist(Percent_diff_SS), na.rm=TRUE)
mean_PdSS
mean_abs_PdSS <- mean(abs(unlist(Percent_diff_SS)), na.rm=TRUE)
mean_abs_PdSS

TD_actual_95 <- lapply(TD_actual, function(i) {
  
  temp <- i[i!=max(i)]
  
  actual_95 <- max(temp, na.rm=TRUE)
  
  return(actual_95)
  
})

TD_actual_95

td95_tempdf <- data.frame(td95 = TD_actual_95, ss = M_sampsize)

td95 <- td95_tempdf %>%
  group_by(ss) %>%
  summarise(td95 = mean(td95, na.rm=TRUE))

td95_lm <- lm(log(td95)~log(ss), td95)
summary(td95_lm)
ss95 <- data.frame(ss = seq(50, 5000, 10))

td95_fit <- data.frame(td95=exp(predict(td95_lm, newdata=ss95, ss=ss95$ss)), ss=ss95$ss)

ggplot(td95_tempdf, aes(x=ss, y=td95))+
  geom_point(size=2)+
  geom_line(data=td95_fit, aes(x=ss, y=td95))+
  scale_y_reverse()
