TD_actual <- vector()
TD_expected <- vector()
Percent_diff_TD <- vector()
Percent_diff_SS <- vector()
SS_expected <- vector()
M_sampsize <- vector()
M_meangr <- vector()
M_sdgr <- vector()
M_tslength <- vector()
counter <- 1
for (i in 14001:14100) {
  
  M_info <- read.csv(file=paste("TestData/", i, "/saved_synth_", i, "_info.csv", sep=""))
  
  M_sampsize[counter] <- M_info$samp_size
  
  M_meangr[counter] <- M_info$mean_gr_raw
  
  M_sdgr[counter] <- M_info$gr_sd_raw
  
  M_tslength[counter] <- M_info$mean_ts_length
  
  TD_actual[counter] <- read.csv(file=paste("TestData/Testing2/fixed4/saved_synth_", i, "_trend_dev_sampled_mean_lambda.csv", sep=""))[,2]
  
 # TD_expected[counter] <- exp((M_sampsize[counter] * (-0.003)) + (M_meangr[counter] * 6.582) + (M_tslength[counter] * (-0.0351)) + (M_sdgr[counter] * 75.31) - 11.64)

 # SS_expected[counter] <- ((M_meangr[counter] * 6.582) + (M_tslength[counter] * (-0.0351)) + (M_sdgr[counter] * 75.31) - log(TD_actual[counter]) - 11.64) / 0.003
  
  TD_expected[counter] <- exp((log(M_sampsize[counter]) * (-0.967)) + (M_meangr[counter] * 6.836) + (M_tslength[counter] * (-0.0346)) + (log(M_sdgr[counter]) * 1.974) + 1.812)
  
  SS_expected[counter] <- exp(((M_meangr[counter] * 6.836) + (M_tslength[counter] * (-0.0346)) + (log(M_sdgr[counter]) * 1.974) - log(TD_actual[counter]) + 1.812) / 0.967)
  
  
  Percent_diff_TD[counter] <- ((TD_actual[counter] - TD_expected[counter]) / ((TD_actual[counter] + TD_expected[counter]) / 2)) * 100
  
  Percent_diff_SS[counter] <- ((M_sampsize[counter] - SS_expected[counter]) / ((M_sampsize[counter] + SS_expected[counter]) / 2)) * 100
  
  counter <- counter + 1
  
}


