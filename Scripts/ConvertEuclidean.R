library(TSdist)
counter <- 1
for (z in 1001:1025) {
  
  grp_final_samp_list <- readRDS(file=paste("saved_synth_", z, "_grpi_final_sampled_list_", k_choice_a[counter], "_lambda.RData", sep=""))
  msi_final_samp_list <- readRDS(file=paste("saved_synth_", z, "_msi_final_sampled_list_", k_choice_a[counter], "_lambda.RData", sep=""))
  grp_real <- readRDS(file=paste("saved_synth_", z, "_group_indices_TrueTrend.RData", sep=""))
  msi_real <- readRDS(file=paste("saved_synth_", z, "_msi_TrueTrend.RData", sep=""))
  grp_final <- readRDS(file=paste("saved_synth_", z, "_group_indices_lambda.RData", sep=""))
  msi_final <- readRDS(file=paste("saved_synth_", z, "_msi_lambda.RData", sep=""))
  
  # test how much the sampled trends deviate from the real trend
  trend.dev.list <- lapply(msi_final_samp_list, function(i) {if(anyNA(i)) {NA} else {EuclideanDistance(as.numeric(i), as.numeric(msi_real))}})
  trend.dev.matrix <- do.call(rbind, trend.dev.list)
  colnames(trend.dev.matrix) <- "MSI"
  write.csv(trend.dev.matrix, file=paste("saved_synth_", z, "_trend_dev2_sampled_list_lambda.csv", sep=""))
  
  mean.trend.dev <- mean(trend.dev.matrix, na.rm=TRUE)
  write.csv(mean.trend.dev, file=paste("saved_synth_", z, "_trend_dev2_sampled_mean_lambda.csv", sep=""))
  
  # test how much the sampled group trends deviate from the real trend
  trend.grp.dev.list <- lapply(grp_final_samp_list, function(i) {lapply(1:length(i), function(j) {
    if (anyNA(i[[j]][,1:c])) {NA} else {EuclideanDistance(as.numeric(i[[j]][,1:c]), as.numeric(grp_real[[j]][,1:c]))}})})
  trend.grp.dev <- as.data.frame(do.call(rbind, trend.grp.dev.list))
  colnames(trend.grp.dev) <- gsub("V", "Grp", colnames(trend.grp.dev))
  write.csv(as.matrix(trend.grp.dev), file=paste("saved_synth_", z, "_trend_dev2_grp_sampled_list_lambda.csv", sep=""))
  
  trend.grp.dev.temp <- as.vector(do.call(rbind, do.call(rbind, trend.grp.dev.list)))
  mean.trend.grp.dev <- vector()
  for (i in 1:length(grp_final)) {
    mean.trend.grp.dev[i] <- mean(trend.grp.dev.temp[((i*bootstrap_size)-bootstrap_size+1):(i*bootstrap_size)], na.rm=TRUE)
  }
  write.csv(mean.trend.grp.dev, file=paste("saved_synth_", z, "_trend_dev2_grp_sampled_mean_lambda.csv", sep=""))
  
  counter <- counter + 1
  
}

