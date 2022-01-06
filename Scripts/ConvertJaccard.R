# recalculate trend deviation values using jaccard distance

library(TSdist)
library(philentropy)

dir_name <- "TestData/Testing2/"

dir_names <- list.dirs(path="TestData", 
                           full.names = TRUE, 
                           recursive = FALSE)

# remove directory Testing2
dir_names <- dir_names[1:300]
dir_names <- paste(dir_names, "/", sep="")

msi_final_samp_files <- vector()
grp_final_samp_files <- vector()
msi_truetrend_files <- vector()
grp_truetrend_files <- vector()
msi_lambda_files <- vector()
grp_lambda_files <- vector()
td_list <- vector()
td_mean <- vector()
td_grp_list <- vector()
td_grp_mean <- vector()

for (i in 1:length(dir_names)) {
  
  msi_final_samp_files[i] <- list.files(paste(dir_names[i], 
                                           sep=""), 
                                     pattern = "msi_final_sampled_list")
  
  grp_final_samp_files[i] <- list.files(paste(dir_names[i], 
                                           sep=""), 
                                     pattern = "grpi_final_sampled_list")
  
  msi_truetrend_files[i] <- list.files(paste(dir_names[i], 
                                          sep=""), 
                                    pattern = "msi_TrueTrend.RData")
  
  grp_truetrend_files[i] <- list.files(paste(dir_names[i], 
                                          sep=""), 
                                    pattern = "group_indices_TrueTrend.RData")
  
  msi_lambda_files[i] <- list.files(paste(dir_names[i], 
                                       sep=""), 
                                 pattern = "msi_lambda.RData")
  
  grp_lambda_files[i] <- list.files(paste(dir_names[i], 
                                       sep=""), 
                                 pattern = "group_indices_lambda.RData")
  
  td_list[i] <- list.files(paste(dir_names[i], 
                              sep=""), 
                        pattern = "trend_dev_sampled_list_lambda")
  
  td_mean[i] <- list.files(paste(dir_names[i], 
                              sep=""), 
                        pattern = "trend_dev_sampled_mean_lambda")
  
  td_grp_list[i] <- list.files(paste(dir_names[i], 
                                  sep=""), 
                            pattern = "trend_dev_grp_sampled_list_lambda")
  
  td_grp_mean[i] <- list.files(paste(dir_names[i], 
                                  sep=""), 
                            pattern = "trend_dev_grp_sampled_mean_lambda")
  
}



#counter <- 1
for (z in 1:length(dir_names)) {
  
  grp_final_samp_list <- readRDS(file=paste(dir_names[z], grp_final_samp_files[z], sep=""))
  msi_final_samp_list <- readRDS(file=paste(dir_names[z], msi_final_samp_files[z], sep=""))
  grp_real <- readRDS(file=paste(dir_names[z], grp_truetrend_files[z], sep=""))
  msi_real <- readRDS(file=paste(dir_names[z], msi_truetrend_files[z], sep=""))
  grp_final <- readRDS(file=paste(dir_names[z], grp_lambda_files[z], sep=""))
  msi_final <- readRDS(file=paste(dir_names[z], msi_lambda_files[z], sep=""))
  
  # test how much the sampled trends deviate from the real trend
  trend.dev.list <- lapply(msi_final_samp_list, function(i) {additive_symm_chi_sq(as.numeric(i), as.numeric(msi_real), testNA=FALSE)})
  trend.dev.matrix <- do.call(rbind, trend.dev.list)
  colnames(trend.dev.matrix) <- "MSI"
  write.csv(trend.dev.matrix, file=paste(dir_name, "fixed/", td_list[z], sep=""))

  mean.trend.dev <- mean(trend.dev.matrix, na.rm=TRUE)
  write.csv(mean.trend.dev, file=paste(dir_name, "fixed/", td_mean[z], sep=""))

  # test how much the sampled group trends deviate from the real trend
  trend.grp.dev.list <- lapply(grp_final_samp_list, function(i) {lapply(1:length(i), function(j) {
    additive_symm_chi_sq(as.numeric(i[[j]][,1:c]), as.numeric(grp_real[[j]][,1:c]), testNA=FALSE)})})
  trend.grp.dev <- as.data.frame(do.call(rbind, trend.grp.dev.list))
  colnames(trend.grp.dev) <- gsub("V", "Grp", colnames(trend.grp.dev))
  write.csv(as.matrix(trend.grp.dev), file=paste(dir_name, "fixed/", td_grp_list[z], sep=""))

  trend.grp.dev.temp <- as.vector(do.call(rbind, do.call(rbind, trend.grp.dev.list)))
  mean.trend.grp.dev <- vector()
  for (i in 1:length(grp_final)) {
    mean.trend.grp.dev[i] <- mean(trend.grp.dev.temp[((i*bootstrap_size)-bootstrap_size+1):(i*bootstrap_size)], na.rm=TRUE)
  }
  write.csv(mean.trend.grp.dev, file=paste(dir_name, "fixed/", td_grp_mean[z], sep=""))

#  counter <- counter + 1
  
}

