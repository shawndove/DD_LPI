# recalculate trend deviation values using jaccard distance

library(TSdist)
library(philentropy)

dir_name <- "TestData/Testing2/"

dir_names <- list.dirs(path="TestData", 
                           full.names = TRUE, 
                           recursive = FALSE)

# remove unwanted directories
dir_names <- dir_names[601:700]
dir_names <- paste(dir_names, "/", sep="")

msi_final_samp_files <- vector()
msi_truetrend_files <- vector()
msi_lambda_files <- vector()
td_list <- vector()
td_mean <- vector()

for (i in 1:length(dir_names)) {
  
  msi_final_samp_files[i] <- list.files(paste(dir_names[i], 
                                           sep=""), 
                                     pattern = "msi_final_sampled_list")
  
  msi_truetrend_files[i] <- list.files(paste(dir_names[i], 
                                          sep=""), 
                                    pattern = "msi_TrueTrend.RData")

  msi_lambda_files[i] <- list.files(paste(dir_names[i], 
                                       sep=""), 
                                 pattern = "msi_lambda.RData")
  
  td_list[i] <- list.files(paste(dir_names[i], 
                              sep=""), 
                        pattern = "trend_dev_sampled_list_lambda")
  
  td_mean[i] <- list.files(paste(dir_names[i], 
                              sep=""), 
                        pattern = "trend_dev_sampled_mean_lambda")
  
}



#counter <- 1
for (z in 1:length(dir_names)) {
  
  msi_final_samp_list <- readRDS(file=paste(dir_names[z], msi_final_samp_files[z], sep=""))
  msi_real <- readRDS(file=paste(dir_names[z], msi_truetrend_files[z], sep=""))
  msi_final <- readRDS(file=paste(dir_names[z], msi_lambda_files[z], sep=""))
  
  # test how much the sampled trends deviate from the real trend
  trend.dev.list <- lapply(msi_final_samp_list, function(i) {jaccard(as.numeric(i), as.numeric(msi_real), testNA=FALSE)})
  trend.dev.matrix <- do.call(rbind, trend.dev.list)
  colnames(trend.dev.matrix) <- "MSI"
  write.csv(trend.dev.matrix, file=paste(dir_name, "fixed2/", td_list[z], sep=""))

  mean.trend.dev <- mean(trend.dev.matrix, na.rm=TRUE)
  write.csv(mean.trend.dev, file=paste(dir_name, "fixed2/", td_mean[z], sep=""))


#  counter <- counter + 1
  
}

