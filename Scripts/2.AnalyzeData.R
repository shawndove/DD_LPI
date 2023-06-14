######################
### Author: Shawn Dove
######################

# This script calculates a trend accuracy model from generated time series data
# and plots results from parameter testing.


## load packages ----

library(tidyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(QuantPsyc)


## load external functions ----

source("Scripts/GatherResultsFunction.R")


## set directory names ----

# create directory for storing .RData files
rd_name <- "RDataFiles" # directory name
if(!dir.exists(paste(rd_name, "/", sep=""))) {dir.create(paste(rd_name, "/", sep=""))} # create directory

# create directory for saving plots
pd_name <- "Plots" # directory name
if(!dir.exists(paste(pd_name, "/", sep=""))) {dir.create(paste(pd_name, "/", sep=""))} # create directory


################
# Accuracy Model
################

# test directory path
dir_name <- "TestData/ModelData" # directory containing the data to analyze

# get results
test_res_full <- gather_results_fn(dir_name)
test_results <- test_res_full[[1]]
test_results_m <- test_res_full[[2]]

# remove datasets in which the mean of any of these parameters falls outside of 
# the LPD parameter ranges
test_results_mc <- test_results_m[test_results_m$SDGRSamp < 0.63 
                             & test_results_m$SDGRSamp > 0.12 
                             & test_results_m$MeanSDSamp > 0.16 
                             & test_results_m$MeanSDSamp < 0.89
                             & test_results_m$MeanGRSamp > -0.19
                             & test_results_m$MeanGRSamp < 0.16
                             & test_results_m$MeanTSLength > 6.0
                             & test_results_m$MeanTSLength < 39,]

saveRDS(test_results_mc, file=paste(rd_name, "/", "test_results_mc.RData", sep=""))

# randomly select data for training the model
model_train_IDs <- sample(test_results_mc$ID, size=0.67*nrow(test_results_mc), replace=FALSE)
model_train_df <- test_results_mc[test_results_mc$ID %in% model_train_IDs,]
saveRDS(model_train_df, file=paste(rd_name, "/", "model_train_df.RData", sep=""))
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
saveRDS(model_test_df_mean, file=paste(rd_name, "/", "model_test_df_mean.RData", sep=""))
saveRDS(model_test_df, file=paste(rd_name, "/", "model_test_df.RData", sep=""))

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
saveRDS(model_pops, file=paste(rd_name, "/", "model_pops.RData", sep=""))

# get beta coefficients
model_pops_beta <- lm.beta(model_pops)
model_pops_beta
saveRDS(model_pops_beta, file=paste(rd_name, "/", "model_pops_beta.RData", sep=""))



##########################################

## build supplementary plots for trend length, number of pops per species,
## clustered vs random distribution of time points, size of data set, and
## amount of observation error


##############
# Trend Length
##############

# test directory path
dir_name <- "TestData/TrendLength" # directory containing the data to analyze

# get results
trendlength <- gather_results_fn(dir_name)[[2]]

trendlength$TotalYears <- factor(trendlength$TotalYears)

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

ggsave(paste(pd_name, "/", "trendlength.tiff", sep=""),
       device = tiff,
       dpi = 1000,
       compression = "lzw",
       width = 5000,
       height = 3000,
       units = "px")


###################
# Observation Error
###################

# test directory path
dir_name <- "TestData/ObsError" # directory containing the data to analyze

# get results
obserror <- gather_results_fn(dir_name)[[2]]

obserror$ObsError <- rep(c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2), each=20)

obserror$ObsError <- factor(obserror$ObsError)

ggplot(obserror, aes(x=ObsError, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
  scale_x_discrete(labels=c("0%", "10%", "20%", "40%", "60%", "80%", "100%", "200%"))+
  ylim(c(0,0.5))+
  ylab("Trend Deviation Value")+
  xlab("Percent Observation Error")+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave(paste(pd_name, "/", "observationerror.tiff", sep=""),
       device = tiff,
       dpi = 1000,
       compression = "lzw",
       width = 5000,
       height = 3000,
       units = "px")


#########################
# Populations Per Species
#########################

# test directory path
dir_name <- "TestData/PPSDist" # directory containing the data to analyze

# get results
pps <- gather_results_fn(dir_name)[[2]]

pps <- pps %>%
  mutate(PopSpec = cut(PopSpec, breaks=c(3,7,12,18,25,40,100),
                       labels=c(5,10,15,20,30,50)))

pps_norm <- pps[pps$TSGenVersion==7,]
pps_exp <- pps[pps$TSGenVersion==8,]
pps_nb <- pps[pps$TSGenVersion==9,]

# normal distribution
pps_norm_plot <- ggplot(pps_norm, aes(x=PopSpec, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
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

# discretized exponential distribution
pps_exp_plot <- ggplot(pps_exp, aes(x=PopSpec, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
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

# zero-truncated negative binomial distribution
pps_nb_plot <- ggplot(pps_nb, aes(x=PopSpec, y=TrendDev))+
  geom_boxplot(show.legend=FALSE, fill="skyblue")+
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

ppsplot <- ggarrange(pps_norm_plot + rremove("ylab") + rremove("xlab"),
                    pps_exp_plot + rremove("ylab") + rremove("xlab"),
                    pps_nb_plot + rremove("ylab") + rremove("xlab"),
                    ncol=3)

ppsplot <- annotate_figure(ppsplot,
                left = textGrob("Trend Deviation Value",
                                rot = 90, vjust = 0.5,
                                gp = gpar(cex = 1.3)),
                bottom = textGrob("Mean Populations Per Species",
                                  hjust = 0.5,
                                  gp = gpar(cex = 1.3)))

ggsave(paste(pd_name, "/", "popsperspec_all.tiff", sep=""),
       plot = ppsplot,
       device = tiff,
       dpi = 1000,
       compression = "lzw",
       width = 5000,
       height = 1500,
       units = "px")


########################
# Timepoint Distribution
########################

# test directory path
dir_name <- "TestData/TPDist" # directory containing the data to analyze

# get results
tpdist <- gather_results_fn(dir_name)[[2]]

tpdist$DegradeType <- factor(tpdist$DegradeType)

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

ggsave(paste(pd_name, "/", "timepointdistribution.tiff", sep=""),
       device = tiff,
       dpi = 1000,
       compression = "lzw",
       width = 5000,
       height = 3000,
       units = "px")


##############################
# Dataset Size and Sample Size
##############################

# test directory path
dir_name <- "TestData/DSSize" # directory containing the data to analyze

# get results
dssize <- gather_results_fn(dir_name)[[2]]

dssize$TotalPops <- factor(dssize$TotalPops)

dssize50 <- dssize[dssize$SampSize==50,]
dssize100 <- dssize[dssize$SampSize==100,]
dssize200 <- dssize[dssize$SampSize==200,]
dssize500 <- dssize[dssize$SampSize==500,]

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

dsplot <- ggarrange(ds50 + rremove("ylab") + rremove("xlab"),
                    ds100 + rremove("ylab") + rremove("xlab"),
                    ds200 + rremove("ylab") + rremove("xlab"),
                    ds500 + rremove("ylab") + rremove("xlab"),
                    ncol=2, nrow=2)

dsplot <- annotate_figure(dsplot,
                left = textGrob("Trend Deviation Value",
                                rot = 90, vjust = 0.5,
                                gp = gpar(cex = 1.3)),
                bottom = textGrob("Total Populations",
                                  hjust = 0.5,
                                  gp = gpar(cex = 1.3)))

ggsave(paste(pd_name, "/", "totalpopsall.tiff", sep=""),
       plot=dsplot,
       device = tiff,
       dpi = 1000,
       compression = "lzw",
       width = 5000,
       height = 3000,
       units = "px")


###########
# Solutions
###########

# test directory path
dir_name <- "TestData/Solutions" # directory containing the data to analyze

# get results
solutions <- gather_results_fn(dir_name)[[2]]

solutions$ep_ratio <- rep(c(1, 0.2, 1, 1), each=20)
solutions$group <- rep(c("Control", "Group A", "Group B", "Group C"), each=20)

solutions$group <- factor(solutions$group, 
                               levels=c("Control", "Group A", "Group B", "Group C"),
                               labels=c("Control", "Group A", "Group B", "Group C"))

ggplot(solutions, aes(x=group, y=TrendDev, fill=factor(group)))+
  geom_boxplot(show.legend=FALSE)+
  ylim(c(0,0.6))+
  ylab("TDV")+
  xlab("")+
  theme_bw()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave(paste(pd_name, "/", "solutions.tiff", sep=""),
       device = tiff,
       dpi = 1000,
       compression = "lzw",
       width = 5000,
       height = 3000,
       units = "px")
