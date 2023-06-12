######################
### Author: Shawn Dove
######################

# This script collates results from LPI trend reliability testing for use in plots

## set directory paths ----

rd_name <- "RDataFiles" # name of directory where .RData files are stored
cd_name <- "CSV" # name of directory where .csv files will be saved

# create csv directory if it does not exist
if(!dir.exists(paste(cd_name, "/", sep=""))) {dir.create(paste(cd_name, "/", sep=""))}

# vector of weightings for taxonomic groups within realms (from McRae et al., 2017)
Weightings_list <- c(0.387205957, 0.197833813, 0.41496023, NA, 
                     0.396527091, 0.172106825, 0.431366084, NA, 
                     0.433535576, 0.249862107, 0.316602317, NA, 
                     0.387661234, 0.127987201, 0.484351565, NA, 
                     0.376366476, 0.249869859, 0.373763665, NA, 
                     NA, NA, NA, NA, 
                     0.192000, 0.009000, 0.207000, 0.590000, 
                     0.176000, 0.008000, 0.321000, 0.493000, 
                     0.211000, 0.015000, 0.179000, 0.592000, 
                     0.107000, 0.010000, 0.298000, 0.584000, 
                     0.203000, 0.013000, 0.217000, 0.565000, 
                     NA, NA, NA, NA, 
                     0.068635, 0.009774, 0.001303, 0.920286, 
                     0.069353, 0.006224, 0.001630, 0.922791, 
                     0.172867, 0.035011, 0.000000, 0.792123, 
                     0.054261, 0.022342, 0.000957, 0.922438, 
                     0.048714, 0.004878, 0.005505, 0.940901, 
                     0.080916, 0.025257, 0.000935, 0.892890)

# vector of weightings for realms within systems (from McRae et al., 2017)
Weightingsr_list <- c(0.189738, 0.292168, 0.116431, 0.321132, 0.061683, NA, 
                      0.211701, 0.225576, 0.123314, 0.365550, 0.060853, NA, 
                      0.146489, 0.214706, 0.014541, 0.099685, 0.456553, 0.068026)


## build data table ----

# import results calculated using 4.CalculateLPITrendReliability.R.
sample_percent <- readRDS(paste(rd_name, "/samp_percent_vec.RData", sep=""))

# create data frame
realms.results.df <- data.frame("SampPercent"=sample_percent)

# add realms
realms.results.df$Realm <- NA
realms.results.df$Realm[1:4] <- "Afrotropic"
realms.results.df$Realm[5:8] <- "Indopacific"
realms.results.df$Realm[9:12] <- "Palearctic"
realms.results.df$Realm[13:16] <- "Neotropic"
realms.results.df$Realm[17:20] <- "Nearctic"
realms.results.df$Realm[21:24] <- "Antarctic"
realms.results.df$Realm[25:28] <- "Afrotropic"
realms.results.df$Realm[29:32] <- "Indopacific"
realms.results.df$Realm[33:36] <- "Palearctic"
realms.results.df$Realm[37:40] <- "Neotropic"
realms.results.df$Realm[41:44] <- "Nearctic"
realms.results.df$Realm[45:48] <- "Antarctic"
realms.results.df$Realm[49:52] <- "AtNoTemp"
realms.results.df$Realm[53:56] <- "AtTrSub"
realms.results.df$Realm[57:60] <- "Arctic"
realms.results.df$Realm[61:64] <- "SoTeAnt"
realms.results.df$Realm[65:68] <- "TroSubIndo"
realms.results.df$Realm[69:72] <- "PaNoTemp"

realms.results.df$Realm <- factor(realms.results.df$Realm, 
                                  levels=c("Afrotropic", 
                                           "Neotropic",
                                           "Indopacific", 
                                           "Palearctic",
                                           "Nearctic",
                                           "Antarctic",
                                           "TroSubIndo",
                                           "AtTrSub",
                                           "AtNoTemp",
                                           "PaNoTemp",
                                           "SoTeAnt",
                                           "Arctic"),
                                  labels=c("Afrotropical", 
                                           "Neotropical",
                                           "IndoPacific", 
                                           "Palearctic",
                                           "Nearctic",
                                           "Antarctica",
                                           "IndoPacific",
                                           "Tropical Atlantic",
                                           "Temperate Atlantic",
                                           "Pacific temperate",
                                           "South temperate",
                                           "Arctic"))

# add systems
realms.results.df$System <- NA
realms.results.df$System[1:24] <- "Terrestrial"
realms.results.df$System[25:48] <- "Freshwater"
realms.results.df$System[49:72] <- "Marine"

realms.results.df$System <- factor(realms.results.df$System, levels=c("Terrestrial",
                                                                      "Freshwater",
                                                                      "Marine"))

# add taxa
realms.results.df$Taxon <- NA
realms.results.df$Taxon <- rep(c("Aves", "Mammalia", "Herps", "Fish"), 18)

realms.results.df$Taxon <- factor(realms.results.df$Taxon, levels = c("Aves", "Mammalia", "Herps", "Fish"), labels=c("Birds", "Mammals", "Reptiles & \n amphibians", "Fishes"))
realms.results.df$Taxon2 <- factor(realms.results.df$Taxon, levels = c("Birds", "Mammals", "Reptiles & \n amphibians", "Fishes"), labels=c("Birds", "Mammals", "Reptiles & amphibians", "Fishes"))

# add current tdv to data frame
realms.results.df$TDV <- readRDS(paste(rd_name, "/tdv_vec.RData", sep=""))

# add current sample size to data frame
realms.results.df$CurrentPops <- readRDS(paste(rd_name, "/pop_size_vec.RData", sep=""))

# add required sample size to data frame
realms.results.df$MinSampSize <- readRDS(paste(rd_name, "/samp_size_vec.RData", sep=""))

# add number of needed extra pops to data frame
realms.results.df$ExtraPopsNeeded <- readRDS(paste(rd_name, "/pops_added_vec.RData", sep=""))

# add taxonomic weightings to data frame
realms.results.df$TaxWeight <- Weightings_list

# add realm weightings to data frame
realms.results.df$RealmWeight <- rep(Weightingsr_list, each = 4)

# calculate relative weightings for each regional taxonomic group and add to data frame
realms.results.df$RelativeWeight <- realms.results.df$TaxWeight * realms.results.df$RealmWeight * 0.33

# add another percent column, this one stopping at 100%
realms.results.df$SampPercent100 <- realms.results.df$SampPercent

realms.results.df$SampPercent100 <- ifelse(realms.results.df$SampPercent100 <= 100, 
                                           realms.results.df$SampPercent100, 
                                           100)

# create results tables for each system
terr.results.df <- realms.results.df[realms.results.df$System=="Terrestrial",]
terr.results.df <- terr.results.df[terr.results.df$Taxon!="Fishes" & terr.results.df$Realm!="Antarctica",]
fw.results.df <- realms.results.df[realms.results.df$System=="Freshwater",]
fw.results.df <- fw.results.df[fw.results.df$Realm!="Antarctica",]
marine.results.df <- realms.results.df[realms.results.df$System=="Marine",]
marine.results.df$Taxon <- factor(marine.results.df$Taxon, levels = c("Birds", "Mammals", "Reptiles & \n amphibians", "Fishes"), labels=c("Birds", "Mammals", "Reptiles", "Fishes"))


# save results tables
saveRDS(realms.results.df, file=paste(rd_name, "/realms_results_df.RData", sep=""))
saveRDS(terr.results.df, file=paste(rd_name, "/terr_results_df.RData", sep=""))
saveRDS(fw.results.df, file=paste(rd_name, "/fw_results_df.RData", sep=""))
saveRDS(marine.results.df, file=paste(rd_name, "/marine_results_df.RData", sep=""))


#########

# Pearsons correlation coefficient testing ----

pcdf <- realms.results.df[!is.na(realms.results.df$SampPercent) & !is.na(realms.results.df$RelativeWeight),]
pcdf.ter <- pcdf[pcdf$System=="Terrestrial",]
pcdf.fw <- pcdf[pcdf$System=="Freshwater",]
pcdf.mar <- pcdf[pcdf$System=="Marine",]
summary(cor(pcdf$SampPercent, pcdf$RelativeWeight, method="pearson"))
cor.test(pcdf$SampPercent, pcdf$RelativeWeight)
cor.test(pcdf.ter$SampPercent, pcdf.ter$RelativeWeight)
cor.test(pcdf.fw$SampPercent, pcdf.fw$RelativeWeight)
cor.test(pcdf.mar$SampPercent, pcdf.mar$RelativeWeight)


#########

## CSV Export ----

# create a new data frame for export to csv and calculation of median minimum sample sizes
SampSizeExportDF <- data.frame(System = realms.results.df$System, 
                               Realm = realms.results.df$Realm, 
                               Taxon = realms.results.df$Taxon,
                               TDV = realms.results.df$TDV,
                               CurrentPops = realms.results.df$CurrentPops,
                               MinSampSize = realms.results.df$MinSampSize, 
                               PopsNeeded = realms.results.df$ExtraPopsNeeded)

# remove NAs
SampSizeExportDF <- SampSizeExportDF[!is.na(SampSizeExportDF$MinSampSize),]
levels(SampSizeExportDF$Taxon) <- c("Birds", "Mammals", "Reptiles & Amphibians", "Fishes")

# remove Antarctica, as it does not contain enough data
SampSizeExportDF <- SampSizeExportDF[SampSizeExportDF$Realm!="Antarctica",]

# save csv file
write.csv(SampSizeExportDF, file=paste(cd_name, "/SampSizeExportTable_3.csv", sep=""))


## Calculate medians ----

# calculate medians for taxa
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Taxon=="Fishes"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Taxon=="Birds"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Taxon=="Mammals"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Taxon=="Reptiles & Amphibians"])

# calculate medians for systems
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$System=="Marine"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$System=="Terrestrial"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$System=="Freshwater"])

# calculate median for global south
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Realm=="Afrotropical" |
                                      SampSizeExportDF$Realm=="IndoPacific" |
                                      SampSizeExportDF$Realm=="Neotropical" |
                                      SampSizeExportDF$Realm=="Tropical Atlantic" |
                                      SampSizeExportDF$Realm=="South temperate"])

# calculate median for global north
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Realm=="Palearctic" |
                                      SampSizeExportDF$Realm=="Nearctic" |
                                      SampSizeExportDF$Realm=="Temperate Atlantic" |
                                      SampSizeExportDF$Realm=="Arctic" |
                                      SampSizeExportDF$Realm=="Pacific temperate"])

