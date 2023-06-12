######################
### Author: Shawn Dove
######################

# This script calculates all regional taxonomic LPI trends, the reliability
# for each trend based on an accuracy model built from test data, and the minimum
# sample size needed to reach a predetermined reliability threshold.


# load packages ----

library(dplyr)
library(mgcv)
library(matrixStats)
library(ggplot2)
library(gridExtra)


# load external functions ----

source("Scripts/CullFunction.R")
source("Scripts/InterpolateFunction.R")
source("Scripts/GrowthRateCalcFunction.R")


# set directory paths ----

rd_name <- "RDataFiles" # name of directory where .RData files are stored
dd_name <- "Data" # name of directory where LPD file is stored
pd_name <- "Plots" # name of directory where plots will be saved

# create plots directory if it does not exist
if(!dir.exists(paste(pd_name, "/", sep=""))) {dir.create(paste(pd_name, "/", sep=""))}

## load publicly available LPD data ----

# read in csv file (change directory and filename if necessary)
LPI_full <- read.csv(file=paste(dd_name, "/", "LPD2022_public.csv", sep=""), sep=",", stringsAsFactors=FALSE)

# fix ID column name
colnames(LPI_full)[1] <- "ID"

# remove duplicates
LPI_full <- LPI_full[LPI_full$Replicate==0,]

# create groupings ----

# taxa (note that this will not work properly with old versions of the LPD,
# as the taxonomy of the fish has changed slightly)
tax_group<- list(Birds=c("Aves"), 
                 Mammals=c("Mammalia"), 
                 Herps=c("Reptilia", "Amphibia"), 
                 Fish=c("Actinopteri", "Elasmobranchii", "Petromyzonti", "Dipneusti", "Holocephali", "Myxini", "Coelacanthi"))

# create an ID key for taxonomic groups
tax_group_key <- lapply(1:length(tax_group), function(i) {
  Group <- names(tax_group[i])
  Class <- tax_group[[i]]
  ID <- i
  c <- data.frame(Group, Class, ID)
  return(c)})

tax_group_key <- do.call(rbind, tax_group_key) # bind list into a data frame

tax_group_IDs <- setNames(tax_group_key$ID, tax_group_key$Class) # convert to vector with named elements

# terrestrial realms
t_realm <- list(Afrotropical=c("Afrotropical"),
              IndoPacific=c("Australasia", "Oceania", "Indo-Malayan"),
              Palearctic=c("Palearctic"),
              Neotropical=c("Neotropical"),
              Nearctic=c("Nearctic"),
              Antarctic=c("Antarctic"))

# create an ID key for terrestrial realms
t_realm_key <- lapply(1:length(t_realm), function(i) {
  Group <- names(t_realm[i])
  Class <- t_realm[[i]]
  ID <- i
  c <- data.frame(Group, Class, ID)
  return(c)})

t_realm_key <- do.call(rbind, t_realm_key) # bind list into a data frame

t_realm_IDs <- setNames(t_realm_key$ID, t_realm_key$Class) # convert to vector with named elements

# freshwater realms
fw_realm <- list(Afrotropical=c("Afrotropical"),
                IndoPacific=c("Australasia", "Oceania", "Indo-Malayan"),
                Palearctic=c("Palearctic"),
                Neotropical=c("Neotropical"),
                Nearctic=c("Nearctic"),
                Antarctic=c("Antarctic"))

# create an ID key for freshwater realms
fw_realm_key <- lapply(1:length(fw_realm), function(i) {
  Group <- names(fw_realm[i])
  Class <- fw_realm[[i]]
  ID <- i
  c <- data.frame(Group, Class, ID)
  return(c)})

fw_realm_key <- do.call(rbind, fw_realm_key) # bind list into a data frame

fw_realm_IDs <- setNames(fw_realm_key$ID, fw_realm_key$Class) + length(unique(t_realm_IDs)) # convert to vector with named elements

# marine realms
m_realm <- list(M_Atlantic_NT=c("Atlantic north temperate"),
                Atlantic_TS=c("Atlantic tropical and subtropical"),
                Arctic=c("Arctic"),
                South_TA=c("South temperate and Antarctic"),
                Tropical_SI=c("Tropical and subtropical Indo-Pacific"),
                Pacific_NT=c("Pacific north temperate"))

# create an ID key for marine realms
m_realm_key <- lapply(1:length(m_realm), function(i) {
  Group <- names(m_realm[i])
  Class <- m_realm[[i]]
  ID <- i
  c <- data.frame(Group, Class, ID)
  return(c)})

m_realm_key <- do.call(rbind, m_realm_key) # bind list into a data frame

m_realm_IDs <- setNames(m_realm_key$ID, m_realm_key$Class)  + length(unique(t_realm_IDs)) + length(unique(fw_realm_IDs)) # convert to vector with named elements

# systems
sys_grp <- list(Terrest=c("Terrestrial"),
                 FW=c("Freshwater"),
                 Marine=c("Marine"))

# create an ID key for systems
sys_key <- lapply(1:length(sys_grp), function(i) {
  Group <- names(sys_grp[i])
  Class <- sys_grp[[i]]
  ID <- i
  c <- data.frame(Group, Class, ID)
  return(c)})

sys_key <- do.call(rbind, sys_key) # bind list into a data frame

sys_IDs <- setNames(sys_key$ID, sys_key$Class) # convert to vector with named elements

# create list of taxonomic groups for index list
tax_list <- rep(unique(tax_group_IDs), length(unique(m_realm_IDs)) + length(unique(t_realm_IDs)) + length(unique(fw_realm_IDs)))

# create list of realms for index list
realm_list <- rep(c(unique(t_realm_IDs), unique(fw_realm_IDs), unique(m_realm_IDs)), each=length(unique(tax_group_IDs)))

# create list of systems for index list
sys_list <- rep(unique(sys_IDs), each=length(unique(tax_group_IDs)) * length(unique(t_realm_IDs)))


####

# create subset of LPD columns with only count data
LPI_trimmed <- LPI_full[,grep("(\\d)", names(LPI_full))]

# convert to numeric data (ignore warnings)
LPI_trimmed <- as.data.frame(sapply(LPI_trimmed, as.numeric))

# create a population ID column by assigning each population a separate ID based on its row number
LPI_trimmed$PopID <- LPI_full$ID

# create a species ID column by assigning each unique species a separate ID number
LPI_trimmed$SpecID <- match(LPI_full$Binomial, unique(LPI_full$Binomial))

# create a taxonomic group ID column by assigning species groups separate ID numbers according to LPI groupings
LPI_trimmed$GrpID <- tax_group_IDs[LPI_full$Class]

# create a taxonomic group ID column by assigning terrestrial realms separate ID numbers according to LPI groupings
LPI_trimmed$TRID <- t_realm_IDs[LPI_full$T_realm]
  
# create a taxonomic group ID column by assigning freshwater realms separate ID numbers according to LPI groupings
LPI_trimmed$FWRID <- fw_realm_IDs[LPI_full$FW_realm]
  
# create a taxonomic group ID column by assigning marine realms separate ID numbers according to LPI groupings
LPI_trimmed$MRID <- m_realm_IDs[LPI_full$M_realm]
  
# create a taxonomic group ID column by assigning systems separate ID numbers according to LPI groupings
LPI_trimmed$SysID <- sys_IDs[LPI_full$System]

pop_list <- list()
# select populations to form each group index
for (i in 1:length(tax_list)) {
  
  temp <- LPI_trimmed$PopID[which(LPI_trimmed$SysID==sys_list[i] & 
                                   LPI_trimmed$GrpID==tax_list[i] & 
                                   (LPI_trimmed$TRID==realm_list[i] | 
                                      LPI_trimmed$FWRID==realm_list[i] | 
                                      LPI_trimmed$MRID==realm_list[i]))]
  
  pop_list[[i]] <- temp
  
}

####

# calculate number of populations (sample sizes) in LPI groups

pop.size.list <- list()

for (i in 1:length(pop_list)) {
  
  pop.size.list[[i]] <- length(pop_list[[i]])
  
}

pop.size.vec <- unlist(pop.size.list)

# save pop size data
saveRDS(pop.size.vec, file=paste(rd_name, "/", "pop_size_vec.RData", sep=""))


# calculate number of species in LPI groups

spec.size.list <- list()

for (i in 1:length(pop_list)) {
  
  spec.size.list[[i]] <- length(unique(LPI_trimmed$SpecID[LPI_trimmed$PopID %in% pop_list[[i]]]))
  
}

spec.size.vec <- unlist(spec.size.list)

# save spec size data
saveRDS(spec.size.vec, file=paste(rd_name, "/", "spec_size_vec.RData", sep=""))


####

## calculate mean growth rate and mean time series length for each actual LPI group

gr.stats.list <- list()

mean.tslength.list <- list()

# make a copy of LPD data
LPI_trimmed2 <- LPI_trimmed

# remove X from year column names
names(LPI_trimmed2) <- gsub("X", "", names(LPI_trimmed2))

# create a vector of years that the LPD dataset covers
year_names <- names(LPI_trimmed2[grep("(\\d)", names(LPI_trimmed2))])

# calculate total number of years that the LPD dataset covers
t_years <- length(year_names)

for (i in 1:length(pop_list)) {
  
  # select all copies of the populations listed in the sample
  group_data <- LPI_trimmed2[LPI_trimmed2$PopID %in% pop_list[[i]],]
  
  # remove all populations with less than 2 data points
  group_data_culled <- cull_fn(group_data, 2, 2, t_years)
  
  # log-linear interpolate all pops
  grp_completed <- interpolate_fn(group_data_culled, t_years, year_names, lessthansix=FALSE, calcsd=TRUE)
  
  # calculate number of non-NA values and divide by total number of time series
  # this gives the mean time series length
  mean.tslength.list[[i]] <- sum(!is.na(as.vector(grp_completed[,1:t_years]))) / pop.size.list[[i]]
  
  if (nrow(grp_completed) >=1) {
    
    # get mean and standard deviation of the mean growth rate
    gr.stats.list[[i]] <- growth_rate_calc_fn(grp_completed, t_years, model=TRUE)
    
  } else {
    
    gr.stats.list[[i]] <- NA
    
  }
  
}

mean.tslength.vec <- unlist(mean.tslength.list)

# save time series length data
saveRDS(mean.tslength.vec, file=paste(rd_name, "/", "mean_ts_length_vec.RData", sep=""))

# save mean growth rate data
saveRDS(gr.stats.list, file=paste(rd_name, "/", "gr_stats_list.RData", sep=""))


####

## find tdv and minimum sample sizes

# load the max_tdv calculated in script 3.FindMaximumTDV.R. This is needed for the model.
max_tdv <- readRDS(paste(rd_name, "/", "max_tdv.RData", sep=""))

tdv.list <- list()

samp.size.list <- list()

for (i in 1:length(gr.stats.list)) {
  
  meangr <- as.numeric(gr.stats.list[[i]][1]) # get mean growth rate for group i
  
  stdev <- as.numeric(gr.stats.list[[i]][2]) # get standard deviation of mean growth rate for group i
  
  tslength <- mean.tslength.list[[i]] # get mean time series length for group i
  
  meansd <- as.numeric(gr.stats.list[[i]][3]) # get mean of population growth rate standard deviations for group i
  
  popsize <- pop.size.list[[i]] # get number of existing populations for group i
  
  if (is.na(stdev)) { # if standard deviation is NA...
    
    tdv.list[[i]] <- NA # set tdv to NA
    
    samp.size.list[[i]] <- NA # set sample size to NA
    
  } else { # otherwise...
    
    tdv.list[[i]] <- exp(
      model_pops$coefficients[1] 
      + (model_pops$coefficients[2] * log(popsize))
      + (model_pops$coefficients[3] * log(stdev))
      + (model_pops$coefficients[4] * meangr)
      + (model_pops$coefficients[5] * meansd)
      + (model_pops$coefficients[6] * tslength)
      )

    samp.size.list[[i]] <- exp(
        (
          model_pops$coefficients[1]
          + (model_pops$coefficients[4] * meangr)
          + (model_pops$coefficients[6] * tslength)
          + (model_pops$coefficients[3] * log(stdev))
          + (model_pops$coefficients[5] * meansd)
          - log(max_tdv)
          )
        / (-model_pops$coefficients[2])
        )
    
  }
  
}

tdv.vec <- unlist(tdv.list)

samp.size.vec <- unlist(samp.size.list)

# save tdv and sample size data
saveRDS(tdv.vec, file=paste(rd_name, "/", "tdv_vec.RData", sep=""))
saveRDS(samp.size.vec, file=paste(rd_name, "/", "samp_size_vec.RData", sep=""))


# calculate percentage of min sample size in data for each group

samp.percent.list <- list()

for (i in 1:length(gr.stats.list)) {
  
  temp.min <- samp.size.list[[i]]
  
  temp.actual <- pop.size.list[[i]]
  
  ratio <- temp.actual/temp.min
  
  if (is.na(ratio)) {
    
    samp.percent.list[[i]] <- NA
    
  } else {
    
    samp.percent.list[[i]] <- ratio * 100
    
  }
  
}

samp.percent.vec <- unlist(samp.percent.list)

# save percentage sample size data
saveRDS(samp.percent.vec, file=paste(rd_name, "/", "samp_percent_vec.RData", sep=""))


# calculate number of populations that need to be added to get below the tdv threshold

pops.added.list <- list()

for (i in 1:length(gr.stats.list)) {
  
  if (is.na(samp.size.list[[i]]) | is.nan(samp.size.list[[i]])) {
    
    pops.added.list[[i]] <- NA
    
  }
  
  else if (pop.size.list[[i]] <= samp.size.list[[i]]) {
    
    pops.added.list[[i]] <- samp.size.list[[i]] - pop.size.list[[i]]
    
  } else {
    
    pops.added.list[[i]] <- 0
    
  }

}

pops.added.vec <- unlist(pops.added.list)

# save data
saveRDS(pops.added.vec, file=paste(rd_name, "/", "pops_added_vec.RData", sep=""))


##########

# build plot of observations per year in the LPD ----

obsyear <- vector()
for (i in 1:(length(LPI_trimmed)-7)) {
  obsyear[i] <- length(which(!is.na(LPI_trimmed[,i])))
}

obsyeardf <- data.frame("Year" = 1950:2020, "Obs" = obsyear)


# generate plot

ggplot(data = obsyeardf, aes(x=Year, y=Obs))+
  geom_point(size=4, colour="skyblue")+
  labs(x="Year", y="Number of Observations")+
  scale_x_continuous(breaks=c(1950,1960,1970,1980,1990,2000,2010,2020))+
  theme_bw()+
  theme(axis.title.x=element_text(size=22),
        axis.title.y=element_text(size=22),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))


# save plot

ggsave(paste(pd_name, "/observations_per_year_update.tiff", sep=""),
       plot = last_plot(),
       device = "tiff",
       width = 12000,
       height = 8000,
       units = "px",
       dpi = 1000,
       compression = "lzw")


################################################################################

# The following section is only for reference. 
# It creates plots of all disaggregated LPI trends.

# load rlpi package
library(rlpi)

ind_name <- "Infiles2" # name of directory to save LPI infiles

# Check if directory to store infiles exists. If not, create it.
if(!dir.exists(paste(ind_name, "/", sep=""))) {dir.create(paste(ind_name, "/", sep=""))} # create directory

# create boolean vectors of populations to include in each trend
T_Afrotropical_Aves <- LPI_trimmed$PopID %in% pop_list[[1]]
T_Afrotropical_Mammalia <- LPI_trimmed$PopID %in% pop_list[[2]]
T_Afrotropical_Herps <- LPI_trimmed$PopID %in% pop_list[[3]]

T_IndoPacific_Aves <- LPI_trimmed$PopID %in% pop_list[[5]]
T_IndoPacific_Mammalia <- LPI_trimmed$PopID %in% pop_list[[6]]
T_IndoPacific_Herps <- LPI_trimmed$PopID %in% pop_list[[7]]

T_Palearctic_Aves <- LPI_trimmed$PopID %in% pop_list[[9]]
T_Palearctic_Mammalia <- LPI_trimmed$PopID %in% pop_list[[10]]
T_Palearctic_Herps <- LPI_trimmed$PopID %in% pop_list[[11]]

T_Neotropical_Aves <- LPI_trimmed$PopID %in% pop_list[[13]]
T_Neotropical_Mammalia <- LPI_trimmed$PopID %in% pop_list[[14]]
T_Neotropical_Herps <- LPI_trimmed$PopID %in% pop_list[[15]]

T_Nearctic_Aves <- LPI_trimmed$PopID %in% pop_list[[17]]
T_Nearctic_Mammalia <- LPI_trimmed$PopID %in% pop_list[[18]]
T_Nearctic_Herps <- LPI_trimmed$PopID %in% pop_list[[19]]

T_Antarctic_Aves <- LPI_trimmed$PopID %in% pop_list[[21]]
T_Antarctic_Mammalia <- LPI_trimmed$PopID %in% pop_list[[22]]

fw_Afrotropical_Aves <- LPI_trimmed$PopID %in% pop_list[[25]]
fw_Afrotropical_Mammalia <- LPI_trimmed$PopID %in% pop_list[[26]]
fw_Afrotropical_Herps <- LPI_trimmed$PopID %in% pop_list[[27]]
fw_Afrotropical_Fish <- LPI_trimmed$PopID %in% pop_list[[28]]

fw_IndoPacific_Aves <- LPI_trimmed$PopID %in% pop_list[[29]]
fw_IndoPacific_Mammalia <- LPI_trimmed$PopID %in% pop_list[[30]]
fw_IndoPacific_Herps <- LPI_trimmed$PopID %in% pop_list[[31]]
fw_IndoPacific_Fish <- LPI_trimmed$PopID %in% pop_list[[32]]

fw_Palearctic_Aves <- LPI_trimmed$PopID %in% pop_list[[33]]
fw_Palearctic_Mammalia <- LPI_trimmed$PopID %in% pop_list[[34]]
fw_Palearctic_Herps <- LPI_trimmed$PopID %in% pop_list[[35]]
fw_Palearctic_Fish <- LPI_trimmed$PopID %in% pop_list[[36]]

fw_Neotropical_Aves <- LPI_trimmed$PopID %in% pop_list[[37]]
fw_Neotropical_Mammalia <- LPI_trimmed$PopID %in% pop_list[[38]]
fw_Neotropical_Herps <- LPI_trimmed$PopID %in% pop_list[[39]]
fw_Neotropical_Fish <- LPI_trimmed$PopID %in% pop_list[[40]]

fw_Nearctic_Aves <- LPI_trimmed$PopID %in% pop_list[[41]]
fw_Nearctic_Mammalia <- LPI_trimmed$PopID %in% pop_list[[42]]
fw_Nearctic_Herps <- LPI_trimmed$PopID %in% pop_list[[43]]
fw_Nearctic_Fish <- LPI_trimmed$PopID %in% pop_list[[44]]

m_AtNoTemp_Aves <- LPI_trimmed$PopID %in% pop_list[[49]]
m_AtNoTemp_Mammalia <- LPI_trimmed$PopID %in% pop_list[[50]]
m_AtNoTemp_Herps <- LPI_trimmed$PopID %in% pop_list[[51]]
m_AtNoTemp_Fish <- LPI_trimmed$PopID %in% pop_list[[52]]

m_AtTrSub_Aves <- LPI_trimmed$PopID %in% pop_list[[53]]
m_AtTrSub_Mammalia <- LPI_trimmed$PopID %in% pop_list[[54]]
m_AtTrSub_Herps <- LPI_trimmed$PopID %in% pop_list[[55]]
m_AtTrSub_Fish <- LPI_trimmed$PopID %in% pop_list[[56]]

m_Arctic_Aves <- LPI_trimmed$PopID %in% pop_list[[57]]
m_Arctic_Mammalia <- LPI_trimmed$PopID %in% pop_list[[58]]
m_Arctic_Fish <- LPI_trimmed$PopID %in% pop_list[[60]]

m_SoTeAnt_Aves <- LPI_trimmed$PopID %in% pop_list[[61]]
m_SoTeAnt_Mammalia <- LPI_trimmed$PopID %in% pop_list[[62]]
m_SoTeAnt_Herps <- LPI_trimmed$PopID %in% pop_list[[63]]
m_SoTeAnt_Fish <- LPI_trimmed$PopID %in% pop_list[[64]]

m_TroSubIndo_Aves <- LPI_trimmed$PopID %in% pop_list[[65]]
m_TroSubIndo_Mammalia <- LPI_trimmed$PopID %in% pop_list[[66]]
m_TroSubIndo_Herps <- LPI_trimmed$PopID %in% pop_list[[67]]
m_TroSubIndo_Fish <- LPI_trimmed$PopID %in% pop_list[[68]]

m_PaNoTemp_Aves <- LPI_trimmed$PopID %in% pop_list[[69]]
m_PaNoTemp_Mammalia <- LPI_trimmed$PopID %in% pop_list[[70]]
m_PaNoTemp_Herps <- LPI_trimmed$PopID %in% pop_list[[71]]
m_PaNoTemp_Fish <- LPI_trimmed$PopID %in% pop_list[[72]]


# create infiles from the boolean vectors.
# (for more info on infiles, see the rlpi documentation)
T_Afrotropical_Aves_infile <- create_infile(LPI_full, index_vector=T_Afrotropical_Aves, name=paste(ind_name, "/T_Afrotropical_Aves", sep=""))
T_Afrotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Afrotropical_Mammalia, name=paste(ind_name, "/T_Afrotropical_Mammalia", sep=""))
T_Afrotropical_Herps_infile <- create_infile(LPI_full, index_vector=T_Afrotropical_Herps, name=paste(ind_name, "/T_Afrotropical_Herps", sep=""))

T_IndoPacific_Aves_infile <- create_infile(LPI_full, index_vector=T_IndoPacific_Aves, name=paste(ind_name, "/T_IndoPacific_Aves", sep=""))
T_IndoPacific_Mammalia_infile <- create_infile(LPI_full, index_vector=T_IndoPacific_Mammalia, name=paste(ind_name, "/T_IndoPacific_Mammalia", sep=""))
T_IndoPacific_Herps_infile <- create_infile(LPI_full, index_vector=T_IndoPacific_Herps, name=paste(ind_name, "/T_IndoPacific_Herps", sep=""))

T_Palearctic_Aves_infile <- create_infile(LPI_full, index_vector=T_Palearctic_Aves, name=paste(ind_name, "/T_Palearctic_Aves", sep=""))
T_Palearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Palearctic_Mammalia, name=paste(ind_name, "/T_Palearctic_Mammalia", sep=""))
T_Palearctic_Herps_infile <- create_infile(LPI_full, index_vector=T_Palearctic_Herps, name=paste(ind_name, "/T_Palearctic_Herps", sep=""))

T_Neotropical_Aves_infile <- create_infile(LPI_full, index_vector=T_Neotropical_Aves, name=paste(ind_name, "/T_Neotropical_Aves", sep=""))
T_Neotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Neotropical_Mammalia, name=paste(ind_name, "/T_Neotropical_Mammalia", sep=""))
T_Neotropical_Herps_infile <- create_infile(LPI_full, index_vector=T_Neotropical_Herps, name=paste(ind_name, "/T_Neotropical_Herps", sep=""))

T_Nearctic_Aves_infile <- create_infile(LPI_full, index_vector=T_Nearctic_Aves, name=paste(ind_name, "/T_Nearctic_Aves", sep=""))
T_Nearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Nearctic_Mammalia, name=paste(ind_name, "/T_Nearctic_Mammalia", sep=""))
T_Nearctic_Herps_infile <- create_infile(LPI_full, index_vector=T_Nearctic_Herps, name=paste(ind_name, "/T_Nearctic_Herps", sep=""))

T_Antarctic_Aves_infile <- create_infile(LPI_full, index_vector=T_Antarctic_Aves, name=paste(ind_name, "/T_Antarctic_Aves", sep=""))
T_Antarctic_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Antarctic_Mammalia, name=paste(ind_name, "/T_Antarctic_Mammalia", sep=""))

fw_Afrotropical_Aves_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Aves, name=paste(ind_name, "/fw_Afrotropical_Aves", sep=""))
fw_Afrotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Mammalia, name=paste(ind_name, "/fw_Afrotropical_Mammalia", sep=""))
fw_Afrotropical_Herps_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Herps, name=paste(ind_name, "/fw_Afrotropical_Herps", sep=""))
fw_Afrotropical_Fish_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Fish, name=paste(ind_name, "/fw_Afrotropical_Fish", sep=""))

fw_IndoPacific_Aves_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Aves, name=paste(ind_name, "/fw_IndoPacific_Aves", sep=""))
fw_IndoPacific_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Mammalia, name=paste(ind_name, "/fw_IndoPacific_Mammalia", sep=""))
fw_IndoPacific_Herps_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Herps, name=paste(ind_name, "/fw_IndoPacific_Herps", sep=""))
fw_IndoPacific_Fish_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Fish, name=paste(ind_name, "/fw_IndoPacific_Fish", sep=""))

fw_Palearctic_Aves_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Aves, name=paste(ind_name, "/fw_Palearctic_Aves", sep=""))
fw_Palearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Mammalia, name=paste(ind_name, "/fw_Palearctic_Mammalia", sep=""))
fw_Palearctic_Herps_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Herps, name=paste(ind_name, "/fw_Palearctic_Herps", sep=""))
fw_Palearctic_Fish_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Fish, name=paste(ind_name, "/fw_Palearctic_Fish", sep=""))

fw_Neotropical_Aves_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Aves, name=paste(ind_name, "/fw_Neotropical_Aves", sep=""))
fw_Neotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Mammalia, name=paste(ind_name, "/fw_Neotropical_Mammalia", sep=""))
fw_Neotropical_Herps_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Herps, name=paste(ind_name, "/fw_Neotropical_Herps", sep=""))
fw_Neotropical_Fish_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Fish, name=paste(ind_name, "/fw_Neotropical_Fish", sep=""))

fw_Nearctic_Aves_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Aves, name=paste(ind_name, "/fw_Nearctic_Aves", sep=""))
fw_Nearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Mammalia, name=paste(ind_name, "/fw_Nearctic_Mammalia", sep=""))
fw_Nearctic_Herps_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Herps, name=paste(ind_name, "/fw_Nearctic_Herps", sep=""))
fw_Nearctic_Fish_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Fish, name=paste(ind_name, "/fw_Nearctic_Fish", sep=""))

m_AtNoTemp_Aves_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Aves, name=paste(ind_name, "/m_AtNoTemp_Aves", sep=""))
m_AtNoTemp_Mammalia_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Mammalia, name=paste(ind_name, "/m_AtNoTemp_Mammalia", sep=""))
m_AtNoTemp_Herps_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Herps, name=paste(ind_name, "/m_AtNoTemp_Herps", sep=""))
m_AtNoTemp_Fish_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Fish, name=paste(ind_name, "/m_AtNoTemp_Fish", sep=""))

m_AtTrSub_Aves_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Aves, name=paste(ind_name, "/m_AtTrSub_Aves", sep=""))
m_AtTrSub_Mammalia_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Mammalia, name=paste(ind_name, "/m_AtTrSub_Mammalia", sep=""))
m_AtTrSub_Herps_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Herps, name=paste(ind_name, "/m_AtTrSub_Herps", sep=""))
m_AtTrSub_Fish_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Fish, name=paste(ind_name, "/m_AtTrSub_Fish", sep=""))

m_Arctic_Aves_infile <- create_infile(LPI_full, index_vector=m_Arctic_Aves, name=paste(ind_name, "/m_Arctic_Aves", sep=""))
m_Arctic_Mammalia_infile <- create_infile(LPI_full, index_vector=m_Arctic_Mammalia, name=paste(ind_name, "/m_Arctic_Mammalia", sep=""))
m_Arctic_Fish_infile <- create_infile(LPI_full, index_vector=m_Arctic_Fish, name=paste(ind_name, "/m_Arctic_Fish", sep=""))

m_SoTeAnt_Aves_infile <- create_infile(LPI_full, index_vector=m_SoTeAnt_Aves, name=paste(ind_name, "/m_SoTeAnt_Aves", sep=""))
m_SoTeAnt_Mammalia_infile <- create_infile(LPI_full, index_vector=m_SoTeAnt_Mammalia, name=paste(ind_name, "/m_SoTeAnt_Mammalia", sep=""))
m_SoTeAnt_Herps_infile <- create_infile(LPI_full, index_vector=m_SoTeAnt_Herps, name=paste(ind_name, "/m_SoTeAnt_Herps", sep=""))
m_SoTeAnt_Fish_infile <- create_infile(LPI_full, index_vector=m_SoTeAnt_Fish, name=paste(ind_name, "/m_SoTeAnt_Fish", sep=""))

m_TroSubIndo_Aves_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Aves, name=paste(ind_name, "/m_TroSubIndo_Aves", sep=""))
m_TroSubIndo_Mammalia_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Mammalia, name=paste(ind_name, "/m_TroSubIndo_Mammalia", sep=""))
m_TroSubIndo_Herps_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Herps, name=paste(ind_name, "/m_TroSubIndo_Herps", sep=""))
m_TroSubIndo_Fish_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Fish, name=paste(ind_name, "/m_TroSubIndo_Fish", sep=""))

m_PaNoTemp_Aves_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Aves, name=paste(ind_name, "/m_PaNoTemp_Aves", sep=""))
m_PaNoTemp_Mammalia_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Mammalia, name=paste(ind_name, "/m_PaNoTemp_Mammalia", sep=""))
m_PaNoTemp_Herps_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Herps, name=paste(ind_name, "/m_PaNoTemp_Herps", sep=""))
m_PaNoTemp_Fish_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Fish, name=paste(ind_name, "/m_PaNoTemp_Fish", sep=""))


# create indices and plot

max_plot_year <- as.integer(max(year_names)) # set maximum year for LPI plots as the final year of dataset

bootstraps <- 100 # increasing this value may improve confidence intervals, but each plot will take longer to generate. 
                  # 100 is the default of the rlpi package

T_Afrotropical_Aves_lpi  <- LPIMain(paste(ind_name, "/T_Afrotropical_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Afrotropical_Mammalia_lpi  <- LPIMain(paste(ind_name, "/T_Afrotropical_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Afrotropical_Herps_lpi  <- LPIMain(paste(ind_name, "/T_Afrotropical_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_IndoPacific_Aves_lpi  <- LPIMain(paste(ind_name, "/T_IndoPacific_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_IndoPacific_Mammalia_lpi  <- LPIMain(paste(ind_name, "/T_IndoPacific_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_IndoPacific_Herps_lpi  <- LPIMain(paste(ind_name, "/T_IndoPacific_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_Palearctic_Aves_lpi  <- LPIMain(paste(ind_name, "/T_Palearctic_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Palearctic_Mammalia_lpi  <- LPIMain(paste(ind_name, "/T_Palearctic_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Palearctic_Herps_lpi  <- LPIMain(paste(ind_name, "/T_Palearctic_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_Neotropical_Aves_lpi  <- LPIMain(paste(ind_name, "/T_Neotropical_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Neotropical_Mammalia_lpi  <- LPIMain(paste(ind_name, "/T_Neotropical_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Neotropical_Herps_lpi  <- LPIMain(paste(ind_name, "/T_Neotropical_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_Nearctic_Aves_lpi  <- LPIMain(paste(ind_name, "/T_Nearctic_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Nearctic_Mammalia_lpi  <- LPIMain(paste(ind_name, "/T_Nearctic_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Nearctic_Herps_lpi  <- LPIMain(paste(ind_name, "/T_Nearctic_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_Antarctic_Aves_lpi <- LPIMain(paste(ind_name, "/T_Antarctic_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Antarctic_Mammalia_lpi <- LPIMain(paste(ind_name, "/T_Antarctic_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_Afrotropical_Aves_lpi  <- LPIMain(paste(ind_name, "/fw_Afrotropical_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Afrotropical_Mammalia_lpi  <- LPIMain(paste(ind_name, "/fw_Afrotropical_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Afrotropical_Herps_lpi  <- LPIMain(paste(ind_name, "/fw_Afrotropical_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Afrotropical_fish_lpi  <- LPIMain(paste(ind_name, "/fw_Afrotropical_Fish_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_IndoPacific_Aves_lpi  <- LPIMain(paste(ind_name, "/fw_IndoPacific_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_IndoPacific_Mammalia_lpi  <- LPIMain(paste(ind_name, "/fw_IndoPacific_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_IndoPacific_Herps_lpi  <- LPIMain(paste(ind_name, "/fw_IndoPacific_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_IndoPacific_fish_lpi  <- LPIMain(paste(ind_name, "/fw_IndoPacific_Fish_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_Palearctic_Aves_lpi  <- LPIMain(paste(ind_name, "/fw_Palearctic_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Palearctic_Mammalia_lpi  <- LPIMain(paste(ind_name, "/fw_Palearctic_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Palearctic_Herps_lpi  <- LPIMain(paste(ind_name, "/fw_Palearctic_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Palearctic_fish_lpi  <- LPIMain(paste(ind_name, "/fw_Palearctic_Fish_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_Neotropical_Aves_lpi  <- LPIMain(paste(ind_name, "/fw_Neotropical_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Neotropical_Mammalia_lpi  <- LPIMain(paste(ind_name, "/fw_Neotropical_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Neotropical_Herps_lpi  <- LPIMain(paste(ind_name, "/fw_Neotropical_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Neotropical_fish_lpi  <- LPIMain(paste(ind_name, "/fw_Neotropical_Fish_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_Nearctic_Aves_lpi  <- LPIMain(paste(ind_name, "/fw_Nearctic_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Nearctic_Mammalia_lpi  <- LPIMain(paste(ind_name, "/fw_Nearctic_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Nearctic_Herps_lpi  <- LPIMain(paste(ind_name, "/fw_Nearctic_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Nearctic_fish_lpi  <- LPIMain(paste(ind_name, "/fw_Nearctic_Fish_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_AtNoTemp_Aves <- LPIMain(paste(ind_name, "/m_AtNoTemp_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtNoTemp_Mammalia <- LPIMain(paste(ind_name, "/m_AtNoTemp_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtNoTemp_Herps <- LPIMain(paste(ind_name, "/m_AtNoTemp_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtNoTemp_Fish <- LPIMain(paste(ind_name, "/m_AtNoTemp_Fish_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_AtTrSub_Aves <- LPIMain(paste(ind_name, "/m_AtTrSub_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtTrSub_Mammalia <- LPIMain(paste(ind_name, "/m_AtTrSub_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtTrSub_Herps <- LPIMain(paste(ind_name, "/m_AtTrSub_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtTrSub_Fish <- LPIMain(paste(ind_name, "/m_AtTrSub_Fish_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_Arctic_Aves <- LPIMain(paste(ind_name, "/m_Arctic_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_Arctic_Mammalia <- LPIMain(paste(ind_name, "/m_Arctic_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_Arctic_Fish <- LPIMain(paste(ind_name, "/m_Arctic_Fish_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_SoTeAnt_Aves <- LPIMain(paste(ind_name, "/m_SoTeAnt_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_SoTeAnt_Mammalia <- LPIMain(paste(ind_name, "/m_SoTeAnt_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_SoTeAnt_Herps <- LPIMain(paste(ind_name, "/m_SoTeAnt_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_SoTeAnt_Fish <- LPIMain(paste(ind_name, "/m_SoTeAnt_Fish_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_TroSubIndo_Aves <- LPIMain(paste(ind_name, "/m_TroSubIndo_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_TroSubIndo_Mammalia <- LPIMain(paste(ind_name, "/m_TroSubIndo_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_TroSubIndo_Herps <- LPIMain(paste(ind_name, "/m_TroSubIndo_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_TroSubIndo_Fish <- LPIMain(paste(ind_name, "/m_TroSubIndo_Fish_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_PaNoTemp_Aves <- LPIMain(paste(ind_name, "/m_PaNoTemp_Aves_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_PaNoTemp_Mammalia <- LPIMain(paste(ind_name, "/m_PaNoTemp_Mammalia_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_PaNoTemp_Herps <- LPIMain(paste(ind_name, "/m_PaNoTemp_Herps_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_PaNoTemp_Fish <- LPIMain(paste(ind_name, "/m_PaNoTemp_Fish_infile.txt", sep=""), REF_YEAR = 1970, PLOT_MAX = max_plot_year, BOOT_STRAP_SIZE = bootstraps, force_recalculation=1, use_weightings=0, use_weightings_B=0)

