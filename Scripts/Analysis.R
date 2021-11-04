# libraries
library(tidyr)
library(dplyr)

# load functions
source("Scripts/Analysis_Functions.R")

# load full LPI dataset
LPI_full <- read.csv(file="Data/LPD_output_20201116.csv", sep=",", stringsAsFactors=FALSE)

# fix ID column name
colnames(LPI_full)[1] <- "ID"

# remove X from years in column names
colnames(LPI_full) <- gsub("X", "", colnames(LPI_full))

# replace "NULL" values with NA
#LPI_full2 <- apply(LPI_full[,65:134], 2, as.numeric)
#LPI_full[,65:134] <- LPI_full2
#remove(LPI_full2)

# get group data
#amphibians <- LPI_full[LPI_full$Class=="Amphibia",]
#reptiles <- LPI_full[LPI_full$Class=="Reptilia",]
#herps <- LPI_full[LPI_full$Class=="Reptilia" | LPI_full$Class=="Amphibia",]
#freshwater <- LPI_full[LPI_full$FWRealm!="NULL",]


# convert to long format
test <- LPI_full %>% 
  pivot_longer('1950':'2019',
               names_to = "year",
               values_to = "value")

# convert years to integer
test$year <- as.integer(test$year)

#convert values to numeric (converts "NULL" to NA)
test$value <- as.numeric(test$value)

test2 <- test %>% 
  group_by(ID) %>% # treat each population separately
  filter(!is.na(value)) %>% # remove NA values
  mutate(start_year = min(year),
         last_year = max(year),
         ts_length = last_year - start_year+1, # calculate length of time series
         n_obs = sum(year > 0)) %>% # count the number of years with values
  select(ID, start_year, last_year, ts_length, n_obs) %>% # remove all extra columns
  distinct() %>% # remove duplicate populations
  ungroup() # remove grouping

# add the new columns to original data
test3 <- left_join(LPI_full, test2, by = "ID")

# select the data you want
amphibians <- test3[test3$Class=="Amphibia",]
reptiles <- test3[test3$Class=="Reptilia",]
herps <- test3[test3$Class=="Reptilia" | test3$Class=="Amphibia",]
freshwater <- test3[test3$FWRealm!="NULL",]



b <- amphibians %>% group_by(Method, Units, Class) %>% summarise(n_methods = n(),
                                                           mean_n_obs = mean(n_obs),
                                                           mean_ts_length = mean(ts_length))

# create subset with only count data
group_vals <- amphibians[,which(colnames(amphibians) %in% (1970:2019))]

m_colnames <- as.character(1970:2019)
collength <- length(m_colnames)

# get rid of time series shorter than 2
group_culled <- cull_fn(group_vals, count_thres = 2, min_ts_length = 2, collength)

# log-linear interpolate populations
group_complete <- complete_time_series(group_culled, collength, m_colnames, calcsd=TRUE)

# get mean and standard deviation of the mean growth rate
gr.stats <- growth_rate_byrow(group_complete[,1:collength], model=TRUE)

max.absolute.gr <- vector()
for (i in 1:length(gr.stats[[1]])) {
  
  max.absolute.gr[i] <- max(abs(c(log(gr.stats[[1]][i]), log(gr.stats[[2]][i]))))
  
}

ts.lengths <- tslength_byrow(group_culled)

ts.numobs <- tsnumobs_byrow(group_culled)

gr_vs_length <- data.frame(Length = ts.lengths, Num_Obs = ts.numobs, Max_GR = gr.stats[[1]], Min_GR = gr.stats[[2]], Mean_GR = gr.stats[[3]], Max_Abs_GR = max.absolute.gr)

amphibians_gr_df <- gr_vs_length

