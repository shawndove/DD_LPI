LPI_pfull <- read.csv(file="Data/LPR2020data_public.csv", sep=",", stringsAsFactors=FALSE)
colnames(LPI_full)[1] <- "ID"
colnames(LPI_full) <- gsub("X", "", colnames(LPI_full))

firstyear <- 1950 # first year of data
startyear <- 1970 # first year of data to use for index
endyear <- 2019 # final year of data
m_colnames <- as.character(startyear:endyear) # index years/column names
m_colnames2 <- as.character(firstyear:endyear) # index years/column names
c <- length(m_colnames) # number of years/columns
c2 <- length(m_colnames2) # number of years/columns

# create subset with only count data
LPI_trimmed <- LPI_full[,which(colnames(LPI_full) %in% (firstyear:endyear))]

# convert to numeric data
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

# select all copies of the populations listed in the sample
group_data <- LPI_trimmed[LPI_trimmed$PopID %in% pop_list[[1]],]

# remove all populations with less than 2 data points
grp_data_culled <- cull_fn(group_data, 2, 2, c2)

# log-linear interpolate populations shorter than 6 years
# populations 6 years or longer will be interpolated with a GAM in the next step
new.grp_data <- complete_time_series(grp_data_culled, c2, m_colnames2, lambda=TRUE)

# GAM the population indices
gam_popmat <- pop_gam_fn(new.grp_data, c2, m_colnames2, n=NA, lambda=TRUE, resample=FALSE, quality=TRUE)

# log-linear interpolate populations that failed gam quality test
new.grp_data2 <- complete_time_series(gam_popmat, c2, m_colnames2, lambda=TRUE)

new.grp_data3 <- new.grp_data2[,21:ncol(new.grp_data2)]

# create species indices from the population indices
full_spec_index <- species_index_fn(new.grp_data3, c, n=NA, n_boot=NA, lambda=TRUE, resample=FALSE)

# create multi species indices from the population indices
msi_full <- group_index_fn(full_spec_index, c, m_colnames, n=NA, n_boot=NA)

# create multi species confidence intervals
msi_ci_full <- ci_fn(full_spec_index, c, m_colnames, lambda=TRUE, grouplevel=2)
