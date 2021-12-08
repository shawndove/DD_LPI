# libraries ----

library(rlpi)
library(dplyr)
library(mgcv)
library(matrixStats)
library(ggplot2)
library(gridExtra)

# load functions ----
source("Scripts/Functions.R")

# create groupings ----

tax_group <- list(Birds=c("Aves"), 
               Mammals=c("Mammalia"), 
               Herps=c("Reptilia", "Amphibia"), 
               Fish=c("Actinopterygii", "Elasmobranchii", "Sarcopterygii", "Cephalaspidomorphi", "Holocephali", "Myxini", "Chondrichthyes"))

# use this alternative for the full LPI (including private), as there are differences in the fish class names
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

sys_grp <- list(Terrest=c("Terrestrial"),
                 FW=c("Freshwater"),
                 Marine=c("Marine"))

# create an ID key for marine realms
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

sys_realm_list <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3)

# create list of weightings for taxonomic groups within realms (copied from McRae et al., 2017)
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

# create list of weightings for realms within systems (copied from McRae et al., 2017)
Weightingsr_list <- c(0.189738, 0.292168, 0.116431, 0.321132, 0.061683, NA, 
                      0.211701, 0.225576, 0.123314, 0.365550, 0.060853, NA, 
                      0.146489, 0.214706, 0.014541, 0.099685, 0.456553, 0.068026)

# load public LPI data ----

# only public data
LPI_full <- read.csv(file="Data/LPR2020data_public.csv", sep=",", stringsAsFactors=FALSE)
LPI_test <- read_excel("Data/LPR2020data_public.csv")

# all data, including private
LPI_full2 <- read.csv(file="Data/LPD_output_20201116.csv", sep=",", stringsAsFactors=FALSE)
# fix ID column name
colnames(LPI_full)[1] <- "ID"
# remove everything added after May 2016
#LPI_full <- LPI_full[LPI_full$ID < 18330,]
#LPI_full <- LPI_full[LPI_full$Confidential==0,]
# remove everything added after March 3, 2015
#LPI_full <- LPI_full[LPI_full$ID < 17527,]

# remove X from years in column names
colnames(LPI_full) <- gsub("X", "", colnames(LPI_full))

firstyear <- 1950 # first year of data
startyear <- 1970 # first year of data to use for index
endyear <- 2015 # final year of data
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
m_SoTeAnt_Fish <- LPI_trimmed$PopID %in% pop_list[[64]]

m_TroSubIndo_Aves <- LPI_trimmed$PopID %in% pop_list[[65]]
m_TroSubIndo_Mammalia <- LPI_trimmed$PopID %in% pop_list[[66]]
m_TroSubIndo_Herps <- LPI_trimmed$PopID %in% pop_list[[67]]
m_TroSubIndo_Fish <- LPI_trimmed$PopID %in% pop_list[[68]]

m_PaNoTemp_Aves <- LPI_trimmed$PopID %in% pop_list[[69]]
m_PaNoTemp_Mammalia <- LPI_trimmed$PopID %in% pop_list[[70]]
m_PaNoTemp_Herps <- LPI_trimmed$PopID %in% pop_list[[71]]
m_PaNoTemp_Fish <- LPI_trimmed$PopID %in% pop_list[[72]]

# add X back to years in column names (to make create_infile work properly)
colnames(LPI_full)[65:134] <- paste("X", colnames(LPI_full)[65:134], sep="")

# create infiles
T_Afrotropical_Aves_infile <- create_infile(LPI_full, index_vector=T_Afrotropical_Aves, name="Infiles/T_Afrotropical_Aves", end_col_name = "X2019")
T_Afrotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Afrotropical_Mammalia, name="Infiles/T_Afrotropical_Mammalia", end_col_name = "X2019")
T_Afrotropical_Herps_infile <- create_infile(LPI_full, index_vector=T_Afrotropical_Herps, name="Infiles/T_Afrotropical_Herps", end_col_name = "X2019")

T_IndoPacific_Aves_infile <- create_infile(LPI_full, index_vector=T_IndoPacific_Aves, name="Infiles/T_IndoPacific_Aves", end_col_name = "X2019")
T_IndoPacific_Mammalia_infile <- create_infile(LPI_full, index_vector=T_IndoPacific_Mammalia, name="Infiles/T_IndoPacific_Mammalia", end_col_name = "X2019")
T_IndoPacific_Herps_infile <- create_infile(LPI_full, index_vector=T_IndoPacific_Herps, name="Infiles/T_IndoPacific_Herps", end_col_name = "X2019")

T_Palearctic_Aves_infile <- create_infile(LPI_full, index_vector=T_Palearctic_Aves, name="Infiles/T_Palearctic_Aves", end_col_name = "X2019")
T_Palearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Palearctic_Mammalia, name="Infiles/T_Palearctic_Mammalia", end_col_name = "X2019")
T_Palearctic_Herps_infile <- create_infile(LPI_full, index_vector=T_Palearctic_Herps, name="Infiles/T_Palearctic_Herps", end_col_name = "X2019")

T_Neotropical_Aves_infile <- create_infile(LPI_full, index_vector=T_Neotropical_Aves, name="Infiles/T_Neotropical_Aves", end_col_name = "X2019")
T_Neotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Neotropical_Mammalia, name="Infiles/T_Neotropical_Mammalia", end_col_name = "X2019")
T_Neotropical_Herps_infile <- create_infile(LPI_full, index_vector=T_Neotropical_Herps, name="Infiles/T_Neotropical_Herps", end_col_name = "X2019")

T_Nearctic_Aves_infile <- create_infile(LPI_full, index_vector=T_Nearctic_Aves, name="Infiles/T_Nearctic_Aves", end_col_name = "X2019")
T_Nearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Nearctic_Mammalia, name="Infiles/T_Nearctic_Mammalia", end_col_name = "X2019")
T_Nearctic_Herps_infile <- create_infile(LPI_full, index_vector=T_Nearctic_Herps, name="Infiles/T_Nearctic_Herps", end_col_name = "X2019")

T_Antarctic_Aves_infile <- create_infile(LPI_full, index_vector=T_Antarctic_Aves, name="T_Antarctic_Aves", end_col_name = "X2019")
T_Antarctic_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Antarctic_Mammalia, name="T_Antarctic_Mammalia", end_col_name = "X2019")

fw_Afrotropical_Aves_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Aves, name="Infiles/fw_Afrotropical_Aves", end_col_name = "X2019")
fw_Afrotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Mammalia, name="Infiles/fw_Afrotropical_Mammalia", end_col_name = "X2019")
fw_Afrotropical_Herps_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Herps, name="Infiles/fw_Afrotropical_Herps", end_col_name = "X2019")
fw_Afrotropical_Fish_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Fish, name="Infiles/fw_Afrotropical_Fish", end_col_name = "X2019")

fw_IndoPacific_Aves_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Aves, name="Infiles/fw_IndoPacific_Aves", end_col_name = "X2019")
fw_IndoPacific_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Mammalia, name="Infiles/fw_IndoPacific_Mammalia", end_col_name = "X2019")
fw_IndoPacific_Herps_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Herps, name="Infiles/fw_IndoPacific_Herps", end_col_name = "X2019")
fw_IndoPacific_Fish_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Fish, name="Infiles/fw_IndoPacific_Fish", end_col_name = "X2019")

fw_Palearctic_Aves_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Aves, name="Infiles/fw_Palearctic_Aves", end_col_name = "X2019")
fw_Palearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Mammalia, name="Infiles/fw_Palearctic_Mammalia", end_col_name = "X2019")
fw_Palearctic_Herps_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Herps, name="Infiles/fw_Palearctic_Herps", end_col_name = "X2019")
fw_Palearctic_Fish_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Fish, name="Infiles/fw_Palearctic_Fish", end_col_name = "X2019")

fw_Neotropical_Aves_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Aves, name="Infiles/fw_Neotropical_Aves", end_col_name = "X2019")
fw_Neotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Mammalia, name="Infiles/fw_Neotropical_Mammalia", end_col_name = "X2019")
fw_Neotropical_Herps_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Herps, name="Infiles/fw_Neotropical_Herps", end_col_name = "X2019")
fw_Neotropical_Fish_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Fish, name="Infiles/fw_Neotropical_Fish", end_col_name = "X2019")

fw_Nearctic_Aves_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Aves, name="Infiles/fw_Nearctic_Aves", end_col_name = "X2019")
fw_Nearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Mammalia, name="Infiles/fw_Nearctic_Mammalia", end_col_name = "X2019")
fw_Nearctic_Herps_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Herps, name="Infiles/fw_Nearctic_Herps", end_col_name = "X2019")
fw_Nearctic_Fish_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Fish, name="Infiles/fw_Nearctic_Fish", end_col_name = "X2019")

m_AtNoTemp_Aves_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Aves, name="Infiles/m_AtNoTemp_Aves", end_col_name = "X2019")
m_AtNoTemp_Mammalia_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Mammalia, name="Infiles/m_AtNoTemp_Mammalia", end_col_name = "X2019")
m_AtNoTemp_Herps_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Herps, name="Infiles/m_AtNoTemp_Herps", end_col_name = "X2019")
m_AtNoTemp_Fish_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Fish, name="Infiles/m_AtNoTemp_Fish", end_col_name = "X2019")

m_AtTrSub_Aves_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Aves, name="Infiles/m_AtTrSub_Aves", end_col_name = "X2019")
m_AtTrSub_Mammalia_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Mammalia, name="Infiles/m_AtTrSub_Mammalia", end_col_name = "X2019")
m_AtTrSub_Herps_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Herps, name="Infiles/m_AtTrSub_Herps", end_col_name = "X2019")
m_AtTrSub_Fish_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Fish, name="Infiles/m_AtTrSub_Fish", end_col_name = "X2019")

m_Arctic_Aves_infile <- create_infile(LPI_full, index_vector=m_Arctic_Aves, name="Infiles/m_Arctic_Aves", end_col_name = "X2019")
m_Arctic_Mammalia_infile <- create_infile(LPI_full, index_vector=m_Arctic_Mammalia, name="Infiles/m_Arctic_Mammalia", end_col_name = "X2019")
m_Arctic_Fish_infile <- create_infile(LPI_full, index_vector=m_Arctic_Fish, name="Infiles/m_Arctic_Fish", end_col_name = "X2019")

m_SoTeAnt_Aves_infile <- create_infile(LPI_full, index_vector=m_SoTeAnt_Aves, name="Infiles/m_SoTeAnt_Aves", end_col_name = "X2019")
m_SoTeAnt_Mammalia_infile <- create_infile(LPI_full, index_vector=m_SoTeAnt_Mammalia, name="Infiles/m_SoTeAnt_Mammalia", end_col_name = "X2019")
m_SoTeAnt_Fish_infile <- create_infile(LPI_full, index_vector=m_SoTeAnt_Fish, name="Infiles/m_SoTeAnt_Fish", end_col_name = "X2019")

m_TroSubIndo_Aves_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Aves, name="Infiles/m_TroSubIndo_Aves", end_col_name = "X2019")
m_TroSubIndo_Mammalia_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Mammalia, name="Infiles/m_TroSubIndo_Mammalia", end_col_name = "X2019")
m_TroSubIndo_Herps_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Herps, name="Infiles/m_TroSubIndo_Herps", end_col_name = "X2019")
m_TroSubIndo_Fish_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Fish, name="Infiles/m_TroSubIndo_Fish", end_col_name = "X2019")

m_PaNoTemp_Aves_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Aves, name="Infiles/m_PaNoTemp_Aves", end_col_name = "X2019")
m_PaNoTemp_Mammalia_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Mammalia, name="Infiles/m_PaNoTemp_Mammalia", end_col_name = "X2019")
m_PaNoTemp_Herps_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Herps, name="Infiles/m_PaNoTemp_Herps", end_col_name = "X2019")
m_PaNoTemp_Fish_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Fish, name="Infiles/m_PaNoTemp_Fish", end_col_name = "X2019")

# create list of taxonomic groups for index list
model_tax_list <- rep(unique(tax_group_IDs), length(unique(sys_IDs)))

# create list of systems for index list
model_sys_list <- rep(unique(sys_IDs), each=length(unique(tax_group_IDs)))

model_pop_list <- list()
# select populations to form each group index
for (i in 1:length(model_tax_list)) {
  
  temp <- LPI_trimmed$PopID[which(LPI_trimmed$SysID==model_sys_list[i] & 
                                   LPI_trimmed$GrpID==model_tax_list[i])]
  
  if (length(temp)==0) {
    
    temp <- NA
    
  }
  
  model_pop_list[[i]] <- temp
  
}

# calculate mean growth rate for each LPI group
gr.stats.list <- list()
# loop over taxonomic groups
for (i in 1:length(model_pop_list)) {
  
  # select all copies of the populations listed in the sample
  group_data <- LPI_trimmed[LPI_trimmed$PopID %in% model_pop_list[[i]],]
  
  # remove all populations with less than 2 data points
  group_data_culled <- cull_fn(group_data, 2, 2, c2)
  
  # log-linear interpolate all pops
  grp_completed <- complete_time_series(group_data_culled, c2, m_colnames2, calcsd=TRUE)
  
  if (nrow(grp_completed) >=1) {
    
    # get mean and standard deviation of the mean growth rate
    gr.stats.list[[i]] <- growth_rate_calc_fn(grp_completed[,1:c2], model=TRUE)
    
  } else {
    
    gr.stats.list[[i]] <- NA
    
  }
  
}

# calculate current sample sizes in LPI groups
pop.size.list <- list()
# loop over groups
for (i in 1:length(model_pop_list)) {
  
  pop.size.list[[i]] <- length(model_pop_list[[i]])
  
}


# put group standard deviations in mean growth rates into the model to find minimum sample sizes
min_acc <- 200 # set minimum accuracy for chosen distance metric
samp.size.list <- list()
# loop over groups
for (i in 1:length(gr.stats.list)) {
  
  stdev <- as.numeric(gr.stats.list[[i]][2]) # get standard deviation of mean growth rate for group i
  
  if (is.na(stdev)) { # if standard deviation is NA...
    
    samp.size.list[[i]] <- NA # ... set sample size to NA
    
  } else { # otherwise...
    
    samp.size.list[[i]] <- exp((10.73 + (0.947*log(stdev)) - log(min_acc)) / 0.637) # calculate sample size
    
  }
  
}

# calculate percentage of min sample size in data for each group
# this is capped at 100% because the min will be used to subsample
samp.percent.list <- list()
# loop over groups
for (i in 1:length(gr.stats.list)) {
  
  temp.min <- samp.size.list[[i]]
  
  temp.actual <- pop.size.list[[i]]
  
  ratio <- temp.actual/temp.min
  
  if (is.na(ratio)) {
    
    samp.percent.list[[i]] <- NA
    
  } else if (ratio > 1) {
    
    samp.percent.list[[i]] <- 100
    
  } else {
    
    samp.percent.list[[i]] <- ratio * 100
    
  }
  
}


# calculate the actual sample size to take from each group
actual.samp.size.list <- list()
# loop over groups
for (i in 1:length(samp.percent.list)) {
  
  if (is.na(samp.percent.list[[i]])) {
    
    actual.samp.size.list[[i]] <- NA 
    
  } else if (samp.percent.list[[i]] == 100) {
    
    actual.samp.size.list[[i]] <- round(samp.size.list[[i]])
    
  } else if (samp.percent.list[[i]] < 100) {
    
    actual.samp.size.list[[i]] <- pop.size.list[[i]]
    
  }
  
}

# calculate expected accuracy for each group
exp.acc.list <- list()
# loop over groups
for (i in 1:length(gr.stats.list)) {
  
  stdev <- as.numeric(gr.stats.list[[i]][2])
  
  ss <- pop.size.list[[i]]
  
  exp.acc.list[[i]] <- exp(10.73 + (0.947*log(stdev)) + (-0.637*log(ss)))
  
}

# calculate percentage of min accuracy achieved for each group
acc.percent.list <- list()
# loop over groups
for (i in 1:length(gr.stats.list)) {
  
  temp <- exp.acc.list[[i]]
  
  ratio <- (1 / (temp / min_acc))
  
  if (is.na(ratio)) {
    
    acc.percent.list[[i]] <- NA
    
  } else if (ratio > 1) {
    
    acc.percent.list[[i]] <- 100
    
  } else {
    
    acc.percent.list[[i]] <- ratio * 100
    
  }
  
}

# calculate percent of total pops used in sampled LPI for each group
used.samp.ratio.list <- list()
# loop over groups
for (i in 1:length(samp.percent.list)) {
  
  used.samp.ratio.list[[i]] <- (actual.samp.size.list[[i]] / pop.size.list[[i]])
  
}

# convert sample percent list to a vector
samp.percent.vec <- do.call(rbind, samp.percent.list)

# convert accuracy percent list to a vector
acc.percent.vec <- do.call(rbind, acc.percent.list)

# calculate penalization weightings
penal.list <- list()
# loop over groups
for (i in unique(model_sys_list)) {
  
  # copy taxonomic weighting data for system
  weights_table <- acc.percent.vec[which(model_sys_list==i)]
  
  # adjust weights so they sum to 1
  penal.list[[i]] <- weights_table * (1 / sum(weights_table, na.rm=TRUE))
  
}

# calculate weightings
birds.ter <- mean(c(0.387, 0.376, 0.387, 0.433, 0.396))

birds.fw <- mean(c(0.192, 0.203, 0.107, 0.211, 0.176))

birds.mar <- mean(c(0.172867, 0.068635, 0.069353, 0.080916, 0.048714, 0.054261))

mammals.ter <- mean(c(0.197, 0.249, 0.127, 0.249, 0.172))

mammals.fw <- mean(c(0.009, 0.013, 0.010, 0.015, 0.008))

mammals.mar <- mean(c(0.035011, 0.009774, 0.006224, 0.025257, 0.004878, 0.022342))

reptiles.ter <- mean(c(0.414, 0.373, 0.484, 0.316, 0.431))

reptiles.fw <- mean(c(0.207, 0.217, 0.298, 0.179, 0.321))

reptiles.mar <- mean(c(0, 0.001303, 0.001630, 0.000935, 0.005505, 0.000957))

fish.ter <- NA

fish.fw <- mean(c(0.590, 0.565, 0.584, 0.592, 0.493))

fish.mar <- mean(c(0.792123, 0.920286, 0.922791, 0.892890, 0.940901, 0.922438))

ter.weight.temp <- c(birds.ter, mammals.ter, reptiles.ter, fish.ter)

fw.weight.temp <- c(birds.fw, mammals.fw, reptiles.fw, fish.fw)

mar.weight.temp <- c(birds.mar, mammals.mar, reptiles.mar, fish.mar)

# adjust weights so they sum to 1
ter.weight.adj <- ter.weight.temp * (1 / sum(ter.weight.temp, na.rm=TRUE))

fw.weight.adj <- fw.weight.temp * (1 / sum(fw.weight.temp, na.rm=TRUE))

mar.weight.adj <- mar.weight.temp * (1 / sum(mar.weight.temp, na.rm=TRUE))

model.sys.weights.adj <- c(ter.weight.adj, fw.weight.adj, mar.weight.adj)


# model infiles
model_sampled_pop_list <- list()
for (i in 1:length(model_pop_list)) {
  
  if (is.na(model_pop_list[[i]])) {

    model_sampled_pop_list[[i]] <- NA
    
    next
    
  }
  
  # select all copies of the populations listed in the sample
  group_data <- LPI_trimmed[LPI_trimmed$PopID %in% model_pop_list[[i]],]
  
  # create list of x randomly sampled populations from the species group, where x is sample size
  sample_pop_id_list <- sample(group_data$PopID, actual.samp.size.list[[i]])
  
  model_sampled_pop_list[[i]] <- sample_pop_id_list
  
}

final_model_weights_list <- list()
# loop over systems
for (i in unique(model_sys_list)) {
  
  # copy penalization weightings for system
  penal_weights <- penal.list[[i]]
  
  # copy LPI weightings for system
  weights_table <- model.sys.weights.adj[which(model_sys_list==i)]
  
  # adjust LPI weightings according to penalization weightings and ensure they sum to 1
  combined_weights <- weights_table * penal_weights
  weights_table_adj <- combined_weights * (1 / sum(combined_weights, na.rm=TRUE))
  
  final_model_weights_list[[i]] <- weights_table_adj
  
}


T_Aves_model <- LPI_trimmed$PopID %in% model_sampled_pop_list[[1]]
T_Mammalia_model <- LPI_trimmed$PopID %in% model_sampled_pop_list[[2]]
T_Herps_model <- LPI_trimmed$PopID %in% model_sampled_pop_list[[3]]
fw_Aves_model <- LPI_trimmed$PopID %in% model_sampled_pop_list[[5]]
fw_Mammalia_model <- LPI_trimmed$PopID %in% model_sampled_pop_list[[6]]
fw_Herps_model <- LPI_trimmed$PopID %in% model_sampled_pop_list[[7]]
fw_Fish_model <- LPI_trimmed$PopID %in% model_sampled_pop_list[[8]]
m_Aves_model <- LPI_trimmed$PopID %in% model_sampled_pop_list[[9]]
m_Mammalia_model <- LPI_trimmed$PopID %in% model_sampled_pop_list[[10]]
m_Herps_model <- LPI_trimmed$PopID %in% model_sampled_pop_list[[11]]
m_Fish_model <- LPI_trimmed$PopID %in% model_sampled_pop_list[[12]]

T_Aves_model_infile <- create_infile(LPI_full, index_vector=T_Aves_model, name="Infiles/T_Aves_model", end_col_name = "X2019")
T_Mammalia_model_infile <- create_infile(LPI_full, index_vector=T_Mammalia_model, name="Infiles/T_Mammalia_model", end_col_name = "X2019")
T_Herps_model_infile <- create_infile(LPI_full, index_vector=T_Herps_model, name="Infiles/T_Herps_model", end_col_name = "X2019")
fw_Aves_model_infile <- create_infile(LPI_full, index_vector=fw_Aves_model, name="Infiles/fw_Aves_model", end_col_name = "X2019")
fw_Mammalia_model_infile <- create_infile(LPI_full, index_vector=fw_Mammalia_model, name="Infiles/fw_Mammalia_model", end_col_name = "X2019")
fw_Herps_model_infile <- create_infile(LPI_full, index_vector=fw_Herps_model, name="Infiles/fw_Herps_model", end_col_name = "X2019")
fw_Fish_model_infile <- create_infile(LPI_full, index_vector=fw_Fish_model, name="Infiles/fw_Fish_model", end_col_name = "X2019")
m_Aves_model_infile <- create_infile(LPI_full, index_vector=m_Aves_model, name="Infiles/m_Aves_model", end_col_name = "X2019")
m_Mammalia_model_infile <- create_infile(LPI_full, index_vector=m_Mammalia_model, name="Infiles/m_Mammalia_model", end_col_name = "X2019")
m_Herps_model_infile <- create_infile(LPI_full, index_vector=m_Herps_model, name="Infiles/m_Herps_model", end_col_name = "X2019")
m_Fish_model_infile <- create_infile(LPI_full, index_vector=m_Fish_model, name="Infiles/m_Fish_model", end_col_name = "X2019")


Amphibian_vec <- LPI_full$Class=="Amphibia"
Reptile_vec <- LPI_full$Class=="Reptilia"
Herps_vec <- LPI_full$Class=="Amphibia" | (LPI_full$Class=="Reptilia" & LPI_full$Family!="Elapidae")
fw_Herps_vec <- (LPI_full$Class=="Amphibia" & LPI_full$FW_realm!="NULL") | (LPI_full$Class=="Reptilia" & LPI_full$FW_realm!="NULL")
T_Herps_vec <- (LPI_full$Class=="Amphibia" & LPI_full$T_realm!="NULL" & LPI_full$ID!=919) | (LPI_full$Class=="Reptilia" & LPI_full$T_realm!="NULL")
m_Herps_vec <- (LPI_full$Class=="Amphibia" & LPI_full$M_realm!="NULL") | (LPI_full$Class=="Reptilia" & LPI_full$M_realm!="NULL")
fw_Amphibian_vec <- (LPI_full$Class=="Amphibia" & LPI_full$FW_realm!="NULL")
fw_Reptile_vec <- (LPI_full$Class=="Reptilia" & LPI_full$FW_realm!="NULL")
T_Amphibian_vec <- (LPI_full$Class=="Amphibia" & LPI_full$T_realm!="NULL")
T_Reptile_vec <- (LPI_full$Class=="Reptilia" & LPI_full$T_realm!="NULL")
m_Amphibian_vec <- (LPI_full$Class=="Amphibia" & LPI_full$M_realm!="NULL")
m_Reptile_vec <- (LPI_full$Class=="Reptilia" & LPI_full$M_realm!="NULL" & LPI_full$Order!="Squamata")

Amphibian_infile <- create_infile(LPI_full, index_vector=Amphibian_vec, name="Infiles/Amphibian", end_col_name = "X2019")
Reptile_infile <- create_infile(LPI_full, index_vector=Reptile_vec, name="Infiles/Reptile", end_col_name = "X2019")
Herps_infile <- create_infile(LPI_full, index_vector=Herps_vec, name="Infiles/Herps", end_col_name = "X2019")
fw_Herps_infile <- create_infile(LPI_full, index_vector=fw_Herps_vec, name="Infiles/fw_Herps", end_col_name = "X2019")
T_Herps_infile <- create_infile(LPI_full, index_vector=T_Herps_vec, name="Infiles/T_Herps", end_col_name = "X2019")
m_Herps_infile <- create_infile(LPI_full, index_vector=m_Herps_vec, name="Infiles/m_Herps", end_col_name = "X2019")
fw_Amphibian_infile <- create_infile(LPI_full, index_vector=fw_Amphibian_vec, name="Infiles/fw_Amphibian", end_col_name = "X2019")
fw_Reptile_infile <- create_infile(LPI_full, index_vector=fw_Reptile_vec, name="Infiles/fw_Reptile", end_col_name = "X2019")
T_Amphibian_infile <- create_infile(LPI_full, index_vector=T_Amphibian_vec, name="Infiles/T_Amphibian", end_col_name = "X2019")
T_Reptile_infile <- create_infile(LPI_full, index_vector=T_Reptile_vec, name="Infiles/T_Reptile", end_col_name = "X2019")
m_Aphibian_infile <- create_infile(LPI_full, index_vector=m_Amphibian_vec, name="Infiles/m_Amphibian", end_col_name = "X2019")
m_Reptile_infile <- create_infile(LPI_full, index_vector=m_Reptile_vec, name="Infiles/m_Reptile", end_col_name = "X2019")

LPI_FW_Afrotropical_Herps2 <- LPI_FW_Afrotropical_Herps[LPI_FW_Afrotropical_Herps$Species!="adspersus" &
                                                          LPI_FW_Afrotropical_Herps$ID!=692 &
                                                          LPI_FW_Afrotropical_Herps$ID!=18238,]
FW_Afrotropical_Herps2_vec <- rownames(LPI_full) %in% rownames(LPI_FW_Afrotropical_Herps2)
FW_Afrotropical_Herps2_infile <- create_infile(LPI_full, index_vector=FW_Afrotropical_Herps2_vec, name="Infiles/FW_Afrotropical_Herps2", end_col_name = "X2019")
FW_Afrotropical_Herps2_lpi <- LPIMain("Infiles/FW_Afrotropical_Herps2_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)



LPI_FW_IndoPacific_Herps <- LPI_full[(LPI_full$Class=="Reptilia" | 
                                        LPI_full$Class=="Amphibia") &
                                       (LPI_full$FW_realm=="Indo-Malayan" |
                                          LPI_full$FW_realm=="Oceania" |
                                          LPI_full$FW_realm=="Australasia"), ]

LPI_FW_IndoPacific_Herps2 <- LPI_FW_IndoPacific_Herps[LPI_FW_IndoPacific_Herps$Common_name!="Northern corroboree frog" & 
                                                        LPI_FW_IndoPacific_Herps$Common_name!="Corroboree frog" & 
                                                        LPI_FW_IndoPacific_Herps$Common_name!="Baw baw frog",]

LPI_FW_IndoPacific_Herps3 <- LPI_FW_IndoPacific_Herps[LPI_FW_IndoPacific_Herps$Species!="australis" & 
                                                        LPI_FW_IndoPacific_Herps$Country!="Nepal",]

LPI_FW_IndoPacific_Herps4 <- LPI_FW_IndoPacific_Herps2[LPI_FW_IndoPacific_Herps2$Species!="australis" & 
                                                        LPI_FW_IndoPacific_Herps2$Country!="Nepal",]

FW_IndoPacific_Herps2_vec <- rownames(LPI_full) %in% rownames(LPI_FW_IndoPacific_Herps2)
FW_IndoPacific_Herps2_infile <- create_infile(LPI_full, index_vector=FW_IndoPacific_Herps2_vec, name="Infiles/FW_IndoPacific_Herps2", end_col_name = "X2019")
FW_IndoPacific_Herps2_lpi <- LPIMain("Infiles/FW_IndoPacific_Herps2_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

FW_IndoPacific_Herps3_vec <- rownames(LPI_full) %in% rownames(LPI_FW_IndoPacific_Herps3)
FW_IndoPacific_Herps3_infile <- create_infile(LPI_full, index_vector=FW_IndoPacific_Herps3_vec, name="Infiles/FW_IndoPacific_Herps3", end_col_name = "X2019")
FW_IndoPacific_Herps3_lpi <- LPIMain("Infiles/FW_IndoPacific_Herps3_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

FW_IndoPacific_Herps4_vec <- rownames(LPI_full) %in% rownames(LPI_FW_IndoPacific_Herps4)
FW_IndoPacific_Herps4_infile <- create_infile(LPI_full, index_vector=FW_IndoPacific_Herps4_vec, name="Infiles/FW_IndoPacific_Herps4", end_col_name = "X2019")
FW_IndoPacific_Herps4_lpi <- LPIMain("Infiles/FW_IndoPacific_Herps4_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)


fw_IndoPacific_Herps_lpi <- fw_IndoPacific_Herps_lpi[complete.cases(fw_IndoPacific_Herps_lpi), ]
FW_IndoPacific_Herps2_lpi <- FW_IndoPacific_Herps2_lpi[complete.cases(FW_IndoPacific_Herps2_lpi), ]
FW_IndoPacific_Herps3_lpi <- FW_IndoPacific_Herps3_lpi[complete.cases(FW_IndoPacific_Herps3_lpi), ]
FW_IndoPacific_Herps4_lpi <- FW_IndoPacific_Herps4_lpi[complete.cases(FW_IndoPacific_Herps4_lpi), ]

ggplot_lpi(fw_IndoPacific_Herps_lpi, ylim=c(0,2))
ggplot_lpi(FW_IndoPacific_Herps2_lpi, ylim=c(0,2))
ggplot_lpi(FW_IndoPacific_Herps3_lpi, ylim=c(0,2))
ggplot_lpi(FW_IndoPacific_Herps4_lpi, ylim=c(0,2))

Herps_analysis_lpi2 <- list(fw_IndoPacific_Herps_lpi, 
                           FW_IndoPacific_Herps2_lpi,
                           FW_IndoPacific_Herps3_lpi,
                           FW_IndoPacific_Herps4_lpi)

Herps_analysis_lpi_plot2 <- ggplot_multi_lpi(Herps_analysis_lpi2, 
                                            title="Freshwater IndoPacific Herps", 
                                            names=c("normal LPI trend",
                                                    "trend with 3 crashing species removed",
                                                    "trend with 2 rising species removed",
                                                    "trend with all 5 species removed"), 
                                            xlims=c(1970, 2018), 
                                            ylims=c(0,2), 
                                            facet=TRUE) +
  theme(legend.position="none", legend.title=element_blank())

ggsave("Plots/IndoPacific_Herps_sensitivity2.jpeg", 
       Herps_analysis_lpi_plot2, 
       height = 2000, 
       width = 8000, 
       units="px")




Amphibian_lpi  <- LPIMain("Infiles/Amphibian_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
Reptile_lpi  <- LPIMain("Infiles/Reptile_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
Herps_lpi  <- LPIMain("Infiles/Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Herps_lpi  <- LPIMain("Infiles/fw_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Herps_lpi  <- LPIMain("Infiles/T_Herps_infile.txt", REF_YEAR = 1976, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_Herps_lpi  <- LPIMain("Infiles/m_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Amphibian_lpi  <- LPIMain("Infiles/fw_Amphibian_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Reptile_lpi  <- LPIMain("Infiles/fw_Reptile_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Amphibian_lpi  <- LPIMain("Infiles/T_Amphibian_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Reptile_lpi  <- LPIMain("Infiles/T_Reptile_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_Amphibian_lpi  <- LPIMain("Infiles/m_Amphibian_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_Reptile_lpi  <- LPIMain("Infiles/m_Reptile_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_TroSubIndo_Herps_lpi  <- LPIMain("Infiles/m_TroSubIndo_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_TroSubIndo_lpi <- LPIMain("Infiles/m_TroSubIndo_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=1, use_weightings_B=0)


# create indices and plot
fw_Afrotropical_fish_lpi  <- LPIMain("Infiles/fw_Afrotropical_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_IndoPacific_fish_lpi  <- LPIMain("Infiles/fw_IndoPacific_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Palearctic_fish_lpi  <- LPIMain("Infiles/fw_Palearctic_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Neotropical_fish_lpi  <- LPIMain("Infiles/fw_Neotropical_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Nearctic_fish_lpi  <- LPIMain("Infiles/fw_Nearctic_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_Afrotropical_Herps_lpi  <- LPIMain("Infiles/T_Afrotropical_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_IndoPacific_Herps_lpi  <- LPIMain("Infiles/T_IndoPacific_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Palearctic_Herps_lpi  <- LPIMain("Infiles/T_Palearctic_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Neotropical_Herps_lpi  <- LPIMain("Infiles/T_Neotropical_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Nearctic_Herps_lpi  <- LPIMain("Infiles/T_Nearctic_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_Afrotropical_Herps_lpi  <- LPIMain("Infiles/fw_Afrotropical_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_IndoPacific_Herps_lpi  <- LPIMain("Infiles/fw_IndoPacific_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Palearctic_Herps_lpi  <- LPIMain("Infiles/fw_Palearctic_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Neotropical_Herps_lpi  <- LPIMain("Infiles/fw_Neotropical_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Nearctic_Herps_lpi  <- LPIMain("Infiles/fw_Nearctic_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_Afrotropical_Aves_lpi  <- LPIMain("Infiles/T_Afrotropical_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_IndoPacific_Aves_lpi  <- LPIMain("Infiles/T_IndoPacific_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Palearctic_Aves_lpi  <- LPIMain("Infiles/T_Palearctic_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Neotropical_Aves_lpi  <- LPIMain("Infiles/T_Neotropical_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Nearctic_Aves_lpi  <- LPIMain("Infiles/T_Nearctic_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_Afrotropical_Aves_lpi  <- LPIMain("Infiles/fw_Afrotropical_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_IndoPacific_Aves_lpi  <- LPIMain("Infiles/fw_IndoPacific_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Palearctic_Aves_lpi  <- LPIMain("Infiles/fw_Palearctic_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Neotropical_Aves_lpi  <- LPIMain("Infiles/fw_Neotropical_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Nearctic_Aves_lpi  <- LPIMain("Infiles/fw_Nearctic_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_Afrotropical_Mammalia_lpi  <- LPIMain("Infiles/T_Afrotropical_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_IndoPacific_Mammalia_lpi  <- LPIMain("Infiles/T_IndoPacific_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Palearctic_Mammalia_lpi  <- LPIMain("Infiles/T_Palearctic_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Neotropical_Mammalia_lpi  <- LPIMain("Infiles/T_Neotropical_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Nearctic_Mammalia_lpi  <- LPIMain("Infiles/T_Nearctic_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_Afrotropical_Mammalia_lpi  <- LPIMain("Infiles/fw_Afrotropical_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_IndoPacific_Mammalia_lpi  <- LPIMain("Infiles/fw_IndoPacific_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Palearctic_Mammalia_lpi  <- LPIMain("Infiles/fw_Palearctic_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Neotropical_Mammalia_lpi  <- LPIMain("Infiles/fw_Neotropical_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Nearctic_Mammalia_lpi  <- LPIMain("Infiles/fw_Nearctic_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_AtNoTemp_Aves <- LPIMain("Infiles/m_AtNoTemp_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtTrSub_Aves <- LPIMain("Infiles/m_AtTrSub_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_Arctic_Aves <- LPIMain("Infiles/m_SoTeAnt_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_SoTeAnt_Aves <- LPIMain("Infiles/m_SoTeAnt_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_TroSubIndo_Aves <- LPIMain("Infiles/m_TroSubIndo_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_PaNoTemp_Aves <- LPIMain("Infiles/m_PaNoTemp_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_AtNoTemp_Mammalia <- LPIMain("Infiles/m_AtNoTemp_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtTrSub_Mammalia <- LPIMain("Infiles/m_AtTrSub_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_Arctic_Mammalia <- LPIMain("Infiles/m_SoTeAnt_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_SoTeAnt_Mammalia <- LPIMain("Infiles/m_SoTeAnt_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_TroSubIndo_Mammalia <- LPIMain("Infiles/m_TroSubIndo_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_PaNoTemp_Mammalia <- LPIMain("Infiles/m_PaNoTemp_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_AtNoTemp_Herps <- LPIMain("Infiles/m_AtNoTemp_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtTrSub_Herps <- LPIMain("Infiles/m_AtTrSub_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_Arctic_Herps <- LPIMain("Infiles/m_SoTeAnt_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_SoTeAnt_Herps <- LPIMain("Infiles/m_SoTeAnt_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_TroSubIndo_Herps <- LPIMain("Infiles/m_TroSubIndo_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_PaNoTemp_Herps <- LPIMain("Infiles/m_PaNoTemp_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_AtNoTemp_Fish <- LPIMain("Infiles/m_AtNoTemp_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtTrSub_Fish <- LPIMain("Infiles/m_AtTrSub_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_Arctic_Fish <- LPIMain("Infiles/m_SoTeAnt_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_SoTeAnt_Fish <- LPIMain("Infiles/m_SoTeAnt_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_TroSubIndo_Fish <- LPIMain("Infiles/m_TroSubIndo_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_PaNoTemp_Fish <- LPIMain("Infiles/m_PaNoTemp_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

# there is an issue with complex infiles being located outside of the project directory...
t_mam_lpi  <- LPIMain("terrestrial_mammals_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=1)
fw_mam_lpi  <- LPIMain("freshwater_mammals_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=1)
fw_herps_lpi  <- LPIMain("freshwater_herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=1)
fw_fish_lpi  <- LPIMain("freshwater_fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=1)
setwd("Infiles")
terrestrial_lpi <- LPIMain("terrestrial_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=1, use_weightings_B=1)
freshwater_lpi <- LPIMain("freshwater_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=1, use_weightings_B=1)
marine_lpi <- LPIMain("marine_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=1, use_weightings_B=1)
full_lpi <- LPIMain("lpi_infile.txt", REF_YEAR = 1970, PLOT_MAX = 100, force_recalculation=1, use_weightings=1, use_weightings_B=1)

terrestrial_model_lpi <- LPIMain("terrestrial_model_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=1, use_weightings_B=0)
freshwater_model_lpi <- LPIMain("freshwater_model_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=1, use_weightings_B=0)
marine_model_lpi <- LPIMain("marine_model_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=1, use_weightings_B=0)
full_model_lpi <- LPIMain("lpi_model_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=1, use_weightings_B=0)

t_mam_lpi <- t_mam_lpi[complete.cases(t_mam_lpi), ]
fw_mam_lpi <- fw_mam_lpi[complete.cases(fw_mam_lpi), ]
terrestrial_lpi <- terrestrial_lpi[complete.cases(terrestrial_lpi), ]
freshwater_lpi <- freshwater_lpi[complete.cases(freshwater_lpi), ]
marine_lpi <- marine_lpi[complete.cases(marine_lpi), ] 
full_lpi <- full_lpi[complete.cases(full_lpi), ] 

T_Afrotropical_Herps_lpi <- T_Afrotropical_Herps_lpi[complete.cases(T_Afrotropical_Herps_lpi), ]
fw_Afrotropical_Herps_lpi <- fw_Afrotropical_Herps_lpi[complete.cases(fw_Afrotropical_Herps_lpi), ]
T_IndoPacific_Herps_lpi <- T_IndoPacific_Herps_lpi[complete.cases(T_IndoPacific_Herps_lpi), ]
fw_IndoPacific_Herps_lpi <- fw_IndoPacific_Herps_lpi[complete.cases(fw_IndoPacific_Herps_lpi), ]

terrestrial_model_lpi <- terrestrial_model_lpi[complete.cases(terrestrial_model_lpi), ]
freshwater_model_lpi <- freshwater_model_lpi[complete.cases(freshwater_model_lpi), ]
marine_model_lpi <- marine_model_lpi[complete.cases(marine_model_lpi), ] 
full_model_lpi <- full_model_lpi[complete.cases(full_model_lpi), ] 
                                               
ggplot_lpi(t_mam_lpi, ylim=c(0,2))
ggplot_lpi(fw_mam_lpi, ylim=c(0,2))
ggplot_lpi(terrestrial_lpi, ylim=c(0,2))
ggplot_lpi(freshwater_lpi, ylim=c(0,2))
ggplot_lpi(marine_lpi, ylim=c(0,2))
ggplot_lpi(full_lpi, ylim=c(0,2))

ggplot_lpi(T_Afrotropical_Herps_lpi, ylim=c(0,5))
ggplot_lpi(T_IndoPacific_Herps_lpi, ylim=c(0,5))
ggplot_lpi(fw_Afrotropical_Herps_lpi, ylim=c(0,1))
ggplot_lpi(fw_IndoPacific_Herps_lpi, ylim=c(0,1))

Herps_analysis_lpi <- list(T_Afrotropical_Herps_lpi, 
                           T_IndoPacific_Herps_lpi,
                           fw_Afrotropical_Herps_lpi, 
                           fw_IndoPacific_Herps_lpi)


Herps_analysis_lpi_plot <- ggplot_multi_lpi(Herps_analysis_lpi, 
                                            title="", 
                                            names=c("Terrestrial Afrotropical Herps", 
                                                    "Terrestrial Indopacific Herps",
                                                    "Freshwater Afrotropical Herps",
                                                    "Freshwater IndoPacific Herps"), 
                                            xlims=c(1970, 2018), 
                                            ylims=c(0,5), 
                                            facet=TRUE) +
  theme(legend.position="none", legend.title=element_blank())

ggsave("Plots/Afro_and_Indo_Herps.jpeg", 
       Herps_analysis_lpi_plot, 
       height = 2000, 
       width = 6000, 
       units="px")


ggplot_lpi(terrestrial_model_lpi, ylim=c(0,2))
ggplot_lpi(freshwater_model_lpi, ylim=c(0,2))
ggplot_lpi(marine_model_lpi, ylim=c(0,2))
ggplot_lpi(full_model_lpi, ylim=c(0,2))

lpi_systems <- list(terrestrial_lpi, freshwater_lpi, marine_lpi)

model_lpi_systems <- list(terrestrial_model_lpi, freshwater_model_lpi, marine_model_lpi)

combined_lpi_terrestrial <- list(terrestrial_lpi, terrestrial_model_lpi)

combined_lpi_freshwater <- list(freshwater_lpi, freshwater_model_lpi)

combined_lpi_marine <- list(marine_lpi, marine_model_lpi)

combined_lpi_systems <- list(combined_lpi_terrestrial, combined_lpi_freshwater, combined_lpi_marine)

combined_lpi_all <- list(full_lpi, full_model_lpi)

ggplot_multi_lpi(lpi_systems, names=c("Terrestrial LPI", "Freshwater LPI", "Marine LPI"), xlims=c(1970, 2018), ylims=c(0,2), facet=TRUE)

ggplot_multi_lpi(model_lpi_systems, names=c("Terrestrial Sampled", "Freshwater Sampled", "Marine Sampled"), xlims=c(1970, 2018), ylims=c(0,2), facet=TRUE)

lpi_plot <- ggplot_multi_lpi(combined_lpi_all, title="Living Planet Index", names=c("Diversity Weighted", "Sampled"), xlims=c(1970, 2018), ylims=c(0,2), facet=FALSE) +
  theme(legend.position=c(0.3,0.9), legend.title=element_blank())

terr_plot <- ggplot_multi_lpi(combined_lpi_terrestrial, title="Terrestrial", names=c("Diversity Weighted", "Sampled"), xlims=c(1970, 2018), ylims=c(0,2), facet=FALSE) +
  theme(legend.position=c(0.3,0.85), legend.title=element_blank())

fw_plot <- ggplot_multi_lpi(combined_lpi_freshwater, title="Freshwater", names=c("Diversity Weighted", "Sampled"), xlims=c(1970, 2018), ylims=c(0,2), facet=FALSE) +
  theme(legend.position=c(0.3,0.85), legend.title=element_blank())

mar_plot <- ggplot_multi_lpi(combined_lpi_marine, title="Marine", names=c("Diversity Weighted", "Sampled"), xlims=c(1970, 2018), ylims=c(0,2), facet=FALSE) +
  theme(legend.position=c(0.3,0.85), legend.title=element_blank())

ggplot_multi_lpi(combined_lpi_systems, names=c("Terrestrial", "Freshwater", "Marine"), xlims=c(1970, 2018), ylims=c(0,2), facet=TRUE)

all_plots <- grid.arrange(terr_plot, fw_plot, mar_plot, ncol=2)

ggsave("Plots/RW_LPI_plots.jpg", all_plots, height = 2300, width = 2600, units="px")


