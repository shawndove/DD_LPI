# libraries ----

library(rlpi)
library(dplyr)
library(mgcv)
library(matrixStats)
library(ggplot2)
library(gridExtra)

# load functions ----
source("Scripts/Functions_NG.R")

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



# load public LPI data ----

# only public data
#LPI_full <- read.csv(file="Data/LPR2020data_public.csv", sep=",", stringsAsFactors=FALSE)
#LPI_test <- read_excel("Data/LPR2020data_public.csv")

# all data, including private
LPI_full <- read.csv(file="Data/LPD_output_20201116.csv", sep=",", stringsAsFactors=FALSE)
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

# add X back to years in column names (to make create_infile work properly)
colnames(LPI_full)[65:134] <- paste("X", colnames(LPI_full)[65:134], sep="") # full
#colnames(LPI_full)[30:98] <- paste("X", colnames(LPI_full)[30:98], sep="") # public

# create infiles
T_Afrotropical_Aves_infile <- create_infile(LPI_full, index_vector=T_Afrotropical_Aves, name="Infiles/T_Afrotropical_Aves", end_col_name = "X2017")
T_Afrotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Afrotropical_Mammalia, name="Infiles/T_Afrotropical_Mammalia", end_col_name = "X2017")
T_Afrotropical_Herps_infile <- create_infile(LPI_full, index_vector=T_Afrotropical_Herps, name="Infiles/T_Afrotropical_Herps", end_col_name = "X2017")

T_IndoPacific_Aves_infile <- create_infile(LPI_full, index_vector=T_IndoPacific_Aves, name="Infiles/T_IndoPacific_Aves", end_col_name = "X2017")
T_IndoPacific_Mammalia_infile <- create_infile(LPI_full, index_vector=T_IndoPacific_Mammalia, name="Infiles/T_IndoPacific_Mammalia", end_col_name = "X2017")
T_IndoPacific_Herps_infile <- create_infile(LPI_full, index_vector=T_IndoPacific_Herps, name="Infiles/T_IndoPacific_Herps", end_col_name = "X2017")

T_Palearctic_Aves_infile <- create_infile(LPI_full, index_vector=T_Palearctic_Aves, name="Infiles/T_Palearctic_Aves", end_col_name = "X2017")
T_Palearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Palearctic_Mammalia, name="Infiles/T_Palearctic_Mammalia", end_col_name = "X2017")
T_Palearctic_Herps_infile <- create_infile(LPI_full, index_vector=T_Palearctic_Herps, name="Infiles/T_Palearctic_Herps", end_col_name = "X2017")

T_Neotropical_Aves_infile <- create_infile(LPI_full, index_vector=T_Neotropical_Aves, name="Infiles/T_Neotropical_Aves", end_col_name = "X2017")
T_Neotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Neotropical_Mammalia, name="Infiles/T_Neotropical_Mammalia", end_col_name = "X2017")
T_Neotropical_Herps_infile <- create_infile(LPI_full, index_vector=T_Neotropical_Herps, name="Infiles/T_Neotropical_Herps", end_col_name = "X2017")

T_Nearctic_Aves_infile <- create_infile(LPI_full, index_vector=T_Nearctic_Aves, name="Infiles/T_Nearctic_Aves", end_col_name = "X2017")
T_Nearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Nearctic_Mammalia, name="Infiles/T_Nearctic_Mammalia", end_col_name = "X2017")
T_Nearctic_Herps_infile <- create_infile(LPI_full, index_vector=T_Nearctic_Herps, name="Infiles/T_Nearctic_Herps", end_col_name = "X2017")

T_Antarctic_Aves_infile <- create_infile(LPI_full, index_vector=T_Antarctic_Aves, name="T_Antarctic_Aves", end_col_name = "X2017")
T_Antarctic_Mammalia_infile <- create_infile(LPI_full, index_vector=T_Antarctic_Mammalia, name="T_Antarctic_Mammalia", end_col_name = "X2017")

fw_Afrotropical_Aves_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Aves, name="Infiles/fw_Afrotropical_Aves", end_col_name = "X2017")
fw_Afrotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Mammalia, name="Infiles/fw_Afrotropical_Mammalia", end_col_name = "X2017")
fw_Afrotropical_Herps_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Herps, name="Infiles/fw_Afrotropical_Herps", end_col_name = "X2017")
fw_Afrotropical_Fish_infile <- create_infile(LPI_full, index_vector=fw_Afrotropical_Fish, name="Infiles/fw_Afrotropical_Fish", end_col_name = "X2017")

fw_IndoPacific_Aves_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Aves, name="Infiles/fw_IndoPacific_Aves", end_col_name = "X2017")
fw_IndoPacific_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Mammalia, name="Infiles/fw_IndoPacific_Mammalia", end_col_name = "X2017")
fw_IndoPacific_Herps_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Herps, name="Infiles/fw_IndoPacific_Herps", end_col_name = "X2017")
fw_IndoPacific_Fish_infile <- create_infile(LPI_full, index_vector=fw_IndoPacific_Fish, name="Infiles/fw_IndoPacific_Fish", end_col_name = "X2017")

fw_Palearctic_Aves_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Aves, name="Infiles/fw_Palearctic_Aves", end_col_name = "X2017")
fw_Palearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Mammalia, name="Infiles/fw_Palearctic_Mammalia", end_col_name = "X2017")
fw_Palearctic_Herps_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Herps, name="Infiles/fw_Palearctic_Herps", end_col_name = "X2017")
fw_Palearctic_Fish_infile <- create_infile(LPI_full, index_vector=fw_Palearctic_Fish, name="Infiles/fw_Palearctic_Fish", end_col_name = "X2017")

fw_Neotropical_Aves_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Aves, name="Infiles/fw_Neotropical_Aves", end_col_name = "X2017")
fw_Neotropical_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Mammalia, name="Infiles/fw_Neotropical_Mammalia", end_col_name = "X2017")
fw_Neotropical_Herps_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Herps, name="Infiles/fw_Neotropical_Herps", end_col_name = "X2017")
fw_Neotropical_Fish_infile <- create_infile(LPI_full, index_vector=fw_Neotropical_Fish, name="Infiles/fw_Neotropical_Fish", end_col_name = "X2017")

fw_Nearctic_Aves_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Aves, name="Infiles/fw_Nearctic_Aves", end_col_name = "X2017")
fw_Nearctic_Mammalia_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Mammalia, name="Infiles/fw_Nearctic_Mammalia", end_col_name = "X2017")
fw_Nearctic_Herps_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Herps, name="Infiles/fw_Nearctic_Herps", end_col_name = "X2017")
fw_Nearctic_Fish_infile <- create_infile(LPI_full, index_vector=fw_Nearctic_Fish, name="Infiles/fw_Nearctic_Fish", end_col_name = "X2017")

m_AtNoTemp_Aves_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Aves, name="Infiles/m_AtNoTemp_Aves", end_col_name = "X2017")
m_AtNoTemp_Mammalia_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Mammalia, name="Infiles/m_AtNoTemp_Mammalia", end_col_name = "X2017")
m_AtNoTemp_Herps_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Herps, name="Infiles/m_AtNoTemp_Herps", end_col_name = "X2017")
m_AtNoTemp_Fish_infile <- create_infile(LPI_full, index_vector=m_AtNoTemp_Fish, name="Infiles/m_AtNoTemp_Fish", end_col_name = "X2017")

m_AtTrSub_Aves_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Aves, name="Infiles/m_AtTrSub_Aves", end_col_name = "X2017")
m_AtTrSub_Mammalia_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Mammalia, name="Infiles/m_AtTrSub_Mammalia", end_col_name = "X2017")
m_AtTrSub_Herps_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Herps, name="Infiles/m_AtTrSub_Herps", end_col_name = "X2017")
m_AtTrSub_Fish_infile <- create_infile(LPI_full, index_vector=m_AtTrSub_Fish, name="Infiles/m_AtTrSub_Fish", end_col_name = "X2017")

m_Arctic_Aves_infile <- create_infile(LPI_full, index_vector=m_Arctic_Aves, name="Infiles/m_Arctic_Aves", end_col_name = "X2017")
m_Arctic_Mammalia_infile <- create_infile(LPI_full, index_vector=m_Arctic_Mammalia, name="Infiles/m_Arctic_Mammalia", end_col_name = "X2017")
m_Arctic_Fish_infile <- create_infile(LPI_full, index_vector=m_Arctic_Fish, name="Infiles/m_Arctic_Fish", end_col_name = "X2017")

m_SoTeAnt_Aves_infile <- create_infile(LPI_full, index_vector=m_SoTeAnt_Aves, name="Infiles/m_SoTeAnt_Aves", end_col_name = "X2017")
m_SoTeAnt_Mammalia_infile <- create_infile(LPI_full, index_vector=m_SoTeAnt_Mammalia, name="Infiles/m_SoTeAnt_Mammalia", end_col_name = "X2017")
m_SoTeAnt_Herps_infile <- create_infile(LPI_full, index_vector=m_SoTeAnt_Herps, name="Infiles/m_SoTeAnt_Herps", end_col_name = "X2017")
m_SoTeAnt_Fish_infile <- create_infile(LPI_full, index_vector=m_SoTeAnt_Fish, name="Infiles/m_SoTeAnt_Fish", end_col_name = "X2017")

m_TroSubIndo_Aves_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Aves, name="Infiles/m_TroSubIndo_Aves", end_col_name = "X2017")
m_TroSubIndo_Mammalia_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Mammalia, name="Infiles/m_TroSubIndo_Mammalia", end_col_name = "X2017")
m_TroSubIndo_Herps_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Herps, name="Infiles/m_TroSubIndo_Herps", end_col_name = "X2017")
m_TroSubIndo_Fish_infile <- create_infile(LPI_full, index_vector=m_TroSubIndo_Fish, name="Infiles/m_TroSubIndo_Fish", end_col_name = "X2017")

m_PaNoTemp_Aves_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Aves, name="Infiles/m_PaNoTemp_Aves", end_col_name = "X2017")
m_PaNoTemp_Mammalia_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Mammalia, name="Infiles/m_PaNoTemp_Mammalia", end_col_name = "X2017")
m_PaNoTemp_Herps_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Herps, name="Infiles/m_PaNoTemp_Herps", end_col_name = "X2017")
m_PaNoTemp_Fish_infile <- create_infile(LPI_full, index_vector=m_PaNoTemp_Fish, name="Infiles/m_PaNoTemp_Fish", end_col_name = "X2017")

# calculate mean growth rate for each actual LPI group
gr.stats.list <- list()
# loop over taxonomic groups
for (i in 1:length(pop_list)) {
  
  # select all copies of the populations listed in the sample
  group_data <- LPI_trimmed[LPI_trimmed$PopID %in% pop_list[[i]],]
  
  # remove all populations with less than 2 data points
  group_data_culled <- cull_fn(group_data, 2, 2, c2)
  
  # log-linear interpolate all pops
  grp_completed <- complete_time_series(group_data_culled, c2, m_colnames2, calcsd=TRUE)
  
  if (nrow(grp_completed) >=1) {
    
    # get mean and standard deviation of the mean growth rate
    gr.stats.list[[i]] <- growth_rate_calc_fn3(grp_completed, c2, model=TRUE)
    
  } else {
    
    gr.stats.list[[i]] <- NA
    
  }
  
}
# save mean growth rate data
saveRDS(gr.stats.list, file="gr_stats_list.RData")
gr.stats.list <- readRDS("gr_stats_list.RData")

# calculate number of populations (sample sizes) in LPI groups
pop.size.list <- list()
# loop over groups
for (i in 1:length(pop_list)) {
  
  pop.size.list[[i]] <- length(pop_list[[i]])
  
}
pop.size.vec <- unlist(pop.size.list)
# save pop size data
saveRDS(pop.size.vec, file="pop_size_vec.RData")

# calculate number of species in LPI groups
spec.size.list <- list()
for (i in 1:length(pop_list)) {
  
  spec.size.list[[i]] <- length(unique(LPI_trimmed$SpecID[LPI_trimmed$PopID %in% pop_list[[i]]]))
  
}
spec.size.vec <- unlist(spec.size.list)
# save spec size data
saveRDS(spec.size.vec, file="spec_size_vec.RData")

#calculate mean time series length in LPI groups
mean.tslength.list <- list()
# loop over groups
for (i in 1:length(pop_list)) {
  
  # select all copies of the populations listed in the sample
  group_data <- LPI_trimmed[LPI_trimmed$PopID %in% pop_list[[i]],]
  
  # remove all populations with less than 2 data points
  group_data_culled <- cull_fn(group_data, 2, 2, c2)
  
  # log-linear interpolate all pops
  grp_completed <- complete_time_series(group_data_culled, c2, m_colnames2, calcsd=TRUE)
  
  # calculate number of non-NA values and divide by total number of time series
  # this gives the mean time series length
  mean.tslength.list[[i]] <- sum(!is.na(as.vector(grp_completed[,1:c2]))) / pop.size.list[[i]]
  
}
mean.tslength.vec <- unlist(mean.tslength.list)
# save time series length data
saveRDS(mean.tslength.vec, file="mean_ts_length_vec.RData")
mean.tslength.vec <- readRDS(file="mean_ts_length_vec.RData")

# find tdv and minimum sample sizes
tdv.list <- list()
samp.size.list <- list()
# loop over groups
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
saveRDS(tdv.vec, file="tdv_vec.RData")
saveRDS(samp.size.vec, file="samp_size_vec.RData")

# calculate percentage of min sample size in data for each group
samp.percent.list <- list()
# loop over groups
for (i in 1:length(gr.stats.list)) {
  
  temp.min <- samp.size.list[[i]]
  
  temp.actual <- pop.size.list[[i]]
  
  ratio <- temp.actual/temp.min
  
  if (is.na(ratio)) {
    
    samp.percent.list[[i]] <- NA
    
#  } else if (ratio > 1) {
    
#    samp.percent.list[[i]] <- 100
    
  } else {
    
    samp.percent.list[[i]] <- ratio * 100
    
  }
  
}
samp.percent.vec <- unlist(samp.percent.list)
# save percentage sample size data
saveRDS(samp.percent.vec, file="samp_percent_vec.RData")
samp.percent.vec <- readRDS("samp_percent_vec.RData")

# calculate number of populations that need to be added to get below the tdv threshold
pops.added.list <- list()
# loop over groups
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
saveRDS(pops.added.vec, file="pops_added_vec.RData")


####


# create indices and plot
T_Afrotropical_Aves_lpi  <- LPIMain("Infiles/T_Afrotropical_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Afrotropical_Mammalia_lpi  <- LPIMain("Infiles/T_Afrotropical_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Afrotropical_Herps_lpi  <- LPIMain("Infiles/T_Afrotropical_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_IndoPacific_Aves_lpi  <- LPIMain("Infiles/T_IndoPacific_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_IndoPacific_Mammalia_lpi  <- LPIMain("Infiles/T_IndoPacific_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_IndoPacific_Herps_lpi  <- LPIMain("Infiles/T_IndoPacific_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_Palearctic_Aves_lpi  <- LPIMain("Infiles/T_Palearctic_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Palearctic_Mammalia_lpi  <- LPIMain("Infiles/T_Palearctic_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Palearctic_Herps_lpi  <- LPIMain("Infiles/T_Palearctic_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_Neotropical_Aves_lpi  <- LPIMain("Infiles/T_Neotropical_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Neotropical_Mammalia_lpi  <- LPIMain("Infiles/T_Neotropical_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Neotropical_Herps_lpi  <- LPIMain("Infiles/T_Neotropical_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_Nearctic_Aves_lpi  <- LPIMain("Infiles/T_Nearctic_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Nearctic_Mammalia_lpi  <- LPIMain("Infiles/T_Nearctic_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Nearctic_Herps_lpi  <- LPIMain("Infiles/T_Nearctic_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

T_Antarctic_Aves_lpi <- LPIMain("Infiles/T_Antarctic_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
T_Antarctic_Mammalia_lpi <- LPIMain("Infiles/T_Antarctic_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_Afrotropical_Aves_lpi  <- LPIMain("Infiles/fw_Afrotropical_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Afrotropical_Mammalia_lpi  <- LPIMain("Infiles/fw_Afrotropical_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Afrotropical_Herps_lpi  <- LPIMain("Infiles/fw_Afrotropical_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Afrotropical_fish_lpi  <- LPIMain("Infiles/fw_Afrotropical_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_IndoPacific_Aves_lpi  <- LPIMain("Infiles/fw_IndoPacific_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_IndoPacific_Mammalia_lpi  <- LPIMain("Infiles/fw_IndoPacific_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_IndoPacific_Herps_lpi  <- LPIMain("Infiles/fw_IndoPacific_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_IndoPacific_fish_lpi  <- LPIMain("Infiles/fw_IndoPacific_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_Palearctic_Aves_lpi  <- LPIMain("Infiles/fw_Palearctic_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Palearctic_Mammalia_lpi  <- LPIMain("Infiles/fw_Palearctic_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Palearctic_Herps_lpi  <- LPIMain("Infiles/fw_Palearctic_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Palearctic_fish_lpi  <- LPIMain("Infiles/fw_Palearctic_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_Neotropical_Aves_lpi  <- LPIMain("Infiles/fw_Neotropical_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Neotropical_Mammalia_lpi  <- LPIMain("Infiles/fw_Neotropical_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Neotropical_Herps_lpi  <- LPIMain("Infiles/fw_Neotropical_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Neotropical_fish_lpi  <- LPIMain("Infiles/fw_Neotropical_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

fw_Nearctic_Aves_lpi  <- LPIMain("Infiles/fw_Nearctic_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Nearctic_Mammalia_lpi  <- LPIMain("Infiles/fw_Nearctic_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Nearctic_Herps_lpi  <- LPIMain("Infiles/fw_Nearctic_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
fw_Nearctic_fish_lpi  <- LPIMain("Infiles/fw_Nearctic_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_AtNoTemp_Aves <- LPIMain("Infiles/m_AtNoTemp_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtNoTemp_Mammalia <- LPIMain("Infiles/m_AtNoTemp_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtNoTemp_Herps <- LPIMain("Infiles/m_AtNoTemp_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtNoTemp_Fish <- LPIMain("Infiles/m_AtNoTemp_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_AtTrSub_Aves <- LPIMain("Infiles/m_AtTrSub_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtTrSub_Mammalia <- LPIMain("Infiles/m_AtTrSub_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtTrSub_Herps <- LPIMain("Infiles/m_AtTrSub_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_AtTrSub_Fish <- LPIMain("Infiles/m_AtTrSub_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_Arctic_Aves <- LPIMain("Infiles/m_Arctic_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_Arctic_Mammalia <- LPIMain("Infiles/m_Arctic_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_Arctic_Fish <- LPIMain("Infiles/m_Arctic_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_SoTeAnt_Aves <- LPIMain("Infiles/m_SoTeAnt_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_SoTeAnt_Mammalia <- LPIMain("Infiles/m_SoTeAnt_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_SoTeAnt_Herps <- LPIMain("Infiles/m_SoTeAnt_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_SoTeAnt_Fish <- LPIMain("Infiles/m_SoTeAnt_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_TroSubIndo_Aves <- LPIMain("Infiles/m_TroSubIndo_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_TroSubIndo_Mammalia <- LPIMain("Infiles/m_TroSubIndo_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_TroSubIndo_Herps <- LPIMain("Infiles/m_TroSubIndo_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_TroSubIndo_Fish <- LPIMain("Infiles/m_TroSubIndo_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

m_PaNoTemp_Aves <- LPIMain("Infiles/m_PaNoTemp_Aves_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_PaNoTemp_Mammalia <- LPIMain("Infiles/m_PaNoTemp_Mammalia_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_PaNoTemp_Herps <- LPIMain("Infiles/m_PaNoTemp_Herps_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)
m_PaNoTemp_Fish <- LPIMain("Infiles/m_PaNoTemp_Fish_infile.txt", REF_YEAR = 1970, PLOT_MAX = 2018, BOOT_STRAP_SIZE = 100, force_recalculation=1, use_weightings=0, use_weightings_B=0)

