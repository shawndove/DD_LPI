## LPI regional taxonomic results

sample_percent <- samp.percent.vec # import results from RWLPI_NG_3.R
face_type <- vector()
for (i in 1:length(sample_percent)) {
  if (is.na(sample_percent[i])) {
    face_type[i] <- NA
  } else if (sample_percent[i] >= 100) {
    face_type[i] <- "smile"
  } else if (sample_percent[i] >= 50 & sample_percent[i] < 100) {
    face_type[i] <- "neutral"
  } else if (sample_percent[i] < 50) {
    face_type[i] <- "frown"
  }
}
face_type <- factor(face_type, levels=c("smile", "neutral", "frown"))

realms.results.df <- data.frame("SampPercent"=sample_percent,
                                "FaceType"=face_type)

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

realms.results.df$System <- NA
realms.results.df$System[1:24] <- "Terrestrial"
realms.results.df$System[25:48] <- "Freshwater"
realms.results.df$System[49:72] <- "Marine"

realms.results.df$System <- factor(realms.results.df$System, levels=c("Terrestrial",
                                                                      "Freshwater",
                                                                      "Marine"))

realms.results.df$Taxon <- NA
realms.results.df$Taxon <- rep(c("Aves", "Mammalia", "Herps", "Fish"), 18)

realms.results.df$Taxon <- factor(realms.results.df$Taxon, levels = c("Aves", "Mammalia", "Herps", "Fish"), labels=c("Birds", "Mammals", "Reptiles & \n amphibians", "Fishes"))
realms.results.df$Taxon2 <- factor(realms.results.df$Taxon, levels = c("Birds", "Mammals", "Reptiles & \n amphibians", "Fishes"), labels=c("Birds", "Mammals", "Reptiles & amphibians", "Fishes"))

# add current tdv to data frame
realms.results.df$TDV <- tdv.vec

# add current sample size to data frame
realms.results.df$CurrentPops <- pop.size.vec

# add required sample size to data frame
realms.results.df$MinSampSize <- samp.size.vec

# add number of needed extra pops to data frame
realms.results.df$ExtraPopsNeeded <- pops.added.vec

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

realms.results.df$MinSampSize <- samp.size.vec

Terr.results.df <- realms.results.df[realms.results.df$System=="Terrestrial",]
Terr.results.df <- Terr.results.df[Terr.results.df$Taxon!="Fishes" & Terr.results.df$Realm!="Antarctica",]
FW.results.df <- realms.results.df[realms.results.df$System=="Freshwater",]
FW.results.df <- FW.results.df[FW.results.df$Realm!="Antarctica",]
Marine.results.df <- realms.results.df[realms.results.df$System=="Marine",]
Marine.results.df$Taxon <- factor(Marine.results.df$Taxon, levels = c("Birds", "Mammals", "Reptiles & \n amphibians", "Fishes"), labels=c("Birds", "Mammals", "Reptiles", "Fishes"))

SampSizeExportDF <- data.frame(System = realms.results.df$System, 
                               Realm = realms.results.df$Realm, 
                               Taxon = realms.results.df$Taxon,
                               TDV = realms.results.df$TDV,
                               CurrentPops = realms.results.df$CurrentPops,
                               MinSampSize = realms.results.df$MinSampSize, 
                               PopsNeeded = realms.results.df$ExtraPopsNeeded)
SampSizeExportDF <- SampSizeExportDF[!is.na(SampSizeExportDF$MinSampSize),]
levels(SampSizeExportDF$Taxon) <- c("Birds", "Mammals", "Reptiles & Amphibians", "Fishes")
SampSizeExportDF <- SampSizeExportDF[SampSizeExportDF$Realm!="Antarctica",]
write.csv(SampSizeExportDF, file="SampSizeExportTable_3.csv")

median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Taxon=="Fishes"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Taxon=="Birds"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Taxon=="Mammals"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Taxon=="Reptiles & Amphibians"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$System=="Marine"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$System=="Terrestrial"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$System=="Freshwater"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Realm=="Afrotropical" |
                                      SampSizeExportDF$Realm=="Afrotropical" |
                                      SampSizeExportDF$Realm=="IndoPacific" |
                                      SampSizeExportDF$Realm=="Neotropical" |
                                      SampSizeExportDF$Realm=="Tropical Atlantic" |
                                      SampSizeExportDF$Realm=="South temperate"])
median(SampSizeExportDF$MinSampSize[SampSizeExportDF$Realm=="Palearctic" |
                                      SampSizeExportDF$Realm=="Nearctic" |
                                      SampSizeExportDF$Realm=="Temperate Atlantic" |
                                      SampSizeExportDF$Realm=="Arctic" |
                                      SampSizeExportDF$Realm=="Pacific temperate"])

