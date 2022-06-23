library(sf)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggsci)
library(svglite)

# make terrestrial realm map
terr_eco_shp <- st_read(dsn="Data/wwf_terr_ecos/wwf_terr_ecos.shp")
terr_eco_shp_robin <- st_transform(terr_eco_shp, crs = "+proj=robin")
terr_eco_shp_robin$REALM[terr_eco_shp_robin$REALM=="AA" | terr_eco_shp_robin$REALM=="IM" | terr_eco_shp_robin$REALM=="OC"] <-"IP" 
terr_eco_shp_robin_realm <- terr_eco_shp_robin %>%
  group_by(REALM) %>% 
  dplyr::summarize()

#world_laea2 = st_transform(terr_eco_shp, crs = "+proj=laea +x_0=0 +y_0=0 +lon_0=-74 +lat_0=40")
#world_laea2 = st_transform(terr_eco_shp, crs = "+proj=mercator")
cbPalette <- c("dark green", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#CC79A7")
ggplot()+
  geom_sf(data=terr_eco_shp_robin_realm, aes(fill=REALM))+
  coord_sf()+
  theme_bw()+
  scale_fill_manual(values=cbPalette, labels=c("Antarctica", 
                                             "Afrotropic", 
                                             "Indopacific", 
                                             "Nearctic",
                                             "Neotropic",
                                             "Palearctic",
                                             "N/A"),
                    guide = "none")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey80"))

ggsave("terr_fw_realms4.tiff",
       plot = last_plot(),
       device = "tiff",
       width = 16000,
       height = 8000,
       units = "px",
       dpi = 1000,
       compression = "lzw")

###########

# make marine realm map
marine_eco_shp <- st_read(dsn="Data/wwf_marine/meow_ecos.shp")
marine_eco_shp_robin <- st_transform(marine_eco_shp, crs = "+proj=robin")
marine_eco_shp_robin$REALM[marine_eco_shp_robin$REALM=="Tropical Eastern Pacific" | 
                             marine_eco_shp_robin$REALM=="Eastern Indo-Pacific" | 
                             marine_eco_shp_robin$REALM=="Central Indo-Pacific" |
                             marine_eco_shp_robin$REALM=="Western Indo-Pacific"] <-"TroSubIndo"
marine_eco_shp_robin$REALM[marine_eco_shp_robin$REALM=="Temperate Southern Africa" | 
                             marine_eco_shp_robin$REALM=="Temperate South America" | 
                             marine_eco_shp_robin$REALM=="Temperate Australasia" |
                             marine_eco_shp_robin$REALM=="Southern Ocean"] <-"SoTeAnt"

#st_is_valid(marine_eco_shp_robin, reason=TRUE)
marine_eco_shp_robin <- marine_eco_shp_robin[-c(69,71),]
marine_eco_shp_robin_realm <- marine_eco_shp_robin %>%
  group_by(REALM) %>% 
  dplyr::summarize()

cbPalette <- c("dark green", "#F0E442", "#0072B2", "#E69F00", "#009E73", "#CC79A7")
ggplot()+
  geom_sf(data=marine_eco_shp_robin_realm, aes(fill=REALM))+
  geom_sf(data=terr_eco_shp_robin_realm, aes())+
  coord_sf()+
  scale_fill_manual(values=cbPalette, labels=c("Arctic", 
                                               "South temperate",
                                               "Temperate Atlantic",
                                               "Pacific temperate",
                                               "Tropical Atlantic",
                                               "IndoPacific"),
                    guide = "none")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey80"))

ggsave("marine_realms4.tiff",
       plot = last_plot(),
       device = "tiff",
       width = 16000,
       height = 8000,
       units = "px",
       dpi = 1000,
       compression = "lzw")
  
########

terr_eco_shp <- st_read(dsn="Data/ecoregions2017/Ecoregions2017.shp")
terr_eco_shp
world_laea2 = st_transform(terr_eco_shp, crs = "+proj=laea +x_0=0 +y_0=0 +lon_0=-74 +lat_0=40")
world_laea2 = st_transform(terr_eco_shp, crs = "+proj=robin")
world_laea2 <- world_laea2[world_laea2$REALM!="Antarctica",]

world_laea3 <- world_laea2 %>%
  group_by(REALM) %>% 
  summarize(across(geometry, ~ st_union(.)), .groups = "keep")
ggplot()+
  geom_sf(data=world_laea2, aes(fill=REALM))+
  coord_sf()




sample_percent <- samp.percent.vec # import results from RWLPI_NG_2.R
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
write.csv(SampSizeExportDF, file="SampSizeExportTable.csv")

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

####################################################
# associate stored images with value categories
emoji_pic <- data.frame(
  FaceType = c("smile", "neutral", "frown"),
  emoji_link = c("Data/images/285740s.png",
                 "Data/images/285741s.png",
                 "Data/images/285742s.png"))

animal_pic <- data.frame(
  Taxon = c("Aves", "Mammalia", "Herps", "Fish"),
  animal_link = c("Data/images/292888s.png",
                  "Data/images/32581s.png",
                  "Data/images/319830s.png",
                  "Data/images/105151s.png")
)

# create a function to link them in ggplot
func_link_to_img <- function(x, size = 16) {
  paste0("<img src='", x, "' width='", size, "'/>")
}

# add image links to results data frame
realms.results.df <- left_join(emoji_pic, realms.results.df, by = c("FaceType"="Facetype"))

# add image links to results data frame
realms.results.df <- left_join(animal_pic, realms.results.df, by = c("Taxon"="Taxon"))
