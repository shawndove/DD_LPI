######################
### Author: Shawn Dove
######################

# This script creates maps of terrestrial and marine LPI realms


# load packages ----

library(sf)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggsci)
library(svglite)


# set directory paths ----

md_name <- "MapData" # name of directory where map files are stored
pd_name <- "Plots" # name of directory where .RData files are stored


# make terrestrial realm map ----

# read the shape file
terr_eco_shp <- st_read(dsn=paste(md_name, "/wwf_terr_ecos/wwf_terr_ecos.shp", sep=""))

# convert to robinson projection
terr_eco_shp_robin <- st_transform(terr_eco_shp, crs = "+proj=robin")

# combine Australasia, IndoMalay, and Oceania into the IndoPacific realm to match LPI realms
terr_eco_shp_robin$REALM[terr_eco_shp_robin$REALM=="AA" | terr_eco_shp_robin$REALM=="IM" | terr_eco_shp_robin$REALM=="OC"] <-"IP" 

# combine objects into realm polygons
terr_eco_shp_robin_realm <- terr_eco_shp_robin %>%
  group_by(REALM) %>% 
  dplyr::summarize()

# set the palette
cbPalette <- c("dark green", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#CC79A7")

# create the map with ggplot
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

# save the map
ggsave(filename=paste(pd_name, "/terr_fw_realms.tiff", sep=""),
       plot = last_plot(),
       device = "tiff",
       width = 16000,
       height = 8000,
       units = "px",
       dpi = 1000,
       compression = "lzw")


###########

# make marine realm map ----

# read shape file
marine_eco_shp <- st_read(dsn=paste(md_name, "/wwf_marine/meow_ecos.shp", sep=""))

# convert to robinson projection
marine_eco_shp_robin <- st_transform(marine_eco_shp, crs = "+proj=robin")

# combine 4 realms into the Tropical and subtropical Indo-Pacific realm, as in the LPI
marine_eco_shp_robin$REALM[marine_eco_shp_robin$REALM=="Tropical Eastern Pacific" | 
                             marine_eco_shp_robin$REALM=="Eastern Indo-Pacific" | 
                             marine_eco_shp_robin$REALM=="Central Indo-Pacific" |
                             marine_eco_shp_robin$REALM=="Western Indo-Pacific"] <-"TroSubIndo"

# combine 4 realms into the South temperate and Antarctic realm, as in the LPI
marine_eco_shp_robin$REALM[marine_eco_shp_robin$REALM=="Temperate Southern Africa" | 
                             marine_eco_shp_robin$REALM=="Temperate South America" | 
                             marine_eco_shp_robin$REALM=="Temperate Australasia" |
                             marine_eco_shp_robin$REALM=="Southern Ocean"] <-"SoTeAnt"

# Remove 2 objects which otherwise cause an error in the next step:
# Fiji Islands and Gilbert/Ellis Islands. While removing objects is not ideal,
# small islands are of little consequence to the overall map.)
marine_eco_shp_robin <- marine_eco_shp_robin[-c(69,71),]

# combine objects into realm polygons
marine_eco_shp_robin_realm <- marine_eco_shp_robin %>%
  group_by(REALM) %>% 
  dplyr::summarize()

#set the palette
cbPalette <- c("dark green", "#F0E442", "#0072B2", "#E69F00", "#009E73", "#CC79A7")

# create the map in ggplot
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

# save the map
ggsave(filename=paste(pd_name, "/marine_realms.tiff", sep=""),
       plot = last_plot(),
       device = "tiff",
       width = 16000,
       height = 8000,
       units = "px",
       dpi = 1000,
       compression = "lzw")
  
