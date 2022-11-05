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

