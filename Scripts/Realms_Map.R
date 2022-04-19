library(sf)
library(ggplot2)
library(dplyr)

terr_eco_shp <- st_read(dsn="Data/wwf_terr_ecos/wwf_terr_ecos.shp")
terr_eco_shp %>%
  group_by(REALM) %>% 
  summarize()
st_geometry_type(terr_eco_shp)
st_crs(terr_eco_shp)
st_bbox(terr_eco_shp)
terr_eco_shp

world_laea2 = st_transform(terr_eco_shp, crs = "+proj=laea +x_0=0 +y_0=0 +lon_0=-74 +lat_0=40")
world_laea2 = st_transform(terr_eco_shp, crs = "+proj=mercator")

ggplot()+
  geom_sf(data=terr_eco_shp, aes(fill=REALM))+
  coord_sf()


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
