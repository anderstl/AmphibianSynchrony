library(tidyverse)
library(sf)
library(maps)

# Read in the shapefile directly as an sf object — no fortify() needed
RASY_shape <- st_read(dsn = "C:/Users/thander/OneDrive - Southern Illinois University Edwardsville/Documents/Projects/RASY/Analysis/iucn_shapefiles/",
                      layer = "data_0")

#load sites
sites<-read.csv('data/coords.csv')

# Make the sites an sf object (x = longitude, y = latitude; WGS84 decimal degrees)
sites_sf <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326)

# Basemap: lower 48 states
Low48_map <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
NAm_map <- st_as_sf(maps::map("world", region = c("USA", "Mexico", "Canada"), plot = FALSE, fill = TRUE))

png('results/fig1.png',res=600,height=3.5,width=3.5,units="in")
ggplot() +
  geom_sf(data = NAm_map, fill = "white", col = "black") +
  geom_sf(data = Low48_map, fill = "white", col = "black") +
  geom_sf(data = RASY_shape, fill = "lightgreen", alpha = 0.8) +
  geom_sf(data = sites_sf, shape = 21, fill = "mediumpurple", col = "white",
          size = 4, stroke = 0.6) +
  theme_void() +
  theme(panel.background = element_rect(fill = "azure"))+
  coord_sf(xlim = c(-95, -69), ylim = c(30, 50))
dev.off()
