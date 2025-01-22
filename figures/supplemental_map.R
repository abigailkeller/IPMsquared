library(ncdf4)
library(tidyverse)
library(reshape2)
library(terra)
library(sf)
library(ggspatial)

ss_raster <- rast("data/DH/bathymetry/SS_Bathymetry.tif_/SS_Bathymetry.tif")

ss_raster_latlon <- project(ss_raster, "+init=epsg:4326")

# cut to Drayton Harbor
dh_ext <- ext(-122.819111, -122.699388, 48.946284, 49.000274)
dh_bath <- crop(ss_raster_latlon, dh_ext)

# read in coastline shapefile
coastline <- vect("data/DH/N45W125")
dh_coast <- crop(coastline, dh_ext)

# Convert SpatRaster to a data frame
dh_bath_df <- as.data.frame(dh_bath, xy = TRUE)

# Convert SpatVector to an sf object
dh_coast_sf <- st_as_sf(dh_coast)

# Plot
map <- ggplot() +
  geom_raster(data = dh_bath_df,
              aes(x = x, y = y, fill = SS_Bathymetry)) +
  geom_sf(data = dh_coast_sf, color = "black",
          fill = NA, linewidth = 1) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(fill = "depth (m)", x = "Longitude", y = "Latitude") +
  annotation_scale(location = "bl", width_hint = 0.25)
ggsave("figures/supplemental_map.png", map,
       dpi = 400, width = 7, height = 3, bg = "white")
