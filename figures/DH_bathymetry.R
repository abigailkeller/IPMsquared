library(ncdf4)
library(tidyverse)
library(reshape2)
library(terra)
library(sf)
library(ggspatial)

ss_raster <- rast('data/DH/bathymetry/SS_Bathymetry.tif_/SS_Bathymetry.tif')
#plot(ss_raster)

ss_raster_latlon <- project(ss_raster,"+init=epsg:4326")

# cut to Drayton Harbor
DH_e <- ext(-122.819111, -122.699388, 48.946284, 49.000274)
DH_bath <- crop(ss_raster_latlon, DH_e)

# read in coastline shapefile
coastline <- vect("data/DH/N45W125")
DH_coast <- crop(coastline, DH_e)

# Convert SpatRaster to a data frame
DH_bath_df <- as.data.frame(DH_bath, xy = TRUE)

# Convert SpatVector to an sf object
DH_coast_sf <- st_as_sf(DH_coast)

# Plot 
map <- ggplot() +
  geom_raster(data = DH_bath_df, 
              aes(x = x, y = y, fill = SS_Bathymetry)) +  
  geom_sf(data = DH_coast_sf, color = "black", 
          fill = NA, linewidth = 1) +            
  scale_fill_viridis_c() +                                     
  theme_minimal() +                                            
  labs(fill = "depth (m)", x = 'Longitude', y = 'Latitude')+
  annotation_scale(location = "bl", width_hint = 0.25)
ggsave('figures/supplemental_map.png',map,
       dpi=400,width=7,height=3,bg='white')

