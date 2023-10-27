# RGB Kew MSc Spatial Analysis Course 2023
# Joseph White, Carolina Tovar, Tarciso Leao, Felix Lim
# Processing Locality Data
# 2023/11/27

### Learning objectives ----
# 1.  Download locality data
# 2.  Clean the coordinate data
# 3.  Apply spatial thinning
# 4.  Create pseudo-absences

# Load in libraries

library(rgbif)
library(DT)
library(tidyverse)
library(CoordinateCleaner)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(mapview)
library(terra)
library(flexsdm)
library(here)

### 1. Download locality data from GBIF ----

# Name your species
spp <- c("Dracaena reflexa")

# download GBIF occurrence data for this species; this takes time if there are many data points!
gbif_download <- occ_data(scientificName = spp, hasCoordinate = TRUE, country = 'MG', limit = 10000)

# select only the data (ignore the other metadata for now)
gbif_data <- gbif_download$data

# visualise
gbif_data %>% DT::datatable(extensions = 'Scroller',
                            options = list(scrollX = TRUE,
                                           scrollY = 200,
                                           scroller = TRUE))

# interactive map
gbif_data %>% 
  st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>% 
  dplyr::select(year) %>% 
  mapview(layer.name = spp)

### 2. Clean GBIF data ----
# using custom filter
gbif_data %>%
  filter(!decimalLatitude == 0 | !decimalLongitude == 0) %>%
  distinct(decimalLongitude,decimalLatitude,speciesKey,datasetKey, .keep_all = TRUE) -> gbif_filt

# using cc_*() functions
gbif_filt %>%
  cc_cen(lat = 'decimalLatitude', lon = 'decimalLongitude', buffer = 2000) %>% # remove country centroids within 2km 
  cc_cap(lat = 'decimalLatitude', lon = 'decimalLongitude', buffer = 2000) %>% # remove capitals centroids within 2km
  cc_inst(lat = 'decimalLatitude', lon = 'decimalLongitude', buffer = 2000) %>% # remove zoo and herbaria within 2km 
  cc_sea(lat = 'decimalLatitude', lon = 'decimalLongitude') -> spp_clean # remove from ocean

# Downloading country borders
# using the rnaturalearth package, download the country border and return it as an sf object
mad <- ne_countries(country = 'Madagascar', scale = 'medium', returnclass = 'sf')

# We can easily plot the spatial layers using ggplot
ggplot() +
  geom_sf(data = mad) +
  geom_point(data = gbif_data, aes(x = decimalLongitude, y = decimalLatitude), col = 'red', size = 3) +
  geom_point(data = spp_clean, aes(x = decimalLongitude, y = decimalLatitude), col = 'blue', size = 3) +
  labs(title = 'Raw vs. cleaned GBIF data') +
  theme_void()

# Clean the dataset by selecting only the columns we need, rename the lon/lat variables
spp_clean %>% 
  dplyr::select(key, year, basisOfRecord, locality, coordinateUncertaintyInMeters, occurrenceRemarks,  recordedBy, lon = decimalLongitude, lat = decimalLatitude) -> spp_sel

# Now we want to convert & project our data to the correct Coordinate Reference System
spp_sf <- st_as_sf(spp_sel, coords = c('lon', 'lat'))
st_crs(spp_sf) # as we converted these directly from GPS points, there is no set CRS

# Let's try plot this
ggplot() +
  geom_sf(data = mad) +
  geom_sf(data = spp_sf)

# provide a crs 
st_crs(spp_sf) = 4326

# Let's try plot this again
ggplot() +
  geom_sf(data = mad) +
  geom_sf(data = spp_sf, size = 3, alpha = 0.5) +
  labs(title = 'Plotting as simple features (sf)') +
  theme_void()

# Interactive maps
spp_sf %>% 
  st_jitter(factor = 0.0001) %>% 
  mapview(zcol = 'year')

### 3. Spatial thinning ----

# Load in environmental variables
# alternative to download worldclim data
# worldclim <- rast(raster::getData("worldclim", var="bio", res=2.5))

# load in worldclim locally
worldclim <- rast(here("data/rast/afr_worldclim.tif"))
names(worldclim)
names(worldclim) <- c("mean_ann_t","mean_diurnal_t_range","isothermality", "t_seas", 'max_t_warm_m','min_t_cold_m',"t_ann_range",'mean_t_wet_q','mean_t_dry_q','mean_t_warm_q','mean_t_cold_q','ann_p', 'p_wet_m','p_dry_m','p_seas','p_wet_q','p_dry_q','p_warm_q','p_cold_q')
plot(worldclim$mean_ann_t, main = 'mean annual temperature (C)')

# Ensure that the worldclim data and Madagascar boundary are in the same project
mad <- st_transform(mad, 'EPSG:4326')

# Crop and mask
wc_mad <- worldclim %>% 
  crop(., vect(mad)) %>% 
  mask(., vect(mad))
plot(wc_mad)

# Thin records by distance
set.seed(42)
spp_filt <- occfilt_geo(
  data = spp_sel %>% rename(x = lon, y = lat),
  x = 'x',
  y = 'y',
  method = c('defined', d = 10), # set a 10 km threshold
  env_layer = wc_mad,
  prj = crs(wc_mad)
)

# join the result of the filter back onto the metadata and convert to sf
spp_filt_sf <- spp_filt %>% 
  st_as_sf(coords = c('x','y'), crs = st_crs(4326))

# identify the localities that were removed
spp_removed_sf <- spp_sel %>% 
  filter(!key %in% spp_filt$key) %>% 
  st_as_sf(coords = c('lon','lat'), crs = st_crs(4326))

# Visualise output
mapview(spp_filt_sf, zcol = NULL, col.regions = 'blue', layer.name = 'localities included') +
  mapview(spp_removed_sf, zcol = NULL, col.regions = 'red', layer.name = 'localities filtered')

### 4. Create pseudo-absences ----
set.seed(42)
pseudo_abs <-  sample_pseudoabs(
  data = spp_filt,
  x = 'x',
  y = 'y',
  n = nrow(spp_filt),
  method = c('geo_const', width = 10000), # threshold of 10 km
  rlayer = wc_mad,
  calibarea = vect(mad)
)
head(pseudo_abs)

# process our presence points to match pseudoabsences
pres <- spp_filt %>% 
  dplyr::select(x, y) %>% 
  mutate(pr_ab = 1)
head(pres)

# combine presences and pseudoabsences
all_pts <- rbind(pres, pseudo_abs)
# convert to spatial format
all_pts_sf <- all_pts %>% 
  mutate(pr_ab = as.factor(pr_ab)) %>% 
  st_as_sf(coords = c('x','y'), crs = st_crs(4326))

# visualise
mapview(all_pts_sf, zcol = 'pr_ab', layer.name = 'all points')

#### Export locality data
write_sf(all_pts_sf, 
         here(paste0('output/species_localities/',gbif_data$genus[1],'_',gbif_data$specificEpithet[1],'_','all_points.shp')))

#### END ####
