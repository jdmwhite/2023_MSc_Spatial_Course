library(tidyverse)
library(rnaturalearth)
library(terra)
library(here)

# # load gbif occurrence records for Madagascar tracheophyta plants ----
# all_gbif <- read_delim(here("data/gbif/gbif_mad_tracheophyta_2023_10_27.csv"),
#                                delim = "\t", escape_double = FALSE,
#                                trim_ws = TRUE)
# # filter out records with no latitude
# all_gbif %>% filter(!is.na(decimalLatitude)) -> all_gbif_coords
# # filter records from Madagascar
# all_gbif_coords %>% filter(countryCode == 'MG') -> mad_gbif_coords
# # filter out records with no species info
# mad_gbif_coords %>% filter(!is.na(species)) -> mad_gbif_coords
# 
# # convert to a terra vect
# mad_gbif_vect <- vect(mad_gbif_coords, geom = c('decimalLongitude','decimalLatitude'))
# 
# # load Madagascar boundary
# mad_vect <- vect(ne_countries(scale = 'medium', country = 'Madagascar', returnclass = 'sf'))
# 
# # check how many data points there are?
# nrow(mad_gbif_vect)
# 
# # mask to Madagascar
# mad_gbif_mask <- mask(mad_gbif_vect, mad_vect)
# # check how many data points there are?
# nrow(mad_gbif_mask)
# 
# # convert to a data frame and save the file
# mad_gbif_df <- as.data.frame(mad_gbif_mask, geom = 'XY')
# mad_gbif_df %>% dplyr::select(family, species, lon = x, lat = y) -> mad_gbif_df
# write_csv(mad_gbif_df, 'data/gbif/gbif_mad_species.csv')

# Load the file back in and count the number of observations per species ----
mad_gbif_df <- read_csv('data/gbif/gbif_mad_species.csv')
names(mad_gbif_df)
mad_gbif_df %>% group_by(species) %>% count() %>% arrange(-n) -> spp_count

spp_count %>% left_join(dplyr::select(mad_gbif_df, species, family), by = 'species') %>% distinct() %>% dplyr::select(family, species, n) -> spp_count

# top 35 species
spp_count %>% head(35)

# remove rice:
spp_count %>% filter(species != 'Oryza sativa') -> spp_count

# top 35 species
spp_count %>% head(35)

# decide on top 35 list
(spp_count_top_35 <- spp_count[1:35,])

# save output
write_csv(spp_count_top_35, 'data/gbif/most_sampled_mad_plants.csv')