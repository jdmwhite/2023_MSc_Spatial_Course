# RGB Kew MSc Spatial Analysis Course 2023
# Joseph White, Carolina Tovar, Tarciso Leao, Felix Lim
# Species distribution model
# 2023/11/31

### Learning objectives ----
# 5.  Process environmental data
# 6.  Prepare extracted data
# 7.  Run & evaluate model
# 8.  Variable insights & predicting habitat suitability

### 5.  Process environmental data ----

#### Install and load packages

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(flexsdm)
library(corrplot)
library(sf)
library(terra)
library(here)
library(SDMtune)
library(virtualspecies)
library(exactextractr)

#### Load in data
# Name your species
spp <- c("Dracaena reflexa")
spp_code <- gsub(' ', '_', spp)

# Read species data 
all_pts <- vect(here(paste0('output/species_localities/',spp_code, '_all_points.shp')))

# using the rnaturalearth package, download the country border and return it as an sf object
mad_sf <- ne_countries(country = 'Madagascar', scale = 'medium', returnclass = 'sf')
mad <- vect(mad_sf)

#### Clean environmental data
# load worldclim data
worldclim <- rast(here("data/rast/afr_worldclim.tif"))

names(worldclim) <- c("mean_ann_t","mean_diurnal_t_range","isothermality", "t_seas", 'max_t_warm_m','min_t_cold_m',"t_ann_range",'mean_t_wet_q','mean_t_dry_q','mean_t_warm_q','mean_t_cold_q','ann_p', 'p_wet_m','p_dry_m','p_seas','p_wet_q','p_dry_q','p_warm_q','p_cold_q')

# Ensure that the worldclim data and Madagascar boundary are in the same project
mad <- project(mad, 'EPSG:4326')

# crop and mask
wc_mad <- worldclim %>% 
  crop(., mad) %>% 
  mask(., mad)

# Re-scale temperature values
wc_mad[[c(1:2,5:11)]] <- wc_mad[[c(1:2,5:11)]]/10
wc_mad[[3:4]] <- wc_mad[[3:4]]/100

#### Visualise raw data 
plot(wc_mad$mean_ann_t)
plot(mad, add = T)
plot(all_pts, add=T)

#### Check projections 
crs(wc_mad) == crs(mad)
crs(wc_mad) == crs(all_pts)
crs(mad) == crs(all_pts)
mad <- project(mad, wc_mad)
crs(wc_mad) == crs(mad)

#### Check for collinearity
# Using Pearson correlation
cov_colin <- correct_colinvar(wc_mad, method = c('pearson', th = "0.7"))
# Take a look at the correlations using corrplot
corrplot(cov_colin$cor_table, type = 'lower', diag = FALSE, order = 'hclust', tl.cex = 0.6, tl.col = 'black', addCoef.col = 'black')

# remove collinearity
set.seed(42)
non_colin <- removeCollinearity(raster::stack(wc_mad), 0.7, method = 'pearson', plot = TRUE, select.variables = TRUE, sample.points = FALSE)
non_colin

# subset variables and run again x1
wc_mad_sel <- wc_mad[[non_colin]]
set.seed(42)
non_colin_check <- removeCollinearity(raster::stack(wc_mad_sel), 0.7, method = 'pearson', plot = TRUE, select.variables = TRUE, sample.points = FALSE)

# subset variables and run again x2
wc_mad_sel <- wc_mad_sel[[non_colin_check]]
set.seed(42)
non_colin_check_2 <- removeCollinearity(raster::stack(wc_mad_sel), 0.7, method = 'pearson', plot = TRUE, select.variables = TRUE, sample.points = FALSE)

### 6.  Prepare extracted data ----

#### Extract data
# Prepare a SWD (Sample with Data), which is a class of data specifically used in the SDMtune package
all_pts_df <- as.data.frame(all_pts, geom = 'XY')

SWDdata <- prepareSWD(
  species = 'Aristida rufescens',
  p = all_pts_df %>% filter(pr_ab == 1) %>% dplyr::select(x, y),
  a = all_pts_df %>% filter(pr_ab == 0) %>% dplyr::select(x, y),
  env = wc_mad_sel
)

# Inspect 10 random rows
sample_n(cbind(SWDdata@pa, SWDdata@data), 10)

#### Train/test split
# Split locations in training (80%) and testing (20%) datasets
split_data <- trainValTest(SWDdata, test = 0.2, seed = 42)
train <- split_data[[1]]
test <- split_data[[2]]

### 7.  Run & evaluate model ----

#### Run a random forest model
set.seed(42)
rf_model <- train(method = c('RF'), ntrees = 500, data = train)
rf_model

# predict probability for training points:
cbind(train@pa, predict(rf_model, data = train)) %>% 
  as.data.frame() %>% 
  rename(true_class = V1, predicted_class = V2) %>% 
  mutate(threshold_0.4 = ifelse(predicted_class >= 0.5, 1, 0),
         threshold_0.5 = ifelse(predicted_class >= 0.5, 1, 0),
         threshold_0.6 = ifelse(predicted_class >= 0.6, 1, 0)) %>%
  DT::datatable()

#### Evaluate model
cm <- confMatrix(rf_model, test = test, th = 0.5)
cm_format <- as.data.frame(rbind(c(cm[,2],cm[,3]), c(cm[,4],cm[,5])))
names(cm_format) <- c('Positive', 'Negative')
rownames(cm_format) <- c('Positive', 'Negative')
cm_format %>% DT::datatable()

# Accuracy
round((cm$tp + cm$tn)/sum(cm$tp, cm$tp, cm$fp, cm$fn)*100, 2)

# True Skill Statistic (TSS)
tss(rf_model, test = test)

# Receiver Operator Curve (ROC)
plotROC(rf_model, test = test)

# Area Under Curve (AUC)
auc(rf_model, test = test)

# Selecting a threshold
(ths <- thresholds(rf_model, test = test))
ths %>% DT::datatable()
cm_custom <- confMatrix(rf_model, test = test, th = 0.71)

# True Positive Rate (Sensitivity)
(TPR <- cm_custom$tp/(cm_custom$tp + cm_custom$fn))

# False Positive Rate (1 - Specificity)
(FPR <- 1 - cm_custom$tn/(cm_custom$tn + cm_custom$fp))

# Receiver Operator Curve (ROC) showing best threshold
plotROC(rf_model, test = test) + 
  geom_point(aes(x = FPR, y = TPR), col = 'red', size = 3) +
  labs(x = 'False Positive Rate (1 - Specificity)', 
       y = 'True Positive Rate (Sensitivity)')

### 8.  Variable insights & predicting habitat suitability ----

#### Variable importance
vi <- varImp(rf_model, permut = 5)
plotVarImp(vi[,1:2])

#### Response curves
SDMtune::plotResponse(rf_model, var = vi$Variable[1], marginal = TRUE, rug = TRUE)
SDMtune::plotResponse(rf_model, var = vi$Variable[2], marginal = TRUE, rug = TRUE)

#### Predict habitat suitability
pred <- predict(rf_model, data = wc_mad_sel)
pred_df <- as.data.frame(pred, xy = TRUE)

ggplot() +
  geom_tile(data = pred_df, aes(x = x, y = y, fill = lyr1)) +
  scale_fill_gradientn(colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"), 
                       na.value = NA, name = 'Habitat\nsuitability') +
  geom_sf(data = mad_sf, fill = NA) +
  geom_point(data = all_pts_df %>% filter(pr_ab == 1), aes(x = x, y = y)) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme_minimal()

#### Predict habitat suitability with threshold
ths <- thresholds(rf_model, test = test)
th = ths[5,2]
# th = 0.5

pred_th <- pred >= th
pred_th_df <- as.data.frame(pred_th, xy = TRUE)

ggplot() +
  geom_tile(data = pred_th_df, aes(x = x, y = y, fill = as.factor(lyr1))) +
  scale_fill_manual(values = c("#2c7bb6", "#d7191c"), na.value = NA,
                    name = paste0('Habitat\nsuitability\nbinary (>= ',th,')'),
                    labels = c('absent','present','')) +
  geom_sf(data = mad_sf, fill = NA) +
  geom_point(data = all_pts_df %>% filter(pr_ab == 1), aes(x = x, y = y)) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme_minimal()

#### Identify overlap with Protected Areas ----
# load in protected area shapefiles. There are 3 files, so we want to load them all in together and then bind them into one file
prot_areas <- list.files(here('data/vect/WDPA_Madagascar'), pattern = '*.shp', full.names = TRUE)
prot_areas_list <- lapply(prot_areas, read_sf)
# bind the 3 files togther
prot_areas_all <- bind_rows(prot_areas_list) %>% filter(MARINE == 0)

#### convert to equal area projection
# convert the protected areas
prot_areas_all %>% 
  st_transform(crs = 'EPSG:29702') -> prot_areas_all_proj

# convert the presence/absence raster
pred_th %>% 
  project(.,vect(prot_areas_all_proj), method = 'near') -> pred_th_proj

# visualise the different projections
par(mfrow=c(1,2))
plot(pred_th)
plot(vect(prot_areas_all), add = TRUE)
plot(pred_th_proj)
plot(vect(prot_areas_all_proj), add = TRUE)

# What is the area of species presences?
# we select and sum only the cells with 1's, then multiply this by the size of the raster cells and lastly divide this by meters to get a result in km2.
pres_area <- (sum(pred_th_proj[] == 1, na.rm = TRUE) * (res(pred_th_proj)[1]*res(pred_th_proj)[2]) / (1000^2))
paste('The area of species presences is',pres_area, 'km2')

# Calculate the area of all cells
all_area <- (sum(!is.na(pred_th_proj[])) * (res(pred_th_proj)[1]*res(pred_th_proj)[2]) / (1000^2))
paste('The area of all cells is',all_area, 'km2')

# And lastly calculate the percentage of coverage of our species across all of Madagascar
paste('The species presences cover',round(pres_area/all_area*100, 2), '% of Madagascar')

#### We now want to work out what % of our species is found within Protected Areas

# create custom function to calculate the proportion of area covered by each Protected Area
sum_cover <- function(x){
  list(x %>%
         group_by(value) %>%
         summarize(total_area = sum(coverage_area)))
}

# extract the amount of area covered 
extract_all <- exact_extract(pred_th_proj, prot_areas_all_proj, coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover)
# add the names of the protected areas back on to our extraction
names(extract_all) <- prot_areas_all_proj$ORIG_NAME

# convert the list to a data frame
extract_df <- bind_rows(extract_all, .id = 'ORIG_NAME')
# take a look at the first 6 rows
head(extract_df)

# we can now sum all of the area that overlaps with the protected areas for presences (i.e. 1's) and divide this by the total area of all presences
area_under_pas <- extract_df %>% 
  filter(value == 1) %>% 
  summarise(sum(total_area)/(1000^2))

paste(round(area_under_pas/pres_area * 100, 2),'% of the predicted presences are found within protected areas')


# Our final step is to join our IUCN protected area categories onto our presence area data.frame. This will provide us with some information on what percentage of our species area is conserved under different categories. This provides important context on both the quality and quantity of protected areas overlapping with our species range:

iucn_cat <- prot_areas_all_proj %>% 
  st_drop_geometry() %>% 
  dplyr::select(ORIG_NAME, IUCN_CAT)

extract_df %>% 
  left_join(iucn_cat, by = 'ORIG_NAME', relationship = 'many-to-many') %>% 
  filter(value == 1) %>%
  group_by(IUCN_CAT) %>%
  summarise(area = sum(total_area)/(1000^2)) %>%
  mutate(perc = round(area/sum(area) * 100, 2))

#### END ####

