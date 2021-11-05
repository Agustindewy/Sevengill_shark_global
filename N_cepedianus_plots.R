
# De Wysiecki et al. - Global ENM projection for the broadnose sevengill shark (Notorynchus cepedianus)

# Additional analyses and figures
# Each figure is self-sufficient, so expect high repetition of lines in script

library(rgdal)
library(rgeos)
library(raster)
library(spThin)
library(ggplot2)
library(kuenm)
library(ellipsenm)
library(rgl)
library(ggspatial)
library(ggsn)
library(viridis)

setwd('SET YOUR WORKING DIRECTORY')

# Global and regional coastlines from the Global Self-consistent, Hierarchical, High-resolution Geography Database
# Dowloaded from https://www.soest.hawaii.edu/pwessel/gshhg/
coast <- readOGR(dsn = 'DOWNLOAD AND READ', layer = 'GSHHS_f_L1_World')
coast <- gBuffer(coast, byid = T, width = 0)
coast1 <- readOGR(dsn = 'DOWNLOAD AND READ', layer = 'GSHHS_f_L1_SouthAmerica')
coast2 <- readOGR(dsn = 'DOWNLOAD AND READ', layer = 'GSHHS_f_L1_Australia')

# Projection
crs <- CRS('+init=epsg:4326')

# Species
Species <- 'Notorynchus cepedianus'

# Function for 10% Minimum Training Presence calculation 
# Taken from https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/
sdm_threshold.10 <- function(sdm, occs, type = 'mtp', binary = FALSE){
  occPredVals <- raster::extract(sdm, occs)
  if(type == 'mtp'){
    thresh <- min(na.omit(occPredVals))
  } else if(type == 'p10'){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
} # MTP 10%

# Functions for response plots
.get_presence <- function(swd) {
  return(swd@data[swd@pa == 1, , drop = F])
}
.get_absence <- function(swd) {
  return(swd@data[swd@pa == 0, , drop = F])
}


#----------------------------------- Figure S1 ---------------------------------------------

# Systematic spatial thinning in SWA

dat <- read.csv('Appendix C.csv', header = T) # mind not all data is available
dat <- subset(dat, Region == 'Southwest Atlantic')
swa_train <- subset(dat, Type == 'calibration')

# Visually find the minimum thinning distance that removes spatial clusters in data
thinned_data <- data.frame()
for(i in seq(10, 50, 10)){
  swa_train_thin <- thin(swa_train, lat.col = 'Latitude', long.col = 'Longitude', spec.col = 'Species',
                         thin.par = i, reps = 1000, locs.thinned.list.return = T, write.files = F, write.log.file = F)
  thinned_data0 <- data.frame(Lon = swa_train_thin[[1]]$Longitude, Lat = swa_train_thin[[1]]$Latitude,
                              Dist = as.factor(paste(i, 'km')))
  thinned_data <- rbind(thinned_data, thinned_data0)
}

# Calibration area (M space)
env <- stack('predictors.tif')
coord <- cbind(dat$Longitude, dat$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point
env.M <- crop(env, occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon

# Erase Pacific areas
lon1 <- c(-84,-69.4,-69.4,-84,-84) 
lat1 <- c(-32,-32,-62,-62,-32)
pol1 <- cbind(lon1, lat1)
Pol1 <- Polygon(pol1, hole = F)
list1 <- list(Pol1)
poly1 <- Polygons(list1, ID = 'p1')
spp1 <- SpatialPolygons(list(poly1), proj4string = crs)
env.M <- mask(env.M, spp1, inverse = T)
Ext1 <- extent(-69.4,-39,-61.4,-22.8)
env.M <- crop(env.M, Ext1)
env.M <- stack(env.M)

empty_raster <- env.M[[1]]
empty_raster[!is.na(empty_raster)] = 1
raster.df <- data.frame(coordinates(empty_raster), values = values(empty_raster))

ggplot(data = thinned_data, aes(x = Lon, y = Lat)) +
  geom_tile(data = raster.df, aes(x = x, y = y, fill = as.factor(values)), inherit.aes = F) +
  geom_point(shape = 21, size = 0.5, color = 'blue') +
  scale_fill_manual(values = 'white', na.value = 'grey95') + coord_equal(expand = F) +
  theme(legend.position = 'none') + labs(x = 'Longitude', y = 'Latitude') +
  facet_wrap(~ Dist, ncol = 5)
ggsave('Figure S1.pdf', width = 30, height = 9, units = 'cm')


#----------------------------------- Figure S2 ---------------------------------------------

# Systematic spatial thinning in AUS

dat <- read.csv('Appendix C.csv', header = T)
dat <- subset(dat, Region == 'Australia')
aus_train <- subset(dat, Type == 'calibration')

# Visually find the minimum thinning distance that removes spatial clusters in data
thinned_data <- data.frame()
for(i in seq(10, 50, 10)){
  aus_train_thin <- thin(aus_train, lat.col = 'Latitude', long.col = 'Longitude', spec.col = 'Species',
                         thin.par = i, reps = 1000, locs.thinned.list.return = T, write.files = F, write.log.file = F)
  thinned_data0 <- data.frame(Lon = aus_train_thin[[1]]$Longitude, Lat = aus_train_thin[[1]]$Latitude,
                              Dist = as.factor(paste(i, 'km')))
  thinned_data <- rbind(thinned_data, thinned_data0)
}

# Calibration area (M space)
env <- stack('predictors.tif')
coord <- cbind(dat$Longitude, dat$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point
env.M <- crop(env, occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon

empty_raster <- env.M[[1]]
empty_raster[!is.na(empty_raster)] = 1
raster.df <- data.frame(coordinates(empty_raster), values = values(empty_raster))

ggplot(data = thinned_data, aes(x = Lon, y = Lat)) +
  geom_tile(data = raster.df, aes(x = x, y = y, fill = as.factor(values)), inherit.aes = F) +
  geom_point(shape = 21, size = 0.5, color = 'blue') +
  scale_fill_manual(values = 'white', na.value = 'grey95') + coord_equal(expand = F) +
  theme(legend.position = 'none') + labs(x = 'Longitude', y = 'Latitude') +
  facet_wrap(~ Dist, ncol = 3)

ggsave('Figure S2.pdf', width = 30, height = 14, units = 'cm')


#----------------------------------- Figure S3 and S4 ---------------------------------------------

# Depth and slope discrepancies between SWA and AUS 

Region1 <- 'SWA'
Region2 <- 'AUS'

# Read predictors
env <- stack('predictors.tif') # create raster in data preparation
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope')
names(env) <- var_names

# Using three variables to simplify analysis
Vars1 <- c('Bathymetry', 'Surface_temperature', 'Kd490')
Vars2 <- c('Kd490', 'Surface_temperature', 'Slope')
Vars3 <- c('Distance_to_coast', 'Surface_temperature', 'Kd490')

# Occurrences
dat1 <- read.csv(paste(Region1, '/N_cepedianus_joint.csv', sep = ''), header = T)
dat2 <- read.csv(paste(Region2, '/N_cepedianus_joint.csv', sep = ''), header = T)

# M space SWA
coord1 <- cbind(dat1$longitude, dat1$latitude)
occ.tot1 <- SpatialPoints(coord1, proj4string = crs) 
occ.buff1 <- buffer(occ.tot1, width = 1000000) 
env.M1 <- crop(env, occ.buff1)
env.M1 <- mask(env.M1, occ.buff1)

# Erase Pacific areas
lon1 <- c(-84,-69.4,-69.4,-84,-84) 
lat1 <- c(-32,-32,-62,-62,-32)
pol1 <- cbind(lon1, lat1)
Pol1 <- Polygon(pol1, hole = F)
list1 <- list(Pol1)
poly1 <- Polygons(list1, ID = 'p1')
spp1 <- SpatialPolygons(list(poly1), proj4string = crs)
env.M1 <- mask(env.M1, spp1, inverse = T)
Ext1 <- extent(-69.4,-39,-61.4,-22.8)
env.M1 <- crop(env.M1, Ext1)
env.M1 <- stack(env.M1)

# M space AUS
coord2 <- cbind(dat2$longitude, dat2$latitude)
occ.tot2 <- SpatialPoints(coord2, proj4string = crs) 
occ.buff2 <- buffer(occ.tot2, width = 1000000) 
env.M2 <- crop(env, occ.buff2)
env.M2 <- mask(env.M2, occ.buff2)
env.M2 <- stack(env.M2)

# Preparing overlap objects to perform analyses
niche1_depth <- overlap_object(data = dat1, species = 'species', longitude = 'longitude', latitude = 'latitude', 
                               method = 'mve1', level = 95, variables = env.M1[[Vars1]])
niche2_depth <- overlap_object(data = dat2, species = 'species', longitude = 'longitude', latitude = 'latitude', 
                               method = 'mve1', level = 95, variables = env.M2[[Vars1]])

niche1_slope <- overlap_object(data = dat1, species = 'species', longitude = 'longitude', latitude = 'latitude', 
                               method = 'mve1', level = 95, variables = env.M1[[Vars2]])
niche2_slope <- overlap_object(data = dat2, species = 'species', longitude = 'longitude', latitude = 'latitude', 
                               method = 'mve1', level = 95, variables = env.M2[[Vars2]])

niche1_dist <- overlap_object(data = dat1, species = 'species', longitude = 'longitude', latitude = 'latitude', 
                               method = 'mve1', level = 95, variables = env.M1[[Vars3]])
niche2_dist <- overlap_object(data = dat2, species = 'species', longitude = 'longitude', latitude = 'latitude', 
                               method = 'mve1', level = 95, variables = env.M2[[Vars3]])

# Overlap
N.test_depth <- ellipsoid_overlap(niche1_depth, niche2_depth, overlap_type = 'back_union', 
                                  significance_test = T, replicates = 1000)

N.test_slope <- ellipsoid_overlap(niche1_slope, niche2_slope, overlap_type = 'back_union', 
                                  significance_test = T, replicates = 1000)

N.test_dist <- ellipsoid_overlap(niche1_dist, niche2_dist, overlap_type = 'back_union', 
                                  significance_test = T, replicates = 1000)

# Plots: ellipsoid overlap and statistical significance test
# Depth
rgl.open(); rgl.bg(color = 'white')
plot_overlap(N.test_depth, niches = c(1, 2), data = T, background = T, proportion = 1, 
             background_type = 'back_union', change_labels = T,
             data_col = c('darkred', 'darkorange'), niche_col = c('darkred', 'darkorange'), 
             background_col = viridis::cividis, legend = F)
decorate3d(aspect = c(1, 1, 1), xlab = '', ylab = '', zlab = '') 
axes3d(xat = c(-4000, -2000, 0, 2000), yat = c(20, 15, 10), zat = c(0, 0.2, 0.4), box = T)
mtext3d('Bathymetry (m)', edge = 'x-+', line = 1.25)
mtext3d('Sea surface temperature (ºC)', edge = 'y-+', line = 0, at = 12)
mtext3d('Diffuse attenuation Kd490 (1/m)', edge = 'z+-', line = 0, at = 0.32)
snapshot3d('Ellipsoid_depth.png', fmt = 'png', top = T, width = 1000, height = 1000, webshot = F)

ggplot() +
  geom_histogram(aes(N.test_depth@significance_results$union_random$Niche_1_vs_2$overlap), 
                 fill = '#00204DFF', bins = 60) +
  geom_vline(xintercept = quantile(N.test_depth@significance_results$union_random$Niche_1_vs_2$overlap, 0.05),
             colour = '#FFEA46FF', size = 0.75, linetype = 'dashed') +  
  geom_vline(xintercept = N.test_depth@union_overlap$overlap, 
             colour = '#7C7B78FF', size = 0.75, linetype = 'dashed') + 
  theme_bw() + xlab('Overlap') + ylab('Frequency') 
ggsave('Test_depth.pdf', dpi = 900, width = 12, height = 10, units = 'cm')

# Slope
rgl.open(); rgl.bg(color = 'white')
plot_overlap(N.test_slope, niches = c(1, 2), data = T, background = T, proportion = 1, 
             background_type = 'back_union', change_labels = T,
             data_col = c('darkred', 'darkorange'), niche_col = c('darkred', 'darkorange'), 
             background_col = viridis::cividis, legend = F)
decorate3d(aspect = c(1, 1, 1), xlab = '', ylab = '', zlab = '') 
axes3d(xat = c(0, 0.2, 0.4), yat = c(20, 15, 10), zat = c(-2, 2, 6), box = T)
mtext3d('Sea surface temperature (ºC)', edge = 'y++', line = 0, at = 12)
mtext3d('Diffuse attenuation Kd490 (1/m)', edge = 'x--', line = 0, at = 0.32)
mtext3d('Slope (º)', edge = 'z+-', line = 1.25)
snapshot3d('Ellipsoid_slope.png', fmt = 'png', top = T, width = 1000, height = 1000, webshot = F)

ggplot() +
  geom_histogram(aes(N.test_slope@significance_results$union_random$Niche_1_vs_2$overlap), 
                 fill = '#00204DFF', bins = 60) +
  geom_vline(xintercept = quantile(N.test_slope@significance_results$union_random$Niche_1_vs_2$overlap, 0.05),
             colour = '#FFEA46FF', size = 0.75, linetype = 'dashed') +  
  geom_vline(xintercept = N.test_slope@union_overlap$overlap, 
             colour = '#7C7B78FF', size = 0.75, linetype = 'dashed') + 
  theme_bw() + xlab('Overlap') + ylab('Frequency') 
ggsave('Test_slope.pdf', dpi = 900, width = 12, height = 10, units = 'cm')

# Distance to coast
rgl.open(); rgl.bg(color = 'white')
plot_overlap(N.test_dist, niches = c(1, 2), data = T, background = T, proportion = 1, 
             background_type = 'back_union', change_labels = T,
             data_col = c('darkred', 'darkorange'), niche_col = c('darkred', 'darkorange'), 
             background_col = viridis::cividis, legend = F)
decorate3d(aspect = c(1, 1, 1), xlab = '', ylab = '', zlab = '') 
axes3d(xat = c(-50, 50, 150, 250), yat = c(20, 15, 10), zat = c(0, 0.2, 0.4), box = T)
mtext3d('Sea surface temperature (ºC)', edge = 'y+-', line = 0, at = 18)
mtext3d('Diffuse attenuation Kd490 (1/m)', edge = 'z--', line = 0, at = 0.09)
mtext3d('Distance to coast (km)', edge = 'x--', line = 1.25)
snapshot3d('Ellipsoid_dist.png', fmt = 'png', top = T, width = 1000, height = 1000, webshot = F)

ggplot() +
  geom_histogram(aes(N.test_dist@significance_results$union_random$Niche_1_vs_2$overlap), 
                 fill = '#00204DFF', bins = 60) +
  geom_vline(xintercept = quantile(N.test_dist@significance_results$union_random$Niche_1_vs_2$overlap, 0.05),
             colour = '#FFEA46FF', size = 0.75, linetype = 'dashed') +  
  geom_vline(xintercept = N.test_dist@union_overlap$overlap, 
             colour = '#7C7B78FF', size = 0.75, linetype = 'dashed') + 
  theme_bw() + xlab('Overlap') + ylab('Frequency') 
ggsave('Test_dist.pdf', dpi = 900, width = 12, height = 10, units = 'cm')


#----------------------------------- Figure 1 ---------------------------------------------

# SWA calibration area 

Region1 <- 'SWA'
Region2 <- 'swa'

final_mod0 <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/calibration_median.tif', sep = '')) # model output
final_mod <- extend(final_mod0, extent(final_mod0) + 1.75)
occ_train <- read.csv(paste(Region1, '/N_cepedianus_joint.csv', sep = ''), header = T) # calibration occurrences
occ_indep <- read.csv(paste(Region1, '/N_cepedianus_ind.csv', sep = ''), header = T) # independent occurrences

# Calibration area
dat <- read.csv('Appendix C.csv', header = T) # mind missing data
dat <- subset(dat, Region == 'Southwest Atlantic')
coord <- cbind(dat$Longitude, dat$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs)
M.buff <- buffer(occ.tot, width = 1000000)
M.buff <- crop(M.buff, final_mod0)
swa_coast <- crop(coast1, final_mod) # crop coastline to calibration area to save time

df <- data.frame(coordinates(final_mod), as.data.frame(final_mod))

ggplot(data = df) +
  geom_tile(aes(x = x, y = y, fill = calibration_median)) + 
  geom_polygon(data = M.buff, aes(x = long, y = lat, group = group), size = 1, color = 'grey5', fill = NA) +
  geom_polygon(data = swa_coast, aes(x = long, y = lat, group = group), color = 'grey50', fill = 'grey50') +
  geom_point(data = occ_train, aes(x = longitude, y = latitude), size = 1.25, color = 'black', fill = 'black') +
  geom_point(data = occ_indep, aes(x = longitude, y = latitude), size = 1.25, color = 'grey95', fill = 'grey95') + 
  scale_fill_viridis_c(option = 'viridis', na.value = 'grey95') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 10),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = 'Latitude (S)', breaks = seq(-30, -60, -10), labels = c('30º', '40º', '50º', '60º')) + 
  scale_x_continuous(name = 'Longitude (W)', breaks = seq(-40, -70, -10), labels = c('40º', '50º', '60º', '70º')) +
  ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T, 
                 dist = 200, st.size = 3, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -40.1, y = -60.3))
ggsave('Figure 1.tiff', dpi = 900, width = 12.8, height = 15, units = 'cm')


#----------------------------------- Figure 2 ---------------------------------------------

# AUS calibration area

Region1 <- 'AUS'
Region2 <- 'aus'

final_mod0 <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/calibration_median.tif', sep = '')) # model output
final_mod <- extend(final_mod0, extent(final_mod0) + 1.75)
occ_train <- read.csv(paste(Region1, '/N_cepedianus_joint.csv', sep = ''), header = T) # calibration occurrences
occ_indep <- read.csv(paste(Region1, '/N_cepedianus_ind.csv', sep = ''), header = T) # independent occurrences

# Calibration area
dat <- read.csv('Appendix C.csv', header = T)
dat <- subset(dat, Region == 'Australia')
coord <- cbind(dat$Longitude, dat$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs)
M.buff <- buffer(occ.tot, width = 1000000)
M.buff <- crop(M.buff, final_mod0)
aus_coast <- crop(coast2, final_mod) # crop coastline to calibration area to save time

df <- data.frame(coordinates(final_mod), as.data.frame(final_mod))

ggplot(data = df) +
  geom_tile(aes(x = x, y = y, fill = calibration_median)) + 
  geom_polygon(data = M.buff, aes(x = long, y = lat, group = group), size = 1, color = 'grey5', fill = NA) +
  geom_polygon(data = aus_coast, aes(x = long, y = lat, group = group), color = 'grey50', fill = 'grey50') +
  geom_point(data = occ_train, aes(x = longitude, y = latitude), size = 1.25, color = 'black', fill = 'black') +
  geom_point(data = occ_indep, aes(x = longitude, y = latitude), size = 1.25, color = 'grey95', fill = 'grey95') + 
  scale_fill_viridis_c(option = 'viridis', na.value = 'grey95') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 10),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = 'Latitude (S)', breaks = seq(-30, -60, -10), labels = c('30º', '40º', '50º', '60º')) + 
  scale_x_continuous(name = 'Longitude (W)', breaks = seq(120, 160, 10), labels = c('120º', '130º', '140º', '150º', '160º')) +
  ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T, 
                 dist = 200, st.size = 3, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = 161.5, y = -52.6))
ggsave('Figure 2.tiff', dpi = 900, width = 18.4, height = 11.7, units = 'cm')

    
#----------------------------------- Figure S5 ---------------------------------------------

# Response plots for SWA model

# Install and call 'SDMtune' package (https://github.com/ConsBiol-unibern/SDMtune)

Region <- 'SWA'

# Read predictors
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope')
names(env) <- var_names

# Subset selected predictors
env <- env[[c('Temperature', 'Surface_temperature', 'Kd490', 'Distance_to_coast')]]

# Occurrences
dat <- read.csv(paste(Region, '/N_cepedianus_joint.csv', sep = ''), header = T)
dat <- dat[, c('longitude', 'latitude')]

# M space
occ.tot <- SpatialPoints(dat, proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env, occ.buff)
env.M <- mask(env.M, occ.buff)

# Erase Pacific areas
lon1 <- c(-84,-69.4,-69.4,-84,-84) 
lat1 <- c(-32,-32,-62,-62,-32)
pol1 <- cbind(lon1, lat1)
Pol1 <- Polygon(pol1, hole = F)
list1 <- list(Pol1)
poly1 <- Polygons(list1, ID = 'p1')
spp1 <- SpatialPolygons(list(poly1), proj4string = crs)
env.M <- mask(env.M, spp1, inverse = T)
Ext1 <- extent(-69.4,-39,-61.4,-22.8)
env.M <- crop(env.M, Ext1)
env.M <- stack(env.M)

p_df = data.frame()
a_df = data.frame()
plot_data_df = data.frame()

for(j in names(env.M)){
  # prepare data
  env = env.M[[j]]
  env <- stack(env)
  Vars = j
  
  # background points
  set.seed(111)
  notna <- which(complete.cases(values(env)))
  samp <- sample(notna, 10000, replace = F)
  samplocs <- as.data.frame(xyFromCell(env, samp))
  
  # SWD object
  data <- prepareSWD(species = Species, p = dat, a = samplocs, env = env)
  
  # run maxent replicates with best parametrization
  folds <- randomFolds(data, k = 10, only_presence = T, seed = 111)
  default_model <- train(method = 'Maxent', data = data, fc = 'lq', reg = 0.3, iter = 1000, folds = folds)
  
  # presences and absences with variable data
  p <- .get_presence(default_model@data)
  p$var <- j
  names(p)[names(p) == j] <- 'values'
  a <- .get_absence(default_model@data)
  a$var <- j
  names(a)[names(a) == j] <- 'values'
  
  pred <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 10))
  for(i in 1:10){
    pred[, i] <- predict(default_model@models[[i]], data = data@data, type = 'logistic')
  }
  
  # plot data
  plot_data <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 4))
  names(plot_data) <- c('mean', 'sd', 'max', 'min')
  plot_data$mean <- rowMeans(pred)
  plot_data$sd <- apply(pred, 1, sd)
  plot_data$max <- plot_data$mean + plot_data$sd
  plot_data$min <- plot_data$mean - plot_data$sd
  plot_data$var <- j
  plot_data$values <- data@data[, j]
  
  # data frames
  p_df <- rbind(p_df, p)
  a_df <- rbind(a_df, a)
  plot_data_df <- rbind(plot_data_df, plot_data)
}

# Re-order
plot_data_df <- transform(plot_data_df, var = factor(var, levels = c('Temperature', 'Distance_to_coast', 
                                                                     'Surface_temperature', 'Kd490')))
a_df <- transform(a_df, var = factor(var, levels = c('Temperature', 'Distance_to_coast', 
                                                     'Surface_temperature', 'Kd490')))
p_df <- transform(p_df, var = factor(var, levels = c('Temperature', 'Distance_to_coast', 
                                                     'Surface_temperature', 'Kd490')))

# Labeller
VAR_names = as_labeller(c(Temperature = 'Temperature~(ºC)', Distance_to_coast = 'Distance~to~coast~(km)',
                          Surface_temperature = 'Sea~Surface~Temperature~(ºC)', Kd490 = 'Kd490~(m^-1)'), 
                          default = label_parsed)

# Plot
ggplot(data = plot_data_df, aes(x = values, y = mean, ymin = min, ymax = max)) + 
  geom_line(colour = '#440154FF') + geom_ribbon(fill = '#440154FF', alpha = 0.2) +
  labs(x = NULL, y = 'Logistic output') + 
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.35),
        panel.grid.major = element_line(size = 0.2, colour = 'grey90'),
        strip.background = element_rect(fill = 'transparent'),
        strip.text = element_text(vjust = -0.5, size = 12),
        panel.spacing.y = unit(-0.3, 'lines'),
        panel.spacing.x = unit(0.7, 'lines'),
        plot.margin = unit(c(-0.2,0.4,0.2,0.2), 'cm'),
        axis.title = element_text(size = 12)) +
  geom_rug(data = p_df, inherit.aes = F, aes(values), sides = 't', color = '#FDE725FF', size = 0.3) + 
  geom_rug(data = a_df, inherit.aes = F, aes(values), sides = 'b', color = '#21908CFF', size = 0.3) + 
  ylim(0, 1) + facet_wrap(~ var, scales = 'free_x', labeller = VAR_names) 

ggsave('Figure S5.pdf', width = 15, height = 12, units = 'cm')


#----------------------------------- Figure S6 ---------------------------------------------

# Response plots for AUS model

# Install and call 'SDMtune' package (https://github.com/ConsBiol-unibern/SDMtune)

Region <- 'AUS'

# Read predictors
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope')
names(env) <- var_names

# Subset selected predictors
env <- env[[c('Temperature', 'Surface_temperature', 'Kd490', 'Distance_to_coast')]]

# Occurrences
dat <- read.csv(paste(Region, '/N_cepedianus_joint.csv', sep = ''), header = T)
dat <- dat[, c('longitude', 'latitude')]

# M space
occ.tot <- SpatialPoints(dat, proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env, occ.buff)
env.M <- mask(env.M, occ.buff)
env.M <- stack(env.M)

p_df = data.frame()
a_df = data.frame()
plot_data_df = data.frame()

for(j in names(env.M)){
  # prepare data
  env = env.M[[j]]
  env <- stack(env)
  Vars = j
  
  # background points
  set.seed(111)
  notna <- which(complete.cases(values(env)))
  samp <- sample(notna, 10000, replace = F)
  samplocs <- as.data.frame(xyFromCell(env, samp))
  
  # SWD object
  data <- prepareSWD(species = Species, p = dat, a = samplocs, env = env)
  
  # run maxent replicates with best parametrization
  folds <- randomFolds(data, k = 10, only_presence = T, seed = 111)
  default_model <- train(method = 'Maxent', data = data, fc = 'lq', reg = 0.3, iter = 1000, folds = folds)
  
  # presences and absences with variable data
  p <- .get_presence(default_model@data)
  p$var <- j
  names(p)[names(p) == j] <- 'values'
  a <- .get_absence(default_model@data)
  a$var <- j
  names(a)[names(a) == j] <- 'values'
  
  pred <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 10))
  for(i in 1:10){
    pred[, i] <- predict(default_model@models[[i]], data = data@data, type = 'logistic')
  }
  
  # plot data
  plot_data <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 4))
  names(plot_data) <- c('mean', 'sd', 'max', 'min')
  plot_data$mean <- rowMeans(pred)
  plot_data$sd <- apply(pred, 1, sd)
  plot_data$max <- plot_data$mean + plot_data$sd
  plot_data$min <- plot_data$mean - plot_data$sd
  plot_data$var <- j
  plot_data$values <- data@data[, j]
  
  # data frames
  p_df <- rbind(p_df, p)
  a_df <- rbind(a_df, a)
  plot_data_df <- rbind(plot_data_df, plot_data)
}

# Re-order
plot_data_df <- transform(plot_data_df, var = factor(var, levels = c('Temperature', 'Distance_to_coast', 
                                                                     'Surface_temperature', 'Kd490')))
a_df <- transform(a_df, var = factor(var, levels = c('Temperature', 'Distance_to_coast', 
                                                     'Surface_temperature', 'Kd490')))
p_df <- transform(p_df, var = factor(var, levels = c('Temperature', 'Distance_to_coast', 
                                                     'Surface_temperature', 'Kd490')))

# Labeller
VAR_names = as_labeller(c(Temperature = 'Temperature~(ºC)', Distance_to_coast = 'Distance~to~coast~(km)',
                          Surface_temperature = 'Sea~Surface~Temperature~(ºC)', Kd490 = 'Kd490~(m^-1)'), 
                        default = label_parsed)

# Plot
ggplot(data = plot_data_df, aes(x = values, y = mean, ymin = min, ymax = max)) + 
  geom_line(colour = '#440154FF') + geom_ribbon(fill = '#440154FF', alpha = 0.2) +
  labs(x = NULL, y = 'Logistic output') + 
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.35),
        panel.grid.major = element_line(size = 0.2, colour = 'grey90'),
        strip.background = element_rect(fill = 'transparent'),
        strip.text = element_text(vjust = -0.5, size = 12),
        panel.spacing.y = unit(-0.3, 'lines'),
        panel.spacing.x = unit(0.7, 'lines'),
        plot.margin = unit(c(-0.2,0.4,0.2,0.2), 'cm'),
        axis.title = element_text(size = 12)) +
  geom_rug(data = p_df, inherit.aes = F, aes(values), sides = 't', color = '#FDE725FF', size = 0.3) + 
  geom_rug(data = a_df, inherit.aes = F, aes(values), sides = 'b', color = '#21908CFF', size = 0.3) + 
  ylim(0, 1) + facet_wrap(~ var, scales = 'free_x', labeller = VAR_names) 

ggsave('Figure S6.pdf', width = 15, height = 12, units = 'cm')


#----------------------------------- Figure 3 ---------------------------------------------

# Niche overlap between SWA and AUS 

Region1 <- 'SWA'
Region2 <- 'AUS'

# Read predictors
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope')
names(env) <- var_names

# Occurrences
dat_swa <- read.csv(paste(Region1, '/N_cepedianus_joint.csv', sep = ''), header = T)
dat_aus <- read.csv(paste(Region2, '/N_cepedianus_joint.csv', sep = ''), header = T)

# M space SWA
coord1 <- cbind(dat_swa$longitude, dat_swa$latitude)
occ.tot1 <- SpatialPoints(coord1, proj4string = crs) 
occ.buff1 <- buffer(occ.tot1, width = 1000000) 
env.M1 <- crop(env, occ.buff1)
env.M1 <- mask(env.M1, occ.buff1)

# Erase Pacific areas
lon1 <- c(-84,-69.4,-69.4,-84,-84) 
lat1 <- c(-32,-32,-62,-62,-32)
pol1 <- cbind(lon1, lat1)
Pol1 <- Polygon(pol1, hole = F)
list1 <- list(Pol1)
poly1 <- Polygons(list1, ID = 'p1')
spp1 <- SpatialPolygons(list(poly1), proj4string = crs)
env.M1 <- mask(env.M1, spp1, inverse = T)
Ext1 <- extent(-69.4,-39,-61.4,-22.8)
env.M1 <- crop(env.M1, Ext1)
env.M_swa <- stack(env.M1)

# M space AUS
coord2 <- cbind(dat_aus$longitude, dat_aus$latitude)
occ.tot2 <- SpatialPoints(coord2, proj4string = crs) 
occ.buff2 <- buffer(occ.tot2, width = 1000000) 
env.M2 <- crop(env, occ.buff2)
env.M2 <- mask(env.M2, occ.buff2)
env.M_aus <- stack(env.M2)

# Ellipsoid plot

# Subset selected predictors minus KD490 not important in models to simplify interpretation
env.M_swa.x <- env.M_swa[[c('Surface_temperature', 'Distance_to_coast', 'Temperature')]]
env.M_aus.x <- env.M_aus[[c('Surface_temperature', 'Distance_to_coast', 'Temperature')]]

# Preparing overlap objects to perform analyses
niche_swa <- overlap_object(data = dat_swa, species = 'species', longitude = 'longitude', latitude = 'latitude', 
                               method = 'mve1', level = 95, variables = env.M_swa.x)
niche_aus <- overlap_object(data = dat_aus, species = 'species', longitude = 'longitude', latitude = 'latitude', 
                               method = 'mve1', level = 95, variables = env.M_aus.x)

# Overlap
N.test <- ellipsoid_overlap(niche_swa, niche_aus, overlap_type = 'back_union', significance_test = F)

# Plots: ellipsoid overlap and statistical significance test
plot_overlap(N.test, niches = c(1, 2), data = T, background = T, proportion = 1, 
             background_type = 'back_union', change_labels = T,
             data_col = c('darkred', 'darkorange'), niche_col = c('darkred', 'darkorange'), 
             background_col = viridis::viridis, legend = F)
rgl.bbox(color = c('#21908CFF', 'black'), xat = c(10, 15, 20), yat = c(0, 100, 200), zat = c(0, 12, 24))
mtext3d('Sea surface temperature (ºC)', edge = 'x++', line = 2.2, at = 16, col = 'black')
mtext3d('Distance to coast (km)', edge = 'y--', line = 0.7, at = -50, col = 'black')
mtext3d('Mean depth temperature (ºC)', edge = 'z-+', line = -0.5, at = 18, col = 'black')
snapshot3d('Ellipsoid.png', fmt = 'png', top = T, width = 1000, height = 1000, webshot = F)
movie3d(spin3d(axis = c(0, 0, 1), rpm = 2), duration = 20, dir = './')

# Test plot

# Subset selected predictors
env.M_swa.x <- env.M_swa[[c('Temperature', 'Distance_to_coast', 'Surface_temperature', 'Kd490')]]
env.M_aus.x <- env.M_aus[[c('Temperature', 'Distance_to_coast', 'Surface_temperature', 'Kd490')]]

# Preparing overlap objects to perform analyses
niche_swa <- overlap_object(data = dat_swa, species = 'species', longitude = 'longitude', latitude = 'latitude', 
                            method = 'mve1', level = 95, variables = env.M_swa.x)
niche_aus <- overlap_object(data = dat_aus, species = 'species', longitude = 'longitude', latitude = 'latitude', 
                            method = 'mve1', level = 95, variables = env.M_aus.x)

# Overlap
N.test <- ellipsoid_overlap(niche_swa, niche_aus, overlap_type = 'back_union', 
                            significance_test = T, replicates = 1000)

ggplot() +
  geom_histogram(aes(N.test@significance_results$union_random$Niche_1_vs_2$overlap), 
                 fill = '#00204DFF', bins = 60) +
  geom_vline(xintercept = quantile(N.test@significance_results$union_random$Niche_1_vs_2$overlap, 0.05),
             colour = '#FFEA46FF', size = 0.75, linetype = 'dashed') +  
  geom_vline(xintercept = N.test@union_overlap$overlap, 
             colour = '#7C7B78FF', size = 0.75, linetype = 'dashed') + 
  theme_bw() + xlab('Overlap') + ylab('Frequency') 
ggsave('Test.pdf', dpi = 900, width = 12, height = 10, units = 'cm')


#----------------------------------- Figure 4 ---------------------------------------------

# Continuos and binary SWA/AUS mosaic global transfers 

# Continuos
Region1 <- 'SWA'
Region2 <- 'AUS'

swa_raw <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
aus_raw <- raster(paste(Region2, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model

# MOP mosaic
swa_mop <- raster(paste(Region1, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
swa_mop[swa_mop > 0] <- NA
swa_mop[swa_mop == 0] <- 1
aus_mop <- raster(paste(Region2, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
aus_mop[aus_mop > 0] <- NA
aus_mop[aus_mop == 0] <- 1
mop_mosaic <- mosaic(swa_mop, aus_mop, fun = sum)
mop_mosaic[mop_mosaic < 1] <- NA
mop_mosaic[mop_mosaic >= 1] <- 2 # areas of strict extrapolation

# SWA/AUS/MOP mosaic
mod_mosaic <- mosaic(swa_raw, aus_raw, fun = median) # final model consensus
mod_mosaic <- mosaic(mod_mosaic, mop_mosaic, fun = max) 
mod_mosaic[mod_mosaic > 1] <- NA # minus areas of strict extrapolation
mod_mosaic <- crop(mod_mosaic, coast)

df <- data.frame(coordinates(mod_mosaic), as.data.frame(mod_mosaic))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = layer)) + 
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), col = 'black', fill = 'grey10', size = 0.1) +
  scale_fill_viridis(na.value = '#440154FF') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(50, 0, -50), labels = c('50ºN', '0º', '50ºS')) + 
  scale_x_continuous(name = NULL, breaks = c(100, 0, -100), labels = c('', '', '')) 
ggsave('Figure 4A.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm')


# Binary
Region1 <- 'SWA'
Region2 <- 'AUS'

swa_raw <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
aus_raw <- raster(paste(Region2, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
bg_binary <- swa_raw
bg_binary[bg_binary >= 0] <- 1 # background raster

# All occurrences used for modelling
swa_occ <- read.csv(paste(Region1, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences
aus_occ <- read.csv(paste(Region2, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences

# 10% Minimum training presence thresholds
swa_mtp10 <- sdm_threshold.10(swa_raw, swa_occ[, c('longitude', 'latitude')], 'p10', binary = F)
swa_bin_10 <- swa_raw
swa_bin_10[swa_bin_10 < minValue(swa_mtp10)] <- NA
swa_bin_10[swa_bin_10 >= minValue(swa_mtp10)] <- 2
aus_mtp10 <- sdm_threshold.10(aus_raw, aus_occ[, c('longitude', 'latitude')], 'p10', binary = F)
aus_bin_10 <- aus_raw
aus_bin_10[aus_bin_10 < minValue(aus_mtp10)] <- NA
aus_bin_10[aus_bin_10 >= minValue(aus_mtp10)] <- 3

# MOP mosaic
swa_mop <- raster(paste(Region1, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
swa_mop[swa_mop > 0] <- NA
swa_mop[swa_mop == 0] <- 8
aus_mop <- raster(paste(Region2, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
aus_mop[aus_mop > 0] <- NA
aus_mop[aus_mop == 0] <- 8
mop_mosaic <- mosaic(swa_mop, aus_mop, fun = sum)
mop_mosaic[mop_mosaic < 8] <- NA
mop_mosaic[mop_mosaic >= 8] <- 8 # areas of strict extrapolation

# SWA mosaic
swa_mosaic <- mosaic(bg_binary, swa_bin_10, mop_mosaic, fun = sum)
swa_mosaic[swa_mosaic > 3] <- 8

# AUS mosaic
aus_mosaic <- mosaic(bg_binary, aus_bin_10, mop_mosaic, fun = sum)
aus_mosaic[aus_mosaic > 4] <- 8

# Final mosaic
mod_mosaic <- mosaic(swa_mosaic, aus_mosaic, fun = sum)
mod_mosaic[mod_mosaic == 2] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 4] <- 2 # swa MTP
mod_mosaic[mod_mosaic == 5] <- 3 # aus MTP
mod_mosaic[mod_mosaic == 7] <- 4 # overlapping MTPs
mod_mosaic[mod_mosaic == 16] <- 5 # strict extrapolation areas
mod_mosaic[mod_mosaic == 9] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 11] <- 5 # strict extrapolation areas
mod_mosaic[mod_mosaic == 12] <- 5 # strict extrapolation areas
mod_mosaic <- crop(mod_mosaic, coast)

# Define bounding box of regions predicted as suitable
swa_ext <- extent(-68.5, -48, -49, -29)        # Southwest Atlantic
aus_ext <- extent(111, 155, -46.5, -24.5)      # Australia
nep_ext <- extent(-135, -107, 21.5, 56)        # Northeast Pacific
nwa_ext <- extent(-78, -56.5, 34, 45)          # Northwest Atlantic
azo_ext <- extent(-34, -22.25, 35, 42)         # Azores Is.
nea_ext <- extent(-20.5, 42.5, 13, 64)         # Northeast Atlantic
nwp_ext <- extent(116, 146.5, 27, 46.5)        # Northwest Pacific
gal_ext <- extent(-93, -87.5, -2.5, 1.25)      # Galápagos Is.
sep_ext <- extent(-84, -69, -52, -0.5)         # Southeast Pacific
tdc_ext <- extent(-15.5, -7, -43, -35)         # Tristan da Cunha Is.
sea_ext <- extent(9, 33, -37.5, -13)           # Southeast Atlantic
ams_ext <- extent(74.5, 80.5, -41.5, -35.25)   # Amsterdam I.
nze1_ext <- extent(163, 179.99, -52, -27)      # New Zealand 1
nze2_ext <- extent(-179.99, -172, -52, -27)    # New Zealand 2

Ids1 <- c('Southwest Atlantic', 'Australia') # local transfers
         
Ids2 <- c('Northeast Pacific', 'Northwest Atlantic', 'Azores Is.', 'Northeast Atlantic', 
         'Northwest Pacific', 'Galápagos Is.', 'Southeast Pacific', 'Tristan da Cunha Is.', 
         'Southeast Atlantic', 'Amsterdam I.', 'New Zealand 1', 'New Zealand 2') # exotic transfers
  
polys1 <- SpatialPolygons(list(
  Polygons(list(Polygon(swa_ext)), Ids1[1]), Polygons(list(Polygon(aus_ext)), Ids1[2])), 
  proj4string = crs) # local transfers

polys2 <- SpatialPolygons(list(
  Polygons(list(Polygon(nep_ext)), Ids2[1]), Polygons(list(Polygon(nwa_ext)), Ids2[2]),
  Polygons(list(Polygon(azo_ext)), Ids2[3]), Polygons(list(Polygon(nea_ext)), Ids2[4]),
  Polygons(list(Polygon(nwp_ext)), Ids2[5]), Polygons(list(Polygon(gal_ext)), Ids2[6]),
  Polygons(list(Polygon(sep_ext)), Ids2[7]), Polygons(list(Polygon(tdc_ext)), Ids2[8]),
  Polygons(list(Polygon(sea_ext)), Ids2[9]), Polygons(list(Polygon(ams_ext)), Ids2[10]),
  Polygons(list(Polygon(nze1_ext)), Ids2[11]), Polygons(list(Polygon(nze2_ext)), Ids2[12])), 
  proj4string = crs) # exotic transfers

polys1 <- SpatialPolygonsDataFrame(polys1, data.frame(ids = Ids1, row.names = names(polys1))) 
polys2 <- SpatialPolygonsDataFrame(polys2, data.frame(ids = Ids2, row.names = names(polys2))) 

# Set colours for each code
Col <- c('grey95', 'grey80', '#FDE725FF', '#21908CFF', '#440154FF', 'grey90')

df <- data.frame(coordinates(mod_mosaic), as.data.frame(mod_mosaic))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  geom_polygon(data = polys1, aes(x = long, y = lat, group = group), col = '#3333CC', fill = NA, size = 0.45) +
  geom_polygon(data = polys2, aes(x = long, y = lat, group = group), col = 'grey40', fill = NA, size = 0.45) +
  coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(50, 0, -50), labels = c('50ºN', '0º', '50ºS')) + 
  scale_x_continuous(name = NULL, breaks = c(100, 0, -100), labels = c('100ºE', '0º', '100ºW')) 
ggsave('Figure 4B.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm')


#----------------------------------- Figure S7 ---------------------------------------------

# Maximum and minimum variability in continuos SWA/AUS mosaic global transfers 

Region1 <- 'SWA'
Region2 <- 'AUS'

swa_raw <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
aus_raw <- raster(paste(Region2, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model

# MOP mosaic
swa_mop <- raster(paste(Region1, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
swa_mop[swa_mop > 0] <- NA
swa_mop[swa_mop == 0] <- 1
aus_mop <- raster(paste(Region2, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
aus_mop[aus_mop > 0] <- NA
aus_mop[aus_mop == 0] <- 1
mop_mosaic <- mosaic(swa_mop, aus_mop, fun = sum)
mop_mosaic[mop_mosaic < 1] <- NA
mop_mosaic[mop_mosaic >= 1] <- 2 # areas of strict extrapolation

# Maximum

# SWA/AUS/MOP mosaic
mod_mosaic <- mosaic(swa_raw, aus_raw, fun = max) # final model consensus
mod_mosaic <- mosaic(mod_mosaic, mop_mosaic, fun = max) 
mod_mosaic[mod_mosaic > 1] <- NA # minus areas of strict extrapolation
mod_mosaic <- crop(mod_mosaic, coast)

df <- data.frame(coordinates(mod_mosaic), as.data.frame(mod_mosaic))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = layer)) + 
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), col = 'black', fill = 'grey10', size = 0.1) +
  scale_fill_viridis(na.value = '#440154FF') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(50, 0, -50), labels = c('50ºN', '0º', '50ºS')) + 
  scale_x_continuous(name = NULL, breaks = c(100, 0, -100), labels = c('', '', '')) 
ggsave('Figure S7A.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm')

# Minimum

# SWA/AUS/MOP mosaic
mod_mosaic <- mosaic(swa_raw, aus_raw, fun = min) # final model consensus
mod_mosaic <- mosaic(mod_mosaic, mop_mosaic, fun = max) 
mod_mosaic[mod_mosaic > 1] <- NA # minus areas of strict extrapolation
mod_mosaic <- crop(mod_mosaic, coast)

df <- data.frame(coordinates(mod_mosaic), as.data.frame(mod_mosaic))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = layer)) + 
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), col = 'black', fill = 'grey10', size = 0.1) +
  scale_fill_viridis(na.value = '#440154FF') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(50, 0, -50), labels = c('50ºN', '0º', '50ºS')) + 
  scale_x_continuous(name = NULL, breaks = c(100, 0, -100), labels = c('100ºE', '0º', '100ºW')) 
ggsave('Figure S7B.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm')


#----------------------------------- Figure 5 ---------------------------------------------

# Binary transfers per continental region with occurrence data and SWA/AUS overlap

Region1 <- 'SWA'
Region2 <- 'AUS'

swa_raw <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
aus_raw <- raster(paste(Region2, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
bg_binary <- swa_raw
bg_binary[bg_binary >= 0] <- 1 # background raster

# All occurrences used for modelling
swa_occ <- read.csv(paste(Region1, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences
aus_occ <- read.csv(paste(Region2, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences

# 10% Minimum training presence thresholds
swa_mtp10 <- sdm_threshold.10(swa_raw, swa_occ[, c('longitude', 'latitude')], 'p10', binary = F)
swa_bin_10 <- swa_raw
swa_bin_10[swa_bin_10 < minValue(swa_mtp10)] <- NA
swa_bin_10[swa_bin_10 >= minValue(swa_mtp10)] <- 2
aus_mtp10 <- sdm_threshold.10(aus_raw, aus_occ[, c('longitude', 'latitude')], 'p10', binary = F)
aus_bin_10 <- aus_raw
aus_bin_10[aus_bin_10 < minValue(aus_mtp10)] <- NA
aus_bin_10[aus_bin_10 >= minValue(aus_mtp10)] <- 3

# MOP mosaic
swa_mop <- raster(paste(Region1, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
swa_mop[swa_mop > 0] <- NA
swa_mop[swa_mop == 0] <- 8
aus_mop <- raster(paste(Region2, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
aus_mop[aus_mop > 0] <- NA
aus_mop[aus_mop == 0] <- 8
mop_mosaic <- mosaic(swa_mop, aus_mop, fun = sum)
mop_mosaic[mop_mosaic < 8] <- NA
mop_mosaic[mop_mosaic >= 8] <- 8 # areas of strict extrapolation

# SWA mosaic
swa_mosaic <- mosaic(bg_binary, swa_bin_10, mop_mosaic, fun = sum)
swa_mosaic[swa_mosaic > 3] <- 8

# AUS mosaic
aus_mosaic <- mosaic(bg_binary, aus_bin_10, mop_mosaic, fun = sum)
aus_mosaic[aus_mosaic > 4] <- 8

# Final mosaic
mod_mosaic <- mosaic(swa_mosaic, aus_mosaic, fun = sum)
mod_mosaic[mod_mosaic == 2] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 4] <- 2 # swa MTP
mod_mosaic[mod_mosaic == 5] <- 3 # aus MTP
mod_mosaic[mod_mosaic == 7] <- 4 # overlapping MTPs
mod_mosaic[mod_mosaic == 16] <- 5 # strict extrapolation areas
mod_mosaic[mod_mosaic == 9] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 11] <- 5 # strict extrapolation areas
mod_mosaic[mod_mosaic == 12] <- 5 # strict extrapolation areas
mod_mosaic <- crop(mod_mosaic, coast)

# Define bounding box of regions predicted as suitable
swa_ext <- extent(-68.5, -48, -49, -29)        # Southwest Atlantic
aus_ext <- extent(111, 155, -46.5, -24.5)      # Australia
nep_ext <- extent(-135, -107, 21.5, 56)        # Northeast Pacific
nwp_ext <- extent(116, 146.5, 27, 46.5)        # Northwest Pacific
sep_ext <- extent(-84, -69, -52, -0.5)         # Southeast Pacific
sea_ext <- extent(9, 33, -37.5, -13)           # Southeast Atlantic
nze1_ext <- extent(163, 179.99, -52, -27)      # New Zealand 1
nze2_ext <- extent(-179.99, -172, -52, -27)    # New Zealand 2

# Read occurrence data
dat <- read.csv('Appendix C.csv', header = T)

## Southwest Atlantic ##

swa_occ <- subset(dat, Region == 'Southwest Atlantic') # subset occurrences
swa_occ1 <- subset(swa_occ, Source1 == 'Published literature')
swa_occ2 <- subset(swa_occ, Source1 == 'GBIF/OBIS')
swa_occ3 <- subset(swa_occ, Source1 %in% c('Social media', 'Grey literature'))

# Check if extent of occurrences is bigger than predictions, if so re-define extent to the broadest area
occ_ext <- extent(SpatialPoints(swa_occ[, c('Longitude', 'Latitude')], crs))
print('occurrences extent'); occ_ext; print('swa extent'); swa_ext
swa_ext <- extent(-69.5, -48, -53, -28.5)
swa_final <- crop(mod_mosaic, swa_ext) # crop to region

swa_df <- data.frame(coordinates(swa_final), as.data.frame(swa_final))

# Set colours for each code
unique(swa_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#FDE725FF', '#21908CFF', '#440154FF', 'grey90')

coast_swa <- crop(coast, swa_final) # crop for faster plotting

ggplot() +
  geom_tile(data = swa_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_swa, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  geom_point(data = swa_occ2, aes(x = Longitude, y = Latitude), shape = 17, size = 1, color = 'darkorange') +
  geom_point(data = swa_occ3, aes(x = Longitude, y = Latitude), shape = 17, size = 1, color = 'red') +
  geom_point(data = swa_occ1, aes(x = Longitude, y = Latitude), size = 2, color = 'white') +
  geom_point(data = swa_occ1, aes(x = Longitude, y = Latitude), size = 1, color = 'blue') +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 6),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(-30, -40, -50), labels = c('-30º', '-40º', '-50º')) + 
  scale_x_continuous(name = NULL, breaks = c(-55, -65), labels = c('-55º', '-65º')) +
  ggsn::scalebar(x.min = min(swa_df$x), x.max = max(swa_df$x), y.min = min(swa_df$y), y.max = max(swa_df$y), transform = T, 
                 dist = 100, st.size = 1.1, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -49.5, y = -51.9), box.fill = 'black')
ggsave('Figure 5_swa.tiff', dpi = 900, width = 5.55, height = 5.57, units = 'cm')


## Australia ##

aus_occ <- subset(dat, Region == 'Australia') # subset occurrences
aus_occ1 <- subset(aus_occ, Source1 == 'Published literature')
aus_occ2 <- subset(aus_occ, Source1 == 'GBIF/OBIS')
aus_occ3 <- subset(aus_occ, Source1 %in% c('Social media', 'Grey literature'))

# Check if extent of occurrences is bigger than predictions, if so re-define extent to the broadest area
occ_ext <- extent(SpatialPoints(aus_occ[, c('Longitude', 'Latitude')], crs))
print('occurrences extent'); occ_ext; print('aus extent'); aus_ext
aus_ext <- aus_ext
aus_final <- crop(mod_mosaic, aus_ext) # crop to region

aus_df <- data.frame(coordinates(aus_final), as.data.frame(aus_final))

# Set colours for each code
unique(aus_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#FDE725FF', '#21908CFF', '#440154FF', 'grey90')

coast_aus <- crop(coast, aus_final) # crop for faster plotting

ggplot() +
  geom_tile(data = aus_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_aus, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  geom_point(data = aus_occ2, aes(x = Longitude, y = Latitude), shape = 17, size = 1.25, color = 'darkorange') +
  geom_point(data = aus_occ3, aes(x = Longitude, y = Latitude), shape = 17, size = 1.25, color = 'red') +
  geom_point(data = aus_occ1, aes(x = Longitude, y = Latitude), size = 2.5, color = 'white') +
  geom_point(data = aus_occ1, aes(x = Longitude, y = Latitude), size = 1.25, color = 'blue') +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 5.5),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(-25, -35, -45), labels = c('-25º', '-35º', '-45º')) + 
  scale_x_continuous(name = NULL, breaks = c(120, 130, 140, 150), labels = c('120º', '130º', '140º', '150º')) +
  ggsn::scalebar(x.min = min(aus_df$x), x.max = max(aus_df$x), y.min = min(aus_df$y), y.max = max(aus_df$y), transform = T, 
                 dist = 200, st.size = 1.4, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = 117, y = -45), box.fill = 'black')
ggsave('Figure 5_aus.tiff', dpi = 900, width = 7.42, height = 3.82, units = 'cm')


## Northeast Pacific ##

nep_occ <- subset(dat, Region == 'Northeast Pacific') # subset occurrences
nep_occ1 <- subset(nep_occ, Source1 == 'Published literature')
nep_occ2 <- subset(nep_occ, Source1 == 'GBIF/OBIS')
nep_occ3 <- subset(nep_occ, Source1 %in% c('Social media', 'Grey literature'))

# Check if extent of occurrences is bigger than predictions, if so re-define extent to the broadest area
occ_ext <- extent(SpatialPoints(nep_occ[, c('Longitude', 'Latitude')], crs))
print('occurrences extent'); occ_ext; print('nep extent'); nep_ext
nep_ext <- nep_ext
nep_final <- crop(mod_mosaic, nep_ext) # crop to region

nep_df <- data.frame(coordinates(nep_final), as.data.frame(nep_final))

# Set colours for each code
unique(nep_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#FDE725FF', '#21908CFF', '#440154FF', 'grey90')

coast_nep <- crop(coast, nep_final) # crop for faster plotting

ggplot() +
  geom_tile(data = nep_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_nep, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  geom_point(data = nep_occ2, aes(x = Longitude, y = Latitude), shape = 17, size = 1.15, color = 'darkorange') +
  geom_point(data = nep_occ3, aes(x = Longitude, y = Latitude), shape = 17, size = 1.15, color = 'red') +
  geom_point(data = nep_occ1, aes(x = Longitude, y = Latitude), size = 2.3, color = 'white') +
  geom_point(data = nep_occ1, aes(x = Longitude, y = Latitude), size = 1.15, color = 'blue') +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 5.5),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(30, 40, 50), labels = c('30º', '40º', '50º')) + 
  scale_x_continuous(name = NULL, breaks = c(-110, -120, -130), labels = c('-110º', '-120º', '-130º')) +
  ggsn::scalebar(x.min = min(nep_df$x), x.max = max(nep_df$x), y.min = min(nep_df$y), y.max = max(nep_df$y), transform = T, 
                 dist = 200, st.size = 1.1, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -130.4, y = 22.9), box.fill = 'black')
ggsave('Figure 5_nep.tiff', dpi = 900, width = 4.59, height = 5.57, units = 'cm')


## Northwest Pacific ##

nwp_occ <- subset(dat, Region == 'Northwest Pacific') # subset occurrences
nwp_occ1 <- subset(nwp_occ, Source1 == 'Published literature')
nwp_occ2 <- subset(nwp_occ, Source1 == 'GBIF/OBIS')
nwp_occ3 <- subset(nwp_occ, Source1 %in% c('Social media', 'Grey literature'))

# Check if extent of occurrences is bigger than predictions, if so re-define extent to the broadest area
occ_ext <- extent(SpatialPoints(nwp_occ[, c('Longitude', 'Latitude')], crs))
print('occurrences extent'); occ_ext; print('nwp extent'); nwp_ext
nwp_ext <- nwp_ext
nwp_final <- crop(mod_mosaic, nwp_ext) # crop to region

nwp_df <- data.frame(coordinates(nwp_final), as.data.frame(nwp_final))

# Set colours for each code
unique(nwp_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#FDE725FF', '#21908CFF', '#440154FF', 'grey90')

coast_nwp <- crop(coast, nwp_final) # crop for faster plotting

ggplot() +
  geom_tile(data = nwp_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_nwp, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  geom_point(data = nwp_occ2, aes(x = Longitude, y = Latitude), shape = 17, size = 1.25, color = 'darkorange') +
  geom_point(data = nwp_occ3, aes(x = Longitude, y = Latitude), shape = 17, size = 1.25, color = 'red') +
  geom_point(data = nwp_occ1, aes(x = Longitude, y = Latitude), size = 2.5, color = 'white') +
  geom_point(data = nwp_occ1, aes(x = Longitude, y = Latitude), size = 1.25, color = 'blue') +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 6),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(30, 40), labels = c('30º', '40º')) + 
  scale_x_continuous(name = NULL, breaks = c(120, 130, 140), labels = c('120º', '130º', '140º')) +
  ggsn::scalebar(x.min = min(nwp_df$x), x.max = max(nwp_df$x), y.min = min(nwp_df$y), y.max = max(nwp_df$y), transform = T, 
                 dist = 200, st.size = 1.2, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = 144.5, y = 28.2), box.fill = 'black')
ggsave('Figure 5_nwp.tiff', dpi = 900, width = 6.5, height = 3.7, units = 'cm')


## Southeast Pacific ##

sep_occ <- subset(dat, Region == 'Southeast Pacific') # subset occurrences
sep_occ1 <- subset(sep_occ, Source1 == 'Published literature')
sep_occ2 <- subset(sep_occ, Source1 == 'GBIF/OBIS')
sep_occ3 <- subset(sep_occ, Source1 %in% c('Social media', 'Grey literature'))

# Check if extent of occurrences is bigger than predictions, if so re-define extent to the broadest area
occ_ext <- extent(SpatialPoints(sep_occ[, c('Longitude', 'Latitude')], crs))
print('occurrences extent'); occ_ext; print('sep extent'); sep_ext
sep_ext <- sep_ext
sep_final <- crop(mod_mosaic, sep_ext) # crop to region

sep_df <- data.frame(coordinates(sep_final), as.data.frame(sep_final))

# Set colours for each code
unique(sep_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#FDE725FF', '#21908CFF', '#440154FF', 'grey90')

coast_sep <- crop(coast, sep_final) # crop for faster plotting

ggplot() +
  geom_tile(data = sep_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_sep, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  geom_point(data = sep_occ2, aes(x = Longitude, y = Latitude), shape = 17, size = 1.25, color = 'darkorange') +
  geom_point(data = sep_occ3, aes(x = Longitude, y = Latitude), shape = 17, size = 1.25, color = 'red') +
  geom_point(data = sep_occ1, aes(x = Longitude, y = Latitude), size = 2.5, color = 'white') +
  geom_point(data = sep_occ1, aes(x = Longitude, y = Latitude), size = 1.25, color = 'blue') +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 6),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(-45, -25, -5), labels = c('-45º', '-25º', '-5º')) + 
  scale_x_continuous(name = NULL, breaks = c(-80, -70), labels = c('-80º', '-70º')) +
  ggsn::scalebar(x.min = min(sep_df$x), x.max = max(sep_df$x), y.min = min(sep_df$y), y.max = max(sep_df$y), transform = T, 
                 dist = 150, st.size = 1.5, height = 0.007, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.3, anchor = c(x = -78, y = -50), box.fill = 'black')
ggsave('Figure 5_sep.tiff', dpi = 900, width = 3.9, height = 10.2, units = 'cm')


## Southeast Atlantic ##

sea_occ <- subset(dat, Region == 'Southeast Atlantic') # subset occurrences
sea_occ1 <- subset(sea_occ, Source1 == 'Published literature')
sea_occ2 <- subset(sea_occ, Source1 == 'GBIF/OBIS')

# Check if extent of occurrences is bigger than predictions, if so re-define extent to the broadest area
occ_ext <- extent(SpatialPoints(sea_occ[, c('Longitude', 'Latitude')], crs))
print('occurrences extent'); occ_ext; print('sea extent'); sea_ext
sea_ext <- sea_ext
sea_final <- crop(mod_mosaic, sea_ext) # crop to region

sea_df <- data.frame(coordinates(sea_final), as.data.frame(sea_final))

# Set colours for each code
unique(sea_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#FDE725FF', '#21908CFF', '#440154FF', 'grey90')

coast_sea <- crop(coast, sea_final) # crop for faster plotting

ggplot() +
  geom_tile(data = sea_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_sea, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  geom_point(data = sea_occ2, aes(x = Longitude, y = Latitude), shape = 17, size = 1.25, color = 'darkorange') +
  geom_point(data = sea_occ1, aes(x = Longitude, y = Latitude), size = 2.5, color = 'white') +
  geom_point(data = sea_occ1, aes(x = Longitude, y = Latitude), size = 1.25, color = 'blue') +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 6),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(-35, -25, -15), labels = c('-35º', '-25º', '-15º')) + 
  scale_x_continuous(name = NULL, breaks = c(10, 20, 30), labels = c('10º', '20º', '30º')) +
  ggsn::scalebar(x.min = min(sea_df$x), x.max = max(sea_df$x), y.min = min(sea_df$y), y.max = max(sea_df$y), transform = T, 
                 dist = 150, st.size = 1.2, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = 13, y = -36.1), box.fill = 'black')
ggsave('Figure 5_sea.tiff', dpi = 900, width = 5.5, height = 5, units = 'cm')


## New Zealand ## - problematic because is located across the 180º meridian

nze_occ <- subset(dat, Region == 'New Zealand') # subset occurrences
nze_occ$Longitude2 <- ifelse(nze_occ$Longitude < 0, (180 + nze_occ$Longitude) + 180, nze_occ$Longitude) # accomodate longitude to 0º - 360º
nze_occ1 <- subset(nze_occ, Source1 == 'Published literature')
nze_occ2 <- subset(nze_occ, Source1 == 'GBIF/OBIS')

# Accomodate raster to 0º - 360º
nze_1 <- crop(mod_mosaic, extent(-180, 0, -68.91667, 83.66667)) # mind pole areas were removed
nze_2 <- crop(mod_mosaic, extent(0, 180, -68.91667, 83.66667))   
extent(nze_1) <- c(180, 360, -68.91667, 83.66667)
mod_mosaic2 <- merge(nze_1, nze_2)
nze1_ext; nze2_ext # create new extent based on previous extents
nze_ext <- extent(163, 188, -52, -27)
nze_final <- crop(mod_mosaic2, nze_ext)
nze_coast <- crop(coast, nze_ext)

# Check if extent of occurrences is bigger than predictions, if so re-define extent to the broadest area
occ_ext <- extent(SpatialPoints(nze_occ[, c('Longitude2', 'Latitude')], crs))
print('occurrences extent'); occ_ext; print('nze extent'); nze_ext
nze_ext <- nze_ext
nze_final <- crop(mod_mosaic2, nze_ext) # crop to region

nze_df <- data.frame(coordinates(nze_final), as.data.frame(nze_final))

# Set colours for each code
unique(nze_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#FDE725FF', '#21908CFF', '#440154FF', 'grey90')

coast_nze <- crop(coast, nze_final) # crop for faster plotting

ggplot() +
  geom_tile(data = nze_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_nze, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  geom_point(data = nze_occ2, aes(x = Longitude2, y = Latitude), shape = 17, size = 1.25, color = 'darkorange') +
  geom_point(data = nze_occ1, aes(x = Longitude2, y = Latitude), size = 2.5, color = 'white') +
  geom_point(data = nze_occ1, aes(x = Longitude2, y = Latitude), size = 1.25, color = 'blue') +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 6),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(-50, -40, -30), labels = c('-50º', '-40º', '-30º')) + 
  scale_x_continuous(name = NULL, breaks = c(170, 180), labels = c('170º', '180º')) +
  ggsn::scalebar(x.min = min(nze_df$x), x.max = max(nze_df$x), y.min = min(nze_df$y), y.max = max(nze_df$y), transform = T, 
                 dist = 100, st.size = 0.9, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = 186, y = -50.5), box.fill = 'black')
ggsave('Figure 5_nze.tiff', dpi = 900, width = 4.6, height = 4.6, units = 'cm')


#----------------------------------- Figure S8 ---------------------------------------------

# Binary transfers per continental region without occurrence data and SWA/AUS overlap

Region1 <- 'SWA'
Region2 <- 'AUS'

swa_raw <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
aus_raw <- raster(paste(Region2, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
bg_binary <- swa_raw
bg_binary[bg_binary >= 0] <- 1 # background raster

# All occurrences used for modelling
swa_occ <- read.csv(paste(Region1, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences
aus_occ <- read.csv(paste(Region2, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences

# 10% Minimum training presence thresholds
swa_mtp10 <- sdm_threshold.10(swa_raw, swa_occ[, c('longitude', 'latitude')], 'p10', binary = F)
swa_bin_10 <- swa_raw
swa_bin_10[swa_bin_10 < minValue(swa_mtp10)] <- NA
swa_bin_10[swa_bin_10 >= minValue(swa_mtp10)] <- 2
aus_mtp10 <- sdm_threshold.10(aus_raw, aus_occ[, c('longitude', 'latitude')], 'p10', binary = F)
aus_bin_10 <- aus_raw
aus_bin_10[aus_bin_10 < minValue(aus_mtp10)] <- NA
aus_bin_10[aus_bin_10 >= minValue(aus_mtp10)] <- 3

# MOP mosaic
swa_mop <- raster(paste(Region1, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
swa_mop[swa_mop > 0] <- NA
swa_mop[swa_mop == 0] <- 8
aus_mop <- raster(paste(Region2, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
aus_mop[aus_mop > 0] <- NA
aus_mop[aus_mop == 0] <- 8
mop_mosaic <- mosaic(swa_mop, aus_mop, fun = sum)
mop_mosaic[mop_mosaic < 8] <- NA
mop_mosaic[mop_mosaic >= 8] <- 8 # areas of strict extrapolation

# SWA mosaic
swa_mosaic <- mosaic(bg_binary, swa_bin_10, mop_mosaic, fun = sum)
swa_mosaic[swa_mosaic > 3] <- 8

# AUS mosaic
aus_mosaic <- mosaic(bg_binary, aus_bin_10, mop_mosaic, fun = sum)
aus_mosaic[aus_mosaic > 4] <- 8

# Final mosaic
mod_mosaic <- mosaic(swa_mosaic, aus_mosaic, fun = sum)
mod_mosaic[mod_mosaic == 2] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 4] <- 2 # swa MTP
mod_mosaic[mod_mosaic == 5] <- 3 # aus MTP
mod_mosaic[mod_mosaic == 7] <- 4 # overlapping MTPs
mod_mosaic[mod_mosaic == 16] <- 5 # strict extrapolation areas
mod_mosaic[mod_mosaic == 9] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 11] <- 5 # strict extrapolation areas
mod_mosaic[mod_mosaic == 12] <- 5 # strict extrapolation areas
mod_mosaic <- crop(mod_mosaic, coast)

# Define bounding box of regions predicted as suitable

nwa_ext <- extent(-78, -56.5, 34, 45)          # Northwest Atlantic
nea_ext <- extent(-20.5, 42.5, 13, 64)         # Northeast Atlantic


## Northwest Atlantic ##

# no occurrences #

nwa_final <- crop(mod_mosaic, nwa_ext) # crop to region

nwa_df <- data.frame(coordinates(nwa_final), as.data.frame(nwa_final))

# Set colours for each code
unique(nwa_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#FDE725FF', '#21908CFF', '#440154FF', 'grey90')

coast_nwa <- crop(coast, nwa_final) # crop for faster plotting

ggplot() +
  geom_tile(data = nwa_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_nwa, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44), labels = c('36º', '40º', '44º')) + 
  scale_x_continuous(name = NULL, breaks = c(-75, -70, -65, -60), labels = c('-75º', '-70º', '-65º', '-60º')) +
  ggsn::scalebar(x.min = min(nwa_df$x), x.max = max(nwa_df$x), y.min = min(nwa_df$y), y.max = max(nwa_df$y), transform = T, 
                 dist = 150, st.size = 1.75, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -57.5, y = 34.7), box.fill = 'black')
ggsave('Figure S8_nwa.tiff', dpi = 900, width = 9.5, height = 5.4, units = 'cm')


## Northeast Atlantic ##

# no occurrences #

nea_final <- crop(mod_mosaic, nea_ext) # crop to region

nea_df <- data.frame(coordinates(nea_final), as.data.frame(nea_final))

# Set colours for each code
unique(nea_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#FDE725FF', '#21908CFF', '#440154FF', 'grey90')

coast_nea <- crop(coast, nea_final) # crop for faster plotting

ggplot() +
  geom_tile(data = nea_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_nea, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(20, 40, 60), labels = c('20º', '40º', '60º')) + 
  scale_x_continuous(name = NULL, breaks = c(-10, 10, 30), labels = c('-10º', '10º', '30º')) +
  ggsn::scalebar(x.min = min(nea_df$x), x.max = max(nea_df$x), y.min = min(nea_df$y), y.max = max(nea_df$y), transform = T, 
                 dist = 400, st.size = 1.75, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -10, y = 15.2), box.fill = 'black')
ggsave('Figure S8_nea.tiff', dpi = 900, width = 9.6, height = 7.9, units = 'cm')


#----------------------------------- Figure S9 ---------------------------------------------

# Binary transfers per offshore islands with and without occurrence data and SWA/AUS overlap

Region1 <- 'SWA'
Region2 <- 'AUS'

swa_raw <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
aus_raw <- raster(paste(Region2, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
bg_binary <- swa_raw
bg_binary[bg_binary >= 0] <- 1 # background raster

# All occurrences used for modelling
swa_occ <- read.csv(paste(Region1, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences
aus_occ <- read.csv(paste(Region2, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences

# 10% Minimum training presence thresholds
swa_mtp10 <- sdm_threshold.10(swa_raw, swa_occ[, c('longitude', 'latitude')], 'p10', binary = F)
swa_bin_10 <- swa_raw
swa_bin_10[swa_bin_10 < minValue(swa_mtp10)] <- NA
swa_bin_10[swa_bin_10 >= minValue(swa_mtp10)] <- 2
aus_mtp10 <- sdm_threshold.10(aus_raw, aus_occ[, c('longitude', 'latitude')], 'p10', binary = F)
aus_bin_10 <- aus_raw
aus_bin_10[aus_bin_10 < minValue(aus_mtp10)] <- NA
aus_bin_10[aus_bin_10 >= minValue(aus_mtp10)] <- 3

# MOP mosaic
swa_mop <- raster(paste(Region1, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
swa_mop[swa_mop > 0] <- NA
swa_mop[swa_mop == 0] <- 8
aus_mop <- raster(paste(Region2, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
aus_mop[aus_mop > 0] <- NA
aus_mop[aus_mop == 0] <- 8
mop_mosaic <- mosaic(swa_mop, aus_mop, fun = sum)
mop_mosaic[mop_mosaic < 8] <- NA
mop_mosaic[mop_mosaic >= 8] <- 8 # areas of strict extrapolation

# SWA mosaic
swa_mosaic <- mosaic(bg_binary, swa_bin_10, mop_mosaic, fun = sum)
swa_mosaic[swa_mosaic > 3] <- 8

# AUS mosaic
aus_mosaic <- mosaic(bg_binary, aus_bin_10, mop_mosaic, fun = sum)
aus_mosaic[aus_mosaic > 4] <- 8

# Final mosaic
mod_mosaic <- mosaic(swa_mosaic, aus_mosaic, fun = sum)
mod_mosaic[mod_mosaic == 2] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 4] <- 2 # swa MTP
mod_mosaic[mod_mosaic == 5] <- 3 # aus MTP
mod_mosaic[mod_mosaic == 7] <- 4 # overlapping MTPs
mod_mosaic[mod_mosaic == 16] <- 5 # strict extrapolation areas
mod_mosaic[mod_mosaic == 9] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 11] <- 5 # strict extrapolation areas
mod_mosaic[mod_mosaic == 12] <- 5 # strict extrapolation areas
mod_mosaic <- crop(mod_mosaic, coast)

# Define bounding box of regions predicted as suitable

gal_ext <- extent(-93, -87.5, -2.5, 1.25)      # Galápagos Is.
azo_ext <- extent(-34, -22.25, 35, 42)         # Azores Is.
tdc_ext <- extent(-15.5, -7, -43, -35)         # Tristan da Cunha Is.
ams_ext <- extent(74.5, 80.5, -41.5, -35.25)   # Amsterdam/Saint-Paul Is.

# Read occurrence data
dat <- read.csv('Appendix C.csv', header = T)


## Azores Is. ## - has no occurrences

azo_final <- crop(mod_mosaic, azo_ext) # crop to region

azo_df <- data.frame(coordinates(azo_final), as.data.frame(azo_final))

# Set colours for each code
unique(azo_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#FDE725FF', '#21908CFF', '#440154FF', 'grey90')

coast_azo <- crop(coast, azo_final) # crop for faster plotting

ggplot() +
  geom_tile(data = azo_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_azo, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(36, 38, 40), labels = c('36º', '38º', '40º')) + 
  scale_x_continuous(name = NULL, breaks = c(-30, -25), labels = c('-30º', '-25º')) +
  ggsn::scalebar(x.min = min(azo_df$x), x.max = max(azo_df$x), y.min = min(azo_df$y), y.max = max(azo_df$y), transform = T, 
                 dist = 50, st.size = 1.5, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -32.5, y = 35.5), box.fill = 'black')
ggsave('Figure S9_azo.tiff', dpi = 900, width = 8.2, height = 5, units = 'cm')


## Galápagos Is. ##

gal_occ <- subset(dat, Region == 'Galápagos Is.') # subset occurrences
gal_occ1 <- subset(gal_occ, Source1 == 'Published literature')

# Check if extent of occurrences is bigger than predictions, if so re-define extent to the broadest area
occ_ext <- extent(SpatialPoints(gal_occ[, c('Longitude', 'Latitude')], crs))
print('occurrences extent'); occ_ext; print('gal extent'); gal_ext
gal_ext <- gal_ext
gal_final <- crop(mod_mosaic, gal_ext) # crop to region

gal_df <- data.frame(coordinates(gal_final), as.data.frame(gal_final))

# Set colours for each code
unique(gal_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#21908CFF', 'grey90')

coast_gal <- crop(coast, gal_final) # crop for faster plotting

ggplot() +
  geom_tile(data = gal_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_gal, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  geom_point(data = gal_occ1, aes(x = Longitude, y = Latitude), size = 2.5, color = 'white') +
  geom_point(data = gal_occ1, aes(x = Longitude, y = Latitude), size = 1.25, color = 'blue') +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 7),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(-2, 0), labels = c('-2º', '0º')) + 
  scale_x_continuous(name = NULL, breaks = c(-92, -90, -88), labels = c('-92º', '-90º', '-88º')) +
  ggsn::scalebar(x.min = min(gal_df$x), x.max = max(gal_df$x), y.min = min(gal_df$y), y.max = max(gal_df$y), transform = T, 
                 dist = 50, st.size = 1.75, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -88, y = -2.2), box.fill = 'black')
ggsave('Figure S9_gal.tiff', dpi = 900, width = 7.2, height = 5, units = 'cm')


## Tristan da Cunha Is. ##

tdc_occ <- subset(dat, Region == 'Tristan da Cunha Is.') # subset occurrences
tdc_occ1 <- subset(tdc_occ, Source1 == 'Published literature')
tdc_occ2 <- subset(tdc_occ, Source1 == 'GBIF/OBIS')

# Check if extent of occurrences is bigger than predictions, if so re-define extent to the broadest area
occ_ext <- extent(SpatialPoints(tdc_occ[, c('Longitude', 'Latitude')], crs))
print('occurrences extent'); occ_ext; print('tdc extent'); tdc_ext
tdc_ext <- tdc_ext
tdc_final <- crop(mod_mosaic, tdc_ext) # crop to region

tdc_df <- data.frame(coordinates(tdc_final), as.data.frame(tdc_final))

# Set colours for each code
unique(tdc_final) # match numbers and colours to be used
Col <- c('grey95', 'grey80', '#FDE725FF', '#440154FF')

coast_tdc <- crop(coast, tdc_final) # crop for faster plotting

ggplot() +
  geom_tile(data = tdc_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_tdc, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  geom_point(data = tdc_occ2, aes(x = Longitude, y = Latitude), shape = 17, size = 1.25, color = 'darkorange') +
  geom_point(data = tdc_occ1, aes(x = Longitude, y = Latitude), size = 2.5, color = 'white') +
  geom_point(data = tdc_occ1, aes(x = Longitude, y = Latitude), size = 1.25, color = 'blue') +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(-40, -36), labels = c('-40º', '-36º')) + 
  scale_x_continuous(name = NULL, breaks = c(-14, -8), labels = c('-14º', '-8º')) +
  ggsn::scalebar(x.min = min(tdc_df$x), x.max = max(tdc_df$x), y.min = min(tdc_df$y), y.max = max(tdc_df$y), transform = T, 
                 dist = 50, st.size = 1.75, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -14, y = -42.5), box.fill = 'black')
ggsave('Figure S9_tdc.tiff', dpi = 900, width = 7.8, height = 6.8, units = 'cm')


## Amsterdam/Saint-Paul Is. ## - has no occurrences

ams_final <- crop(mod_mosaic, ams_ext) # crop to region

ams_df <- data.frame(coordinates(ams_final), as.data.frame(ams_final))

# Set colours for each code
unique(ams_final) # match numbers and colours to be used
Col <- c('grey80', '#FDE725FF', '#440154FF', 'grey90')

coast_ams <- crop(coast, ams_final) # crop for faster plotting

ggplot() +
  geom_tile(data = ams_df, aes(x = x, y = y, fill = as.factor(layer))) + scale_fill_manual(values = Col) +
  geom_polygon(data = coast_ams, aes(x = long, y = lat, group = group), col = 'grey30', fill = 'grey95', size = 0.1) +
  coord_equal(expand = 0) +  
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 7),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(-40, -38, -36), labels = c('-40º', '-38º', '-36º')) + 
  scale_x_continuous(name = NULL, breaks = c(75, 80), labels = c('75º', '80º')) +
  ggsn::scalebar(x.min = min(ams_df$x), x.max = max(ams_df$x), y.min = min(ams_df$y), y.max = max(ams_df$y), transform = T, 
                 dist = 50, st.size = 1.75, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = 80, y = -41.1), box.fill = 'black')
ggsave('Figure S9_ams.tiff', dpi = 900, width = 6.4, height = 6.1, units = 'cm')


#----------------------------------- Table 2 ---------------------------------------------

# Evaluation scores

dat <- read.csv('Appendix C.csv', header = T)

swa_occ <- subset(dat, Region == 'Southwest Atlantic')
swa_occ1 <- subset(swa_occ, Source1 == 'Published literature')
swa_occ2 <- subset(swa_occ, Source1 != 'Published literature')
swa_occ1$Classifier <- extract(swa_final, SpatialPoints(swa_occ1[, c('Longitude', 'Latitude')], crs)) # extract values
swa_occ2$Classifier <- extract(swa_final, SpatialPoints(swa_occ2[, c('Longitude', 'Latitude')], crs)) # extract values
table(swa_occ1$Classifier) # calculate %
table(swa_occ2$Classifier)

aus_occ <- subset(dat, Region == 'Australia')
aus_occ1 <- subset(aus_occ, Source1 == 'Published literature')
aus_occ2 <- subset(aus_occ, Source1 != 'Published literature')
aus_occ1$Classifier <- extract(aus_final, SpatialPoints(aus_occ1[, c('Longitude', 'Latitude')], crs)) 
aus_occ2$Classifier <- extract(aus_final, SpatialPoints(aus_occ2[, c('Longitude', 'Latitude')], crs)) 
table(aus_occ1$Classifier)
table(aus_occ2$Classifier)

nep_occ <- subset(dat, Region == 'Northeast Pacific')
nep_occ1 <- subset(nep_occ, Source1 == 'Published literature')
nep_occ2 <- subset(nep_occ, Source1 != 'Published literature')
nep_occ1$Classifier <- extract(nep_final, SpatialPoints(nep_occ1[, c('Longitude', 'Latitude')], crs)) 
nep_occ2$Classifier <- extract(nep_final, SpatialPoints(nep_occ2[, c('Longitude', 'Latitude')], crs)) 
table(nep_occ1$Classifier)
table(nep_occ2$Classifier)

sea_occ <- subset(dat, Region == 'Southeast Atlantic')
sea_occ1 <- subset(sea_occ, Source1 == 'Published literature')
sea_occ2 <- subset(sea_occ, Source1 != 'Published literature')
sea_occ1$Classifier <- extract(sea_final, SpatialPoints(sea_occ1[, c('Longitude', 'Latitude')], crs)) 
sea_occ2$Classifier <- extract(sea_final, SpatialPoints(sea_occ2[, c('Longitude', 'Latitude')], crs)) 
table(sea_occ1$Classifier)
table(sea_occ2$Classifier)

nze_occ <- subset(dat, Region == 'New Zealand')
nze_occ1 <- subset(nze_occ, Source1 == 'Published literature')
nze_occ2 <- subset(nze_occ, Source1 != 'Published literature')
nze_occ1$Classifier <- extract(nze_final, SpatialPoints(nze_occ1[, c('Longitude', 'Latitude')], crs)) 
nze_occ2$Classifier <- extract(nze_final, SpatialPoints(nze_occ2[, c('Longitude', 'Latitude')], crs)) 
table(nze_occ1$Classifier)
table(nze_occ2$Classifier)

sep_occ <- subset(dat, Region == 'Southeast Pacific')
sep_occ1 <- subset(sep_occ, Source1 == 'Published literature')
sep_occ2 <- subset(sep_occ, Source1 != 'Published literature')
sep_occ1$Classifier <- extract(sep_final, SpatialPoints(sep_occ1[, c('Longitude', 'Latitude')], crs)) 
sep_occ2$Classifier <- extract(sep_final, SpatialPoints(sep_occ2[, c('Longitude', 'Latitude')], crs)) 
table(sep_occ1$Classifier)
table(sep_occ2$Classifier)

nwp_occ <- subset(dat, Region == 'Northwest Pacific')
nwp_occ1 <- subset(nwp_occ, Source1 == 'Published literature')
nwp_occ2 <- subset(nwp_occ, Source1 != 'Published literature')
nwp_occ1$Classifier <- extract(nwp_final, SpatialPoints(nwp_occ1[, c('Longitude', 'Latitude')], crs)) 
nwp_occ2$Classifier <- extract(nwp_final, SpatialPoints(nwp_occ2[, c('Longitude', 'Latitude')], crs)) 
table(nwp_occ1$Classifier)
table(nwp_occ2$Classifier)

tdc_occ <- subset(dat, Region == 'Tristan da Cunha Is.')
tdc_occ1 <- subset(tdc_occ, Source1 == 'Published literature')
tdc_occ2 <- subset(tdc_occ, Source1 != 'Published literature')
tdc_occ1$Classifier <- extract(tdc_final, SpatialPoints(tdc_occ1[, c('Longitude', 'Latitude')], crs)) 
tdc_occ2$Classifier <- extract(tdc_final, SpatialPoints(tdc_occ2[, c('Longitude', 'Latitude')], crs)) 
table(tdc_occ1$Classifier)
table(tdc_occ2$Classifier)

gal_occ <- subset(dat, Region == 'Galápagos Is.')
gal_occ1 <- subset(gal_occ, Source1 == 'Published literature')
gal_occ1$Classifier <- extract(gal_final, SpatialPoints(gal_occ1[, c('Longitude', 'Latitude')], crs)) 
table(gal_occ1$Classifier)


#----------------------------------- Table S1 ---------------------------------------------

# Global and regional overlap scores in binary transfers per continental region with occurrence data

Region1 <- 'SWA'
Region2 <- 'AUS'

swa_raw <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
aus_raw <- raster(paste(Region2, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
bg_binary <- swa_raw
bg_binary[bg_binary >= 0] <- 1 # background raster

# All occurrences used for modelling
swa_occ <- read.csv(paste(Region1, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences
aus_occ <- read.csv(paste(Region2, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences

# 10% Minimum training presence thresholds
swa_mtp10 <- sdm_threshold.10(swa_raw, swa_occ[, c('longitude', 'latitude')], 'p10', binary = F)
swa_bin_10 <- swa_raw
swa_bin_10[swa_bin_10 < minValue(swa_mtp10)] <- NA
swa_bin_10[swa_bin_10 >= minValue(swa_mtp10)] <- 2
aus_mtp10 <- sdm_threshold.10(aus_raw, aus_occ[, c('longitude', 'latitude')], 'p10', binary = F)
aus_bin_10 <- aus_raw
aus_bin_10[aus_bin_10 < minValue(aus_mtp10)] <- NA
aus_bin_10[aus_bin_10 >= minValue(aus_mtp10)] <- 3

# MOP mosaic
swa_mop <- raster(paste(Region1, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
swa_mop[swa_mop > 0] <- NA
swa_mop[swa_mop == 0] <- 8
aus_mop <- raster(paste(Region2, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
aus_mop[aus_mop > 0] <- NA
aus_mop[aus_mop == 0] <- 8
mop_mosaic <- mosaic(swa_mop, aus_mop, fun = sum)
mop_mosaic[mop_mosaic < 8] <- NA
mop_mosaic[mop_mosaic >= 8] <- 8 # areas of strict extrapolation

# SWA mosaic
swa_mosaic <- mosaic(bg_binary, swa_bin_10, mop_mosaic, fun = sum)
swa_mosaic[swa_mosaic > 3] <- 8

# AUS mosaic
aus_mosaic <- mosaic(bg_binary, aus_bin_10, mop_mosaic, fun = sum)
aus_mosaic[aus_mosaic > 4] <- 8

# Final mosaic
mod_mosaic <- mosaic(swa_mosaic, aus_mosaic, fun = sum)
mod_mosaic[mod_mosaic == 2] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 4] <- 2 # swa MTP
mod_mosaic[mod_mosaic == 5] <- 3 # aus MTP
mod_mosaic[mod_mosaic == 7] <- 4 # overlap between swa MTP and aus MTP
mod_mosaic[mod_mosaic == 16] <- 5 # strict extrapolation areas
mod_mosaic[mod_mosaic == 9] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 11] <- 5 # strict extrapolation areas
mod_mosaic[mod_mosaic == 12] <- 5 # strict extrapolation areas
mod_mosaic <- crop(mod_mosaic, coast)

# Bounding box of regions predicted as suitable, crop to region and % overlap calculation 

# Global
length(mod_mosaic[mod_mosaic %in% c(2, 4)]) * 100 / length(mod_mosaic[mod_mosaic %in% 2:4]) # % of suitable areas predicted by swa
length(mod_mosaic[mod_mosaic %in% c(3, 4)]) * 100 / length(mod_mosaic[mod_mosaic %in% 2:4]) # % of suitable areas predicted by aus
length(mod_mosaic[mod_mosaic == 4]) * 100 / length(mod_mosaic[mod_mosaic %in% 2:4])         # % of suitable areas predicted by both
length(mod_mosaic[mod_mosaic == 4]) * 100 / length(mod_mosaic[mod_mosaic %in% c(2, 4)])     # % of suitable areas predicted by swa also predicted by aus
length(mod_mosaic[mod_mosaic == 4]) * 100 / length(mod_mosaic[mod_mosaic %in% c(3, 4)])     # % of suitable areas predicted by aus also predicted by swa
length(mod_mosaic[mod_mosaic %in% c(3, 4)]) / length(mod_mosaic[mod_mosaic %in% c(2, 4)])   # swa/aus suitability ratio

# Southwest Atlantic
swa_ext <- extent(-68.5, -48, -49, -29)        
swa_final <- crop(mod_mosaic, swa_ext) 
length(swa_final[swa_final %in% c(2, 4)]) * 100 / length(swa_final[swa_final %in% 2:4]) 
length(swa_final[swa_final %in% c(3, 4)]) * 100 / length(swa_final[swa_final %in% 2:4]) 
length(swa_final[swa_final == 4]) * 100 / length(swa_final[swa_final %in% 2:4])         
length(swa_final[swa_final == 4]) * 100 / length(swa_final[swa_final %in% c(2, 4)])     
length(swa_final[swa_final == 4]) * 100 / length(swa_final[swa_final %in% c(3, 4)])     
length(swa_final[swa_final %in% c(3, 4)]) / length(swa_final[swa_final %in% c(2, 4)])   

# Australia
aus_ext <- extent(111, 155, -46.5, -24.5)
aus_final <- crop(mod_mosaic, aus_ext) 
length(aus_final[aus_final %in% c(2, 4)]) * 100 / length(aus_final[aus_final %in% 2:4]) 
length(aus_final[aus_final %in% c(3, 4)]) * 100 / length(aus_final[aus_final %in% 2:4]) 
length(aus_final[aus_final == 4]) * 100 / length(aus_final[aus_final %in% 2:4])         
length(aus_final[aus_final == 4]) * 100 / length(aus_final[aus_final %in% c(2, 4)])     
length(aus_final[aus_final == 4]) * 100 / length(aus_final[aus_final %in% c(3, 4)])     
length(aus_final[aus_final %in% c(3, 4)]) / length(aus_final[aus_final %in% c(2, 4)])   

# Northeast Pacific
nep_ext <- extent(-135, -107, 21.5, 56)        
nep_final <- crop(mod_mosaic, nep_ext) 
length(nep_final[nep_final %in% c(2, 4)]) * 100 / length(nep_final[nep_final %in% 2:4]) 
length(nep_final[nep_final %in% c(3, 4)]) * 100 / length(nep_final[nep_final %in% 2:4]) 
length(nep_final[nep_final == 4]) * 100 / length(nep_final[nep_final %in% 2:4])         
length(nep_final[nep_final == 4]) * 100 / length(nep_final[nep_final %in% c(2, 4)])     
length(nep_final[nep_final == 4]) * 100 / length(nep_final[nep_final %in% c(3, 4)])    
length(nep_final[nep_final %in% c(3, 4)]) / length(nep_final[nep_final %in% c(2, 4)])  

# Northwest Pacific
nwp_ext <- extent(116, 146.5, 27, 46.5)        
nwp_final <- crop(mod_mosaic, nwp_ext) 
length(nwp_final[nwp_final %in% c(2, 4)]) * 100 / length(nwp_final[nwp_final %in% 2:4]) 
length(nwp_final[nwp_final %in% c(3, 4)]) * 100 / length(nwp_final[nwp_final %in% 2:4]) 
length(nwp_final[nwp_final == 4]) * 100 / length(nwp_final[nwp_final %in% 2:4])         
length(nwp_final[nwp_final == 4]) * 100 / length(nwp_final[nwp_final %in% c(2, 4)])    
length(nwp_final[nwp_final == 4]) * 100 / length(nwp_final[nwp_final %in% c(3, 4)])    
length(nwp_final[nwp_final %in% c(3, 4)]) / length(nwp_final[nwp_final %in% c(2, 4)])   

# Southeast Pacific
sep_ext <- extent(-84, -69, -52, -0.5)         
sep_final <- crop(mod_mosaic, sep_ext) 
length(sep_final[sep_final %in% c(2, 4)]) * 100 / length(sep_final[sep_final %in% 2:4]) 
length(sep_final[sep_final %in% c(3, 4)]) * 100 / length(sep_final[sep_final %in% 2:4]) 
length(sep_final[sep_final == 4]) * 100 / length(sep_final[sep_final %in% 2:4])         
length(sep_final[sep_final == 4]) * 100 / length(sep_final[sep_final %in% c(2, 4)])     
length(sep_final[sep_final == 4]) * 100 / length(sep_final[sep_final %in% c(3, 4)])     
length(sep_final[sep_final %in% c(3, 4)]) / length(sep_final[sep_final %in% c(2, 4)])   

# Southeast Atlantic
sea_ext <- extent(9, 33, -37.5, -13)           
sea_final <- crop(mod_mosaic, sea_ext) 
length(sea_final[sea_final %in% c(2, 4)]) * 100 / length(sea_final[sea_final %in% 2:4]) 
length(sea_final[sea_final %in% c(3, 4)]) * 100 / length(sea_final[sea_final %in% 2:4]) 
length(sea_final[sea_final == 4]) * 100 / length(sea_final[sea_final %in% 2:4])         
length(sea_final[sea_final == 4]) * 100 / length(sea_final[sea_final %in% c(2, 4)])    
length(sea_final[sea_final == 4]) * 100 / length(sea_final[sea_final %in% c(3, 4)])     
length(sea_final[sea_final %in% c(3, 4)]) / length(sea_final[sea_final %in% c(2, 4)])  

# New Zealand
nze1_ext <- extent(163, 179.99, -52, -27)      
nze2_ext <- extent(-179.99, -172, -52, -27)    
nze1_final <- crop(mod_mosaic, nze1_ext) 
nze2_final <- crop(mod_mosaic, nze2_ext) 
(length(nze1_final[nze1_final %in% c(2, 4)]) + length(nze2_final[nze2_final %in% c(2, 4)])) * 100 / (length(nze1_final[nze1_final %in% 2:4]) + length(nze2_final[nze2_final %in% 2:4]))     
(length(nze1_final[nze1_final %in% c(3, 4)]) + length(nze2_final[nze2_final %in% c(3, 4)])) * 100 / (length(nze1_final[nze1_final %in% 2:4]) + length(nze2_final[nze2_final %in% 2:4]))     
(length(nze1_final[nze1_final == 4]) + length(nze2_final[nze2_final == 4])) * 100 / (length(nze1_final[nze1_final %in% 2:4]) + length(nze2_final[nze2_final %in% 2:4]))                     
(length(nze1_final[nze1_final == 4]) + length(nze2_final[nze2_final == 4])) * 100 / (length(nze1_final[nze1_final %in% c(2, 4)]) + length(nze2_final[nze2_final %in% c(2, 4)]))             
(length(nze1_final[nze1_final == 4]) + length(nze2_final[nze2_final == 4])) * 100 / (length(nze1_final[nze1_final %in% c(3, 4)]) + length(nze2_final[nze2_final %in% c(3, 4)]))             
(length(nze1_final[nze1_final %in% c(3, 4)]) + length(nze2_final[nze2_final %in% c(3, 4)])) / (length(nze1_final[nze1_final %in% c(2, 4)]) + length(nze2_final[nze2_final %in% c(2, 4)]))   

# Northwest Atlantic
nwa_ext <- extent(-78, -56.5, 34, 45)          
nwa_final <- crop(mod_mosaic, nwa_ext) 
length(nwa_final[nwa_final %in% c(2, 4)]) * 100 / length(nwa_final[nwa_final %in% 2:4]) 
length(nwa_final[nwa_final %in% c(3, 4)]) * 100 / length(nwa_final[nwa_final %in% 2:4])
length(nwa_final[nwa_final == 4]) * 100 / length(nwa_final[nwa_final %in% 2:4])        
length(nwa_final[nwa_final == 4]) * 100 / length(nwa_final[nwa_final %in% c(2, 4)])    
length(nwa_final[nwa_final == 4]) * 100 / length(nwa_final[nwa_final %in% c(3, 4)])    
length(nwa_final[nwa_final %in% c(3, 4)]) / length(nwa_final[nwa_final %in% c(2, 4)])   

# Northeast Atlantic
nea_ext <- extent(-20.5, 42.5, 13, 64)         
nea_final <- crop(mod_mosaic, nea_ext) 
length(nea_final[nea_final %in% c(2, 4)]) * 100 / length(nea_final[nea_final %in% 2:4]) 
length(nea_final[nea_final %in% c(3, 4)]) * 100 / length(nea_final[nea_final %in% 2:4])
length(nea_final[nea_final == 4]) * 100 / length(nea_final[nea_final %in% 2:4])        
length(nea_final[nea_final == 4]) * 100 / length(nea_final[nea_final %in% c(2, 4)])    
length(nea_final[nea_final == 4]) * 100 / length(nea_final[nea_final %in% c(3, 4)])    
length(nea_final[nea_final %in% c(3, 4)]) / length(nea_final[nea_final %in% c(2, 4)])   

# Galápagos Is.
gal_ext <- extent(-93, -87.5, -2.5, 1.25)      
gal_final <- crop(mod_mosaic, gal_ext) 
length(gal_final[gal_final %in% c(2, 4)]) * 100 / length(gal_final[gal_final %in% 2:4])
length(gal_final[gal_final %in% c(3, 4)]) * 100 / length(gal_final[gal_final %in% 2:4])
length(gal_final[gal_final == 4]) * 100 / length(gal_final[gal_final %in% 2:4])        

# Azores Is.
azo_ext <- extent(-34, -22.25, 35, 42)         
azo_final <- crop(mod_mosaic, azo_ext) 
length(azo_final[azo_final %in% c(2, 4)]) * 100 / length(azo_final[azo_final %in% 2:4]) 
length(azo_final[azo_final %in% c(3, 4)]) * 100 / length(azo_final[azo_final %in% 2:4])
length(azo_final[azo_final == 4]) * 100 / length(azo_final[azo_final %in% 2:4])        
length(azo_final[azo_final == 4]) * 100 / length(azo_final[azo_final %in% c(2, 4)])    
length(azo_final[azo_final == 4]) * 100 / length(azo_final[azo_final %in% c(3, 4)])    
length(azo_final[azo_final %in% c(3, 4)]) / length(azo_final[azo_final %in% c(2, 4)])   

# Tristan da Cunha Is.
tdc_ext <- extent(-15.5, -7, -43, -35)         
tdc_final <- crop(mod_mosaic, tdc_ext) 
length(tdc_final[tdc_final %in% c(2, 4)]) * 100 / length(tdc_final[tdc_final %in% 2:4]) 
length(tdc_final[tdc_final %in% c(3, 4)]) * 100 / length(tdc_final[tdc_final %in% 2:4])
length(tdc_final[tdc_final == 4]) * 100 / length(tdc_final[tdc_final %in% 2:4])        
length(tdc_final[tdc_final == 4]) * 100 / length(tdc_final[tdc_final %in% c(2, 4)])    
length(tdc_final[tdc_final == 4]) * 100 / length(tdc_final[tdc_final %in% c(3, 4)])    
length(tdc_final[tdc_final %in% c(3, 4)]) / length(tdc_final[tdc_final %in% c(2, 4)])   

# Amsterdam I.
ams_ext <- extent(74.5, 80.5, -41.5, -35.25)   
ams_final <- crop(mod_mosaic, ams_ext) 
length(ams_final[ams_final %in% c(2, 4)]) * 100 / length(ams_final[ams_final %in% 2:4]) 
length(ams_final[ams_final %in% c(3, 4)]) * 100 / length(ams_final[ams_final %in% 2:4])
length(ams_final[ams_final == 4]) * 100 / length(ams_final[ams_final %in% 2:4])        
length(ams_final[ams_final == 4]) * 100 / length(ams_final[ams_final %in% c(2, 4)])    
length(ams_final[ams_final == 4]) * 100 / length(ams_final[ams_final %in% c(3, 4)])    
length(ams_final[ams_final %in% c(3, 4)]) / length(ams_final[ams_final %in% c(2, 4)])   



#----------------------------------- Figure 6 ---------------------------------------------

# Dispersal and connectivity between regional suitable patch areas
# remotes::install_github("luismurao/bam")

library(bam)
library(matrixStats)
library(animation)
library(furrr)
library(magrittr)

Region1 <- 'SWA'
Region2 <- 'AUS'

set.seed(111)

# Functions from Nuñez-Penichet et al. (2021) available at https://github.com/townpeterson/vespa
source('Functions_Nuñez_Penichet_2021.R')

# Preparing data for simulations

# Read calibration and global presence data
swa_occ <- read.csv(paste(Region1, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences
aus_occ <- read.csv(paste(Region2, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences
calib_occ <- rbind(swa_occ, aus_occ)
global_occ <- read.csv('Appendix C.csv')

# Read SDM outputs
swa_raw <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
aus_raw <- raster(paste(Region2, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
mod_mosaic <- mosaic(swa_raw, aus_raw, fun = median) # final model consensus
mod_mosaic <- crop(mod_mosaic, coast)

# Add seamount/knoll data expected suitability based on inverse distance weighting from centroid 
# By relevant regions because computation is limiting

# Load seamount/knoll polygon areas as suitable stepping-stones-like dispersal corridors
# Data downloaded from Yesson et al. (2011)
Seamounts <- readOGR(dsn = 'D:/D/ENV SDM/Seamounts/01_Data/SeamountsBaseArea', layer = 'SeamountsBaseArea')
Seamounts <- Seamounts[which(as.numeric(Seamounts$DEPTH) > -601), ] # ecologically meaningful
Knolls <- readOGR(dsn = 'D:/D/ENV SDM/Seamounts/01_Data/KnollsBaseArea', layer = 'KnollsBaseArea')
Knolls <- Knolls[which(-as.numeric(Knolls$DEPTH) > -601), ] # ecologically meaningful
seam_knoll <- rbind(Seamounts, Knolls)
seam_knoll <- crop(seam_knoll, mod_mosaic)

# East Pacific corridor
mod_mosaic_EP <- crop(mod_mosaic, extent(-120, -77, -10, 32))
crs(mod_mosaic_EP) <- CRS('+init=epsg:4326')
seam_knoll_EP <- crop(seam_knoll, extent(-120, -77, -10, 32))
centroids_EP <- gCentroid(as(seam_knoll_EP, 'SpatialPolygons'), byid = T)
seam_knoll_EP_dis <- distanceFromPoints(mod_mosaic_EP, centroids_EP) # calculate distance to nearest seamount
seam_knoll_EP_dis <- mask(seam_knoll_EP_dis/1000, mod_mosaic_EP) # mask to study area
centroids_EP_buf <- buffer(centroids_EP, width = 160000) # create buffers based on dispersal ability of species
seam_knoll_EP_dis <- mask(seam_knoll_EP_dis, centroids_EP_buf)
seam_knoll_EP_dis <- (seam_knoll_EP_dis / (maxValue(seam_knoll_EP_dis) / maxValue(mod_mosaic_EP))) / 2
mod_mosaic_EP <- extend(seam_knoll_EP_dis, extent(mod_mosaic))
mod_mosaic_EP <- mod_mosaic_EP * -1 + maxValue(mod_mosaic_EP) + minValue(mod_mosaic_EP) # turn around suitability values

# Southwest Atlantic corridor
mod_mosaic_SWA <- crop(mod_mosaic, extent(-54, -10, -39, -26))
crs(mod_mosaic_SWA) <- CRS('+init=epsg:4326')
seam_knoll_SWA <- crop(seam_knoll, extent(-54, -10, -39, -26))
centroids_SWA <- gCentroid(as(seam_knoll_SWA, 'SpatialPolygons'), byid = T)
seam_knoll_SWA_dis <- distanceFromPoints(mod_mosaic_SWA, centroids_SWA) 
seam_knoll_SWA_dis <- mask(seam_knoll_SWA_dis/1000, mod_mosaic_SWA) 
centroids_SWA_buf <- buffer(centroids_SWA, width = 160000) 
seam_knoll_SWA_dis <- mask(seam_knoll_SWA_dis, centroids_SWA_buf)
seam_knoll_SWA_dis <- (seam_knoll_SWA_dis / (maxValue(seam_knoll_SWA_dis) / maxValue(mod_mosaic_SWA))) / 2
mod_mosaic_SWA <- extend(seam_knoll_SWA_dis, extent(mod_mosaic))
mod_mosaic_SWA <- mod_mosaic_SWA * -1 + maxValue(mod_mosaic_SWA) + minValue(mod_mosaic_SWA) 

# Southeast Atlantic corridor
mod_mosaic_SEA <- crop(mod_mosaic, extent(-15, 15, -39, -5))
crs(mod_mosaic_SEA) <- CRS('+init=epsg:4326')
seam_knoll_SEA <- crop(seam_knoll, extent(-15, 15, -39, -5))
centroids_SEA <- gCentroid(as(seam_knoll_SEA, 'SpatialPolygons'), byid = T)
seam_knoll_SEA_dis <- distanceFromPoints(mod_mosaic_SEA, centroids_SEA) 
seam_knoll_SEA_dis <- mask(seam_knoll_SEA_dis/1000, mod_mosaic_SEA) 
centroids_SEA_buf <- buffer(centroids_SEA, width = 160000) 
seam_knoll_SEA_dis <- mask(seam_knoll_SEA_dis, centroids_SEA_buf)
seam_knoll_SEA_dis <- (seam_knoll_SEA_dis / (maxValue(seam_knoll_SEA_dis) / maxValue(mod_mosaic_SEA))) / 2
mod_mosaic_SEA <- extend(seam_knoll_SEA_dis, extent(mod_mosaic))
mod_mosaic_SEA <- mod_mosaic_SEA * -1 + maxValue(mod_mosaic_SEA) + minValue(mod_mosaic_SEA) 

# Oceania corridor
mod_mosaic_OCE <- crop(mod_mosaic, extent(148, 177, -40, -19))
crs(mod_mosaic_OCE) <- CRS('+init=epsg:4326')
seam_knoll_OCE <- crop(seam_knoll, extent(148, 177, -40, -19))
centroids_OCE <- gCentroid(as(seam_knoll_OCE, 'SpatialPolygons'), byid = T)
seam_knoll_OCE_dis <- distanceFromPoints(mod_mosaic_OCE, centroids_OCE) 
seam_knoll_OCE_dis <- mask(seam_knoll_OCE_dis/1000, mod_mosaic_OCE) 
centroids_OCE_buf <- buffer(centroids_OCE, width = 160000) 
seam_knoll_OCE_dis <- mask(seam_knoll_OCE_dis, centroids_OCE_buf)
seam_knoll_OCE_dis <- (seam_knoll_OCE_dis / (maxValue(seam_knoll_OCE_dis) / maxValue(mod_mosaic_OCE))) / 2
mod_mosaic_OCE <- extend(seam_knoll_OCE_dis, extent(mod_mosaic))
mod_mosaic_OCE <- mod_mosaic_OCE * -1 + maxValue(mod_mosaic_OCE) + minValue(mod_mosaic_OCE)

# Northwest Pacific corridor
mod_mosaic_NWP <- crop(mod_mosaic, extent(132, 158, 9, 38))
crs(mod_mosaic_NWP) <- CRS('+init=epsg:4326')
seam_knoll_NWP <- crop(seam_knoll, extent(132, 158, 9, 38))
centroids_NWP <- gCentroid(as(seam_knoll_NWP, 'SpatialPolygons'), byid = T)
seam_knoll_NWP_dis <- distanceFromPoints(mod_mosaic_NWP, centroids_NWP) 
seam_knoll_NWP_dis <- mask(seam_knoll_NWP_dis/1000, mod_mosaic_NWP) 
centroids_NWP_buf <- buffer(centroids_NWP, width = 160000) 
seam_knoll_NWP_dis <- mask(seam_knoll_NWP_dis, centroids_NWP_buf)
seam_knoll_NWP_dis <- (seam_knoll_NWP_dis / (maxValue(seam_knoll_NWP_dis) / maxValue(mod_mosaic_NWP))) / 2
mod_mosaic_NWP <- extend(seam_knoll_NWP_dis, extent(mod_mosaic))
mod_mosaic_NWP <- mod_mosaic_NWP * -1 + maxValue(mod_mosaic_NWP) + minValue(mod_mosaic_NWP)

# Chilean Islands corridor
mod_mosaic_CHI <- crop(mod_mosaic, extent(-86, -69, -39, -23))
crs(mod_mosaic_CHI) <- CRS('+init=epsg:4326')
seam_knoll_CHI <- crop(seam_knoll, extent(-86, -69, -39, -23))
centroids_CHI <- gCentroid(as(seam_knoll_CHI, 'SpatialPolygons'), byid = T)
seam_knoll_CHI_dis <- distanceFromPoints(mod_mosaic_CHI, centroids_CHI) 
seam_knoll_CHI_dis <- mask(seam_knoll_CHI_dis/1000, mod_mosaic_CHI) 
centroids_CHI_buf <- buffer(centroids_CHI, width = 160000) 
seam_knoll_CHI_dis <- mask(seam_knoll_CHI_dis, centroids_CHI_buf)
seam_knoll_CHI_dis <- (seam_knoll_CHI_dis / (maxValue(seam_knoll_CHI_dis) / maxValue(mod_mosaic_CHI))) / 2
mod_mosaic_CHI <- extend(seam_knoll_CHI_dis, extent(mod_mosaic))
mod_mosaic_CHI <- mod_mosaic_CHI * -1 + maxValue(mod_mosaic_CHI) + minValue(mod_mosaic_CHI)

# Add expected suitability at seamounts/knolls to ENM final model
mod_mosaic <- mosaic(mod_mosaic, mod_mosaic_EP, mod_mosaic_SWA, mod_mosaic_SEA, 
                    mod_mosaic_OCE, mod_mosaic_NWP, mod_mosaic_CHI, fun = max) 

# Extract suitability values at presence points
nat_suits <- extract(mod_mosaic, SpatialPoints(calib_occ[, 2:3], CRS('+init=epsg:4326')))
inv_suits <- extract(mod_mosaic, SpatialPoints(global_occ[, 3:4], CRS('+init=epsg:4326')))

# Remove NA values
inv_suits <- na.omit(inv_suits)
nat_suits <- na.omit(nat_suits)

# suitability values of invasion and native points
all_suits <- data.frame(model = c(nat_suits, inv_suits))

# Create a raster stack from SDMs outputs
models_amp <- mod_mosaic
nlayers <- nlayers(models_amp)

# Define vector of threshold values
th_vec <- c(0.01, 0.10)

# Apply threshold function for each threshold value and model
rst_vals <- lapply(th_vec, function(x){
  mods_ras <- lapply(1:nlayers, function(y){
    r <- thfunc(median_mod = models_amp[[y]],
                suits = all_suits[, y], percent = x)
  })
  return(mods_ras)
})

# Suitability threshold values for 3%
suits_minimas <- unlist(rst_vals[1])

# Suitability threshold values for 10%
suits_maximas <- unlist(rst_vals[2])

# All suitability values in a matrix of 
# 2 (thresholds for 1% and 10%) X 1 (one and only model)
smin_max <-rbind(suits_minimas, suits_maximas)

# Binarize models using thresholds in the smin_max matrix 
maps_x <- lapply(1:ncol(smin_max), function(x){
  # A sequence of ten suitability values between 1% and 10% threshold
  umbrales <- seq(smin_max[1, x],
                  smin_max[2, x],
                  length.out = 10)
  
  # Binarize models
  r1 <- lapply(seq_along(umbrales), function(y){
    r0 <-models_amp[[x]] >= umbrales[y]
  })
  names(r1) <- paste0('th_', round(umbrales, 2))
  r1 <- stack(r1)
  return(r1)
})
maps_x <- maps_x[[1]]


## East Pacific corridor ##

maps_x_EP <- crop(maps_x, extent(-120, -77, -10, 32))
crs(maps_x_EP) <- CRS('+init=epsg:4326')
maps_x_EP <- stack(maps_x_EP)

# Dispersal simulations

# Connectivity matrices

# First, we convert the model into a diagonal sparse matrix
sparse_temp <-  bam::model2sparse(maps_x_EP[[1]])

# and, we convert the occurrence points into a sparse matrix
# Use all data in each region as seed for simulations 
sub1 <- subset(global_occ, Region == 'Northeast Pacific')
sub2 <- subset(global_occ, Region == 'Galápagos Is.')
sub3 <- subset(global_occ, Region == 'Southeast Pacific')
invasion_occ1 <- SpatialPoints(sub1[,c('Longitude', 'Latitude')], CRS('+init=epsg:4326'))
invasion_occ1 <- crop(invasion_occ1, maps_x_EP)
invasion_occ2 <- SpatialPoints(sub2[,c('Longitude', 'Latitude')], CRS('+init=epsg:4326'))
invasion_occ2 <- crop(invasion_occ2, maps_x_EP)
invasion_occ3 <- SpatialPoints(sub3[,c('Longitude', 'Latitude')], CRS('+init=epsg:4326'))
invasion_occ3 <- crop(invasion_occ3, maps_x_EP)
invasion_occ <- rbind(invasion_occ1, invasion_occ2, invasion_occ3)

occs_sparse <- bam::occs2sparse(modelsparse = sparse_temp, 
                                occs = data.frame(longitude = invasion_occ@coords[,1], 
                                                  latitude = invasion_occ@coords[,2]))

# Second, we create adjacency matrices for different distance values 
Values <- 1:10
pasos <- 200
for(i in Values){
  ad_L <- lapply(i, function(x){
    adj_matrices <- bam::adj_mat(sparse_temp, ngbs = x)
    print(x)
    return(adj_matrices)
  })
  
  # Name of the matrix 
  names(ad_L) <- i
  
  # Apply the simulation function
  simul_x_EP <- sim_disperal(adL = ad_L, bin_thmodels = maps_x_EP,
                            base_name = 'EastPacific dispersal', pasos = pasos)
}

# Calculate the sum of all the resulting maps for the different distance values
ras_EP <- stack()
for(i in Values){
  ras <- raster(paste('EastPacific dispersal/Sim_EastPacific dispersal_D_', i, '.tif', sep = ''))
  ras_EP <- stack(ras_EP, ras)
}
simul_x_EP_sum <- calc(ras_EP, sum)
writeRaster(simul_x_EP_sum, paste('EastPacific dispersal_', pasos, '.tif', sep = ''), format = 'GTiff', overwrite = T)

EPdis <- raster('EastPacific dispersal_200.tif')
EPdis[EPdis == 0] <- NA
df <- data.frame(coordinates(EPdis), as.data.frame(EPdis))
coast0 <- crop(coast, EPdis)
seam_knoll0 <- crop(seam_knoll, EPdis)
invasion_occ <- data.frame(coordinates(invasion_occ))
df_mexican_gulf <- data.frame(y = c( 14.989, 26.859, 27.872, 6.718, 7.613, 9.282, 8.323, 14.989),
                              x = c(-89.958, -90.1, -77, -77, -77.557, -78.977, -81.695, -89.958))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = EastPacific_dispersal_200)) + 
  geom_polygon(data = df_mexican_gulf, aes(x = x, y = y), color = 'grey95', fill = 'grey95', size = 0.1) +
  geom_polygon(data = seam_knoll0, aes(x = long, y = lat, group = group), col = 'darkorange', 
               fill = 'darkorange', size = 0.15) +
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), col = 'black', fill = 'black', size = 0.1) +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 2, color = 'white') +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 1, color = 'blue') +
  scale_fill_viridis(na.value = 'grey95') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 6),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) +
  scale_y_continuous(name = NULL, breaks = c(0, 20), labels = c('0º', '20º')) + 
  scale_x_continuous(name = NULL, breaks = c(-120, -100, -80), labels = c('-120º', '-100º', '-80º')) +
  ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T, 
                 dist = 400, st.size = 1.1, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -111.5, y = -8), box.fill = 'black')
ggsave('Figure 6_EP.tiff', dpi = 900, width = 5.7, height = 5.4, units = 'cm')


## Southwest Atlantic corridor ##

maps_x_SWA <- crop(maps_x, extent(-54, -10, -39, -26))
crs(maps_x_SWA) <- CRS('+init=epsg:4326')
maps_x_SWA <- stack(maps_x_SWA)

# Disperal simulations
sparse_temp <-  bam::model2sparse(maps_x_SWA[[1]])
sub1 <- subset(global_occ, Region == 'Southwest Atlantic')
sub2 <- subset(global_occ, Region == 'Tristan da Cunha Is.')
invasion_occ1 <- SpatialPoints(sub1[,c('Longitude', 'Latitude')], CRS('+init=epsg:4326'))
invasion_occ1 <- crop(invasion_occ1, maps_x_SWA)
invasion_occ2 <- SpatialPoints(sub2[,c('Longitude', 'Latitude')], CRS('+init=epsg:4326'))
invasion_occ2 <- crop(invasion_occ2, maps_x_SWA)
invasion_occ <- rbind(invasion_occ1, invasion_occ2)
occs_sparse <- bam::occs2sparse(modelsparse = sparse_temp, 
                                occs = data.frame(longitude = invasion_occ@coords[,1], 
                                                  latitude = invasion_occ@coords[,2]))
Values <- 1:10
pasos <- 200
for(i in Values){
  ad_L <- lapply(i, function(x){
    adj_matrices <- bam::adj_mat(sparse_temp, ngbs = x)
    print(x)
    return(adj_matrices)
  })
  
  # Name of the matrx 
  names(ad_L) <- i
  
  # Apply the simulation function to each scenario
  simul_x_SWA <- sim_disperal(adL = ad_L, bin_thmodels = maps_x_SWA,
                             base_name = 'SouthwestAtlantic dispersal', pasos = pasos)
}
ras_SWA <- stack()
for(i in Values){
  ras <- raster(paste('SouthwestAtlantic dispersal/Sim_SouthwestAtlantic dispersal_D_', i, '.tif', sep = ''))
  ras_SWA <- stack(ras_SWA, ras)
}
simul_x_SWA_sum <- calc(ras_SWA, sum)
writeRaster(simul_x_SWA_sum, paste('SouthwestAtlantic dispersal_', pasos, '.tif', sep = ''), format = 'GTiff', overwrite = T)

SWAdis <- raster('SouthwestAtlantic dispersal_200.tif')
SWAdis[SWAdis == 0] <- NA
df <- data.frame(coordinates(SWAdis), as.data.frame(SWAdis))
coast0 <- crop(coast, SWAdis)
seam_knoll0 <- crop(seam_knoll, SWAdis)
invasion_occ <- data.frame(coordinates(invasion_occ))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = SouthwestAtlantic_dispersal_200)) + 
  geom_polygon(data = seam_knoll0, aes(x = long, y = lat, group = group), col = 'darkorange', 
               fill = 'darkorange', size = 0.15) +
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), col = 'black', fill = 'black', size = 0.1) +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 2, color = 'white') +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 1, color = 'blue') +
  scale_fill_viridis(na.value = 'grey95') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 6),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) +
  scale_y_continuous(name = NULL, breaks = c(-35, -30), labels = c('-35º', '-30º')) + 
  scale_x_continuous(name = NULL, breaks = c(-50, -40, -30, -20), labels = c('-50º', '-40º', '-30º', '-20º')) +
  ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T, 
                 dist = 200, st.size = 1.1, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -48, y = -38.25), box.fill = 'black')
ggsave('Figure 6_SWA.tiff', dpi = 900, width = 11.4, height = 3.5, units = 'cm')


## Southeast Atlantic corridor ##

maps_x_SEA <- crop(maps_x, extent(-15, 15, -39, -5))
crs(maps_x_SEA) <- CRS('+init=epsg:4326')
maps_x_SEA <- stack(maps_x_SEA)

# Dispersal simulations
sparse_temp <-  bam::model2sparse(maps_x_SEA[[1]])
sub1 <- subset(global_occ, Region == 'Southeast Atlantic')
sub2 <- subset(global_occ, Region == 'Tristan da Cunha Is.')
invasion_occ1 <- SpatialPoints(sub1[,c('Longitude', 'Latitude')], CRS('+init=epsg:4326'))
invasion_occ1 <- crop(invasion_occ1, maps_x_SEA)
invasion_occ2 <- SpatialPoints(sub2[,c('Longitude', 'Latitude')], CRS('+init=epsg:4326'))
invasion_occ2 <- crop(invasion_occ2, maps_x_SEA)
invasion_occ <- rbind(invasion_occ1, invasion_occ2)
occs_sparse <- bam::occs2sparse(modelsparse = sparse_temp, 
                                occs = data.frame(longitude = invasion_occ@coords[,1], 
                                                  latitude = invasion_occ@coords[,2]))
Values <- 1:10
pasos <- 200
for(i in Values){
  ad_L <- lapply(i, function(x){
    adj_matrices <- bam::adj_mat(sparse_temp, ngbs = x)
    print(x)
    return(adj_matrices)
  })
  
  # Name of the matrx 
  names(ad_L) <- i
  
  # Apply the simulation function to each scenario
  simul_x_SEA <- sim_disperal(adL = ad_L, bin_thmodels = maps_x_SEA,
                              base_name = 'SoutheastAtlantic dispersal', pasos = pasos)
}
ras_SEA <- stack()
for(i in Values){
  ras <- raster(paste('SoutheastAtlantic dispersal/Sim_SoutheastAtlantic dispersal_D_', i, '.tif', sep = ''))
  ras_SEA <- stack(ras_SEA, ras)
}
simul_x_SEA_sum <- calc(ras_SEA, sum)
writeRaster(simul_x_SEA_sum, paste('SoutheastAtlantic dispersal_', pasos, '.tif', sep = ''), format = 'GTiff', overwrite = T)

SEAdis <- raster('SoutheastAtlantic dispersal_200.tif')
SEAdis <- crop(SEAdis, extent(-15, 15, -39, -8))
SEAdis[SEAdis == 0] <- NA
df <- data.frame(coordinates(SEAdis), as.data.frame(SEAdis))
coast0 <- crop(coast, SEAdis)
seam_knoll0 <- crop(seam_knoll, SEAdis)
invasion_occ <- data.frame(coordinates(invasion_occ))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = SoutheastAtlantic_dispersal_200)) + 
  geom_polygon(data = seam_knoll0, aes(x = long, y = lat, group = group), col = 'darkorange', 
               fill = 'darkorange', size = 0.15) +
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), col = 'black', fill = 'black', size = 0.1) +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 2, color = 'white') +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 1, color = 'blue') +
  scale_fill_viridis(na.value = 'grey95') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 6),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) +
  scale_y_continuous(name = NULL, breaks = c(-35, -25, -15), labels = c('-35º', '-25º', '-15º')) + 
  scale_x_continuous(name = NULL, breaks = c(-10, 0, 10), labels = c('-10º', '0º', '10º')) +
  ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T, 
                 dist = 200, st.size = 1.1, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = 13.5, y = -37.75), box.fill = 'black')
ggsave('Figure 6_SEA.tiff', dpi = 900, width = 4.8, height = 4.85, units = 'cm')


## Oceania corridor ##

maps_x_OCE <- crop(maps_x, extent(148, 177, -40, -19))
crs(maps_x_OCE) <- CRS('+init=epsg:4326')
maps_x_OCE <- stack(maps_x_OCE)

# Dispersal simulations
sparse_temp <-  bam::model2sparse(maps_x_OCE[[1]])
sub1 <- subset(global_occ, Region == 'Australia')
sub2 <- subset(global_occ, Region == 'New Zealand')
invasion_occ1 <- SpatialPoints(sub1[,c('Longitude', 'Latitude')], CRS('+init=epsg:4326'))
invasion_occ1 <- crop(invasion_occ1, maps_x_OCE)
invasion_occ2 <- SpatialPoints(sub2[,c('Longitude', 'Latitude')], CRS('+init=epsg:4326'))
invasion_occ2 <- crop(invasion_occ2, maps_x_OCE)
invasion_occ <- rbind(invasion_occ1, invasion_occ2)
occs_sparse <- bam::occs2sparse(modelsparse = sparse_temp, 
                                occs = data.frame(longitude = invasion_occ@coords[,1], 
                                                  latitude = invasion_occ@coords[,2]))
Values <- 1:10
pasos <- 200
for(i in Values){
  ad_L <- lapply(i, function(x){
    adj_matrices <- bam::adj_mat(sparse_temp, ngbs = x)
    print(x)
    return(adj_matrices)
  })
  
  # Name of the matrx 
  names(ad_L) <- i
  
  # Apply the simulation function to each scenario
  simul_x_OCE <- sim_disperal(adL = ad_L, bin_thmodels = maps_x_OCE,
                              base_name = 'Oceania dispersal', pasos = pasos)
}
ras_OCE <- stack()
for(i in Values){
  ras <- raster(paste('Oceania dispersal/Sim_Oceania dispersal_D_', i, '.tif', sep = ''))
  ras_OCE <- stack(ras_OCE, ras)
}
simul_x_OCE_sum <- calc(ras_OCE, sum)
writeRaster(simul_x_OCE_sum, paste('Oceania dispersal_', pasos, '.tif', sep = ''), format = 'GTiff', overwrite = T)

OCEdis <- raster('Oceania dispersal_200.tif')
OCEdis <- crop(OCEdis, extent(148, 177, -38, -19))
OCEdis[OCEdis == 0] <- NA
df <- data.frame(coordinates(OCEdis), as.data.frame(OCEdis))
coast0 <- crop(coast, OCEdis)
seam_knoll0 <- crop(seam_knoll, OCEdis)
invasion_occ <- data.frame(coordinates(invasion_occ))
invasion_occ <- subset(invasion_occ, Latitude >= -38)

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = Oceania_dispersal_200)) + 
  geom_polygon(data = seam_knoll0, aes(x = long, y = lat, group = group), col = 'darkorange', 
               fill = 'darkorange', size = 0.15) +
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), col = 'black', fill = 'black', size = 0.1) +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 1.5, color = 'white') +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 0.75, color = 'blue') +
  scale_fill_viridis(na.value = 'grey95') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 5),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) +
  scale_y_continuous(name = NULL, breaks = c(-35, -25), labels = c('-35º', '-25º')) + 
  scale_x_continuous(name = NULL, breaks = c(150, 160, 170), labels = c('150º', '160º', '170º')) +
  ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T, 
                 dist = 200, st.size = 1.1, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = 158, y = -37), box.fill = 'black')
ggsave('Figure 6_OCE.tiff', dpi = 900, width = 5.9, height = 3.9, units = 'cm')


## Northwest Pacific corridor ##

maps_x_NWP <- crop(maps_x, extent(132, 158, 9, 38))
crs(maps_x_NWP) <- CRS('+init=epsg:4326')
maps_x_NWP <- stack(maps_x_NWP)

# Dispersal simulations
sparse_temp <-  bam::model2sparse(maps_x_NWP[[1]])
sub1 <- subset(global_occ, Region == 'Northwest Pacific')
invasion_occ <- SpatialPoints(sub1[,c('Longitude', 'Latitude')], CRS('+init=epsg:4326'))
invasion_occ <- crop(invasion_occ, maps_x_NWP)
occs_sparse <- bam::occs2sparse(modelsparse = sparse_temp, 
                                occs = data.frame(longitude = invasion_occ@coords[,1], 
                                                  latitude = invasion_occ@coords[,2]))
Values <- 1:10
pasos <- 200
for(i in Values){
  ad_L <- lapply(i, function(x){
    adj_matrices <- bam::adj_mat(sparse_temp, ngbs = x)
    print(x)
    return(adj_matrices)
  })
  
  # Name of the matrx 
  names(ad_L) <- i
  
  # Apply the simulation function to each scenario
  simul_x_NWP <- sim_disperal(adL = ad_L, bin_thmodels = maps_x_NWP,
                              base_name = 'Northwest Pacific dispersal', pasos = pasos)
}
ras_NWP <- stack()
for(i in Values){
  ras <- raster(paste('Northwest Pacific dispersal/Sim_Northwest Pacific dispersal_D_', i, '.tif', sep = ''))
  ras_NWP <- stack(ras_NWP, ras)
}
simul_x_NWP_sum <- calc(ras_NWP, sum)
writeRaster(simul_x_NWP_sum, paste('Northwest Pacific dispersal_', pasos, '.tif', sep = ''), format = 'GTiff', overwrite = T)

NWPdis <- raster('Northwest Pacific dispersal_200.tif')
NWPdis[NWPdis == 0] <- NA
df <- data.frame(coordinates(NWPdis), as.data.frame(NWPdis))
coast0 <- crop(coast, NWPdis)
seam_knoll0 <- crop(seam_knoll, NWPdis)
invasion_occ <- data.frame(coordinates(invasion_occ))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = Northwest_Pacific_dispersal_200)) + 
  geom_polygon(data = seam_knoll0, aes(x = long, y = lat, group = group), col = 'darkorange', 
               fill = 'darkorange', size = 0.15) +
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), col = 'black', fill = 'black', size = 0.1) +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 2, color = 'white') +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 1, color = 'blue') +
  scale_fill_viridis(na.value = 'grey95') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 6),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) +
  scale_y_continuous(name = NULL, breaks = c(35, 25, 15), labels = c('35º', '25º', '15º')) + 
  scale_x_continuous(name = NULL, breaks = c(135, 145, 155), labels = c('135º', '145º', '155º')) +
  ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T, 
                 dist = 200, st.size = 1.1, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = 156.5, y = 10.2), box.fill = 'black')
ggsave('Figure 6_NWP.tiff', dpi = 900, width = 4.8, height = 5.2, units = 'cm')


## Chilean Islands corridor ##

maps_x_CHI <- crop(maps_x, extent(-86, -69, -39, -24))
crs(maps_x_CHI) <- CRS('+init=epsg:4326')
maps_x_CHI <- stack(maps_x_CHI)

# Dispersal simulations
sparse_temp <-  bam::model2sparse(maps_x_CHI[[1]])
sub1 <- subset(global_occ, Region == 'Southeast Pacific')
invasion_occ <- SpatialPoints(sub1[,c('Longitude', 'Latitude')], CRS('+init=epsg:4326'))
invasion_occ <- crop(invasion_occ, maps_x_CHI)
occs_sparse <- bam::occs2sparse(modelsparse = sparse_temp, 
                                occs = data.frame(longitude = invasion_occ@coords[,1], 
                                                  latitude = invasion_occ@coords[,2]))
Values <- 1:10
pasos <- 200
for(i in Values){
  ad_L <- lapply(i, function(x){
    adj_matrices <- bam::adj_mat(sparse_temp, ngbs = x)
    print(x)
    return(adj_matrices)
  })
  
  # Name of the matrx 
  names(ad_L) <- i
  
  # Apply the simulation function to each scenario
  simul_x_CHI <- sim_disperal(adL = ad_L, bin_thmodels = maps_x_CHI,
                              base_name = 'Chilean Islands dispersal', pasos = pasos)
}
ras_CHI <- stack()
for(i in Values){
  ras <- raster(paste('Chilean Islands dispersal/Sim_Chilean Islands dispersal_D_', i, '.tif', sep = ''))
  ras_CHI <- stack(ras_CHI, ras)
}
simul_x_CHI_sum <- calc(ras_CHI, sum)
writeRaster(simul_x_CHI_sum, paste('Chilean Islands dispersal_', pasos, '.tif', sep = ''), format = 'GTiff', overwrite = T)

CHIdis <- raster('Chilean Islands dispersal_200.tif')
CHIdis <- crop(CHIdis, extent(-84, -71, -37, -24))
CHIdis[CHIdis == 0] <- NA
df <- data.frame(coordinates(CHIdis), as.data.frame(CHIdis))
coast0 <- crop(coast, CHIdis)
seam_knoll0 <- crop(seam_knoll, CHIdis)
invasion_occ <- data.frame(coordinates(invasion_occ))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = Chilean_Islands_dispersal_200)) + 
  geom_polygon(data = seam_knoll0, aes(x = long, y = lat, group = group), col = 'darkorange', 
               fill = 'darkorange', size = 0.15) +
  geom_polygon(data = coast0, aes(x = long, y = lat, group = group), col = 'black', fill = 'black', size = 0.1) +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 2, color = 'white') +
  geom_point(data = invasion_occ, aes(x = Longitude, y = Latitude), size = 1, color = 'blue') +
  scale_fill_viridis(na.value = 'grey95') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 6),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) +
  scale_y_continuous(name = NULL, breaks = c(-35, -30, -25), labels = c('-35º', '-30º', '-25º')) + 
  scale_x_continuous(name = NULL, breaks = c(-80, -75), labels = c('-80º', '-75º')) +
  ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T, 
                 dist = 100, st.size = 1.1, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -81.5, y = -36.5), box.fill = 'black')
ggsave('Figure 6_CHI.tiff', dpi = 900, width = 5.6, height = 5.4, units = 'cm')


## World ##

ep_ext <- extent(-120, -77, -10, 32)        
swa_ext <- extent(-54, -10, -39, -26)      
sea_ext <- extent(-15, 15, -39, -8)        
oce_ext <- extent(148, 177, -38, -19) 
nwp_ext <- extent(132, 158, 9, 38)         
chi_ext <- extent(-84, -71, -37, -24)         

Ids <- c('EP', 'SWA', 'SEA', 'OCE', 'NWP', 'CHI')

polys <- SpatialPolygons(list(
  Polygons(list(Polygon(ep_ext)), Ids[1]), Polygons(list(Polygon(swa_ext)), Ids[2]),
  Polygons(list(Polygon(sea_ext)), Ids[3]), Polygons(list(Polygon(oce_ext)), Ids[4]),
  Polygons(list(Polygon(nwp_ext)), Ids[5]), Polygons(list(Polygon(chi_ext)), Ids[6])),
  proj4string = crs)

polys <- SpatialPolygonsDataFrame(polys, data.frame(ids = Ids, row.names = names(polys))) 

ggplot() +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), col = 'black', fill = 'black', size = 0.1) +
  geom_polygon(data = polys, aes(x = long, y = lat, group = group), col = 'grey40', fill = NA, size = 0.5) +
  coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 5),
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) +
  scale_y_continuous(name = NULL, breaks = c(50, 0, -50), labels = c('50ºN', '0º', '50ºS')) + 
  scale_x_continuous(name = NULL, breaks = c(100, 0, -100), labels = c('100ºE', '0º', '100ºW')) 
ggsave('Figure 6_world.tiff', dpi = 900, width = 10, height = 4.3, units = 'cm')


#----------------------------------- Figure 7 ---------------------------------------------

# Continuos SWA/AUS mosaic global transfers + dispersal routes

Region1 <- 'SWA'
Region2 <- 'AUS'

# Model mosaic
swa_raw <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
aus_raw <- raster(paste(Region2, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
mod_mosaic <- mosaic(swa_raw, aus_raw, fun = median) # final model consensus
mod_mosaic <- crop(mod_mosaic, coast)

# MOP mosaic
swa_mop <- raster(paste(Region1, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
swa_mop[swa_mop > 0] <- NA
swa_mop[swa_mop == 0] <- 1
aus_mop <- raster(paste(Region2, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
aus_mop[aus_mop > 0] <- NA
aus_mop[aus_mop == 0] <- 1
mop_mosaic <- mosaic(swa_mop, aus_mop, fun = sum)
mop_mosaic[mop_mosaic < 1] <- NA
mop_mosaic[mop_mosaic >= 1] <- 2 # areas of strict extrapolation
mop_mosaic <- crop(mop_mosaic, mod_mosaic)

# Final mosaic
mod_mosaic <- mosaic(mod_mosaic, mop_mosaic, fun = max) 
mod_mosaic[mod_mosaic > 1] <- NA # minus areas of strict extrapolation

# Load seamount/knoll polygon areas as suitable stepping-stones-like dispersal corridors
# Data downloaded from Yesson et al. (2011)
Seamounts <- readOGR(dsn = 'D:/D/ENV SDM/Seamounts/01_Data/SeamountsBaseArea', layer = 'SeamountsBaseArea')
Seamounts <- Seamounts[which(as.numeric(Seamounts$DEPTH) > -1001), ] # ecologically meaningful
Knolls <- readOGR(dsn = 'D:/D/ENV SDM/Seamounts/01_Data/KnollsBaseArea', layer = 'KnollsBaseArea')
Knolls <- Knolls[which(-as.numeric(Knolls$DEPTH) > -1001), ] # ecologically meaningful
seam_knoll <- rbind(Seamounts, Knolls)
seam_knoll <- crop(seam_knoll, mod_mosaic)

df <- data.frame(coordinates(mod_mosaic), as.data.frame(mod_mosaic))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = layer)) + 
  geom_polygon(data = seam_knoll, aes(x = long, y = lat, group = group), col = 'darkorange', 
               fill = 'darkorange', size = 0.15) +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), col = 'black', fill = 'grey10', size = 0.1) +
  scale_fill_viridis(na.value = '#440154FF') + coord_equal(expand = 0) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(50, 0, -50), labels = c('50ºN', '0º', '50ºS')) + 
  scale_x_continuous(name = NULL, breaks = c(100, 0, -100), labels = c('100ºW', '0º', '100ºE')) 
ggsave('Figure 7.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm')


#----------------------------------- Figure S10 ---------------------------------------------

# Binary SWA/AUS mosaic global transfers + dispersal routes

Region1 <- 'SWA'
Region2 <- 'AUS'

swa_raw <- raster(paste(Region1, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
aus_raw <- raster(paste(Region2, '/Final_Model_Stats/Statistics_E/current_median.tif', sep = '')) # best model
bg_binary <- swa_raw
bg_binary[bg_binary >= 0] <- 1 # background raster

# All occurrences used for modelling
swa_occ <- read.csv(paste(Region1, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences
aus_occ <- read.csv(paste(Region2, '/N_cepedianus_joint.csv', sep = ''), header = T) # occurrences

# 10% Minimum training presence thresholds
swa_mtp10 <- sdm_threshold.10(swa_raw, swa_occ[, c('longitude', 'latitude')], 'p10', binary = F)
swa_bin_10 <- swa_raw
swa_bin_10[swa_bin_10 < minValue(swa_mtp10)] <- NA
swa_bin_10[swa_bin_10 >= minValue(swa_mtp10)] <- 2
aus_mtp10 <- sdm_threshold.10(aus_raw, aus_occ[, c('longitude', 'latitude')], 'p10', binary = F)
aus_bin_10 <- aus_raw
aus_bin_10[aus_bin_10 < minValue(aus_mtp10)] <- NA
aus_bin_10[aus_bin_10 >= minValue(aus_mtp10)] <- 3

# MOP mosaic
swa_mop <- raster(paste(Region1, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
swa_mop[swa_mop > 0] <- NA
swa_mop[swa_mop == 0] <- 8
aus_mop <- raster(paste(Region2, '/MOP_Results/Set46/MOP_10%_current.tif', sep = ''))
aus_mop[aus_mop > 0] <- NA
aus_mop[aus_mop == 0] <- 8
mop_mosaic <- mosaic(swa_mop, aus_mop, fun = sum)
mop_mosaic[mop_mosaic < 8] <- NA
mop_mosaic[mop_mosaic >= 8] <- 8 # areas of strict extrapolation

# SWA mosaic
swa_mosaic <- mosaic(bg_binary, swa_bin_10, mop_mosaic, fun = sum)
swa_mosaic[swa_mosaic > 3] <- 8

# AUS mosaic
aus_mosaic <- mosaic(bg_binary, aus_bin_10, mop_mosaic, fun = sum)
aus_mosaic[aus_mosaic > 4] <- 8

# Final mosaic
mod_mosaic <- mosaic(swa_mosaic, aus_mosaic, fun = sum)
mod_mosaic[mod_mosaic < 3] <- 1 # unsuitable areas
mod_mosaic[mod_mosaic == 4] <- 2 # swa MTP
mod_mosaic[mod_mosaic == 5] <- 2 # aus MTP
mod_mosaic[mod_mosaic == 7] <- 2 # overlapping MTPs
mod_mosaic[mod_mosaic > 7] <- 1 # strict extrapolation areas
mod_mosaic <- crop(mod_mosaic, coast)

# Load seamount/knoll polygon areas as suitable stepping-stones-like dispersal corridors
# Data downloaded from Yesson et al. (2011)
Seamounts <- readOGR(dsn = 'D:/D/ENV SDM/Seamounts/01_Data/SeamountsBaseArea', layer = 'SeamountsBaseArea')
Seamounts <- Seamounts[which(as.numeric(Seamounts$DEPTH) > -1001), ] # ecologically meaningful
Knolls <- readOGR(dsn = 'D:/D/ENV SDM/Seamounts/01_Data/KnollsBaseArea', layer = 'KnollsBaseArea')
Knolls <- Knolls[which(-as.numeric(Knolls$DEPTH) > -1001), ] # ecologically meaningful
seam_knoll <- rbind(Seamounts, Knolls)
seam_knoll <- crop(seam_knoll, mod_mosaic)

df <- data.frame(coordinates(mod_mosaic), as.data.frame(mod_mosaic))

ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = as.factor(layer))) + 
  geom_polygon(data = seam_knoll, aes(x = long, y = lat, group = group), col = 'darkorange', 
               fill = 'darkorange', size = 0.15) +
  scale_fill_manual(values = c('grey95', '#440154FF')) +  coord_equal(expand = 0) +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), col = 'grey80', fill = 'white', size = 0.1) +
  theme(panel.background = element_rect(fill = 'transparent'), panel.grid = element_blank(), 
        legend.position = 'none', axis.text = element_text(size = 8),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
  scale_y_continuous(name = NULL, breaks = c(50, 0, -50), labels = c('50ºN', '0º', '50ºS')) + 
  scale_x_continuous(name = NULL, breaks = c(100, 0, -100), labels = c('100ºW', '0º', '100ºE')) 
ggsave('Figure S10.tiff', dpi = 900, width = 20, height = 8.5, units = 'cm')


# ----------------------------------- END -------------------------------------
