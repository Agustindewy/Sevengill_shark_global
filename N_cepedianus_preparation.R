
# De Wysiecki et al. - World ENM projection for the broadnose sevengill shark (Notorynchus cepedianus)

# Preparation data for modelling

library(maps)
library(rgdal)
library(dplyr)
library(spThin)
library(raster)
library(rgeos)
library(sdmpredictors)
library(ggplot2)
library(combinat)

setwd('SET YOUR WORKING DIRECTORY')

# Projection
crs <- CRS('+init=epsg:4326')

# Species
Species <- 'Notorynchus cepedianus'

#----------------------------------------- Predictors ------------------------------------------------

# Bio-Oracle
list_layers(c('Bio-ORACLE'))[, 2:4] # explore layers

# Check large-scale correlation
layers_correlation(c('BO2_tempmean_bdmean', 'BO_sstmean', 'BO2_ppmean_bdmean', 'BO_damean', 'BO_chlomean',
                     'BO2_salinitymean_bdmean', 'BO_bathymean'))
# BO_chlomean vs BO_damean corr > 0.8, we dropped BO_chlomean

bio <- load_layers(c('BO2_tempmean_bdmean', 'BO_sstmean', 'BO2_ppmean_bdmean', 'BO_damean', 
                     'BO2_salinitymean_bdmean', 'BO_bathymean'))

# Distance to the coast from the Global Self-consistent, Hierarchical, High-resolution Geography Database
# Dowloaded from https://www.soest.hawaii.edu/pwessel/gshhg/
dis <- raster('DOWNLOAD AND READ THE RASTER.tif') #read tif
dis <- aggregate(dis, fact = 2, fun = mean) # reduce resolution to match bio-oracle predictors

# Slope (as 14-9-2021 MARSPEC web site not working, so create slope from bathymetry raster)
slo <- terrain(bio[['BO_bathymean']], opt = 'slope', unit = 'degrees', neighbors = 4)
slo <- projectRaster(slo, bio)

# Stack final predictors 
env <- stack(bio, dis, slo)

# Check large-scale correlation again now we added additional predictors
layerStats(env, 'pearson', na.rm = T) # all good

# Check if environmental values in each variable make sense
env # bathymetry has some positive values, NA them
env[[6]][env[[6]] > -1] <- NA

writeRaster(env, filename = 'predictors.tif') #write as tif

#------------------------------------------ G space --------------------------------------------

# Read predictors and name them
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope')
names(env) <- var_names

# Discard bathymetry and slope because discrepancies at shelf margin (see analysis for Figure S3 and S4) 
var_set <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity', 
             'Distance_to_coast')

#------------------------------------ Southwest Atlantic (SWA) ------------------------------------

# IMPORTANT!! occurrence data for SWA is incomplete as some data has confidenciality agreements in place 
# and could not be released, so unfortunately you won't get the same results

Region <- 'SWA'

# Seed
set.seed(111)

# Read occurrence data
dat <- read.csv('Appendix C.csv', header = T)
dat <- subset(dat, Region == 'Southwest Atlantic')

# Calibration area (M space)
coord <- cbind(dat$Longitude, dat$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point
env.M_swa <- crop(env, occ.buff) # crop raster stack to buffer polygon
env.M_swa <- mask(env.M_swa, occ.buff) # mask raster stack to buffer polygon

# Erase Pacific areas
lon1 <- c(-84,-69.4,-69.4,-84,-84) 
lat1 <- c(-32,-32,-62,-62,-32)
pol1 <- cbind(lon1, lat1)
Pol1 <- Polygon(pol1, hole = F)
list1 <- list(Pol1)
poly1 <- Polygons(list1, ID = 'p1')
spp1 <- SpatialPolygons(list(poly1), proj4string = crs)
env.M_swa <- mask(env.M_swa, spp1, inverse = T)
Ext1 <- extent(-69.4,-39,-61.4,-22.8)
env.M_swa <- crop(env.M_swa, Ext1)
env.M_swa <- stack(env.M_swa)

# Occurrence subsets for traning, testing and independent evaluation
swa_train <- subset(dat, Type == 'calibration')

# 50 km is the minimum distance to remove clusters (see script for Figure S1)
swa_train_thin <- thin(swa_train, lat.col = 'Latitude', long.col = 'Longitude', spec.col = 'Species',
                       thin.par = 50, reps = 1000, locs.thinned.list.return = T, write.files = F, write.log.file = F)
swa_train <- data.frame(species = Species, longitude = round(swa_train_thin[[1]]$Longitude, 2), 
                        latitude = round(swa_train_thin[[1]]$Latitude, 2))
swa_train$Type <- 'calibration'

# Independent data
swa_indep <- subset(dat, Type == 'independent')
swa_indep <- data.frame(species = Species, longitude = swa_indep$Longitude, latitude = swa_indep$Latitude, Type = 'independent')

# Final data for modelling
swa <- rbind(swa_train, swa_indep)

# Save any points to nearest pixel with data that may be falling outside a layer of the M space raster 
library(rSDM) #package in GitHub (Pakillo/rSDM), install and call 'rSDM'
for(i in 1:length(env.M_swa@layers)){
  spp <- SpatialPoints(swa[, c('longitude', 'latitude')], crs)
  spp_corrected <- points2nearestcell(locs = spp, ras = env.M_swa, layer = i, move = T) 
  swa$longitude <- round(spp_corrected@coords[, 1], 2)  # replace coordinates including those new if any
  swa$latitude <- round(spp_corrected@coords[, 2], 2)
} # plots will appear if any correction applies

# Create csv files for training and testing (75% and 25%)
dat <- subset(swa, Type == 'calibration')
dat$check <- paste(dat[, 'longitude'], dat[, 'latitude'], sep = '_')
train <- dat[sample(nrow(dat), round((length(dat[, 1]) / 10 * 7.5))), ]
test <- dat[!dat[, 5] %in% train[, 5], ]
dat$check <- NULL; train$check <- NULL; test$check <- NULL
dir.create(Region)
write.csv(dat[, 1:3], paste(Region, '/N_cepedianus_joint.csv', sep = ''), row.names = F)
write.csv(train[, 1:3], paste(Region, '/N_cepedianus_train.csv', sep = ''), row.names = F)
write.csv(test[, 1:3], paste(Region, '/N_cepedianus_test.csv', sep = ''), row.names = F)

# Create csv files of independent data
dat <- subset(swa, Type == 'independent')
write.csv(dat[, 1:3], paste(Region, '/N_cepedianus_ind.csv', sep = ''), row.names = F)

# Check correlation in calibration areas
layerStats(env.M_swa[[var_set]], 'pearson', na.rm = T) # all good (0.8 threshold)
var_set_swa <- var_set

# Variable selection
var_subsets_swa <- unlist(lapply(1:length(var_set_swa), combinat::combn, x = var_set_swa, simplify = F), recursive = F)

# Subset list to those with at least Surface_temperature, Distance_to_coast and Kd490 variables
var_subsets_swa <- var_subsets_swa[which(sapply(var_subsets_swa, function(x) is.element(c('Surface_temperature'), x)))]
var_subsets_swa <- var_subsets_swa[which(sapply(var_subsets_swa, function(x) is.element(c('Distance_to_coast'), x)))]
var_subsets_swa <- var_subsets_swa[which(sapply(var_subsets_swa, function(x) is.element(c('Kd490'), x)))]

# M variables for variable selection with 'kuenm' package
dir.create(paste(Region, '/Var_selection', sep = ''))
dir.create(paste(Region, '/Var_selection/M_variables', sep = ''))
for(i in 1:length(var_subsets_swa)){
  dir.create(paste(Region, '/Var_selection/M_variables/Set', i, sep = ''))
  writeRaster(env.M_swa[[var_subsets_swa[[i]]]], filename = paste(Region, '/Var_selection/M_variables/Set', i,'/env.asc', sep = ''),
              format = 'ascii', bylayer = T, suffix = var_subsets_swa[[i]])
}


#---------------------------------------- Australia (AUS) -----------------------------------------

# Full data available, full reproducibility

Region <- 'AUS'

# Seed
set.seed(111)

# Read occurrence data
dat <- read.csv('Appendix C.csv', header = T)
dat <- subset(dat, Region == 'Australia')

# Calibration area (M space)
coord <- cbind(dat$Longitude, dat$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point
env.M_aus <- crop(env, occ.buff) # crop raster stack to buffer polygon
env.M_aus <- mask(env.M_aus, occ.buff) # mask raster stack to buffer polygon
env.M_aus <- stack(env.M_aus)

# Occurrence subsets for traning, testing and independent evaluation
swa_train <- subset(dat, Type == 'calibration')

# 50 km is the minimum distance to remove clusters (see script for Figure S2)
swa_train_thin <- thin(swa_train, lat.col = 'Latitude', long.col = 'Longitude', spec.col = 'Species',
                       thin.par = 50, reps = 1000, locs.thinned.list.return = T, write.files = F, write.log.file = F)
swa_train <- data.frame(species = Species, longitude = round(swa_train_thin[[1]]$Longitude, 2), 
                        latitude = round(swa_train_thin[[1]]$Latitude, 2))
swa_train$Type <- 'calibration'

# Independent data
swa_indep <- subset(dat, Type == 'independent')
swa_indep <- data.frame(species = Species, longitude = swa_indep$Longitude, latitude = swa_indep$Latitude, Type = 'independent')

# Final data for modelling
swa <- rbind(swa_train, swa_indep)

# Save any points to nearest pixel with data that may be falling outside a layer of the M space raster 
library(rSDM) #package in GitHub (Pakillo/rSDM), install and call 'rSDM'
for(i in 1:length(env.M_aus@layers)){
  spp <- SpatialPoints(swa[, c('longitude', 'latitude')], crs)
  spp_corrected <- points2nearestcell(locs = spp, ras = env.M_aus, layer = i, move = T) 
  swa$longitude <- round(spp_corrected@coords[, 1], 2)  # replace coordinates including those new if any
  swa$latitude <- round(spp_corrected@coords[, 2], 2)
} # plots will appear if any correction applies

# Create csv files for training and testing (75% and 25%)
dat <- subset(swa, Type == 'calibration')
dat$check <- paste(dat[, 'longitude'], dat[, 'latitude'], sep = '_')
train <- dat[sample(nrow(dat), round((length(dat[, 1]) / 10 * 7.5))), ]
test <- dat[!dat[, 5] %in% train[, 5], ]
dat$check <- NULL; train$check <- NULL; test$check <- NULL
dir.create(Region)
write.csv(dat[, 1:3], paste(Region, '/N_cepedianus_joint.csv', sep = ''), row.names = F)
write.csv(train[, 1:3], paste(Region, '/N_cepedianus_train.csv', sep = ''), row.names = F)
write.csv(test[, 1:3], paste(Region, '/N_cepedianus_test.csv', sep = ''), row.names = F)

# Create csv files of independent data
dat <- subset(swa, Type == 'independent')
write.csv(dat[, 1:3], paste(Region, '/N_cepedianus_ind.csv', sep = ''), row.names = F)

# Check correlation in calibration areas
layerStats(env.M_aus[[var_set]], 'pearson', na.rm = T) # Salinity highly correlated with Temperature, drop salinity
var_set_aus <- var_set[! var_set == 'Salinity']

# Variable selection
var_subsets_aus <- unlist(lapply(1:length(var_set_aus), combinat::combn, x = var_set_aus, simplify = F), recursive = F)

# Subset list to those with at least Surface_temperature and Distance_to_coast variables
var_subsets_aus <- var_subsets_aus[which(sapply(var_subsets_aus, function(x) is.element(c('Surface_temperature'), x)))]
var_subsets_aus <- var_subsets_aus[which(sapply(var_subsets_aus, function(x) is.element(c('Distance_to_coast'), x)))]
var_subsets_aus <- var_subsets_aus[which(sapply(var_subsets_aus, function(x) is.element(c('Kd490'), x)))]

# M variables for variable selection with 'kuenm' package
dir.create(paste(Region, '/Var_selection', sep = ''))
dir.create(paste(Region, '/Var_selection/M_variables', sep = ''))
for(i in 1:length(var_subsets_aus)){
  dir.create(paste(Region, '/Var_selection/M_variables/Set', i, sep = ''))
  writeRaster(env.M_aus[[var_subsets_aus[[i]]]], filename = paste(Region, '/Var_selection/M_variables/Set', i,'/env.asc', sep = ''),
              format = 'ascii', bylayer = T, suffix = var_subsets_aus[[i]])
}

#----------------------------- M and G variables for final modelling ---------------------------------

# IMPORTANT!! run this last chunk once you determine best set of variables for both regions 
# Go to modelling - var selection first then come back here

var_subsets <- unlist(lapply(1:length(var_set), combinat::combn, x = var_set, simplify = F), recursive = F)
best_set_swa <- 2 # must be an integer
best_set_aus <- 2 # must be an integer
best_set0 <- unique(c(var_subsets_swa[[best_set_swa]], var_subsets_aus[[best_set_aus]]))

SameElements <- function(a, b) return(identical(sort(a), sort(b)))
for(i in 1:length(var_subsets)){
  resul <- SameElements(var_subsets[[i]], best_set0)
  if(resul == T) {best_set <- i}
}

Region <- 'SWA' 
dir.create(paste(Region, '/M_variables', sep = ''))
dir.create(paste(Region, '/M_variables/Set', best_set, sep = ''))
writeRaster(env.M_swa[[var_subsets[[best_set]]]], filename = paste(Region, '/M_variables/Set', best_set,'/env.asc', sep = ''),
            format = 'ascii', bylayer = T, suffix = var_subsets[[best_set]])
dir.create(paste(Region, '/G_variables', sep = ''))
dir.create(paste(Region, '/G_variables/Set', best_set, sep = ''))
dir.create(paste(Region, '/G_variables/Set', best_set,'/current', sep = ''))
writeRaster(env[[var_subsets[[best_set]]]], filename = paste(Region, '/G_variables/Set', best_set,'/current/env.asc', sep = ''),
            format = 'ascii', bylayer = T, suffix = var_subsets[[best_set]])

Region <- 'AUS' 
dir.create(paste(Region, '/M_variables', sep = ''))
dir.create(paste(Region, '/M_variables/Set', best_set, sep = ''))
writeRaster(env.M_aus[[var_subsets[[best_set]]]], filename = paste(Region, '/M_variables/Set', best_set,'/env.asc', sep = ''),
            format = 'ascii', bylayer = T, suffix = var_subsets[[best_set]])
dir.create(paste(Region, '/G_variables', sep = ''))
dir.create(paste(Region, '/G_variables/Set', best_set, sep = ''))
dir.create(paste(Region, '/G_variables/Set', best_set,'/current', sep = ''))
writeRaster(env[[var_subsets[[best_set]]]], filename = paste(Region, '/G_variables/Set', best_set,'/current/env.asc', sep = ''),
            format = 'ascii', bylayer = T, suffix = var_subsets[[best_set]])


# -------------------------------------------- END ------------------------------------------------
