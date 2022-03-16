
# De Wysiecki et al. - Global ENM projection for the broadnose sevengill shark (Notorynchus cepedianus)

# Preparation data for modelling
library(maps)
library(rgdal)
library(dplyr)
library(spThin)
library(rgeos)
library(sdmpredictors)
library(ggplot2)
library(combinat)
library(ENMeval) # find it on GitHub (jamiemkass/ENMeval)
library(raster)
library(sqldf)
library(blockCV)
library(dismo)
library(kuenm) # find it on GitHub (marlonecobos/kuenm)
library(ntbox)
library(rSDM) # find it on GitHub (Pakillo/rSDM)

setwd('SET YOUR WORKING DIRECTORY')

# Projection
crs <- CRS('+init=epsg:4326')

# Species
Species <- 'Notorynchus cepedianus'

# Function to include partial ROC statistic in ENM evaluation
proc <- function(vars) {
  proc <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
  out <- data.frame(proc_auc_ratio = proc$pROC_summary[1], 
                    proc_pval = proc$pROC_summary[2], row.names = NULL)
  return(out)
} 

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

writeRaster(env, filename = 'predictors.tif') # write as tif

#------------------------------------------ G space --------------------------------------------

# Read predictors and name them
env <- stack('predictors.tif')
var_names <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity',
               'Bathymetry', 'Distance_to_coast', 'Slope')
names(env) <- var_names

# Discard bathymetry and slope because discrepancies at shelf margin (see analysis for Figure S3 and S4) 
var_set <- c('Temperature', 'Surface_temperature', 'Primary_productivity', 'Kd490', 'Salinity', 
             'Distance_to_coast')

# G space
env.G <- env[[var_set]]

#----------------------------------------- Global fit ---------------------------------------

Fit <- 'GLOBAL'
dir.create(Fit)

set.seed(111) # seed

#--- Occurrence data processing -------------------------------------------------------------

# Read occurrence data
dat <- read.csv('Appendix C.csv') # note some confidential data is not available

# Get rid of duplicates
dups <- duplicated(dat[c('Longitude', 'Latitude')])
occ_cal <- dat[!dups, ]

# Spatial thinning: 50 km is the minimum distance to remove spatial clusters (see script for Figure S1)
train_thin <- thin(occ_cal, lat.col = 'Latitude', long.col = 'Longitude', spec.col = 'Species',
                       thin.par = 50, reps = 10, locs.thinned.list.return = T, write.files = F, write.log.file = F)
occ_cal <- data.frame(Species = Species, Longitude = round(train_thin[[1]]$Longitude, 2), 
                        Latitude = round(train_thin[[1]]$Latitude, 2))
occ_cal$Type <- 'calibration'

# Initial M space for environmental filtering
# some data in New Zealand need buffers across 180º meridian so first recenter points to 110ºW longitude to avoid this
occ_cal$Longitude2 <- occ_cal$Longitude + 110
occ_cal$Longitude2 <- ifelse(occ_cal$Longitude2 > 180, (occ_cal$Longitude2 - 180) - 180, occ_cal$Longitude2)
coord <- cbind(occ_cal$Longitude2, occ_cal$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point
# now recenter raster stack to 110ºW longitude
env1 <- crop(env.G, extent(-180, 70, -90, 90))
env2 <- crop(env.G, extent(70, 180, -90, 90))   
extent(env1) <- c(-70, 180, -90, 90)
extent(env2) <- c(-180, -70, -90, 90)
env0 <- merge(env1, env2)
names(env0) <- var_set
env.M <- crop(env0, occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon
env.M <- stack(env.M)

# More occurrences processing
# Save any points to nearest pixel with data that may be falling outside a layer of the M space raster 
for(i in 1:length(env.M@layers)){
  spp <- SpatialPoints(occ_cal[, c('Longitude2', 'Latitude')], crs)
  spp_corrected <- points2nearestcell(locs = spp, ras = env.M, layer = i, move = T) 
  occ_cal$Longitude2 <- round(spp_corrected@coords[, 1], 3)  # replace coordinates including those new if any
  occ_cal$Latitude <- round(spp_corrected@coords[, 2], 3)
} # plots will appear if any correction applies

# Get rid of duplicates again in case some corrected points overlap
dups <- duplicated(occ_cal[c('Longitude2', 'Latitude')])
occ_cal <- occ_cal[!dups, ]

# Environmental filtering
source('envSample.R')
coords <- occ_cal[, c('Longitude2', 'Latitude')]
env.data <- extract(env.M, coords)
env.data <- as.data.frame(env.data)
coords <- envSample(coords, filters = list(env.data$Surface_temperature, env.data$Distance_to_coast),
                    res = list(0.5, 1), do.plot = T)
occ_cal <- merge(occ_cal, coords, by.x = c('Longitude2', 'Latitude'), by.y = c('lon', 'lat'))

# Detecting environmental outliers in species occurrences (Cobos et al., 2018)
# Getting data from the variables
variables_values <- na.omit(values(env.M)) # for the region of interest
occ_variables <- na.omit(cbind(as.numeric(row.names(occ_cal)), extract(env.M, occ_cal[, c('Longitude2', 'Latitude')])))
occ_variables <- as.data.frame(occ_variables)
if(dim(variables_values)[1] > 10000) { # sample of 10000 values if more pixels exist
  variables_values <- variables_values[sample(1:dim(variables_values)[1], 10000), ]
}

# Plot for searching for potential environmental outliers
# plot all combinations, then scroll plots and manually detect outliers and write them down in code below
for(i in 1:dim(env.M)[3]){
  for(j in 1:dim(env.M)[3]){
    par(mfcol = c(1, 2), mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
    plot(variables_values[, i], variables_values[, j], col = 'grey65', pch = 1,
         xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    points(occ_variables[, i + 1], occ_variables[, j + 1], col = 'blue', pch = 19, cex = 1.5)
    legend('bottomright', legend = c('Region of interest', 'Occurrences'),
           pch = c(1, 19), col = c('gray65', 'black'), bty = 'n')
    plot(variables_values[, i], variables_values[, j], col = 'gray65',
         pch = 1, xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    text(occ_variables[, i + 1], occ_variables[, j + 1], occ_variables[, 1], cex = 1, col = 'blue')
    legend('bottomright', legend = 'Occurrence ID', bty = 'n')
  }
}

# Remove outliers from calibration records
occ_variables <- subset(occ_variables, !V1 %in% c()) # PUT THE INTEGER POINT CODES TO REMOVE WITHIN THE c()
occ_cal <- occ_cal[which(row.names(occ_cal) %in% as.character(occ_variables$V1)), ]

# Final calibration points
occ_cal <- rbind(occ_cal[, c(3, 4, 2, 5)])
dir.create(Fit)
write.csv(occ_cal, paste(Fit, '/', Fit, '_calibration_points.csv', sep = ''), row.names = F)


#--- Calibration areas -----------------------------------------------------------------

occ_cal <- read.csv(paste(Fit, '/', Fit, '_calibration_points.csv', sep = ''))
  
# Final M space (calibration areas)
# some data in New Zealand need buffers across 180º meridian so first recenter points to 110ºW longitude to avoid this
occ_cal$Longitude2 <- occ_cal$Longitude + 110
occ_cal$Longitude2 <- ifelse(occ_cal$Longitude2 > 180, (occ_cal$Longitude2 - 180) - 180, occ_cal$Longitude2)
coord <- cbind(occ_cal$Longitude2, occ_cal$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point

# now recenter raster stack to 110ºW longitude
env1 <- crop(env.G, extent(-180, 70, -90, 90))
env2 <- crop(env.G, extent(70, 180, -90, 90))   
extent(env1) <- c(-70, 180, -90, 90)
extent(env2) <- c(-180, -70, -90, 90)
env0 <- merge(env1, env2)
names(env0) <- var_set
env.M <- crop(env0, occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon

# recenter raster stack back to 0º longitude
env.M <- extend(env.M, extent(-180, 180, -61.33333, 55.58333)) # extend back to -180º, 180º
env.M1 <- crop(env.M, extent(-180, -70, -61.33333, 55.58333))
env.M2 <- crop(env.M, extent(-70, 180, -61.33333, 55.58333))   
extent(env.M1) <- c(70, 180, -61.33333, 55.58333)
extent(env.M2) <- c(-180, 70, -61.33333, 55.58333)
env.M <- merge(env.M2, env.M1)
names(env.M) <- var_set
env.M <- stack(env.M)

# Check correlation in M space
layerStats(env.M, 'pearson', na.rm = T) # all < 0.8, all good

writeRaster(env.M, filename = paste(Fit, '/', Fit, '_calibration_areas.tiff', sep = '')) # write as tif


#--- Variable selection -------------------------------------------------------------

# Based on kuenm package from Cobos et al. (2019)

# Calibration points
occ_cal <- read.csv(paste(Fit, '/', Fit, '_calibration_points.csv', sep = ''))
occ_cal <- occ_cal[, -4]
names(occ_cal) <- c('species', 'longitude', 'latitude')

# Create training and testing sets (75% and 25%)
occ_cal$check <- paste(occ_cal[, 'longitude'], occ_cal[, 'latitude'], sep = '_')
train <- occ_cal[sample(nrow(occ_cal), round((length(occ_cal[, 1]) / 10 * 7.5))), ]
test <- occ_cal[!occ_cal[, 4] %in% train[, 4], ]
dir.create(paste(Fit, '/Var_selection', sep = ''))
write.csv(occ_cal[, -4], paste(Fit, '/Var_selection/', Fit, '_occ_joint.csv', sep = ''), row.names = F)
write.csv(train[, -4], paste(Fit, '/Var_selection/', Fit, '_occ_train.csv', sep = ''), row.names = F)
write.csv(test[, -4], paste(Fit, '/Var_selection/', Fit, '_occ_test.csv', sep = ''), row.names = F)

# Calibration areas
env.M <- stack(paste(Fit, '/', Fit, '_calibration_areas.tiff', sep = ''))
names(env.M) <- var_set

# Variable combinations
var_subsets <- unlist(lapply(1:length(var_set), combinat::combn, x = var_set, simplify = F), recursive = F)

# Subset combinations to those with at least Surface_temperature, Distance_to_coast and Kd490 variables
var_subsets <- var_subsets[which(sapply(var_subsets, function(x) is.element(c('Surface_temperature'), x)))]
var_subsets <- var_subsets[which(sapply(var_subsets, function(x) is.element(c('Distance_to_coast'), x)))]
var_subsets <- var_subsets[which(sapply(var_subsets, function(x) is.element(c('Kd490'), x)))]

# Create sets of variables for selection with 'kuenm' package
dir.create(paste(Fit, '/Var_selection/M_variables', sep = ''))
for(i in 1:length(var_subsets)){
  dir.create(paste(Fit, '/Var_selection/M_variables/Set', i, sep = ''))
  writeRaster(env.M[[var_subsets[[i]]]], filename = paste(Fit, '/Var_selection/M_variables/Set', i,'/env.asc', sep = ''),
              format = 'ascii', bylayer = T, suffix = var_subsets[[i]])
}

# Candidate models
occ_joint <- paste(Fit, '/Var_selection/', Fit, '_occ_joint.csv', sep = '')
occ_tra <- paste(Fit, '/Var_selection/', Fit, '_occ_train.csv', sep = '')
M_var_dir <- paste(Fit, '/Var_selection/M_variables', sep = '')
batch_cal <- paste(Fit, '/Var_selection/Candidate_Models', sep = '')
out_dir <- paste(Fit, '/Var_selection/Candidate_Models', sep = '')
reg_mult <- c(0.1, 0.5, 1, 2.5, 5) # coarse set of values
f_clas <- c('lq', 'lp', 'lqp')
args <- NULL
maxent_path <- 'PATH TO MAXENT JAR FILE' 
wait <- TRUE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal, 
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args, maxent.path = maxent_path, 
          wait = wait, run = run)

# Evaluation and selection of best models
occ_test <- paste(Fit, '/Var_selection/', Fit, '_occ_test.csv', sep = '')
out_eval <- paste(Fit, '/Var_selection/Calibration_Results', sep = '')
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- FALSE
selection <- 'OR_AICc'
paral_proc <- FALSE 

kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, 
            batch = batch_cal, out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, 
            iterations = iterations, kept = kept, selection = selection, parallel.proc = paral_proc)

# Best set number 8 = all predictors


#--- Model calibration, evaluation and selection -----------------------------------------------

# Calibration points
occ_cal <- read.csv(paste(Fit, '/', Fit, '_calibration_points.csv', sep = ''))
occ_cal <- occ_cal[, c('Longitude', 'Latitude')]
colnames(occ_cal) <- c('longitude', 'latitude')

# Calibration areas
env.M <- stack(paste(Fit, '/', Fit, '_calibration_areas.tiff', sep = ''))
names(env.M) <- var_set

# Background points
bg_points <- as.data.frame(randomPoints(env.M[[1]], n = 10000))
colnames(bg_points) <- c('longitude', 'latitude')

# Presence-background points
all_pts <- rbind(occ_cal, bg_points)
all_pts <- SpatialPoints(all_pts, proj4string = crs) # transform to spatial points

# Check the effective range of spatial autocorrelation
eff.range <- spatialAutoRange(rasterLayer = env.M, sampleNumber = 10000, doParallel = T, showPlots = T)
median_range <- 830000 # actual range of 833717 round up to 830000

# Spatial block design based on 'blockCV' by Valavi et al. (2019)
sp_block <- spatialBlock(speciesData = all_pts, rasterLayer = env.M[[1]], theRange = median_range, k = 5, selection = 'random')

# Get block information for modelling
user.grp <- list(occs.grp = sp_block$foldID[1:nrow(occ_cal)],
                 bg.grp = sp_block$foldID[(nrow(occ_cal) + 1):length(sp_block$foldID)])

# Running and evaluating candidate models
mod_global <- ENMevaluate(occs = occ_cal, envs = env.G, bg = bg_points, user.grp = user.grp, 
                          algorithm = 'maxnet', partitions = 'user', user.eval = proc, doClamp = F, 
                          tune.args = list(fc = c('LQ', 'LP', 'LQP'), rm = c(seq(0.1, 2, 0.1), 2.5, 3, 4, 5, 7, 10)))

# Model selection
res <- eval.results(mod_global) # overall results
opt.seq <- res %>% 
  filter(proc_pval.avg == min(proc_pval.avg)) %>%
  filter(or.10p.avg == min(or.10p.avg)) # 3-step criteria in Cobos et al. (2019)

# Optimal model settings
mod.seq <- eval.models(mod_global)[[opt.seq$tune.args]]

# Model predictions
pred.seq <- eval.predictions(mod_global)[[opt.seq$tune.args]]
# pred.seq <- mosaic(pred.seq[[1]], pred.seq[[2]], fun = median) # in case two or more model are selected as best models
writeRaster(pred.seq, filename = paste(Fit, '/', Fit, '_predictions.tiff', sep = '')) # write as tif


#------------------------------------- Regional merged fit ----------------------------------

Fit <- 'REGIONAL_MERGED' # SWA and AUS as a whole, i.e. well-represented regions
dir.create(Fit)

set.seed(111) # seed

#--- Occurrence data processing -------------------------------------------------------------

# Read occurrence data
dat <- read.csv('Appendix C.csv') # note some confidential data is not available
dat <- subset(dat, Region %in% c('Southwest Atlantic', 'Australia'))

# Get rid of duplicates
dups <- duplicated(dat[c('Longitude', 'Latitude')])
occ_cal <- dat[!dups, ]

# Spatial thinning: 50 km is the minimum distance to remove spatial clusters (see script for Figure S1)
train_thin <- thin(occ_cal, lat.col = 'Latitude', long.col = 'Longitude', spec.col = 'Species',
                   thin.par = 50, reps = 10, locs.thinned.list.return = T, write.files = F, write.log.file = F)
occ_cal <- data.frame(Species = Species, Longitude = round(train_thin[[1]]$Longitude, 2), 
                      Latitude = round(train_thin[[1]]$Latitude, 2))
occ_cal$Type <- 'calibration'

# Initial M space for environmental filtering
coord <- cbind(occ_cal$Longitude, occ_cal$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point
env.M <- crop(env.G, occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon
env.M <- stack(env.M)

# Erase Pacific areas from SWA as they inclusion is an artifact of creating the spatial buffers
lon1 <- c(-84, -69.4, -69.4, -84, -84) 
lat1 <- c(-32, -32, -62, -62, -32)
pol1 <- cbind(lon1, lat1)
Pol1 <- Polygon(pol1, hole = F)
list1 <- list(Pol1)
poly1 <- Polygons(list1, ID = 'p1')
spp1 <- SpatialPolygons(list(poly1), proj4string = crs)
env.M <- mask(env.M, spp1, inverse = T)
Ext1 <- extent(-69.4, 162.5, -61.4, -22.8)
env.M <- crop(env.M, Ext1)
env.M <- stack(env.M)

# More occurrences processing
# Save any points to nearest pixel with data that may be falling outside a layer of the M space raster 
for(i in 1:length(env.M@layers)){
  spp <- SpatialPoints(occ_cal[, c('Longitude', 'Latitude')], crs)
  spp_corrected <- points2nearestcell(locs = spp, ras = env.M, layer = i, move = T) 
  occ_cal$Longitude <- round(spp_corrected@coords[, 1], 3)  # replace coordinates including those new if any
  occ_cal$Latitude <- round(spp_corrected@coords[, 2], 3)
} # plots will appear if any correction applies

# Get rid of duplicates again in case some corrected points overlap
dups <- duplicated(occ_cal[c('Longitude', 'Latitude')])
occ_cal <- occ_cal[!dups, ]

# Environmental filtering
source('envSample.R')
coords <- occ_cal[, c('Longitude', 'Latitude')]
env.data <- extract(env.M, coords)
env.data <- as.data.frame(env.data)
coords <- envSample(coords, filters = list(env.data$Surface_temperature, env.data$Distance_to_coast),
                    res = list(0.5, 1), do.plot = T)
occ_cal <- merge(occ_cal, coords, by.x = c('Longitude', 'Latitude'), by.y = c('lon', 'lat'))

# Detecting environmental outliers in species occurrences (Cobos et al., 2018)
# Getting data from the variables
variables_values <- na.omit(values(env.M)) # for the region of interest
occ_variables <- na.omit(cbind(as.numeric(row.names(occ_cal)), extract(env.M, occ_cal[, c('Longitude', 'Latitude')])))
occ_variables <- as.data.frame(occ_variables)
if(dim(variables_values)[1] > 10000) { # sample of 10000 values if more pixels exist
  variables_values <- variables_values[sample(1:dim(variables_values)[1], 10000), ]
}

# Plot for searching for potential environmental outliers
# plot all combinations, then scroll plots and manually detect outliers and write them down in code below
for(i in 1:dim(env.M)[3]){
  for(j in 1:dim(env.M)[3]){
    par(mfcol = c(1, 2), mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
    plot(variables_values[, i], variables_values[, j], col = 'grey65', pch = 1,
         xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    points(occ_variables[, i + 1], occ_variables[, j + 1], col = 'blue', pch = 19, cex = 1.5)
    legend('bottomright', legend = c('Region of interest', 'Occurrences'),
           pch = c(1, 19), col = c('gray65', 'black'), bty = 'n')
    plot(variables_values[, i], variables_values[, j], col = 'gray65',
         pch = 1, xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    text(occ_variables[, i + 1], occ_variables[, j + 1], occ_variables[, 1], cex = 1, col = 'blue')
    legend('bottomright', legend = 'Occurrence ID', bty = 'n')
  }
}

# Remove outliers from calibration records
occ_variables <- subset(occ_variables, !V1 %in% c()) # PUT THE INTEGER POINT CODES TO REMOVE WITHIN THE c() 
occ_cal <- occ_cal[which(row.names(occ_cal) %in% as.character(occ_variables$V1)), ]

# Final calibration points
occ_cal <- rbind(occ_cal[, c(3, 1, 2, 4)])
dir.create(Fit)
write.csv(occ_cal, paste(Fit, '/', Fit, '_calibration_points.csv', sep = ''), row.names = F)


#--- Calibration areas -----------------------------------------------------------------

occ_cal <- read.csv(paste(Fit, '/', Fit, '_calibration_points.csv', sep = ''))

# Final M space (calibration areas)
coord <- cbind(occ_cal$Longitude, occ_cal$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point
env.M <- crop(env.G, occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon

# Erase Pacific areas from SWA as they inclusion is an artifact of creating the spatial buffers
lon1 <- c(-84, -69.4, -69.4, -84, -84) 
lat1 <- c(-32, -32, -62, -62, -32)
pol1 <- cbind(lon1, lat1)
Pol1 <- Polygon(pol1, hole = F)
list1 <- list(Pol1)
poly1 <- Polygons(list1, ID = 'p1')
spp1 <- SpatialPolygons(list(poly1), proj4string = crs)
env.M <- mask(env.M, spp1, inverse = T)
Ext1 <- extent(-69.4, 162.5, -61.4, -22.8)
env.M <- crop(env.M, Ext1)
env.M <- stack(env.M)

# Check correlation in M space
layerStats(env.M, 'pearson', na.rm = T) # all < 0.8, all good

writeRaster(env.M, filename = paste(Fit, '/', Fit, '_calibration_areas.tiff', sep = '')) # write as tif


#--- Model calibration, evaluation and selection -----------------------------------------------

# Calibration points
occ_cal <- read.csv(paste(Fit, '/', Fit, '_calibration_points.csv', sep = ''))
occ_cal <- occ_cal[, c('Longitude', 'Latitude')]
colnames(occ_cal) <- c('longitude', 'latitude')

# Calibration areas
env.M <- stack(paste(Fit, '/', Fit, '_calibration_areas.tiff', sep = ''))
names(env.M) <- var_set

# Background points
bg_points <- as.data.frame(randomPoints(env.M[[1]], n = 10000))
colnames(bg_points) <- c('longitude', 'latitude')

# Presence-background points
all_pts <- rbind(occ_cal, bg_points)
all_pts <- SpatialPoints(all_pts, proj4string = crs) # transform to spatial points

# Check the effective range of spatial autocorrelation
eff.range <- spatialAutoRange(rasterLayer = env.M, sampleNumber = 10000, doParallel = T, showPlots = T)
median_range <- 1200000 # actual range of 1191747 round up to 1200000

# Spatial block design based on 'blockCV' by Valavi et al. (2019)
sp_block <- spatialBlock(speciesData = all_pts, rasterLayer = env.M[[1]], theRange = median_range, k = 5, selection = 'random')

# Get block information for modelling
user.grp <- list(occs.grp = sp_block$foldID[1:nrow(occ_cal)],
                 bg.grp = sp_block$foldID[(nrow(occ_cal) + 1):length(sp_block$foldID)])

# Running and evaluating candidate models
mod_regional_merged <- ENMevaluate(occs = occ_cal, envs = env.G, bg = bg_points, user.grp = user.grp, 
                                   algorithm = 'maxnet', partitions = 'user', user.eval = proc, doClamp = F, 
                                   tune.args = list(fc = c('LQ', 'LP', 'LQP'), rm = c(seq(0.1, 2, 0.1), 2.5, 3, 4, 5, 7, 10)))

# Model selection
res <- eval.results(mod_regional_merged) # overall results
opt.seq <- res %>% 
  filter(proc_pval.avg == min(proc_pval.avg)) %>%
  filter(or.10p.avg == min(or.10p.avg)) # 3-step criteria in Cobos et al. (2019)

# Optimal model settings
mod.seq <- eval.models(mod_regional_merged)[[opt.seq$tune.args]]

# Model predictions
pred.seq <- eval.predictions(mod_regional_merged)[[opt.seq$tune.args]]
# pred.seq <- mosaic(pred.seq[[1]], pred.seq[[2]], fun = median) # in case two or more model are selected as best models
writeRaster(pred.seq, filename = paste(Fit, '/', Fit, '_predictions.tiff', sep = '')) # write as tif


#----------------------------------- Regional independent fit ---------------------------------

Fit <- 'REGIONAL_INDEPENDENT'
dir.create(Fit)

#--- Southwest Atlantic (SWA) ---------------------------------------------------------------

Region <- 'SWA'

set.seed(111) # seed

#--- Occurrence data processing -------------------------------------------------------------

# Read occurrence data
dat <- read.csv('Appendix C.csv') # note some confidential data is not available
dat <- subset(dat, Region == 'Southwest Atlantic')

# Get rid of duplicates
dups <- duplicated(dat[c('Longitude', 'Latitude')])
occ_cal <- dat[!dups, ]

# Spatial thinning: 50 km is the minimum distance to remove spatial clusters (see script for Figure S1)
train_thin <- thin(occ_cal, lat.col = 'Latitude', long.col = 'Longitude', spec.col = 'Species',
                   thin.par = 50, reps = 10, locs.thinned.list.return = T, write.files = F, write.log.file = F)
occ_cal <- data.frame(Species = Species, Longitude = round(train_thin[[1]]$Longitude, 2), 
                      Latitude = round(train_thin[[1]]$Latitude, 2))
occ_cal$Type <- 'calibration'

# Initial M space for environmental filtering
coord <- cbind(occ_cal$Longitude, occ_cal$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point
env.M <- crop(env.G, occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon
env.M <- stack(env.M)

# Erase Pacific areas from SWA as they inclusion is an artifact of creating the spatial buffers
lon1 <- c(-84, -69.4, -69.4, -84, -84) 
lat1 <- c(-32, -32, -62, -62, -32)
pol1 <- cbind(lon1, lat1)
Pol1 <- Polygon(pol1, hole = F)
list1 <- list(Pol1)
poly1 <- Polygons(list1, ID = 'p1')
spp1 <- SpatialPolygons(list(poly1), proj4string = crs)
env.M <- mask(env.M, spp1, inverse = T)
Ext1 <- extent(-69.4, -39, -61.4, -22.8)
env.M <- crop(env.M, Ext1)
env.M <- stack(env.M)

# More occurrences processing
# Save any points to nearest pixel with data that may be falling outside a layer of the M space raster 
for(i in 1:length(env.M@layers)){
  spp <- SpatialPoints(occ_cal[, c('Longitude', 'Latitude')], crs)
  spp_corrected <- points2nearestcell(locs = spp, ras = env.M, layer = i, move = T) 
  occ_cal$Longitude <- round(spp_corrected@coords[, 1], 3)  # replace coordinates including those new if any
  occ_cal$Latitude <- round(spp_corrected@coords[, 2], 3)
} # plots will appear if any correction applies

# Get rid of duplicates again in case some corrected points overlap
dups <- duplicated(occ_cal[c('Longitude', 'Latitude')])
occ_cal <- occ_cal[!dups, ]

# Environmental filtering
source('envSample.R')
coords <- occ_cal[, c('Longitude', 'Latitude')]
env.data <- extract(env.M, coords)
env.data <- as.data.frame(env.data)
coords <- envSample(coords, filters = list(env.data$Surface_temperature, env.data$Distance_to_coast),
                    res = list(0.5, 1), do.plot = T)
occ_cal <- merge(occ_cal, coords, by.x = c('Longitude', 'Latitude'), by.y = c('lon', 'lat'))

# Detecting environmental outliers in species occurrences (Cobos et al., 2018)
# Getting data from the variables
variables_values <- na.omit(values(env.M)) # for the region of interest
occ_variables <- na.omit(cbind(as.numeric(row.names(occ_cal)), extract(env.M, occ_cal[, c('Longitude', 'Latitude')])))
occ_variables <- as.data.frame(occ_variables)
if(dim(variables_values)[1] > 10000) { # sample of 10000 values if more pixels exist
  variables_values <- variables_values[sample(1:dim(variables_values)[1], 10000), ]
}

# Plot for searching for potential environmental outliers
# plot all combinations, then scroll plots and manually detect outliers and write them down in code below
for(i in 1:dim(env.M)[3]){
  for(j in 1:dim(env.M)[3]){
    par(mfcol = c(1, 2), mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
    plot(variables_values[, i], variables_values[, j], col = 'grey65', pch = 1,
         xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    points(occ_variables[, i + 1], occ_variables[, j + 1], col = 'blue', pch = 19, cex = 1.5)
    legend('bottomright', legend = c('Region of interest', 'Occurrences'),
           pch = c(1, 19), col = c('gray65', 'black'), bty = 'n')
    plot(variables_values[, i], variables_values[, j], col = 'gray65',
         pch = 1, xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    text(occ_variables[, i + 1], occ_variables[, j + 1], occ_variables[, 1], cex = 1, col = 'blue')
    legend('bottomright', legend = 'Occurrence ID', bty = 'n')
  }
}

# Remove outliers from calibration records
occ_variables <- subset(occ_variables, !V1 %in% c(6)) # PUT THE INTEGER POINT CODES TO REMOVE WITHIN THE c() 
occ_cal <- occ_cal[which(row.names(occ_cal) %in% as.character(occ_variables$V1)), ]

# Final calibration points
occ_cal <- rbind(occ_cal[, c(3, 1, 2, 4)])
dir.create(Fit)
dir.create(paste(Fit, '/', Region, sep = ''))
write.csv(occ_cal, paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_points.csv', sep = ''), row.names = F)


#--- Calibration areas -----------------------------------------------------------------

occ_cal <- read.csv(paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_points.csv', sep = ''))

# Final M space (calibration areas)
coord <- cbind(occ_cal$Longitude, occ_cal$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point
env.M <- crop(env.G, occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon

# Erase Pacific areas from SWA as they inclusion is an artifact of creating the spatial buffers
lon1 <- c(-84, -69.4, -69.4, -84, -84) 
lat1 <- c(-32, -32, -62, -62, -32)
pol1 <- cbind(lon1, lat1)
Pol1 <- Polygon(pol1, hole = F)
list1 <- list(Pol1)
poly1 <- Polygons(list1, ID = 'p1')
spp1 <- SpatialPolygons(list(poly1), proj4string = crs)
env.M <- mask(env.M, spp1, inverse = T)
Ext1 <- extent(-69.4, -39, -61.4, -22.8)
env.M <- crop(env.M, Ext1)
env.M <- stack(env.M)

# Check correlation in M space
layerStats(env.M, 'pearson', na.rm = T) # all < 0.8, all good

writeRaster(env.M, filename = paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_areas.tiff', sep = '')) # write as tif


#--- Model calibration, evaluation and selection -----------------------------------------------

# Calibration points
occ_cal <- read.csv(paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_points.csv', sep = ''))
occ_cal <- occ_cal[, c('Longitude', 'Latitude')]
colnames(occ_cal) <- c('longitude', 'latitude')

# Calibration areas
env.M <- stack(paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_areas.tiff', sep = ''))
names(env.M) <- var_set

# Background points
bg_points <- randomPoints(env.M[[1]], n = 10000)
colnames(bg_points) <- c('longitude', 'latitude')

# Presence-background points
all_pts <- rbind(occ_cal, bg_points)
all_pts <- SpatialPoints(all_pts, proj4string = crs) # transform to spatial points

# Check the effective range of spatial autocorrelation
eff.range <- spatialAutoRange(rasterLayer = env.M, sampleNumber = 10000, doParallel = T, showPlots = T)
median_range <- 1400000 # actual range of 1379327 round up to 1400000

# Spatial block design based on 'blockCV' by Valavi et al. (2019)
sp_block <- spatialBlock(speciesData = all_pts, rasterLayer = env.M[[1]], theRange = median_range, k = 3, selection = 'random')

# Get block information for modelling
user.grp <- list(occs.grp = sp_block$foldID[1:nrow(occ_cal)],
                 bg.grp = sp_block$foldID[(nrow(occ_cal) + 1):length(sp_block$foldID)])

# Running and evaluating candidate models
mod_regional_swa <- ENMevaluate(occs = occ_cal, envs = env.G, bg = bg_points, user.grp = user.grp, 
                                algorithm = 'maxnet', partitions = 'user', user.eval = proc, doClamp = F, 
                                tune.args = list(fc = c('LQ', 'LP', 'LQP'), rm = c(seq(0.1, 2, 0.1), 2.5, 3, 4, 5, 7, 10)))

# Model selection
res <- eval.results(mod_regional_swa) # overall results
opt.seq <- res %>% 
  filter(proc_pval.avg == min(proc_pval.avg)) %>%
  filter(or.10p.avg == min(or.10p.avg)) # 3-step criteria in Cobos et al. (2019)

# Optimal model settings
mod.seq <- eval.models(mod_regional_swa)[[opt.seq$tune.args]]

# Model predictions
pred.seq <- eval.predictions(mod_regional_swa)[[opt.seq$tune.args]]
# pred.seq <- mosaic(pred.seq[[1]], pred.seq[[2]], fun = median) # in case two or more model are selected as best models
writeRaster(pred.seq, filename = paste(Fit, '/', Region, '/', Fit, '_', Region, '_predictions.tiff', sep = '')) # write as tif


#----------------------------------- Regional independent fit ---------------------------------

Fit <- 'REGIONAL_INDEPENDENT'
dir.create(Fit)

#--- Australia (AUS) ---------------------------------------------------------------

Region <- 'AUS'

set.seed(111) # seed

#--- Occurrence data processing -------------------------------------------------------------

# Read occurrence data
dat <- read.csv('Appendix C.csv')
dat <- subset(dat, Region == 'Australia')

# Get rid of duplicates
dups <- duplicated(dat[c('Longitude', 'Latitude')])
occ_cal <- dat[!dups, ]

# Spatial thinning: 50 km is the minimum distance to remove spatial clusters (see script for Figure S1)
train_thin <- thin(occ_cal, lat.col = 'Latitude', long.col = 'Longitude', spec.col = 'Species',
                   thin.par = 50, reps = 10, locs.thinned.list.return = T, write.files = F, write.log.file = F)
occ_cal <- data.frame(Species = Species, Longitude = round(train_thin[[1]]$Longitude, 2), 
                      Latitude = round(train_thin[[1]]$Latitude, 2))
occ_cal$Type <- 'calibration'

# Initial M space for environmental filtering
coord <- cbind(occ_cal$Longitude, occ_cal$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point
env.M <- crop(env.G, occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon
env.M <- stack(env.M)

# More occurrences processing
# Save any points to nearest pixel with data that may be falling outside a layer of the M space raster 
for(i in 1:length(env.M@layers)){
  spp <- SpatialPoints(occ_cal[, c('Longitude', 'Latitude')], crs)
  spp_corrected <- points2nearestcell(locs = spp, ras = env.M, layer = i, move = T) 
  occ_cal$Longitude <- round(spp_corrected@coords[, 1], 3)  # replace coordinates including those new if any
  occ_cal$Latitude <- round(spp_corrected@coords[, 2], 3)
} # plots will appear if any correction applies

# Get rid of duplicates again in case some corrected points overlap
dups <- duplicated(occ_cal[c('Longitude', 'Latitude')])
occ_cal <- occ_cal[!dups, ]

# Environmental filtering
source('envSample.R')
coords <- occ_cal[, c('Longitude', 'Latitude')]
env.data <- extract(env.M, coords)
env.data <- as.data.frame(env.data)
coords <- envSample(coords, filters = list(env.data$Surface_temperature, env.data$Distance_to_coast),
                    res = list(0.5, 1), do.plot = T)
occ_cal <- merge(occ_cal, coords, by.x = c('Longitude', 'Latitude'), by.y = c('lon', 'lat'))

# Detecting environmental outliers in species occurrences (Cobos et al., 2018)
# Getting data from the variables
variables_values <- na.omit(values(env.M)) # for the region of interest
occ_variables <- na.omit(cbind(as.numeric(row.names(occ_cal)), extract(env.M, occ_cal[, c('Longitude', 'Latitude')])))
occ_variables <- as.data.frame(occ_variables)
if(dim(variables_values)[1] > 10000) { # sample of 10000 values if more pixels exist
  variables_values <- variables_values[sample(1:dim(variables_values)[1], 10000), ]
}

# Plot for searching for potential environmental outliers
# plot all combinations, then scroll plots and manually detect outliers and write them down in code below
for(i in 1:dim(env.M)[3]){
  for(j in 1:dim(env.M)[3]){
    par(mfcol = c(1, 2), mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
    plot(variables_values[, i], variables_values[, j], col = 'grey65', pch = 1,
         xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    points(occ_variables[, i + 1], occ_variables[, j + 1], col = 'blue', pch = 19, cex = 1.5)
    legend('bottomright', legend = c('Region of interest', 'Occurrences'),
           pch = c(1, 19), col = c('gray65', 'black'), bty = 'n')
    plot(variables_values[, i], variables_values[, j], col = 'gray65',
         pch = 1, xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
    text(occ_variables[, i + 1], occ_variables[, j + 1], occ_variables[, 1], cex = 1, col = 'blue')
    legend('bottomright', legend = 'Occurrence ID', bty = 'n')
  }
}

# Remove outliers from calibration records
occ_variables <- subset(occ_variables, !V1 %in% c(51)) # PUT IN c() INTEGER POINT CODES TO REMOVE 
occ_cal <- occ_cal[which(row.names(occ_cal) %in% as.character(occ_variables$V1)), ]

# Final calibration points
occ_cal <- rbind(occ_cal[, c(3, 1, 2, 4)])
dir.create(paste(Fit, '/', Region, sep = ''))
write.csv(occ_cal, paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_points.csv', sep = ''), row.names = F)


#--- Calibration areas -----------------------------------------------------------------

occ_cal <- read.csv(paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_points.csv', sep = ''))

# Final M space (calibration areas)
coord <- cbind(occ_cal$Longitude, occ_cal$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # sampling bias: 1000 km buffer around each point
env.M <- crop(env.G, occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon

# Check correlation in M space
layerStats(env.M, 'pearson', na.rm = T) # all < 0.8, all good

writeRaster(env.M, filename = paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_areas.tiff', sep = '')) # write as tif


#--- Model calibration, evaluation and selection -----------------------------------------------

# Calibration points
occ_cal <- read.csv(paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_points.csv', sep = ''))
occ_cal <- occ_cal[, c('Longitude', 'Latitude')]
colnames(occ_cal) <- c('longitude', 'latitude')

# Calibration areas
env.M <- stack(paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_areas.tiff', sep = ''))
names(env.M) <- var_set

# Background points
bg_points <- randomPoints(env.M[[1]], n = 10000)
colnames(bg_points) <- c('longitude', 'latitude')

# Presence-background points
all_pts <- rbind(occ_cal, bg_points)
all_pts <- SpatialPoints(all_pts, proj4string = crs) # transform to spatial points

# Check the effective range of spatial autocorrelation
eff.range <- spatialAutoRange(rasterLayer = env.M, sampleNumber = 10000, doParallel = T, showPlots = T)
median_range <- 860000 # actual range of 858016 round up to 860000

# Spatial block design based on 'blockCV' by Valavi et al. (2019)
sp_block <- spatialBlock(speciesData = all_pts, rasterLayer = env.M[[1]], theRange = median_range, k = 3, selection = 'random')

# Get block information for modelling
user.grp <- list(occs.grp = sp_block$foldID[1:nrow(occ_cal)],
                 bg.grp = sp_block$foldID[(nrow(occ_cal) + 1):length(sp_block$foldID)])

# Running and evaluating candidate models
mod_regional_aus <- ENMevaluate(occs = occ_cal, envs = env.G, bg = bg_points, user.grp = user.grp, 
                                algorithm = 'maxnet', partitions = 'user', user.eval = proc, doClamp = F, 
                                tune.args = list(fc = c('LQ', 'LP', 'LQP'), rm = c(seq(0.1, 2, 0.1), 2.5, 3, 4, 5, 7, 10)))

# Model selection
res <- eval.results(mod_regional_aus) # overall results
opt.seq <- res %>% 
  filter(proc_pval.avg == min(proc_pval.avg)) %>%
  filter(or.10p.avg == min(or.10p.avg)) # 3-step criteria in Cobos et al. (2019)

# Optimal model settings
mod.seq <- eval.models(mod_regional_aus)[[opt.seq$tune.args]]

# Model predictions
pred.seq <- eval.predictions(mod_regional_aus)[[opt.seq$tune.args]]
# pred.seq <- mosaic(pred.seq[[1]], pred.seq[[2]], fun = median) # in case two or more model are selected as best models
writeRaster(pred.seq, filename = paste(Fit, '/', Region, '/', Fit, '_', Region, '_predictions.tiff', sep = '')) # write as tif


#----------------------------------- Extrapolation risk (MOP) ---------------------------------

# Based on Mobility Oriented Parity (MOP) from Owens et al. (2013) 

# Global fit
Fit <- 'GLOBAL'
env.M <- stack(paste(Fit, '/', Fit, '_calibration_areas.tiff', sep = ''))
names(env.M) <- var_set
mop_analysis <- mop(M_stack = env.M, G_stack = env.G, percent = 10, comp_each = 2000, parallel = T)
writeRaster(mop_analysis, filename = paste(Fit, '/', Fit, '_mop.tiff', sep = '')) # write as tif

# Regional merged fit
Fit <- 'REGIONAL_MERGED'
env.M <- stack(paste(Fit, '/', Fit, '_calibration_areas.tiff', sep = ''))
names(env.M) <- var_set
mop_analysis <- mop(M_stack = env.M, G_stack = env.G, percent = 10, comp_each = 2000, parallel = T)
writeRaster(mop_analysis, filename = paste(Fit, '/', Fit, '_mop.tiff', sep = '')) # write as tif

# Regional independent fit
Fit <- 'REGIONAL_INDEPENDENT'
Region <- 'SWA'
env.M <- stack(paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_areas.tiff', sep = ''))
names(env.M) <- var_set
mop_analysis <- mop(M_stack = env.M, G_stack = env.G, percent = 10, comp_each = 2000, parallel = T)
writeRaster(mop_analysis, filename = paste(Fit, '/', Region, '/', Fit, '_', Region, '_mop.tiff', sep = '')) # write as tif

Region <- 'AUS'
env.M <- stack(paste(Fit, '/', Region, '/', Fit, '_', Region, '_calibration_areas.tiff', sep = ''))
names(env.M) <- var_set
mop_analysis <- mop(M_stack = env.M, G_stack = env.G, percent = 10, comp_each = 2000, parallel = T)
writeRaster(mop_analysis, filename = paste(Fit, '/', Region, '/', Fit, '_', Region, '_mop.tiff', sep = '')) # write as tif


# -------------------------------------------- END ------------------------------------------------
