
thfunc <- function(median_mod,suits,percent){
  # sort suitability values into ascending order and get rid of NAs
  suits <- na.omit(sort(suits))
  npts <- length(suits)
  # case 1:  all the suitability values are considered
  if(percent==0){
    threshold <- suits[1]
    #bin_mod <- median_mod>=threshold
    return(threshold)
  }
  # case 2: only a percentage of the suitability values are considered
  # define the smallest integer not less than the corresponding percentage of presence points
  pth <- ceiling(npts*percent)
  threshold <- suits[pth]
  
  #bin_mod <- median_mod>threshold
  return(threshold)
}


sim_disperal <- function(adL,bin_thmodels,base_name,pasos){
  # Run the process in parallel
  plan(multiprocess)
  # Increase the available memory to run process
  options(future.globals.maxSize = 3000 * 1024^2)
  # Run simulations for each dispersal scenario
  results_sims <- seq_along(adL) %>% furrr::future_map(function(x){
    # Number of theresholded maps trimed by the MOP (10)
    nlay <- raster::nlayers(bin_thmodels)
    # Adjacency matrix
    ad_matrix <- adL[[x]]
    # apply bam simulation process using occurrence points as the initial conditions
    sim_bin_thmodels <- lapply(1:nlay, function(y){
      sparse_mod <- bam::model2sparse(bin_thmodels[[y]])
      # Run the simulation 
      sim <- bam::sdm_sim(set_A = sparse_mod,
                          set_M =  ad_matrix,
                          initial_points = occs_sparse,
                          nsteps = pasos)
      # convert simulation results to raster
      sim_final <- bam::sim2Raster(sim,pasos)
      
    })
    sim_bin_thmodels <- raster::stack(sim_bin_thmodels)
    # Estimate concensus cells (cells that were
    # invaded in the dispersal scenarios)
    sim_bin_thmodels_sum <- raster::calc(sim_bin_thmodels,sum)
    
    return(sim_bin_thmodels_sum)
  },.progress = TRUE)
  
  # Create an animation with the maps obtained with the different distance values
  plan(sequential)
  names(results_sims) <- names(adL)
  if(!dir.exists(base_name)) dir.create(base_name)
  
  base_path <- normalizePath(base_name)
  base_path <- gsub("[\\]","/",base_path)
  # Animation path 
  anipath <- file.path( base_path,
                        paste0("animation_sims_",
                               base_name,".gif"))
  
  #Pathos to the raster files 
  rasters_path <- file.path( base_path,
                             paste0("Sim_",
                                    base_name,"_D_",
                                    names(adL),".tif"))
  # Save the animations 
  
  animation::saveGIF({
    for(i in seq_along(results_sims)){
      plot(results_sims[[i]], 
           main= paste("Dispersal distance ",names(results_sims[i])))
    }
  },movie.name = anipath,interval = 0.8,
  ani.width = 1200, ani.height = 1200, ani.res = 300)
  
  res_d <- lapply(seq_along(rasters_path), function(x){
    writeRaster(results_sims[[x]],  rasters_path[x],overwrite=T)
  })
  return(res_d)
}
