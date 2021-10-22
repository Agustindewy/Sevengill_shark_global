
# De Wysiecki et al. - World ENM projection for the broadnose sevengill shark (Notorynchus cepedianus)

# Modelling using automated protocol in 'kuenm' package (Cobos et al., 2019)

library(kuenm)
library(rgdal)
library(raster)
library(LSD)
library(rgeos)

setwd('SET YOUR WORKING DIRECTORY')

# Species
Species1 <- 'Notorynchus cepedianus'
Species2 <- 'Notorynchus_cepedianus'

# Set region you want to model, either 'SWA' or 'AUS'
Region <- 'SWA'  

#----------------------------------- Variable selection ---------------------------------------------

# Seed
set.seed(111)

#------------------------------------ Candidate models ----------------------------------------------

occ_joint <- paste(Region, '/N_cepedianus_joint.csv', sep = '')
occ_tra <- paste(Region, '/N_cepedianus_train.csv', sep = '')
M_var_dir <- paste(Region, '/Var_selection/M_variables', sep = '')
batch_cal <- paste(Region, '/Var_selection/Candidate_Models', sep = '')
out_dir <- paste(Region, '/Var_selection/Candidate_Models', sep = '')
reg_mult <- c(0.1, 0.5, 1, 2.5, 5) # coarse set of values
f_clas <- c('lq', 'lp', 'lqp')
args <- NULL
maxent_path <- 'PATH TO MAXENT JAR FILE'
wait <- TRUE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal, 
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args, maxent.path = maxent_path, 
          wait = wait, run = run)

#-------------------------- Evaluation and selection of best models ----------------------------

occ_test <- paste(Region, '/N_cepedianus_test.csv', sep = '')
out_eval <- paste(Region, '/Var_selection/Calibration_Results', sep = '')
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- FALSE
selection <- 'AICc'
paral_proc <- FALSE 

kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, 
            batch = batch_cal, out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, 
            iterations = iterations, kept = kept, selection = selection, parallel.proc = paral_proc)

#---------------------------------- Final model creation --------------------------------

batch_fin <- paste(Region, '/Var_selection/Final_Models', sep = '')
mod_dir <- paste(Region, '/Var_selection/Final_Models', sep = '')
rep_n <- 10
rep_type <- 'Bootstrap'
jackknife <- FALSE
out_format <- 'logistic'
project <- FALSE
G_var_dir <- NULL
ext_type <- 'no_ext'
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- TRUE
run1 <- TRUE
args <- NULL

kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, 
          rep.n = rep_n, rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, 
          out.format = out_format, project = project, G.var.dir = G_var_dir, ext.type = ext_type, 
          write.mess = write_mess, write.clamp = write_clamp, maxent.path = maxent_path,args = args, 
          wait = wait1, run = run1)

#------------------------------------- END variable selection -----------------------------------------

# Once you created M and G final files come back here

#------------------------------------- Final modelling ---------------------------------------------

# Set region to model, either 'SWA' or 'AUS'
Region <- 'SWA'  

# Seed
set.seed(111)

#------------------------------------ Candidate models -----------------------------------------

occ_joint <- paste(Region, '/N_cepedianus_joint.csv', sep = '')
occ_tra <- paste(Region, '/N_cepedianus_train.csv', sep = '')
M_var_dir <- paste(Region, '/M_variables', sep = '')
batch_cal <- paste(Region, '/Candidate_Models', sep = '')
out_dir <- paste(Region, '/Candidate_Models', sep = '')
reg_mult <- c(seq(0.1, 2, 0.1), 2.5, 3, 4, 5, 7, 10) # finer set of values
f_clas <- c('lq', 'lp', 'lqp')
args <- NULL
maxent_path <- 'PATH TO MAXENT JAR FILE'
wait <- TRUE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal, 
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args, maxent.path = maxent_path, 
          wait = wait, run = run)

#-------------------------- Evaluation and selection of best models ----------------------------

occ_test <- paste(Region, '/N_cepedianus_test.csv', sep = '')
out_eval <- paste(Region, '/Calibration_Results', sep = '')
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- FALSE
selection <- 'AICc'
paral_proc <- FALSE 

kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, 
            batch = batch_cal, out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, 
            iterations = iterations, kept = kept, selection = selection, parallel.proc = paral_proc)

#---------------------------------- Final model creation --------------------------------

batch_fin <- paste(Region, '/Final_Models', sep = '')
mod_dir <- paste(Region, '/Final_Models', sep = '')
rep_n <- 10
rep_type <- 'Bootstrap'
jackknife <- FALSE
out_format <- 'logistic'
project <- TRUE
G_var_dir <- paste(Region, '/G_variables', sep = '')
ext_type <- 'ext'
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- TRUE
run1 <- TRUE
args <- NULL

kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, 
          rep.n = rep_n, rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, 
          out.format = out_format, project = project, G.var.dir = G_var_dir, ext.type = ext_type, 
          write.mess = write_mess, write.clamp = write_clamp, maxent.path = maxent_path,args = args, 
          wait = wait1, run = run1)

#----------------------------- Evaluation with independent data -----------------------------------

occ_ind <- paste(Region, '/N_cepedianus_ind.csv', sep = '')
replicates <- TRUE
out_feval <- paste(Region, '/Final_Evaluation', sep = '')

fin_eval <- kuenm_feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
                        out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
                        iterations = iterations, parallel.proc = paral_proc)

#--------------------- Uncertainty: median and range of best selected models ------------------------

sp_name <- Species2
format <- 'asc'
project <- TRUE
stats <- c('median', 'range')
rep <- TRUE
scenarios <- 'current'
ext_type <- 'E'
out_dir <- paste(Region, '/Final_Model_Stats', sep = '')

kuenm_modstats(sp.name = sp_name, fmod.dir = mod_dir, format = format, project = project,
               statistics = stats, replicated = rep, proj.scenarios = scenarios,
               ext.type = ext_type, out.dir = out_dir)

#----------------------------------- Extrapolation risk -------------------------------------------

# Based on MOP (Mobility Oriented Parity) from Owens et al. (2013) 

G_var_dir <- paste(Region, '/G_variables', sep = '')
M_var_dir <- paste(Region, '/M_variables', sep = '')
sets_var <- 'SetXX' # indicate best selected set of variables
out_mop <- paste(Region, '/MOP_Results', sep = '')
percent <- 10
paral <- FALSE
is.swd <- FALSE

kuenm_mmop(G.var.dir = G_var_dir, M.var.dir = M_var_dir, sets.var = sets_var,
           out.mop = out_mop, percent = percent, parallel = paral)

#------------------------------------------- END -----------------------------------------------
