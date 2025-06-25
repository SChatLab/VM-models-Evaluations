#### Load Generic Libraries ####

packages <- c('pkgmaker', 'stringr', 'foreach', 'doParallel', 'MASS', 'data.table',
              'matrixStats', 'missMethyl', 'evora', 'caret', 'pROC', 'jlst', 'qvalue', 
              'tidyverse', 'readr', 'peakRAM', 'onewaytests', 'dglm', 'pETM')
for (i in packages){
  print(i)
  print(packageVersion(i))
  suppressWarnings(
    suppressPackageStartupMessages(library(i, character.only = TRUE)))
}

#### Set Path ####

dir <- 'path_to_directory/Normal Simulations'

#### Load functions ####

pkgmaker::source_files(paste(dir, "Codes", "fit_functions", sep = '/'),'*.R')

#### Generate Simulated data with 100 replicates (in-parallel) ####

ncores <- 10
nreps <- 100
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)
exports <- c('fit.Bartletttest', 'fit.BrownForsythetest', 'fit.DGLM', 'fit.DiffVar', 'fit.EVORA', 
             'fit.iEVORA', 'iEVORA', 'fit.JLSsc', 'fit.JLSp', 'fit.OLS', 'perfMetrics', 'permMetrics',
             'simvar', 'doDV', 'doTT', 'fit.pETM')
restmp <- foreach(r = 1:nreps, .combine = rbind, .packages = packages, .export = exports) %dopar% {
  seedR <- 1233 + r
  simdata <- normal.simvar(Ncpgs = 10000,
                           nsamples = nsamples,
                           var.effect = var.effect,
                           grp.split.pc = grp.split.pc,
                           pDM = pDM,
                           outpc = outpc,
                           seed = seedR)
  
  #### Gathering Simulated Data ####
  
  mvals <- simdata$mval
  metadata <- simdata$metadata
  dm.cpgs <- simdata$dm.cpgs
  CpG.Info <- simdata$CpG.Info
  expVar <- 'Exposure'
  coVars <- NULL
  
  #### Fit models ####
  
  fun <- paste(paste('fit', Method, sep = '.'),"(mvals,metadata,expVar,coVars)", sep = '')
  fit <- eval(parse(text = fun))
  
  #### Collect performance metrics ####
  
  if(pDM == 0){
    performance <- permMetrics(fit,
                               simulation = r, 
                               alpha = 0.05)
  }else{
    performance <- perfMetrics(fit, 
                               dm.cpgs, 
                               simulation = r, 
                               alpha = 0.05, 
                               all.metrics = TRUE)
  }
  return(performance)
}
stopCluster(cl)

#### Save Results ####

sim.scene <- paste(Method, nsamples, var.effect, grp.split.pc, pDM, sep = '_')
savenam <- paste(sim.scene, ".RData", sep = '')

if(pDM == 0){
  save(restmp, file = paste(dir, "Results/Null", savenam, sep = '/'))
}else{
  save(restmp, file = paste(dir, "Results/Power", savenam, sep = '/'))
}
