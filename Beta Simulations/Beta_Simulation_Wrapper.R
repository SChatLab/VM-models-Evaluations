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

dir <- 'path_to_directory/Beta Simulations'

#### Load functions ####

pkgmaker::source_files(paste(dir, "Codes", "fit_functions", sep = '/'),'*.R')

#### Load Template Data ####

load(paste(dir, 'Data', 'Preeclampsia_450K.RData', sep = '/'))
betavals.input <- data$beta_vals

#### Generate Simulated data with 100 replicates (in-parallel) ####

ncores <- 10
nreps <- 100
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)
exports <- c('fit.Bartletttest', 'fit.BrownForsythetest', 'fit.DGLM', 'fit.DiffVar', 'fit.EVORA', 
             'fit.iEVORA', 'iEVORA', 'fit.JLSsc', 'fit.JLSp', 'fit.OLS', 'perfMetrics', 'permMetrics',
             'new.bvalsimvar.sitebased', 'doDV', 'doTT', 'fit.pETM', 'estBetaParams')
restmp <- foreach(r = 1:nreps, .combine = rbind, .packages = packages, .export = exports) %dopar% {
  seedR <- 1233 + r
  simdata <- beta.simvar(ncpgs = 10000,
                         betavals.input = betavals.input,
                         nsamples = nsamples,
                         mu.eff = mu.eff,
                         var.eff = var.eff,
                         pM = pM,
                         gprop = gprop,
                         outpc = outpc,
                         noutpc = 0.02,
                         seed = seedR)
  
  #### Gathering Simulated Data ####
  
  beta_vals <- simdata$beta_vals
  metadata <- simdata$metadata
  dM.cpgs <- simdata$dm.cpgs
  CpG.Info <- simdata$CpG.Info
  expVar <- 'Exposure'
  coVars <- NULL
  
  #### Fit models ####
  
  fun <- paste(paste('fit', Method, sep = '.'),"(beta_vals,metadata,expVar,coVars)", sep = '')
  fit <- eval(parse(text = fun))
  
  #### Collect performance metrics ####
  
  if(pM == 0 | var.eff == 1){
    performance <- permMetrics(fit,
                               simulation = r, 
                               alpha = 0.05)
  }else{
    performance <- perfMetrics(fit, 
                               dM.cpgs, 
                               simulation = r, 
                               alpha = 0.05, 
                               all.metrics = TRUE)
  }
  return(performance)
}
stopCluster(cl)

#### Save Results ####

sim.scene <- paste(Method, nsamples, mu.eff, var.eff, gprop, outpc, pM, sep = '_')
savenam <- paste(sim.scene, ".RData", sep = '')

if(pM == 0 | var.eff == 1){
  save(restmp, file = paste(dir, "Results/Null", savenam, sep = '/'))
}else{
  save(restmp, file = paste(dir, "Results/Power", savenam, sep = '/'))
}
