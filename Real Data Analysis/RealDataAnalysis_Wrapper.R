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

dir <- 'path_to_directory/Real Data Analysis'

#### Load Data ####

load(paste(dir, 'Data', paste(datnam, 'RData', sep = '.'), sep = '/'))

beta_vals <- data$beta_vals
metadata <- data$metadata
expVar <- 'Exposure'
coVars <- NULL
beta_vals[beta_vals >= 1] <- 0.99
beta_vals[beta_vals <= 0] <- 0.001

#### Load functions ####

pkgmaker::source_files(paste(dir, "Codes", "fit_functions", sep = '/'),'*.R')

#### Fit Real Data ####

fun <- paste(paste('fit', Method, sep = '.'),"(beta_vals,metadata,expVar,coVars)", sep = '')
fit <- eval(parse(text = fun))

#### Save fit results ####

savenam <- paste(Method, paste(datnam, 'RData', sep = '.'), sep = '_')
save(fit, file = paste(dir, "Results", savenam, sep = '/'))
