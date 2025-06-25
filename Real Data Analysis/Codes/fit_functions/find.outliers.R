###################################################
#### Finding Influential observations per CpG ####
##################################################

find.outliers <- function(betas, metadata) {
  n <- nrow(betas)
  outlier_list <- vector("list", n)
  for (i in seq_len(n)) {
    regData <- data.frame(y = as.numeric(betas[i, ]), Exposure = metadata$Exposure)
    model <- lm(y ~ Exposure, data = regData)
    cooksD <- cooks.distance(model)
    threshold <- 3 * mean(cooksD, na.rm = TRUE)
    influential <- sum(cooksD > threshold, na.rm = TRUE)
    outlier_list[[i]] <- data.frame(
      CpG = rownames(betas)[i],
      outlier.pc = round(100 * influential / ncol(betas))
    )
  }
  full.out.df <- do.call(rbind, outlier_list)
  return(full.out.df)
}

find.outliers.RLM <- function(betas, metadata) {
  n <- nrow(betas)
  outlier_list <- vector("list", n)
  for (i in seq_len(n)) {
    regData <- data.frame(y = as.numeric(betas[i, ]), Exposure = metadata$Exposure)
    model <- lmrob(y ~ Exposure, data = regData)
    fit.summary <- summary(model)
    out.thresh <- fit.summary$control$eps.outlier
    out.liers <- fit.summary$rweights[which(fit.summary$rweights <= out.thresh)]
    outlier_list[[i]] <- data.frame(
      CpG = rownames(betas)[i],
      outlier.pc = round(100 * length(out.liers) / ncol(betas))
    )
  }
  full.out.df <- do.call(rbind, outlier_list)
  return(full.out.df)
}

dir <- 'path_to_directory/Real Data Analysis'
datnam <- c("Preeclampsia_450K.RData", "BreastCancer_850K.RData")

for(d in datnam){
  load(paste(dir, 'Data', d, sep = '/'))
  
  beta_vals <- data$beta_vals
  metadata <- data$metadata
  
  beta_vals[beta_vals >= 1] <- 0.99
  beta_vals[beta_vals <= 0] <- 0.001
  
  out <- find.outliers(beta_vals, metadata)
  #out <- find.outliers.RLM(beta_vals, metadata)
  savenam1 <- paste('outpc', paste0(d, '.RData'), sep = '_')
  #savenam2 <- paste('outpcRLM', paste0(d, '.RData'), sep = '_')
  save(out, file = paste(dir, 'Results', savenam1, sep = '/'))
  #save(b, file = paste(dir, 'Results', savenam2, sep = '/'))
}
