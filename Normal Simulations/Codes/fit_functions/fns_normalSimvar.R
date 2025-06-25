#### Normal Simulator Function ####

normal.simvar <- function(Ncpgs = 10000,
                          nsamples = 50,
                          var.effect = 1.5,
                          grp.split.pc = 70,
                          pDM = 0.1,
                          outpc = 0.1,
                          seed = seedR){
  
  set.seed(seed)
  
  group1_size <- round(nsamples * (grp.split.pc / 100))
  group2_size <- nsamples - group1_size
  
  # Set up prior parameters
  d0 <- 20 #degrees of freedom
  s0 <- 0.8 #first group sd
  
  # Number of differentially variable CpGs
  ndm <- Ncpgs * pDM
  dm <- c(rep(1, ndm), rep(0, Ncpgs - ndm))
  
  # Set up empty vectors and matrices to store results
  sigma2 <- matrix(0, ncol = 2, nrow = Ncpgs)
  
  # Sample the variances from a scaled inverse chisquare distribution
  sigma2[, 1] <- d0 * s0^2/(rchisq(Ncpgs, d0))
  sigma2[, 2] <- sigma2[, 1]
  
  # Sample the variances for the differentially variable CpGs to be roughly double in the second group compared to the first group
  if(ndm != 0){
    sigma2[1:ndm, 2] <- var.effect * d0/(rchisq(ndm, d0))
  }
  
  # Generate M values from a mixture of two normal distributions
  y <- matrix(0, ncol = nsamples, nrow = Ncpgs)
  index <- sample(1:Ncpgs, Ncpgs/2)
  
  y[index, ] <- matrix(rnorm(Ncpgs * nsamples/2, mean = -2), Ncpgs/2, nsamples)
  y[-index, ] <- matrix(rnorm(Ncpgs * nsamples/2, mean = 2), Ncpgs/2, nsamples)
  
  # Add noise into the M value signal from generated variances  
  y[, 1:group1_size] <- sqrt(sigma2[, 1]) * y[, 1:group1_size]
  y[, (group1_size + 1):nsamples] <- sqrt(sigma2[, 2]) * y[, (group1_size + 1):nsamples]
  
  # Add Outliers
  if(outpc > 0){
    outliers <- sample((ndm + 1):Ncpgs, Ncpgs * outpc)
    y[outliers, (group1_size + 1)] <- max(y)
    is.out <- rep('no', Ncpgs)
    is.out[outliers] <- 'yes'
  }else{
    is.out <- rep('no', Ncpgs)
  }
  
  # Summarize Output Data
  ids = rep(1:nsamples)
  group = as.factor(c(rep(0, group1_size), rep(1, group2_size)))
  metadata <- data.frame(Sample = ids,
                         Exposure = group)
  mvals <- as.data.frame(y)
  
  if(pDM == 0){
    rownames(mvals) <- c(paste('cpg', (ndm + 1):Ncpgs, sep = ''))
    dm.cpgs <- 'Type-1 simulation'
  }else{
    rownames(mvals) <- c(paste('cpg', 1:ndm, '_TP', sep = ''), paste('cpg', (ndm + 1):Ncpgs, sep = ''))
    dm.cpgs <- rownames(mvals[grep('_TP', rownames(mvals)),])
  }
  sim.mu.eff <- rowMeans(mvals[, (group1_size + 1):nsamples]) - rowMeans(mvals[, 1:group1_size])
  sim.var.eff <- rowVars(as.matrix(mvals[, (group1_size + 1):nsamples]))/rowVars(as.matrix(mvals[, 1:group1_size]))
  CpG.Info <- data.frame(CpGs = rownames(mvals), Mean.Effect = sim.mu.eff, Var.Effect = sim.var.eff, Outlier = is.out)
  return(list(mval = mvals,
              metadata = metadata,
              CpG.Info = CpG.Info,
              dm.cpgs = dm.cpgs))
}