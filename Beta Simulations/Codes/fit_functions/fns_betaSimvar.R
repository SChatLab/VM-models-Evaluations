##################################################################
#### Site-based simulator from variable dispersion beta model ####
##################################################################

beta.simvar <- function(ncpgs = 10000,
                        betavals.input = betavals.input,
                        nsamples = 50,
                        mu.eff = 0.01,
                        var.eff = 2,
                        pM = 0.1,
                        gprop = 50,
                        outpc = 0,
                        noutpc = 0,
                        seed = 1234){
  
  #### Fix seed ####
  
  set.seed(seed)
  
  #### Assigning sample sizes ####
  
  n1 <- round(nsamples * (gprop/100))
  n2 <- nsamples - n1
  
  #### Generate Exposure and Covariates ####
  
  Exposure <- as.factor(rep(0:1, times = c(n1, n2)))
  
  #### Get empirical mean and variances ####
  
  betavals.input <- betavals.input[complete.cases(betavals.input), ]
  muvar <- data.frame(mu = rowMeans(betavals.input),
                      var = rowVars(as.matrix(betavals.input)))
  muvar <- muvar[muvar$mu < 0.5, ]
  muvar <- muvar[sample(nrow(muvar), ncpgs), ]
  
  #### Convert mean and variance to shape parameters ####
  
  params1 <- estBetaParams(muvar$mu, muvar$var)
  rownames(params1) <- rownames(muvar)
  params1 <- subset(params1,  apply(params1, 1, function(x) all(x >= 0))) 
  
  #### Spiking parameters ####
  
  if(pM == 0 | var.eff == 1){
    params2 <- estBetaParams(muvar$mu, muvar$var)
    rownames(params2) <- rownames(muvar)
    params2$spiked <- rep(2, ncpgs)
    params2 <- subset(params2,  apply(params2, 1, function(x) all(x >= 0)))
  }else{
    params21 <- estBetaParams((mu.eff + muvar$mu[1:(pM*ncpgs)]), 
                              (var.eff * muvar$var[1:(pM*ncpgs)]))
    params22 <- estBetaParams((muvar$mu[((pM*ncpgs) + 1):ncpgs]), 
                              (muvar$var[((pM*ncpgs) + 1):ncpgs]))
    params2 <- data.frame(rbind(params21, params22))
    rownames(params2) <- rownames(muvar)
    params2$spiked <- rep(c(1, 2), c(pM*ncpgs, ncpgs - (pM*ncpgs)))
    params2 <- subset(params2,  apply(params2, 1, function(x) all(x >= 0)))
  }
  
  #### Gathering parameters ####
  
  params <- transform(merge(params1, params2, by = 'row.names', all = TRUE), 
                      row.names=Row.names, Row.names = NULL)
  params <- params[complete.cases(params), ]
  names(params)[1:4] <- c('alpha1', 'beta1', 'alpha2', 'beta2')
  
  #### Generate beta variates as vector operation ####
  
  beta_samples1 <- t(mapply(function(a, b) rbeta(n1, a, b), params$alpha1, params$beta1))
  beta_samples2 <- t(mapply(function(a, b) rbeta(n2, a, b), params$alpha2, params$beta2))
  grp1 <- matrix(beta_samples1, nrow = dim(beta_samples1)[1], ncol = n1)
  grp2 <- matrix(beta_samples2, nrow = dim(beta_samples1)[1], ncol = n2)
  
  #### Organize simulated data ####
  
  beta_vals <- data.frame(cbind(grp1, grp2))
  beta_vals[beta_vals >= 1] <- 0.99
  beta_vals[beta_vals <= 0] <- 0.001
  rownames(beta_vals) <- rownames(params)
  
  #### Adding Outliers ####
  
  if(outpc > 0){
    get.indices <- rownames(params[params$spiked == 2, ])
    get.indices <- sample(get.indices, noutpc*ncpgs)
    out.beta_vals <- beta_vals[rownames(beta_vals) %in% get.indices, ]
    out.params <- params[rownames(params) %in% get.indices, ]
    df.betaOuts <- NULL
    for(k in 1:length(get.indices)){
      tmp.gen <- as.numeric(out.beta_vals[rownames(out.beta_vals) %in% get.indices[k], ])
      names(tmp.gen) <- paste('S', 1:length(tmp.gen), sep = '')
      target.samps <- ceiling(outpc * (n1 + n2))
      tmp.gen.target <- names(sample(tmp.gen, target.samps))
      tmp.gen <- tmp.gen[!names(tmp.gen) %in% tmp.gen.target]
      q <- stats::qbeta(0.999, shape1 = out.params[rownames(out.params) %in% get.indices[k], 'alpha1'], 
                        shape2 = out.params[rownames(out.params) %in% get.indices[k], 'beta1'])
      out.add <- stats::runif(target.samps, q, 0.99)
      names(out.add) <- tmp.gen.target
      tmp.gen <- c(tmp.gen, out.add)
      tmp.gen <- as.numeric(tmp.gen)
      df.betaOuts <- data.frame(cbind(df.betaOuts, tmp.gen))
    }
    
    #### Inserting created outlier cpgs inside simulated data ####
    
    df.betaOuts <- data.frame(t(df.betaOuts))
    rownames(df.betaOuts) <- get.indices
    beta_vals <- beta_vals[!rownames(beta_vals) %in% get.indices, ]
    beta_vals <- data.frame(rbind(beta_vals, df.betaOuts))
    beta_vals[beta_vals >= 1] <- 0.99
    beta_vals[beta_vals <= 0] <- 0.001
    Is.outlier <- ifelse(rownames(params) %in% get.indices, 'yes', 'no')
  }else{
    Is.outlier <- 'no'
  }
  
  #### Cosmetic Organization ####
  
  metadata <- data.frame(Sample = paste0("Sample", 1:nsamples), Exposure = Exposure)
  names(beta_vals) <- metadata$Sample
  sim.mu.eff <- rowMeans(beta_vals[, (n1 + 1):nsamples]) - rowMeans(beta_vals[, 1:n1])
  sim.var.eff <- rowVars(as.matrix(beta_vals[, (n1 + 1):nsamples]))/rowVars(as.matrix(beta_vals[, 1:n1]))
  CpG.Info <- data.frame(params, sim.mu.eff, sim.var.eff, Is.outlier)
  if(pM == 0 | var.eff == 1){
    dm.cpgs <- 'Type-1 simulation'
  }else{
    dm.cpgs <- rownames(params[params$spiked == 1, ])
  }
  simData <- list(beta_vals = beta_vals, metadata = metadata, 
                  CpG.Info = CpG.Info, dm.cpgs = dm.cpgs)
  return(simData)
}