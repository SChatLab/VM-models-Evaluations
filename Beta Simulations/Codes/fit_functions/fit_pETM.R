#### Fit pETM To A Dataset ####

fit.pETM <- function(beta_vals, metadata, expVar, coVars){
  
  #### Tracking time starts ####
  
  start.time <- Sys.time()
  
  #### Standard pETM pipeline for 2 groups ####
  mem <- peakRAM::peakRAM({
    conds <- unique(metadata[, expVar])
    if(is.data.frame(conds)){
      conds <- conds$Exposure
    }
    group0indexes <- which(metadata[, expVar] == conds[1])
    group1indexes <- which(metadata[, expVar] == conds[2])
    grpvar <- c(rep(0, length(group0indexes)), rep(1, length(group1indexes)))
    beta_vals_na <- na.omit(beta_vals)
    beta_vals_na <- as.data.frame(t(beta_vals_na))
    gr <- NULL
    res <- pETM::pETM(x = beta_vals_na, y = grpvar, cx = NULL, alpha = 0.1, maxit = 100000, 
                      thre = 1e-6, group = gr, lambda = NULL, type = "ring", etm = "beta", 
                      psub = 0.5, nlam = 10, kb = 10, K = 100)
    
    #### Extracting Results ####
    
    paras <- data.frame(res$selprob)
    names(paras)[1] <- 'selprob'
    paras$cpgs <- row.names(paras)
    paras$adjPval <- ifelse(paras$selprob > 0.6, 0.01, 1)
    paras$Pval <- paras$adjPval
    paras <- paras[order(paras$adjPval, decreasing = FALSE),]
    
  })
  stop.time <- Sys.time()
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  finres = list(Method = 'pETM',
                res = paras,
                Time.min = Time.min,
                peak.memory.gb = peak.memory.gb)
  return(finres)
}