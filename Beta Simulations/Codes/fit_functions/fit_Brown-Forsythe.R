#### Fit Brown-Forsythe test To A Dataset ####

fit.BrownForsythetest <- function(beta_vals, metadata, expVar, coVars){
  
  #### Tracking time starts ####
  
  start.time <- Sys.time()
  
  #### Standard Brown-Forsythe test pipeline for 2 groups ####
  mem <- peakRAM::peakRAM({
    spe <- function(x, metadata, expVar, coVars){
      if(is.null(coVars)){
        regData <- metadata[, c(expVar), drop = FALSE]
      }else{
        regData <- metadata[, c(expVar, coVars)]
      }
      regData.new <- data.frame(expr = x, regData)
      formula <- as.formula(paste("expr ~ ", paste(colnames(regData), collapse = "+")))
      fit <- onewaytests::bf.test(formula = formula, data = regData.new, alpha = 0.05, na.rm = TRUE, verbose = FALSE)
      resTab <- data.frame(statistic = as.numeric(fit$statistic), Pval = fit$p.value)
      return(resTab)
    }
    
    #### Per-cpg Fitting ####
    
    beta_vals_na <- na.omit(beta_vals)
    res <- apply(beta_vals_na, 1, function(x) spe(x = x, metadata = metadata, expVar = expVar, coVars = coVars))
    res <- do.call('rbind', res)
    
    #### Extracting Results ####
    
    paras <- data.frame(cpgs = rownames(res),
                        Statistic = res$statistic,
                        Pval = res$Pval)
    paras$adjPval <- p.adjust(paras$Pval, method = 'BH')
    paras <- paras[order(paras$adjPval, decreasing = FALSE),]
  })
  stop.time <- Sys.time()
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  cpgs <- data.frame(cpgs = row.names(beta_vals_na))
  paras <- merge(cpgs, paras, id = 'cpgs', all = TRUE)
  finres = list(Method = 'Brown-Forsythe Test',
                res = paras,
                Time.min = Time.min,
                peak.memory.gb = peak.memory.gb)
  return(finres)
}