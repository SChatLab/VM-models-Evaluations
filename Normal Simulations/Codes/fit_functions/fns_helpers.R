
#### Function to calculate evaluation metrics ####

perfMetrics <- function(fit, dm.cpgs, simulation, alpha, all.metrics = TRUE){
  Method <- fit$Method
  resi <- fit$res
  resi[is.na(resi[, "adjPval"]), "adjPval"] <- 1
  adjPval <- resi$adjPval
  names(adjPval) <- resi$cpgs
  response <- rep(FALSE, length(adjPval))
  names(response) <- resi$cpgs
  response[names(response) %in% dm.cpgs] <- TRUE
  response <- factor(response, levels = c("TRUE", "FALSE"))
  predictor <- factor(adjPval <= 0.05, levels = c("TRUE","FALSE"))
  xtab <- table(predictor, response)
  cm <- caret::confusionMatrix(xtab)
  
  #### Getting Metrics ####
  
  if(all.metrics){
    Sensitivity <- format(cm$byClass['Sensitivity'], scientific = FALSE)
    Sensitivity <- as.numeric(substr(as.character(Sensitivity), 1, 4))
    Specificity <- format(cm$byClass['Specificity'], scientific = FALSE)
    Specificity <- as.numeric(substr(as.character(Specificity), 1, 4))
    FDR <- format(1 - cm$byClass['Precision'], scientific = FALSE)
    FDR <- as.numeric(substr(as.character(FDR), 1, 4))
    F1 <- format(cm$byClass['F1'], scientific = FALSE)
    F1 <- as.numeric(substr(as.character(F1), 1, 4))
    MCC <- format(round(mltools::mcc(preds = predictor, actuals = response), 2), scientific = FALSE)
    MCC <- as.numeric(substr(as.character(MCC), 1, 4))
    
    #### ROC ####
    
    response <- rep(FALSE, length(adjPval))
    names(response) <- resi$cpgs
    response[names(response) %in% dm.cpgs] <- TRUE
    response <- factor(response)
    predictor <- adjPval
    proc_out1 <- pROC::roc(response = response, predictor = predictor, direction = ">", quiet = TRUE)
    ROC <- as.numeric(proc_out1$auc)
    
    #### pAUROC ####
    
    proc_out2 <- roc(response = response, predictor = predictor,
                     direction = ">",quiet = TRUE, partial.auc = c(1, 0.85), 
                     partial.auc.correct = TRUE)
    pAUROC <- as.numeric(proc_out2$auc)
    
    #### Summarizing Metrics ####
    
    Measures <- data.frame(Method = Method,
                           Power = Sensitivity,
                           FPR = 1 - Specificity,
                           FDR = FDR,
                           F1 = F1,
                           MCC = MCC,
                           ROC = ROC,
                           pAUROC = pAUROC,
                           alpha = alpha,
                           replicate = simulation,
                           Time = fit$Time.min,
                           PeakMemory_GB = fit$peak.memory.gb)
  }else{
    Sensitivity <- round(as.numeric(cm$byClass['Sensitivity']), 2)
    Specificity <- round(as.numeric(cm$byClass['Specificity']), 2)
    FDR <- round(as.numeric(1 - cm$byClass['Precision']), 2)
    F1 <- round(as.numeric(cm$byClass['F1']), 2)
    MCC <- round(mltools::mcc(preds = predictor, actuals = response), 2)
    
    #### Summarizing Metrics ####
    
    Measures <- data.frame(Method = Method,
                           Power = Sensitivity,
                           FPR = 1 - Specificity,
                           FDR = FDR,
                           F1 = F1,
                           MCC = MCC,
                           alpha = alpha,
                           replicate = simulation,
                           Time = fit$Time.min,
                           PeakMemory_GB = fit$peak.memory.gb)
  }
  return(Measures)
}

#### Function to Benchmark Models for Null Analysis ####

permMetrics <- function(fit, simulation, alpha){
  Method <- fit$Method
  resi <- fit$res
  Time <- fit$Time.min
  resi[is.na(resi[, "Pval"]), "Pval"] <- 1
  resi[is.na(resi[, "adjPval"]), "adjPval"] <- 1
  DMcpgs.p <- resi[resi$Pval < alpha, ]
  DMcpgs.fdr <- resi[resi$adjPval < alpha, ]
  if(dim(DMcpgs.p)[1] == 0 & dim(DMcpgs.fdr)[1] == 0){
    permfit <- data.frame(Method = Method,
                          nfeatures = nrow(resi),
                          nDMcpgs.p = 0,
                          nDMcpgs.fdr = 0,
                          pDMcpgs.p = 0,
                          pDMcpgs.fdr = 0,
                          simulation = simulation,
                          Time = Time,
                          PeakMemory_GB = fit$peak.memory.gb)
  }else{
    permfit <- data.frame(Method = Method,
                          nfeatures = nrow(resi),
                          nDMcpgs.p = nrow(DMcpgs.p),
                          nDMcpgs.fdr = nrow(DMcpgs.fdr),
                          pDMcpgs.p = nrow(DMcpgs.p)/nrow(resi),
                          pDMcpgs.fdr = nrow(DMcpgs.fdr)/nrow(resi),
                          simulation = simulation,
                          Time = Time,
                          PeakMemory_GB = fit$peak.memory.gb)
  }
  return(permfit)
}