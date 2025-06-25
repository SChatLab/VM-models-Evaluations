#### Variance test between 2 Groups Function ####

fit.heteroTest <- function(features, metadata, test, ncores){
  
  #### Assign variance test ####
  
  geneVar <- function(x, test){
    if(!is.factor(metadata$Exposure)){
      group <- as.factor(metadata$Exposure)
    }else{
      group <- metadata$Exposure
    }
    if(test == 'Fligner'){
      fit <- fligner.test(x ~ group)
      pval <- fit$p.value
    }else if(test == 'MBF'){
      fit <- MBF(data = x, group = group)
      pval <- fit[1, 4]
    }else if(test == 'BF'){
      fit <- bf.test(x ~ group, data = data.frame(x = x, group = group))
      pval <- fit$p.value
    }else if(test == 'Levene'){
      fit <- leveneTest(x ~ group)
      pval <- fit[1, 3]
    }else if(test == 'F-test'){
      fit <- var.test(x ~ group, data = data.frame(x = x, group = group))
      pval <- fit$p.value
    }else if(test == 'ncv'){
      mod <- lm(x ~ group)
      res <- as.numeric(residuals(mod))
      resnorm <- ks.test(res,"pnorm")
      resP <- resnorm$p.value
      fit <- ncvTest(mod)
      pval <- fit$p
    }else if(test == 'Koenkar'){
      mod <- lm(x ~ group)
      fit <- breusch_pagan(mod)
      pval <- fit$p.value
    }else if(test == 'OBrien'){
      fit <- lawstat::levene.test(y = x, group = group, location = 'median',
                                  correction.method = "correction.factor",
                                  bootstrap = FALSE, kruskal.test = TRUE)
      pval <- fit$p.value
    }else if(test == 'whiteboot'){
      mod <- lm(x ~ group)
      fit <- white_test_boot(mod, bootstraps = 1000)
      pval <- fit$p_value
    }
    return(pval)
  }
  
  #### Testing Equal Variance for all genes ####
  
  cl <- parallel::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)
  packages <- c('edgeR', 'pkgmaker', 'stringr', 'foreach', 'doParallel', 'pbapply', 'car',
                'lawstat', 'lmtest', 'MASS', 'skedastic', 'doex', 'whitestrap',
                'onewaytests')
  pb <- txtProgressBar(max = nrow(features), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  p <- foreach(j = 1:nrow(features), .combine = 'c', .packages = packages, 
               .options.snow = opts) %dopar% {
                 x <- as.numeric(features[j, ])
                 p <- geneVar(x, test = test)
                 return(p)
               }
  close(pb)
  stopCluster(cl)
  fit <- data.frame(Genes = rownames(features), Pval = p)
  fit$adj.Pval <- p.adjust(fit$Pval, method = 'BH')
  return(fit)
}