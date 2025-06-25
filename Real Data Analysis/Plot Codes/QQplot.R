#### Load Libraries ####

packages <- c('doParallel', 'snow', 'robustbase', 'skedastic', 'doex', 'whitestrap',
              'ggplot2', 'qqplotr', 'ggpubr', 'randomcoloR')
for (p in packages){
  print(p)
  print(packageVersion(p))
  suppressWarnings(
    suppressPackageStartupMessages(library(p, character.only = TRUE)))
}

#### Set Path ####

dir <- 'path_to_directory/Real Data Analysis'

#### Load Plot Data ####

all.files <- list.files(paste(dir, 'Data', sep = '/'), pattern = '*.RData')
new.file <- NULL
datname <- c('Preeclampsia_450K', 'BreastCancer_850K')
for(j in 1:length(all.files)){
  
  #### Load Data ####
  
  load(paste(dir, 'Data', all.files[j], sep = '/'))
  
  #### Extract Data ####
  
  datnam <- tools::file_path_sans_ext(all.files[j])
  beta_vals <- data$beta_vals
  metadata <- data$metadata
  new.file[j] <- c(paste(ncol(beta_vals), datnam, sep = ''))
}
new.files <- gtools::mixedsort(new.file)

#### Plot ####

hetlist <- list()
for(p in 1:length(new.files)){
  
  #### Extract Data ####
  
  nsamples <- as.numeric(gsub("([0-9]+).*$", "\\1", new.files[p]))
  datnam <- datname[p]
  
  #### Load Data ####
  
  tmp.name <- paste(datnam, '.RData', sep = '')
  load(paste(dir, 'Data', tmp.name, sep = '/'))
  
  #### Filter Data ####
  
  message("Data ", datnam, " running")
  beta_vals <- data$beta_vals
  metadata <- data$metadata
  
  #### Unequal Variance Test ####
  
  het_test <- fit.heteroTest(beta_vals, metadata, test = 'Levene', ncores = 4)
  
  #### Plot Data ####
  
  title <- paste(LETTERS[p], ') ', datnam, ' (N = ', nsamples, ')', sep = '')
  ci <- 0.95
  ps = as.numeric(het_test$Pval)
  n  <- length(ps)
  df <- data.frame(observed = -log10(sort(ps)),
                   expected = -log10(ppoints(n)),
                   clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
                   cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1)))
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  hetplot <- ggplot(df, aes(x = reorder(dataset, nsamples))) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      fill = 'lightgrey',
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3, fill = 'steelblue') +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_bw(base_size = 24) +
    theme(axis.ticks = element_line(size = 1), panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16)) +
    ggtitle(title)
  hetlist[[p]] <- hetplot
}

#### Multi-panel plot ####

png(paste(dir, "Plots", 'QQplot.png', sep = '/'), units = 'in',
    width = 20, height = 20, res = 300)
ggarrange(hetlist[[1]], hetlist[[2]], nrow = 2, ncol = 1)
dev.off()
