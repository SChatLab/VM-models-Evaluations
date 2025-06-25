#### Load Required Libraries ####

reqpkg = c("ComplexUpset", "dplyr", "purrr", "ggplot2", "ggpubr")
for (i in reqpkg){
  print(i)
  print(packageVersion(i))
  suppressPackageStartupMessages(library(i, character.only = TRUE)) 
}

#### Set Path ####

dir <- 'path_to_directory/Real Data Analysis'

#### Load Files ####

datnam <- c('Preeclampsia_450K', 'BreastCancer_850K')
files <- list.files(paste(dir, 'Results', sep = '/'),
                    pattern = paste('.*', 'RData', sep = ''))

for(d in 1:length(datnam)){
  files.sub <- files[grep(datnam[d], files)]
  fulldf <- data.frame()
  for(i in 1:length(files.sub)){
    load(paste(dir, 'Results', files.sub[i], sep = '/'))
    Method.i <- strsplit(files.sub[i], "_")[[1]][1]
    fitres <- fit[["res"]]
    fitres$binary <- ifelse(fitres$adjPval < 0.05, 1, 0)
    if(nrow(fulldf) == 0){
      fulldf <- data.frame(CpG = fitres$cpgs)
      fulldf <- data.frame(fulldf, fitres$binary)
      colnames(fulldf)[[i + 1]] <- Method.i
    }else{
      fulldf <- data.frame(fulldf, fitres$binary)
      colnames(fulldf)[[i + 1]] <- Method.i
    }
  }
  
  fulldf_filtered <- fulldf[apply(fulldf[, 2:ncol(fulldf)] != 0, 1, any), ]
  models <- names(fulldf_filtered)[2:ncol(fulldf_filtered)]
  
  ncpgs_per_method <- colSums(fulldf_filtered[, 2:ncol(fulldf_filtered)] == 1)
  new.labels <- setNames(
    paste0(names(ncpgs_per_method), " (", ncpgs_per_method, ")"),
    names(ncpgs_per_method)
  )
  
  upPlot <- upset(fulldf_filtered,
                  models,
                  queries = list(
                    upset_query(set = 'DiffVar', fill = '#e7298a'),
                    upset_query(set = 'iEVORA', fill = '#66a61e'),
                    upset_query(set = 'JLSp', fill = '#e6ab02'),
                    upset_query(set = 'JLSsc', fill = '#7570b3')),
                  base_annotations = list(
                    'Intersection size' = (
                      intersection_size(
                        bar_number_threshold = 1,
                        width = 0.8,
                        text = list(size = 5),)
                      + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
                      + theme(
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = 'black'),
                        axis.title = element_text(size = 16, face = 'bold'),
                        axis.text = element_text(size = 14, face = 'bold')))),
                  stripes = upset_stripes(
                    geom = geom_segment(size = 12),
                    colors = c('grey95', 'white')
                  ),
                  matrix =( intersection_matrix(
                    geom = geom_point(
                      shape ='circle filled',
                      size = 5,
                      stroke = 0.25))
                  + theme(
                    axis.title = element_text(size = 16, face = 'bold'),
                    axis.text = element_text(size = 14, face = 'bold')
                  )),
                  set_sizes=FALSE,
                  labeller = ggplot2::as_labeller(new.labels),
                  sort_sets = 'ascending',
                  sort_intersections = 'descending',
                  themes=upset_default_themes(axis.title = element_text(size = 16, face = 'bold'),
                                              axis.text=element_text(size = 16, face = 'bold')))
  
  
  #### Saving all Plots ####
  
  plotnam <- paste("UpsetPlot", paste(datnam[d], "png", sep = '.'), sep = '_')
  ggsave(paste(dir, 'Plots', plotnam, sep = '/'), width = 15, height = 10, dpi = 500)
}
