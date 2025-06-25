packages <- c('ggplot2', 'stringr', "reshape2", "writexl",
              'stringr', 'gaston', 'gridExtra', 'ggpubr', 'dplyr')
for (p in packages){
  print(p)
  print(packageVersion(p))
  suppressWarnings(
    suppressPackageStartupMessages(library(p, character.only = TRUE)))
}

dir <- 'path_to_directory'

files <- list.files(paste(dir, 'Results/Power', sep = '/'),
                    pattern = paste('.*', 'RData', sep = ''))

findf <- data.frame()
for(i in 1:length(files)){
  load(paste(dir, 'Results/Power', files[i], sep = '/'))
  files.i <- strsplit(files[i], "_")[[1]]
  method.i <- files.i[1]
  nsamples.i <- as.numeric(files.i[2])
  mean.effect.i <- as.numeric(files.i[3])
  var.effect.i <- as.numeric(files.i[4])
  grp.split.pc.i <- as.numeric(files.i[5])
  outpc.i <- as.numeric(files.i[6])
  pDM.i <- as.numeric(strsplit(files.i[7], ".RD")[[1]][1])
  
  if(nrow(findf) == 0){
    findf <- data.frame(restmp,
                        nsamples = nsamples.i,
                        mean.effect = mean.effect.i,
                        var.effect = var.effect.i,
                        grp.split.pc = grp.split.pc.i,
                        outpc = outpc.i,
                        pDM = pDM.i)
  }else{
    findf.sub <- data.frame(restmp,
                            nsamples = nsamples.i,
                            mean.effect = mean.effect.i,
                            var.effect = var.effect.i,
                            grp.split.pc = grp.split.pc.i,
                            outpc = outpc.i,
                            pDM = pDM.i)
    findf <- rbind(findf, findf.sub)
  }
}

summary <- findf %>%
  group_by(Method) %>%
  summarise(
    median_FDR = median(FDR, na.rm = TRUE),
    median_Power = median(Power, na.rm = TRUE),
    median_MCC = median(MCC, na.rm = TRUE),
    median_Time = median(Time, na.rm = TRUE)
  )

summary$Method <- str_remove_all(summary$Method, regex("-?Test", ignore_case = TRUE))

assign_symbol <- function(values, higher_is_better = TRUE) {
  q <- quantile(values, probs = c(1/3, 2/3), na.rm = TRUE)
  
  if (!higher_is_better) {
    return(case_when(
      values <= q[1] ~ "●",
      values <= q[2] ~ "◐",
      TRUE ~ "○"
    ))
  } else {
    return(case_when(
      values <= q[1] ~ "○",
      values <= q[2] ~ "◐",
      TRUE ~ "●"
    ))
  }
}

findf.params <- unique(findf[-c(1:12)])

summary <- summary %>%
  mutate(
    FDR_symbol = assign_symbol(median_FDR, higher_is_better = FALSE),
    Power_symbol = assign_symbol(median_Power, higher_is_better = TRUE),
    MCC_symbol = assign_symbol(median_MCC, higher_is_better = TRUE),
    Time_symbol = assign_symbol(median_Time, higher_is_better = FALSE)
  )

write_xlsx(summary, path = "meta_summary.xlsx")
