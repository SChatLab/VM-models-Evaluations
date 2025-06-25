packages <- c('ggplot2', 'stringr', "reshape2", "forcats",
              'stringr', 'gaston', 'gridExtra', 'ggpubr', 'dplyr')
for (p in packages){
  print(p)
  print(packageVersion(p))
  suppressWarnings(
    suppressPackageStartupMessages(library(p, character.only = TRUE)))
}

dir <- 'path_to_directory'

files <- list.files(paste(dir, 'Results/Null', sep = '/'),
                    pattern = paste('.*', 'RData', sep = ''))

findf <- data.frame()
for(i in 1:length(files)){
  load(paste(dir, 'Results/Null', files[i], sep = '/'))
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

## Plots
FPR_list <- list()

var_eff <- 1
grp_split <- c(50, 70)
out_pc_list <- c(0, 0.1)
i <- 1

for(k in 1:length(out_pc_list)){
  out_pc <- out_pc_list[k]
  
  for(j in 1:length(grp_split)){
    grp_split_pc = grp_split[j]
    res <- findf %>%
      filter(var.effect == var_eff & grp.split.pc == grp_split_pc & outpc == out_pc)
    
    if(grp_split_pc == 50){
      split <- 'Equal Group Split'
    }else{
      split <- 'Unequal Group Split'
    }
    
    if(out_pc == 0){
      out <- 'Without Outliers'
    }else{
      out <- 'With Outliers'
    }
    
    title <- paste0(LETTERS[i], ') ', split, ', ', out)
    
    # Plot FPR
    
    summary_FPR <- res %>%
      group_by(nsamples, Method) %>%
      summarise(median_FPR = median(pDMcpgs.p), .groups = 'drop')
    
    summary_FPR$Method <- str_remove_all(summary_FPR$Method, regex("-?Test", ignore_case = TRUE))
    
    summary_FPR <- summary_FPR %>%
      mutate(across(everything(), ~ifelse(is.na(.), 0, .)))
    
    write.xlsx(summary_FPR, file = paste(dir, "Plots/Main/Excel", 
                                         paste0(paste("Type1error",
                                                      'grp_split_pc', grp_split_pc, 
                                                      'out_pc', out_pc, sep = '_'), ".xlsx"), sep = '/'))
    
    y_cap <- 0.1
    
    # Cap the FPR values
    summary_FPR <- summary_FPR %>%
      mutate(median_FPR_capped = pmin(median_FPR, y_cap))
    
    summary_FPR2 <- summary_FPR[, -3]
    
    summary_FPR2 <- summary_FPR2 %>%
      group_by(nsamples) %>% 
      arrange(median_FPR_capped, .by_group = TRUE)
    
    summary_FPR2 <- summary_FPR2 %>%
      group_by(nsamples) %>%
      mutate(Method_ordered = factor(Method, levels = Method)) %>%
      ungroup()
    
    FPR_list[[i]] <- ggplot(summary_FPR2, aes(x = factor(nsamples), y = median_FPR_capped, fill = Method_ordered)) +
      geom_bar(
        aes(group = interaction(nsamples, Method_ordered)),
        stat = "identity",
        position = position_dodge(width = 0.8)
      ) +
      geom_hline(yintercept = 0.05, color = "red", linetype = "dashed", size = 1) +
      labs(x = "No. of Samples", y = "Type-I Error", title = title) +
      scale_y_continuous(limits = c(0, y_cap)) +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "top",
        legend.background = element_rect(color = "black", size = 0.5),
        legend.box.background = element_rect(color = "black", size = 0.5)
      ) +
      scale_fill_manual(
        name = "Method",
        values = c(
          "Bartlett" = "#1b9e77",        
          "Brown-Forsythe " = "#d95f02", 
          "DGLM" = "#7570b3",            
          "DiffVar" = "#e7298a",         
          "iEVORA" = "#66a61e",          
          "JLSp" = "#e6ab02",            
          "JLSsc" = "#a6761d",           
          "pETM" = "#666666"             
        )
      )
    
    i <- i + 1
  }
}

batch_size <- 4
row <- 2
col <- 2
batches_FPR <- split(FPR_list, ceiling(seq_along(FPR_list) / batch_size))

# FPR Plots
for (i in seq_along(batches_FPR)) {
  output_file <- paste(dir, "Plots/Main", paste0(paste("Type1error", sep = '_'), ".png"), sep = '/')
  png(output_file, units = 'in', width = 20, height = 20, res = 300)
  print(do.call(ggarrange, c(batches_FPR[[i]], nrow = row, ncol = col)))
  dev.off()
}
