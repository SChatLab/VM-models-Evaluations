packages <- c('ggplot2', 'stringr', "reshape2",
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
  var.effect.i <- as.numeric(files.i[3]) / 0.64
  grp.split.pc.i <- as.numeric(files.i[4])
  pDM.i <- as.numeric(files.i[5])
  outpc.i <- as.numeric(strsplit(files.i[6], ".RD")[[1]][1])
  
  if(nrow(findf) == 0){
    findf <- data.frame(restmp,
                        nsamples = nsamples.i,
                        var.effect = var.effect.i,
                        grp.split.pc = grp.split.pc.i,
                        pDM = pDM.i,
                        outpc = outpc.i)
  }else{
    findf.sub <- data.frame(restmp,
                            nsamples = nsamples.i,
                            var.effect = var.effect.i,
                            grp.split.pc = grp.split.pc.i,
                            pDM = pDM.i,
                            outpc = outpc.i)
    findf <- rbind(findf, findf.sub)
  }
}

## Plots
Power_list <- list()
MCC_list <- list()
AUROC_list <- list()
F1_list <- list()

var_eff <- 1.5
grp_split <- c(50, 70)
out_pc_list <- c(0, 0.1)
p_DM <- 0.1
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
    
    title <- str_wrap(paste0(LETTERS[i], ') ' , split, ', ', out), 40)
    
    # Plot Power
    summary_Power <- res %>%
      group_by(nsamples, Method) %>%
      summarise(median_power = median(Power), .groups = 'drop')
    
    summary_Power$Method <- str_remove_all(summary_Power$Method, regex("-?Test", ignore_case = TRUE))
    
    summary_Power <- summary_Power %>%
      mutate(across(everything(), ~ifelse(is.na(.), 0, .)))
    
    Power_list[[i]] <- ggplot(summary_Power, aes(x = nsamples, y = median_power, color = Method, group = Method)) +
      geom_line(size = 1) +
      geom_point(shape = 16, size = 3) +
      labs(
        title = title,
        x = "No. of Samples",
        y = "Power",
        color = "Methods"
      ) +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.background = element_rect(color = "black", size = 0.5, fill = "white"),
        text = element_text(size = 12),
        legend.position = "top",
        plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16)
      ) +
      scale_x_continuous(limits = c(50, 150), breaks = c(50, 75, 100, 125, 150)) +  
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_manual(
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
    
    # Plot MCC
    summary_MCC <- res %>%
      group_by(nsamples, Method) %>%
      summarise(median_MCC = median(MCC), .groups = 'drop')
    
    summary_MCC$Method <- str_remove_all(summary_MCC$Method, regex("-?Test", ignore_case = TRUE))
    
    summary_MCC <- summary_MCC %>%
      mutate(across(everything(), ~ifelse(is.na(.), 0, .)))
    
    MCC_list[[i]] <- ggplot(summary_MCC, aes(x = nsamples, y = median_MCC, color = Method, group = Method)) +
      geom_line(size = 1) +
      geom_point(shape = 16, size = 3) +
      labs(
        title = title,
        x = "No. of Samples",
        y = "MCC Score",
        color = "Methods"
      ) +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.background = element_rect(color = "black", size = 0.5, fill = "white"),
        text = element_text(size = 12),
        legend.position = "top",
        plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16)
      ) +
      scale_x_continuous(limits = c(50, 150), breaks = c(50, 75, 100, 125, 150)) +  
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_manual(
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
    
    # Plot F1
    summary_F1 <- res %>%
      group_by(nsamples, Method) %>%
      summarise(median_F1 = median(F1), .groups = 'drop')
    
    summary_F1$Method <- str_remove_all(summary_F1$Method, regex("-?Test", ignore_case = TRUE))
    
    summary_F1 <- summary_F1 %>%
      mutate(across(everything(), ~ifelse(is.na(.), 0, .)))
    
    F1_list[[i]] <- ggplot(summary_F1, aes(x = nsamples, y = median_F1, color = Method, group = Method)) +
      geom_line(size = 1) +
      geom_point(shape = 16, size = 3) +
      labs(
        title = title,
        x = "No. of Samples",
        y = "F1 Score",
        color = "Methods"
      ) +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.background = element_rect(color = "black", size = 0.5, fill = "white"),
        text = element_text(size = 12),
        legend.position = "top",
        plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16)
      ) +
      scale_x_continuous(limits = c(50, 150), breaks = c(50, 75, 100, 125, 150)) +  
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_manual(
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
    
    # Plot AUROC
    summary_ROC <- res %>%
      group_by(nsamples, Method) %>%
      summarise(median_ROC = median(ROC), .groups = 'drop')
    
    summary_ROC$Method <- str_remove_all(summary_ROC$Method, regex("-?Test", ignore_case = TRUE))
    
    summary_ROC <- summary_ROC %>%
      mutate(across(everything(), ~ifelse(is.na(.), 0, .)))
    
    AUROC_list[[i]] <- ggplot(summary_ROC, aes(x = nsamples, y = median_ROC, color = Method, group = Method)) +
      geom_line(size = 1) +
      geom_point(shape = 16, size = 3) +
      labs(
        title = title,
        x = "No. of Samples",
        y = "AuROC Score",
        color = "Methods"
      ) +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.background = element_rect(color = "black", size = 0.5, fill = "white"),
        text = element_text(size = 12),
        legend.position = "top",
        plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16)
      ) +
      scale_x_continuous(limits = c(50, 150), breaks = c(50, 75, 100, 125, 150)) +  
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_manual(
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
batches_power <- split(Power_list, ceiling(seq_along(Power_list) / batch_size))
batches_mcc <- split(MCC_list, ceiling(seq_along(MCC_list) / batch_size))
batches_auc <- split(AUROC_list, ceiling(seq_along(AUROC_list) / batch_size))
batches_f1 <- split(F1_list, ceiling(seq_along(F1_list) / batch_size))


# Power Plots
for (i in seq_along(batches_power)) {
  output_file <- paste(dir, "Plots/Main", paste0(paste("Power_vareff", var_eff, "pDM", p_DM, sep = '_'), ".png"), sep = '/')
  png(output_file, units = 'in', width = 20, height = 20, res = 300)
  print(do.call(ggarrange, c(batches_power[[i]], nrow = row, ncol = col)))
  dev.off()
}

# MCC Plots
for (i in seq_along(batches_mcc)) {
  output_file <- paste(dir, "Plots/Main", paste0(paste("MCC_vareff", var_eff, "pDM", p_DM, sep = '_'), ".png"), sep = '/')
  png(output_file, units = 'in', width = 20, height = 20, res = 300)
  print(do.call(ggarrange, c(batches_mcc[[i]], nrow = row, ncol = col)))
  dev.off()
}

# F1 Plots
for (i in seq_along(batches_f1)) {
  output_file <- paste(dir, "Plots/Main", paste0(paste("F1_vareff", var_eff, "pDM", p_DM, sep = '_'), ".png"), sep = '/')
  png(output_file, units = 'in', width = 20, height = 20, res = 300)
  print(do.call(ggarrange, c(batches_f1[[i]], nrow = row, ncol = col)))
  dev.off()
}

# AUROC Plots
for (i in seq_along(batches_auc)) {
  output_file <- paste(dir, "Plots/Main", paste0(paste("AUROC_vareff", var_eff, 'pDM', p_DM, sep = '_'), ".png"), sep = '/')
  png(output_file, units = 'in', width = 20, height = 20, res = 300)
  print(do.call(ggarrange, c(batches_auc[[i]], nrow = row, ncol = col)))
  dev.off()
}
