
## clear environment
rm(list=ls())


# install.packages("ggplot2")
library(ggplot2)
# install.packages("patchwork")
library(patchwork)
# install.packages("readr")
library(readr)
# install.packages("purrr")
library(purrr)


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 1A - CORE CLOCK TIMESERIES, JUST PROTECTED


## clear environment
rm(list=ls())


make_timeseries_plots_PROTECTED_ONLY <- function(gene_ids_to_plot, text_labels=TRUE, output_path="plots/") {
  
  
  data_label_PROTECTED <- "PROTECTED"
  times_PROTECTED <- readRDS(paste("processed_data/times_", data_label_PROTECTED, ".RDS", sep=""))
  ids_PROTECTED <- readRDS(paste("processed_data/ids_", data_label_PROTECTED, ".RDS", sep=""))
  id_factors_PROTECTED <- readRDS(paste("processed_data/id_factors_", data_label_PROTECTED, ".RDS", sep=""))
  data_PROTECTED <- readRDS(paste("processed_data/gene_data_", data_label_PROTECTED, ".RDS", sep=""))
  
  
  filename <- paste0("results/RDS_files/id_to_name.RDS")
  id_to_name <- readRDS(filename)
  
  for (gene_to_find in gene_ids_to_plot){
    
    
    ########################################################################################################################################################################
    # PROTECTED
    
    expression_PROTECTED <- unlist(as.vector(data_PROTECTED[gene_to_find,]))
    
    trig_model_PROTECTED <- lm(expression_PROTECTED ~ cos(times_PROTECTED * pi / 12) + sin(times_PROTECTED * pi / 12)  + id_factors_PROTECTED) 
    
    cos_fit_PROTECTED <- unname(trig_model_PROTECTED$coefficients[2])
    sin_fit_PROTECTED <- unname(trig_model_PROTECTED$coefficients[3])
    
    amp_PROTECTED <- sqrt(cos_fit_PROTECTED**2 + sin_fit_PROTECTED**2)
    phase_rad_fit_PROTECTED <- atan2(sin_fit_PROTECTED, cos_fit_PROTECTED)
    phase_hour_fit_PROTECTED <- unname(phase_rad_fit_PROTECTED * (12/pi))
    
    if (phase_hour_fit_PROTECTED < 0) {
      phase_hour_fit_PROTECTED <- 24 + phase_hour_fit_PROTECTED
    }
    
    c_PROTECTED <- trig_model_PROTECTED$coefficients[1]
    subject_factors_PROTECTED <- trig_model_PROTECTED$coefficients[4:length(trig_model_PROTECTED$coefficients)]
    subject_mesors_PROTECTED <- subject_factors_PROTECTED + c_PROTECTED
    subject_mesors_PROTECTED <- append(c_PROTECTED, subject_mesors_PROTECTED)
    average_mesor_PROTECTED <- mean(subject_mesors_PROTECTED)
    
    adjusted_expression_PROTECTED <- expression_PROTECTED
    
    for (sub in 1:20){
      adjusted_expression_PROTECTED[ids_PROTECTED == sub] <- adjusted_expression_PROTECTED[ids_PROTECTED == sub] - subject_mesors_PROTECTED[sub]
    }
    
    y_min <- min(adjusted_expression_PROTECTED)
    y_max <- max(adjusted_expression_PROTECTED)
    
    to_plot_PROTECTED <- data.frame(times = times_PROTECTED, adjusted_expression = adjusted_expression_PROTECTED)
    
    dup_to_plot_PROTECTED <- to_plot_PROTECTED
    dup_to_plot_PROTECTED$times <- dup_to_plot_PROTECTED$times+24
    
    fit_x <- seq(0, 48, by = 1)
    fit_line_PROTECTED <- data.frame(
      fit_x = fit_x,
      fit_y = cos_fit_PROTECTED * cos(fit_x * (2 * pi / 24)) + sin_fit_PROTECTED * sin(fit_x * (2 * pi / 24)) 
    )
    
    fit_line_PROTECTED_solid <- subset(fit_line_PROTECTED, fit_x <= 18)
    fit_line_PROTECTED_dotted <- subset(fit_line_PROTECTED, fit_x > 18)
    

    ########################################################################################################################################################################
    # plot
    
    
    if(text_labels){
      
      
      
      plot_PROTECTED <- ggplot() +
        geom_point(data = to_plot_PROTECTED, aes(x = times, y = adjusted_expression), color = "black" , size   = 2.5) +
        geom_point(data = dup_to_plot_PROTECTED, aes(x = times, y = adjusted_expression), color = "gray", alpha = 0.9 , size   = 2.5) +
        geom_line(data = fit_line_PROTECTED_solid, aes(x = fit_x, y = fit_y), color = "blue", linetype = "solid", size   = 1.2) +
        geom_line(data = fit_line_PROTECTED_dotted, aes(x = fit_x, y = fit_y), color = "blue", linetype = "dotted", size   = 1.2) +
        scale_x_continuous(limits = c(0, 48), breaks = seq(0, 48, by = 12)) +
        scale_y_continuous(limits = c(y_min - 0.5, y_max + 0.5)) +
        labs(x = "Time (hours)", y = "MESOR-Corrected Expression", title = "PROTECTED")+
        theme_minimal() +
        theme(
          axis.text.y  = element_text(size = 16),
          axis.ticks.y = element_blank(),
          axis.text.x  = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18)
        ) +
        geom_vline(xintercept = phase_hour_fit_PROTECTED, linetype = "dashed", color = "red") +
        annotate("text", 
                 x = phase_hour_fit_PROTECTED, 
                 y = max(to_plot_PROTECTED$adjusted_expression, na.rm = TRUE), 
                 label = round(phase_hour_fit_PROTECTED, 2), 
                 vjust = -0.5, 
                 hjust = -0.1, 
                 color = "red")
      
      
      gene_name <- id_to_name$gene_name[id_to_name$gene_id==gene_to_find]
      
      if(gene_name != "na"){
        gene_name_title <- paste0(gene_to_find," (",gene_name,")")
      }else{
        gene_name_title <- gene_to_find
        gene_name <- ""
      }
      
      
      
      plot_PROTECTED <- plot_PROTECTED +
        plot_annotation(
          title = gene_name_title,
          theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
        )      
      
      
      ggsave(paste0(output_path,"timeseries_",gene_name,"_", gene_to_find, ".png"), plot_PROTECTED, width = 7.5, height = 4.5, units = "in", dpi=300)
      
      
      
      
      
      
      
    }else{
      
      plot_PROTECTED <- ggplot() +
        geom_point(data = to_plot_PROTECTED, aes(x = times, y = adjusted_expression), color = "black" , size   = 2.5) +
        geom_point(data = dup_to_plot_PROTECTED, aes(x = times, y = adjusted_expression), color = "gray", alpha = 0.9 , size   = 2.5) +
        geom_line(data = fit_line_PROTECTED_solid, aes(x = fit_x, y = fit_y), color = "blue", linetype = "solid", size   = 1.2) +
        geom_line(data = fit_line_PROTECTED_dotted, aes(x = fit_x, y = fit_y), color = "blue", linetype = "dotted", size   = 1.2) +
        scale_x_continuous(limits = c(0, 48), breaks = seq(0, 48, by = 12)) +
        scale_y_continuous(limits = c(y_min - 0.5, y_max + 0.5)) +
        labs(x = NULL, y = NULL, title = NULL) +
        theme_minimal() +
        theme(
          axis.text.y  = element_text(size = 22),
          axis.ticks.y = element_blank(),
          axis.text.x  = element_text(size = 22)
        )    
      
      
      gene_name <- id_to_name$gene_name[id_to_name$gene_id==gene_to_find]
      
      if(gene_name != "na"){
        gene_name_title <- paste0(gene_to_find," (",gene_name,")")
      }else{
        gene_name_title <- gene_to_find
        gene_name <- ""
      }
      
      
      
      plot_PROTECTED <- plot_PROTECTED +
        plot_annotation(
          title = gene_name_title,
          theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
        )      
      
      
      
      
      ggsave(paste0(output_path,"timeseries_",gene_name,"_", gene_to_find, ".png"), plot_PROTECTED, width = 7.5, height = 4.5, units = "in", dpi=300)
      
      
      
    }
    
  }
  
  
}


##########################################################################################################################################
## core clock genes

coreclock_genes <- c( "NR1D2", "PER1", "PER2", "CRY2", "RORC", "BHLHE40", "PER3", "CRY1", "NPAS2", "NFIL3", "CIART", "TEF", "BMAL1", "CLOCK", "NR1D1", "DBP", "HLF", "RORA", "BHLHE41", "CSNK1D")
# 20 genes to check (not all are cycling)


filename <- paste0("results/RDS_files/first_pass_results_PROTECTED.RDS")
results_PROTECTED <- readRDS(filename)


rel_amp_thres <- 0.2
bhq_thres <- 0.05


core_clock_ids <- results_PROTECTED[results_PROTECTED$bh_q <= bhq_thres & results_PROTECTED$rel_amp >= rel_amp_thres & results_PROTECTED$gene_name %in% coreclock_genes,]$gene_id


make_timeseries_plots_PROTECTED_ONLY(core_clock_ids, output_path="plots/fig1_NOLABELS/timeseries_coreclock_PROTECTED/", text_labels=FALSE)
make_timeseries_plots_PROTECTED_ONLY(core_clock_ids, output_path="plots/fig1/timeseries_coreclock_PROTECTED/", text_labels=TRUE)










##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 1B - ROSE PLOT, JUST PROTECTED


## clear environment
rm(list=ls()) 



data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)




rel_amp_thres <- 0.2
bhq_thres <- 0.05


sig_results_PROTECTED <- results_PROTECTED[results_PROTECTED$bh_q <= bhq_thres & results_PROTECTED$rel_amp >= rel_amp_thres,]


peaktimes <- sig_results_PROTECTED$phase_hour




num_bins <- 48
breaks <- seq(0, 24, length.out = num_bins + 1)

bins <- cut(peaktimes, breaks = breaks, include.lowest = TRUE, right = FALSE)



bin_counts <- as.data.frame(table(bins))
colnames(bin_counts) <- c("bin", "count")
bin_counts$dataset <- data_label



bin_midpoints <- (breaks[-1] + breaks[-length(breaks)]) / 2
bin_angles <- (bin_midpoints / 24) * 2 * pi

bins_df <- data.frame(
  bin = levels(bins),
  bin_midpoint = bin_midpoints,
  bin_rad = bin_angles
)

data_counts <- merge(bins_df, bin_counts, by = "bin")


fill_colors <- c("#0000FF")
names(fill_colors) <- c( data_label)


hour_labels <- 0:23
break_angles <- (hour_labels / 24) * 2 * pi



rose_plot <- ggplot(data_counts, aes(x = bin_rad, y = count, fill = dataset)) +
  geom_bar(
    stat = "identity",
    color = NA,
    width = (2 * pi) / num_bins,
    alpha = 0.5,
    position = "identity"
  ) +
  coord_polar(start = 4*pi / 2) + 
  scale_x_continuous(
    limits = c(0, 2*pi + 1e-5),
    breaks = break_angles,
    labels = paste0(hour_labels, "h")
  ) +
  scale_fill_manual(values = fill_colors) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Dataset",
    title = paste0("Peak Times - PROTECTED")
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 14, hjust = 0.5)
  )



ggsave(filename = "plots/fig1/rose_plot_PROTECTED.png", plot = rose_plot, width = 8, height = 6, dpi = 300)








########################################################################################
# WITHOUT LABELS  
rose_plot <- ggplot(data_counts, aes(x = bin_rad, y = count, fill = dataset)) +
  geom_bar(
    stat = "identity",
    color = NA,
    width = (2 * pi) / num_bins,
    alpha = 0.5,
    position = "identity"
  ) +
  coord_polar(start = 4*pi / 2) + 
  scale_x_continuous(
    limits = c(0, 2*pi + 1e-5),
    breaks = break_angles,
    labels = paste0(hour_labels, "h")
  ) +
  scale_fill_manual(values = fill_colors) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL,
    title = NULL
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 16),  # keep x tick labels (e.g., 0h, 6h...)
    plot.title = element_blank(),
    legend.position = "none"
  )



ggsave(filename = "plots/fig1_NOLABELS/rose_plot_PROTECTED.png", plot = rose_plot, width = 8, height = 6, dpi = 300)
########################################################################################


