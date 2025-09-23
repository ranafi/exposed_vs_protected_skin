
## clear environment
rm(list=ls()) 


# install.packages("ggplot2")
library(ggplot2)
# install.packages("dplyr")
library(dplyr)
# install.packages("twosamples")
library(twosamples)

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 2A,B - RELATIVE AMPLITUDE


## clear environment
rm(list=ls()) 



data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)

data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)



filename <- paste0("results/RDS_files/overlap_ids.RDS")
cycling_in_both <- readRDS(filename)


# keep genes cycling in both

results_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_id %in% cycling_in_both,]
results_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_id %in% cycling_in_both,]

# check matching order - should be 0
# sum(results_EXPOSED$gene_id != results_PROTECTED$gene_id)



df1 <- data.frame(value = results_EXPOSED$rel_amp, group = "EXPOSED")
df2 <- data.frame(value = results_PROTECTED$rel_amp, group = "PROTECTED")

df_combined <- bind_rows(df1, df2)


p <- ggplot(df_combined, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50, color = "black") +
  labs(title = paste0("Relative Amplitude Distributions - Cycling in Both (",length(results_EXPOSED$rel_amp)," Transcripts)" ),
       x = "Relative Amplitude",
       y = "Frequency",
       fill = "Group") +
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1.5, "lines"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 13)

  )



p <- p + scale_fill_manual(values = c("EXPOSED" = "#FF0000", "PROTECTED" = "#0000FF"))


filename <- "plots/fig3/amp_dists_CYCLING_IN_BOTH.png"
ggsave(filename, plot = p, width = 8, height = 6, units = "in", dpi = 300)



########################################################################################
# WITHOUT LABELS    
p <- ggplot(df_combined, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50, color = "black") +
  labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
  )



p <- p + scale_fill_manual(values = c("EXPOSED" = "#FF0000", "PROTECTED" = "#0000FF"))


filename <- "plots/fig3_NOLABELS/amp_dists_CYCLING_IN_BOTH.png"
ggsave(filename, plot = p, width = 4, height = 4, units = "in", dpi = 300)
########################################################################################



amp_ratios <- results_EXPOSED$rel_amp / results_PROTECTED$rel_amp


df <- data.frame(value = amp_ratios)



  

p <- ggplot(df, aes(x = log2(value))) +
  geom_histogram(bins = 40, fill = "#D9A6F9", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", , size = 1.2) +
  labs(title = paste0( "Relative Amplitude Ratio - Cycling in Both (",length(results_EXPOSED$rel_amp)," Transcripts)"),
       x = "log2(EXPOSED amp / PROTECTED amp)",
       y = "Frequency")+
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1.5, "lines"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 20)
    
  )


filename <- "plots/fig3/amp_ratios_CYCLING_IN_BOTH.png"
ggsave(filename, plot = p, width = 8, height = 6, units = "in", dpi = 300)





########################################################################################
# WITHOUT LABELS   
p <- ggplot(df, aes(x = log2(value))) +
  geom_histogram(bins = 40, fill = "#D9A6F9", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1.2) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_blank()
  )



filename <- "plots/fig3_NOLABELS/amp_ratios_CYCLING_IN_BOTH.png"
ggsave(filename, plot = p, width = 4, height = 4, units = "in", dpi = 300)
########################################################################################







log2_test <- t.test(log2(amp_ratios))
print(log2_test)

# log2 dist mean
round(mean(log2(amp_ratios)),3)



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 2C ACROPHASE ROSE PLOT


## clear environment
rm(list=ls()) 


data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)

data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)


filename <- paste0("results/RDS_files/overlap_ids.RDS")
cycling_in_both <- readRDS(filename)


# keep genes cycling in both

results_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_id %in% cycling_in_both,]
results_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_id %in% cycling_in_both,]




peaktimes1 <- results_EXPOSED$phase_hour
peaktimes2 <- results_PROTECTED$phase_hour

data_label1 <- "EXPOSED"
data_label2 <- "PROTECTED"



num_bins <- 48
breaks <- seq(0, 24, length.out = num_bins + 1)

bins1 <- cut(peaktimes1, breaks = breaks, include.lowest = TRUE, right = FALSE)
bins2 <- cut(peaktimes2, breaks = breaks, include.lowest = TRUE, right = FALSE)

bin_counts1 <- as.data.frame(table(bins1))
colnames(bin_counts1) <- c("bin", "count")
bin_counts1$dataset <- data_label1

bin_counts2 <- as.data.frame(table(bins2))
colnames(bin_counts2) <- c("bin", "count")
bin_counts2$dataset <- data_label2



bin_midpoints <- (breaks[-1] + breaks[-length(breaks)]) / 2
bin_angles <- (bin_midpoints / 24) * 2 * pi

bins_df <- data.frame(
  bin = levels(bins1),
  bin_midpoint = bin_midpoints,
  bin_rad = bin_angles
)

data_counts1 <- merge(bins_df, bin_counts1, by = "bin")
data_counts2 <- merge(bins_df, bin_counts2, by = "bin")

data_counts_combined <- rbind(data_counts1, data_counts2)

fill_colors <- c("#FF0000", "#0000FF")
names(fill_colors) <- c(data_label1, data_label2)


hour_labels <- 0:23
break_angles <- (hour_labels / 24) * 2 * pi





rose_plot <- ggplot(data_counts_combined, aes(x = bin_rad, y = count, fill = dataset)) +
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
    title = paste0("Peak Time Comparison - Cycling in Both (",length(peaktimes1)," Transcripts)")
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



ggsave(filename = "plots/fig3/rose_plot_CYCLING_IN_BOTH.png", plot = rose_plot, width = 8, height = 6, dpi = 300)








########################################################################################
# WITHOUT LABELS  
rose_plot <- ggplot(data_counts_combined, aes(x = bin_rad, y = count, fill = dataset)) +
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



ggsave(filename = "plots/fig3_NOLABELS/rose_plot_CYCLING_IN_BOTH.png", plot = rose_plot, width = 8, height = 6, dpi = 300)
########################################################################################







# convert to radians
phase_EXPOSED_rads <- results_EXPOSED$phase_hour / 24 * 2 * pi
phase_PROTECTED_rads <- results_PROTECTED$phase_hour / 24 * 2 * pi


# kuiper test
# https://search.r-project.org/CRAN/refmans/twosamples/html/kuiper_test.html


kuiper_test(phase_EXPOSED_rads,phase_PROTECTED_rads) # radians
kuiper_test(results_EXPOSED$phase_hour,results_PROTECTED$phase_hour) # hours
# documentation doesn't specifiy which unit to use, but should be same for both





##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 2D - ACROPHASE DIFFERENCE PLOT


## clear environment
rm(list=ls()) 


data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)

data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)



filename <- paste0("results/RDS_files/overlap_ids.RDS")
cycling_in_both <- readRDS(filename)


# keep genes cycling in both

results_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_id %in% cycling_in_both,]
results_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_id %in% cycling_in_both,]





phase_difference <- results_EXPOSED$phase_hour - results_PROTECTED$phase_hour



phase_difference[phase_difference > 12] <- phase_difference[phase_difference > 12] - 24
phase_difference[phase_difference < -12] <- phase_difference[phase_difference < -12] + 24





df <- data.frame(value = phase_difference)





p <- ggplot(df, aes(x = value)) +
  geom_histogram(bins = 40, fill = "#D9A6F9", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", , size = 1.2) +
  labs(title = paste0( "Phase Difference - (EXPOSED - PROTECTED) - Cycling in Both (",length(phase_difference)," Transcripts)"),
       x = "Phase Difference - (EXPOSED - PROTECTED)",
       y = "Frequency")+
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1.5, "lines"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 14)
    
  )


filename <- paste0("plots/fig3/phase_difference_CYCLING_IN_BOTH.png")
ggsave(filename, plot = p, width = 8, height = 6, units = "in", dpi = 300)


  
########################################################################################
# WITHOUT LABELS  
p <- ggplot(df, aes(x = value)) +
  geom_histogram(bins = 40, fill = "#D9A6F9", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1.2) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
  )




filename <- paste0("plots/fig3_NOLABELS/phase_difference_CYCLING_IN_BOTH.png")
ggsave(filename, plot = p, width = 4, height = 4, units = "in", dpi = 300)
########################################################################################

  

  

## calculate circular mean of differences

diffs_hours <- results_EXPOSED$phase_hour  - results_PROTECTED$phase_hour 
diffs_radians <- 2 * pi * diffs_hours / 24

x <- mean(cos(diffs_radians))
y <- mean(sin(diffs_radians))

circular_mean_radians <- atan2(y, x)

circular_mean_hours <- (circular_mean_radians / (2 * pi)) * 24

if (circular_mean_hours < 0) {
  circular_mean_hours <- circular_mean_hours + 24
}


# 25th and 75th percentiles
q25 <- quantile(diffs_hours, 0.25)
q75 <- quantile(diffs_hours, 0.75)



