
## clear environment
rm(list=ls()) 


# install.packages("ggplot2")
library(ggplot2)
# install.packages("dplyr")
library(dplyr)
# install.packages("pheatmap")
library(pheatmap)
# install.packages("gridExtra")
library(gridExtra)
library(grid)





##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

# FIG 3C - HEATMAPS

## clear environment
rm(list=ls()) 


DR_ids <- readRDS("results/RDS_files/DR_ids.RDS")


data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)
results_PROTECTED_DR <- results_PROTECTED[results_PROTECTED$gene_id %in% DR_ids,]

# acrophase order based on PROTECTED results
acrophase_order <- order(results_PROTECTED_DR$phase_hour)
row_order <- DR_ids[acrophase_order]





times_EXPOSED <- readRDS("processed_data/times_EXPOSED.RDS")
ids_EXPOSED   <- readRDS("processed_data/ids_EXPOSED.RDS")
id_factors_EXPOSED <- readRDS("processed_data/id_factors_EXPOSED.RDS")
data_EXPOSED  <- readRDS("processed_data/gene_data_EXPOSED.RDS")
data_EXPOSED <- data_EXPOSED[rownames(data_EXPOSED) %in% DR_ids, ]
data_EXPOSED <- data_EXPOSED[row_order, ]

times_PROTECTED <- readRDS("processed_data/times_PROTECTED.RDS")
ids_PROTECTED   <- readRDS("processed_data/ids_PROTECTED.RDS")
id_factors_PROTECTED <- readRDS("processed_data/id_factors_PROTECTED.RDS")
data_PROTECTED  <- readRDS("processed_data/gene_data_PROTECTED.RDS")

data_PROTECTED <- data_PROTECTED[rownames(data_PROTECTED) %in% DR_ids, ]
data_PROTECTED <- data_PROTECTED[row_order, ]




data_combined <- cbind(data_PROTECTED, data_EXPOSED)

data_combined_labels <- c(rep("PROTECTED",length(times_PROTECTED)), rep("EXPOSED",length(times_EXPOSED)))


# convert to z-scores
for (i in 1:dim(data_combined)[1]) {
  row_vals <- data_combined[i, ]
  row_mean <- mean(row_vals)
  row_sd   <- sd(row_vals)
  
  if (row_sd > 0) {
    data_combined[i, ] <- (row_vals - row_mean) / row_sd
  } else {
    data_combined[i, ] <- 0
  }
}



data_PROTECTED_Z <- data_combined[, data_combined_labels == "PROTECTED"]
data_EXPOSED_Z   <- data_combined[, data_combined_labels == "EXPOSED"]




###########################################################################################
# PROTECTED heatmap

col_order_PROTECTED <- order(times_PROTECTED) # ordering cols by sample time

ordered_data_PROTECTED_Z <- data_PROTECTED_Z[, col_order_PROTECTED]
ordered_times_PROTECTED <- times_PROTECTED[col_order_PROTECTED]

unique_times_PROTECTED <- sort(unique(times_PROTECTED))


average_z_PROTECTED <- matrix(
  NA, 
  nrow = nrow(ordered_data_PROTECTED_Z), 
  ncol = length(unique_times_PROTECTED),
  dimnames = list(rownames(ordered_data_PROTECTED_Z), unique_times_PROTECTED)
)


for (i in 1: dim(ordered_data_PROTECTED_Z)[1]) {
  for (j in 1:length(unique_times_PROTECTED)) {
    unique_time <- unique_times_PROTECTED[j]
    average_z_PROTECTED[i, j] <- mean(ordered_data_PROTECTED_Z[i, ordered_times_PROTECTED == unique_time])
  }
}



zlim <- 3
colors <- colorRampPalette(c("navy", "white", "firebrick"))(50)
breaks_z  <- seq(-zlim, zlim, length.out = length(colors) + 1)


annotation_colors <- list(
  Time = c('0' = "grey", '6' = "lightblue", '12' = "blue", '18' = "beige")
)


times_factor_avg <- factor(colnames(average_z_PROTECTED), levels = colnames(average_z_PROTECTED))
avg_annotation_col <- data.frame(Time = times_factor_avg)
rownames(avg_annotation_col) <- colnames(average_z_PROTECTED)


avgerage_heatmap_PROTECTED <- pheatmap(
  average_z_PROTECTED,
  scale = "none",
  cluster_cols = FALSE, 
  cluster_rows = FALSE,
  annotation_col = avg_annotation_col,
  annotation_colors = annotation_colors,
  color = colors,
  breaks = breaks_z,
  fontsize_row = 1,
  fontsize_col = 8,
  legend = TRUE,
  silent = TRUE
)


times_factor_full <- factor(ordered_times_PROTECTED, levels = sort(unique(ordered_times_PROTECTED)))
full_annotation_col <- data.frame(Time = times_factor_full)
rownames(full_annotation_col) <- colnames(ordered_data_PROTECTED_Z)



full_heatmap_PROTECTED <- pheatmap(
  ordered_data_PROTECTED_Z,
  scale = "none",
  cluster_cols = FALSE, 
  cluster_rows = FALSE,
  annotation_col = full_annotation_col,
  annotation_colors = annotation_colors,
  color = colors,
  breaks = breaks_z,
  fontsize_row = 1,
  fontsize_col = 4,
  legend = TRUE,
  silent = TRUE
)


combined_plot_PROTECTED <- grid.arrange(
  grobs = list(avgerage_heatmap_PROTECTED$gtable, full_heatmap_PROTECTED$gtable),
  ncol = 2,
  widths = c(1, 3),
  padding = unit(1, "line")
)


png(
  filename = "plots/fig4/heatmap_PROTECTED.png",
  width = 15,
  height = 8,
  units = "in",
  res = 300
)
grid.draw(combined_plot_PROTECTED)
dev.off()




###########################################################################################
# EXPOSED heatmap

col_order_EXPOSED <- order(times_EXPOSED) # ordering cols by sample time

ordered_data_EXPOSED_Z <- data_EXPOSED_Z[, col_order_EXPOSED]
ordered_times_EXPOSED <- times_EXPOSED[col_order_EXPOSED]

unique_times_EXPOSED <- sort(unique(times_EXPOSED))


average_z_EXPOSED <- matrix(
  NA, 
  nrow = nrow(ordered_data_EXPOSED_Z), 
  ncol = length(unique_times_EXPOSED),
  dimnames = list(rownames(ordered_data_EXPOSED_Z), unique_times_EXPOSED)
)


for (i in 1: dim(ordered_data_EXPOSED_Z)[1]) {
  for (j in 1:length(unique_times_EXPOSED)) {
    unique_time <- unique_times_EXPOSED[j]
    average_z_EXPOSED[i, j] <- mean(ordered_data_EXPOSED_Z[i, ordered_times_EXPOSED == unique_time])
  }
}



zlim <- 3
colors <- colorRampPalette(c("navy", "white", "firebrick"))(50)
breaks_z  <- seq(-zlim, zlim, length.out = length(colors) + 1)


annotation_colors <- list(
  Time = c('0' = "grey", '6' = "lightblue", '12' = "blue", '18' = "beige")
)


times_factor_avg <- factor(colnames(average_z_EXPOSED), levels = colnames(average_z_EXPOSED))
avg_annotation_col <- data.frame(Time = times_factor_avg)
rownames(avg_annotation_col) <- colnames(average_z_EXPOSED)


avgerage_heatmap_EXPOSED <- pheatmap(
  average_z_EXPOSED,
  scale = "none",
  cluster_cols = FALSE, 
  cluster_rows = FALSE,
  annotation_col = avg_annotation_col,
  annotation_colors = annotation_colors,
  color = colors,
  breaks = breaks_z,
  fontsize_row = 1,
  fontsize_col = 8,
  legend = TRUE,
  silent = TRUE
)


times_factor_full <- factor(ordered_times_EXPOSED, levels = sort(unique(ordered_times_EXPOSED)))
full_annotation_col <- data.frame(Time = times_factor_full)
rownames(full_annotation_col) <- colnames(ordered_data_EXPOSED_Z)



full_heatmap_EXPOSED <- pheatmap(
  ordered_data_EXPOSED_Z,
  scale = "none",
  cluster_cols = FALSE, 
  cluster_rows = FALSE,
  annotation_col = full_annotation_col,
  annotation_colors = annotation_colors,
  color = colors,
  breaks = breaks_z,
  fontsize_row = 1,
  fontsize_col = 4,
  legend = TRUE,
  silent = TRUE
)


combined_plot_EXPOSED <- grid.arrange(
  grobs = list(avgerage_heatmap_EXPOSED$gtable, full_heatmap_EXPOSED$gtable),
  ncol = 2,
  widths = c(1, 3),
  padding = unit(1, "line")
)


png(
  filename = "plots/fig4/heatmap_EXPOSED.png",
  width = 15,
  height = 8,
  units = "in",
  res = 300
)
grid.draw(combined_plot_EXPOSED)
dev.off()



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 3D - RELATIVE AMPLITUDE 


## clear environment
rm(list=ls()) 

data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)

data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)



filename <- paste0("results/RDS_files/DR_ids.RDS")
DR_ids <- readRDS(filename)




# keep differentially rhythmic genes

results_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_id %in% DR_ids,]
results_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_id %in% DR_ids,]

# check matching order - should be 0
# sum(results_EXPOSED$gene_id != results_PROTECTED$gene_id)



df1 <- data.frame(value = results_EXPOSED$rel_amp, group = "EXPOSED")
df2 <- data.frame(value = results_PROTECTED$rel_amp, group = "PROTECTED")

df_combined <- bind_rows(df1, df2)





p <- ggplot(df_combined, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50, color = "black") +
  labs(title = paste0("Relative Amplitude Distributions - Differentially Rhythmic (",length(results_EXPOSED$rel_amp)," Transcripts)" ),
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
    plot.title = element_text(size = 12)
    
  )


p <- p + scale_fill_manual(values = c("EXPOSED" = "#FF0000", "PROTECTED" = "#0000FF"))


filename <- "plots/fig4/amp_dists_DIFFERENTIALLY_RHYTHMIC.png"
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


filename <- "plots/fig4_NOLABELS/amp_dists_DIFFERENTIALLY_RHYTHMIC.png"
ggsave(filename, plot = p, width = 4, height = 4, units = "in", dpi = 300)
########################################################################################



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 3E - ACROPHASE ROSE PLOT


## clear environment
rm(list=ls()) 


data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)

data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)


filename <- paste0("results/RDS_files/DR_ids.RDS")
DR_ids <- readRDS(filename)


# keep differentially rhythmic genes

results_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_id %in% DR_ids,]
results_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_id %in% DR_ids,]



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
    title = paste0("Peak Time Comparison - Differentially Rhythmic (",length(peaktimes1)," Transcripts)")
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


ggsave(filename = "plots/fig4/rose_plot_DIFFERENTIALLY_RHYTHMIC.png", plot = rose_plot, width = 8, height = 6, dpi = 300)


  
  
  
  
  



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
    axis.text.x = element_text(size = 16),
    plot.title = element_blank(),
    legend.position = "none"
  )



ggsave(filename = "plots/fig4_NOLABELS/rose_plot_DIFFERENTIALLY_RHYTHMIC.png", plot = rose_plot, width = 8, height = 6, dpi = 300)
########################################################################################

