
## clear environment
rm(list=ls()) 

# install.packages("ggplot2")
library(ggplot2)
# install.packages("dplyr")
library(dplyr)
# # install.packages("eulerr")
# library(eulerr)
# install.packages("ggrepel")
library(ggrepel)
# install.packages("ggnewscale")
library(ggnewscale)
# install.packages("scales")
library(scales)
# install.packages("readr")
library(readr)
# install.packages("purrr")
library(purrr)

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG S1B,C - ABSOLUTE AMPLITUDE


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



df1 <- data.frame(value = results_EXPOSED$abs_amp, group = "EXPOSED")
df2 <- data.frame(value = results_PROTECTED$abs_amp, group = "PROTECTED")

df_combined <- bind_rows(df1, df2)








p <- ggplot(df_combined, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50, color = "black") +
  scale_x_log10() +  # apply log scale transformation
  labs(title = paste0("Absolute Amplitude Distributions - Cycling in Both (", 
                      length(results_EXPOSED$rel_amp)," Transcripts)"),
       x = "Absolute Amplitude (log10 scale)",
       y = "Frequency",
       fill = "Group") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1.5, "lines"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 13)
  ) +
  scale_fill_manual(values = c("EXPOSED" = "#FF0000", "PROTECTED" = "#0000FF"))




p <- p + scale_fill_manual(values = c("EXPOSED" = "#FF0000", "PROTECTED" = "#0000FF"))


filename <- "plots/figS1/absolute_amp_dists_CYCLING_IN_BOTH.png"
ggsave(filename, plot = p, width = 8, height = 6, units = "in", dpi = 300)







########################################################################################
# WITHOUT LABELS   
p <- ggplot(df_combined, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50, color = "black") +
  scale_x_log10() +  # apply log scale transformation
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


filename <- "plots/figS1_NOLABELS/absolute_amp_dists_CYCLING_IN_BOTH.png"
ggsave(filename, plot = p, width = 4, height = 4, units = "in", dpi = 300)
########################################################################################




amp_ratios <- results_EXPOSED$abs_amp / results_PROTECTED$abs_amp


df <- data.frame(value = amp_ratios)




p <- ggplot(df, aes(x = log2(value))) +
  geom_histogram(bins = 40, fill = "#D9A6F9", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", , size = 1.2) +
  labs(title = paste0( "Absolute Amplitude Ratio - Cycling in Both (",length(results_EXPOSED$rel_amp)," Transcripts)"),
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


filename <- "plots/figS1/absolute_amp_ratios_CYCLING_IN_BOTH.png"
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



filename <- "plots/figS1_NOLABELS/absolute_amp_ratios_CYCLING_IN_BOTH.png"
ggsave(filename, plot = p, width = 4, height = 4, units = "in", dpi = 300)
########################################################################################


log2_test <- t.test(log2(amp_ratios))
print(log2_test)


# log2 dist mean
round(mean(log2(amp_ratios)),3)


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG S1D - ROSE PLOT, RHYTHMIC TRANSCRIPTS IN EACH SITE



## clear environment
rm(list=ls()) 


data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)

data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)




rel_amp_thres <- 0.2
bhq_thres <- 0.05


sig_results_EXPOSED <- results_EXPOSED[results_EXPOSED$bh_q <= bhq_thres & results_EXPOSED$rel_amp >= rel_amp_thres,]
sig_results_PROTECTED <- results_PROTECTED[results_PROTECTED$bh_q <= bhq_thres & results_PROTECTED$rel_amp >= rel_amp_thres,]


peaktimes1 <- sig_results_EXPOSED$phase_hour
peaktimes2 <- sig_results_PROTECTED$phase_hour

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
    fill = "Site",
    title = paste0("Peak Time Distributions for Rhythmic Transcripts in Each Site")
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



ggsave(filename = "plots/figS1/rose_plot_rhythmic_in_each.png", plot = rose_plot, width = 8, height = 6, dpi = 300)



#### percentage checks

# about 64% of EXPOSED peaks were between 10pm and 6am
sum(sig_results_EXPOSED$phase_hour >= 22 | sig_results_EXPOSED$phase_hour <= 6) / length(sig_results_EXPOSED$phase_hour)

# about 51% of EXPOSED peaks were between 10pm and 6am
sum(sig_results_PROTECTED$phase_hour >= 22 | sig_results_PROTECTED$phase_hour <= 6) / length(sig_results_PROTECTED$phase_hour)




##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG S1E - DNA repair ring plot


## clear environment
rm(list=ls())



read_gmt <- function(file_path) {
  gmt_data <- read_lines(file_path)
  
  gene_sets <- map(gmt_data, function(line) {
    tokens <- strsplit(line, "\t")[[1]]
    list(name = tokens[1], description = tokens[2], genes = tokens[-(1:2)])
  })
  
  return(gene_sets)
}



id_to_name <- readRDS("results/RDS_files/id_to_name.RDS")
rhythmic_in_either_ids <- readRDS("results/RDS_files/rhythmic_in_either_ids.RDS")
id_to_name_RHYTHMIC <- id_to_name[id_to_name$gene_id %in% rhythmic_in_either_ids,]


gene_sets_HALLMARK <- read_gmt("pathway_analysis_inputs/h.all.v2024.1.Hs.symbols.gmt")

term_name <- "HALLMARK_DNA_REPAIR"
term <- gene_sets_HALLMARK %>% keep(~ .x$name == term_name)
term_genes <- term[[1]]$genes
genes_to_plot <- id_to_name_RHYTHMIC$gene_name[id_to_name_RHYTHMIC$gene_name %in% term_genes]








rel_amp_thres <- 0.2
bhq_thres <- 0.05


data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)

sig_results_EXPOSED <- results_EXPOSED[
  results_EXPOSED$bh_q <= bhq_thres & 
    results_EXPOSED$rel_amp >= rel_amp_thres &
    results_EXPOSED$gene_name %in% genes_to_plot,
]

data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)

sig_results_PROTECTED <- results_PROTECTED[
  results_PROTECTED$bh_q <= bhq_thres & 
    results_PROTECTED$rel_amp >= rel_amp_thres &
    results_PROTECTED$gene_name %in% genes_to_plot,
]

shared_genes <- intersect(sig_results_EXPOSED$gene_name, sig_results_PROTECTED$gene_name)

shared_colors <- hue_pal()(length(shared_genes))
names(shared_colors) <- shared_genes

assign_color <- function(gene) {
  if (gene %in% names(shared_colors)) {
    return(shared_colors[gene])
  } else {
    return("black")
  }
}

outer_df <- data.frame(
  Gene      = sig_results_PROTECTED$gene_name,
  Acrophase = sig_results_PROTECTED$phase_hour,
  Condition = "PROTECTED",
  Radius    = 2
)

inner_df <- data.frame(
  Gene      = sig_results_EXPOSED$gene_name,
  Acrophase = sig_results_EXPOSED$phase_hour,
  Condition = "EXPOSED",
  Radius    = 1
)

acrophase_data <- rbind(outer_df, inner_df)

acrophase_data$Angle <- 2 * pi * acrophase_data$Acrophase / 24

acrophase_data$Color <- sapply(acrophase_data$Gene, assign_color)

p <- ggplot(acrophase_data, aes(x = Angle, y = Radius, label = Gene)) +
  geom_hline(yintercept = 1, color = "grey50") +
  geom_hline(yintercept = 2, color = "grey50") +
  geom_point(aes(color = Color), size = 3) +
  scale_x_continuous(
    limits = c(0, 2 * pi),
    breaks = seq(0, 2 * pi, length.out = 13),
    labels = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
  ) +
  scale_y_continuous(
    limits = c(0, 2.3),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_polar(start = -0 * pi / 2) +
  scale_color_identity() +
  labs(title = "", x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )



p <- p + geom_text(x = -pi/4, y = 1, label = "EXPOSED", color = "black",
                   size = 4, vjust = -1) +
  geom_text(x = -pi/4, y = 2, label = "PROTECTED", color = "black",
            size = 4, vjust = -1)




ggsave(
  filename = "plots/figS1_NOLABELS/ring_plot_dna_repair.png",
  plot     = p,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 500,
  bg       = "white"
)


p <- p +
  geom_text_repel(
    size           = 2.5,
    box.padding    = 0.5,
    point.padding  = 0.3,
    segment.color  = "black",
    segment.size   = 0.5,
    min.segment.length = 0,
    force = 2
  )



ggsave(
  filename = "plots/figS1/ring_plot_dna_repair.png",
  plot     = p,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 500,
  bg       = "white"
)

