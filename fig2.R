
## clear environment
rm(list=ls()) 



# install.packages("eulerr")
library(eulerr)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("ggrepel")
library(ggrepel)
# install.packages("ggnewscale")
library(ggnewscale)
# install.packages("scales")
library(scales)


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 1B - CORECLOCK TWO-RING PLOT

## clear environment
rm(list = ls())


genes_to_plot <- c(
  "PER1", "PER2", "PER3",
  "CRY1", "CRY2", "NR1D2",
  "NPAS2", "BMAL1", "NR1D1", 
  "TEF", "CIART", "DBP", "HLF", "RORC",
  "NFIL3", "RORA", "BHLHE40"
)




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



p <- p + geom_text(x = pi/4, y = 1, label = "EXPOSED", color = "black",
                   size = 4, vjust = -1) +
  geom_text(x = pi/4, y = 2, label = "PROTECTED", color = "black",
            size = 4, vjust = -1)




ggsave(
  filename = "plots/fig2_NOLABELS/ring_plot_core_clock.png",
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
  filename = "plots/fig2/ring_plot_core_clock.png",
  plot     = p,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 500,
  bg       = "white"
)





##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 1C - PATHWAYS: GSEA EXPOSED vs GSEA PROTECTED


## clear environment
rm(list=ls()) 


results <- read.table(paste0('results/pathway_analysis/gsea_results/gsea_results_EXPOSED/gsea_report_for_na_pos_1747671560554.tsv'), header = TRUE, sep ='\t', row.names=1)
gene_set_sizes_EXPOSED <- results$SIZE
gsea_ES_EXPOSED <- results$ES
names(gsea_ES_EXPOSED) <- rownames(results)
names(gene_set_sizes_EXPOSED) <- rownames(results)




results <- read.table(paste0('results/pathway_analysis/gsea_results/gsea_results_PROTECTED/gsea_report_for_na_pos_1747672538069.tsv'), header = TRUE, sep ='\t', row.names=1)
gene_set_sizes_PROTECTED <- results$SIZE
gsea_ES_PROTECTED <- results$ES
names(gsea_ES_PROTECTED) <- rownames(results)
names(gene_set_sizes_PROTECTED) <- rownames(results)


# get PROTECTED in same order as EXPOSED

gsea_ES_PROTECTED <- gsea_ES_PROTECTED[names(gsea_ES_EXPOSED)]

gene_set_sizes <- gene_set_sizes_EXPOSED


gsea_ES_MEAN <- (gsea_ES_EXPOSED + gsea_ES_PROTECTED) / 2

###################
labels <- names(gsea_ES_EXPOSED)
labels <- gsub("HALLMARK_", "", labels)
labels <- gsub("_", " ", labels)

names(labels) <- names(gsea_ES_EXPOSED)

############################################################################################################
# names(gene_set_sizes) <- names(neg_log_psea)

plot_data <- data.frame(
  x = gsea_ES_PROTECTED,
  y = gsea_ES_EXPOSED,
  size_variable = gene_set_sizes,
  label = labels,
  ES_mean = gsea_ES_MEAN
)

plot_data$scaled_x <- (plot_data$x - min(plot_data$x)) / (max(plot_data$x) - min(plot_data$x))
plot_data$scaled_y <- (plot_data$y - min(plot_data$y)) / (max(plot_data$y) - min(plot_data$y))

plot_data$label_score <- plot_data$scaled_x + plot_data$scaled_y

plot_data$rank <- rank(-plot_data$label_score, ties.method = "first")

label_size_min <- 1.5 
label_size_max <- 4  

plot_data$label_size <- label_size_max / sqrt(plot_data$rank)

plot_data$label_size <- pmax(plot_data$label_size, label_size_min)

xlabel <- "GSEA Enrichment Score - PROTECTED" 
ylabel <- "GSEA Enrichment Score - EXPOSED"



plot_data$label[plot_data$label=="MYC TARGETS V2"] <- "MYC\nTARGETS V2"
plot_data$label[plot_data$label=="MYC TARGETS V1"] <- "MYC\nTARGETS V1"

plot_data$label[plot_data$label=="UNFOLDED PROTEIN RESPONSE"] <- "UNFOLDED PROTEIN\nRESPONSE"






scatter_plot <- ggplot(plot_data, aes(x = x, y = y, color = ES_mean)) +
  geom_point(aes(size = size_variable), alpha = 0.6) +
  scale_size_continuous(name = "Gene Set Size", range = c(1, 7)) +
  new_scale("size") +
  geom_text_repel(
    data = subset(plot_data, x > 0.5 & y > 0.2),
    aes(label = label),
    size = 6,
    segment.color = 'grey50',
    segment.size = 0.1,
    arrow = arrow(length = unit(0.005, 'npc')),
    force = 1,
    max.overlaps = Inf,
    box.padding = 0.75,
    point.padding = 0.5,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(label_size_min, label_size_max), guide = "none") +
  
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  theme_minimal(base_size = 14) +
  xlab(xlabel) +
  ylab(ylabel) +
  ggtitle(paste0("GSEA Enriched Terms: EXPOSED vs PROTECTED")) +
  theme(legend.position = "bottom") +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  scale_color_gradient(low = "black", high = "#D9A6F9", name = "Mean Enrichment Score")

filename <- paste0("plots/fig2/pathways_gsea_vs_gsea.png")
ggsave(filename, plot = scatter_plot, width = 12, height = 8, bg = "white",, dpi = 300)







########################################################################################
# WITHOUT LABELS  

scatter_plot <- ggplot(plot_data, aes(x = x, y = y, color = ES_mean)) +
  geom_point(aes(size = size_variable), alpha = 0.6, show.legend = TRUE) +
  # new_scale("size") +
  geom_text_repel(
    data = subset(plot_data, x > 0.6 & y > 0.5),
    aes(label = label),
    size = 9,
    segment.color = 'grey50',
    segment.size = 0.1,
    arrow = arrow(length = unit(0.005, 'npc')),
    force = 1,
    max.overlaps = Inf,
    box.padding = 0.75,
    point.padding = 0.5,
    nudge_x = 0.05,
    nudge_y = 0.05,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = subset(plot_data, x < 0.6 & y > 0.51),
    aes(label = label),
    size = 9,
    segment.color = 'grey50',
    segment.size = 0.1,
    arrow = arrow(length = unit(0.005, 'npc')),
    force = 1,
    max.overlaps = Inf,
    box.padding = 0.75,
    point.padding = 0.5,
    nudge_x = -0.05,
    nudge_y = 0.05,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = subset(plot_data, x > 0.51 & y < 0.51),
    aes(label = label),
    size = 7,
    segment.color = 'grey50',
    segment.size = 0.1,
    arrow = arrow(length = unit(0.005, 'npc')),
    force = 1,
    max.overlaps = Inf,
    box.padding = 0.75,
    point.padding = 0.5,
    nudge_x = 0.07,
    nudge_y = -0.01,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = subset(plot_data, x > 0.5 & y < 0.5),
    aes(label = label),
    size = 7,
    segment.color = 'grey50',
    segment.size = 0.1,
    arrow = arrow(length = unit(0.005, 'npc')),
    force = 1,
    max.overlaps = Inf,
    box.padding = 0.75,
    point.padding = 0.5,
    nudge_x = 0.05,
    nudge_y = -0.03,
    show.legend = FALSE
  ) +
  scale_size_continuous(name = NULL, range = c(3, 10)) +
  scale_color_gradient(low = "black", high = "#D9A6F9", name = NULL) + 
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  theme_minimal(base_size = 14) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 14),
    plot.title = element_blank(),
    legend.title = element_blank()
  )


filename <- paste0("plots/fig2_NOLABELS/pathways_gsea_vs_gsea.png")
ggsave(filename, plot = scatter_plot, width = 14, height = 8, bg = "white",, dpi = 300)
########################################################################################




##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 1D - VENN DIAGRAM


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


cycling_gene_ids_EXPOSED <- results_EXPOSED[results_EXPOSED$bh_q <= bhq_thres & results_EXPOSED$rel_amp >= rel_amp_thres,]$gene_id
cycling_gene_ids_PROTECTED <- results_PROTECTED[results_PROTECTED$bh_q <= bhq_thres & results_PROTECTED$rel_amp >= rel_amp_thres,]$gene_id


sets <- list(
  "EXPOSED" = cycling_gene_ids_EXPOSED,
  "PROTECTED" = cycling_gene_ids_PROTECTED
)



filename <- "plots/fig2/venn_diagram.png"

png(filename, 
    width = 7, height = 7, units = "in", res = 300) 

venn <- euler(sets)



plot(
  venn,
  quantities = TRUE,
  fills = c(
    "EXPOSED"            = "#FF0000",  # red
    "PROTECTED"          = "#0000FF",  # blue
    "EXPOSED&PROTECTED"  = "#D9A6F9"   # purple
  ),
  labels = list(font = 4),
  main   = "Cycling Genes",
  alpha=0.4
)




dev.off()

# ## checks
# length(cycling_gene_ids_EXPOSED)
# length(cycling_gene_ids_EXPOSED) / dim(results_EXPOSED)[1]
# 
# length(cycling_gene_ids_PROTECTED)
# length(cycling_gene_ids_PROTECTED) / dim(results_PROTECTED)[1]




coreclock_genes <- c( "NR1D2", "PER1", "PER2", "CRY2", "RORC", "BHLHE40", "PER3", "CRY1", "NPAS2", "NFIL3", "CIART", "TEF", "BMAL1", "CLOCK", "NR1D1", "DBP", "HLF", "RORA", "BHLHE41", "CSNK1D")


cycling_gene_names_EXPOSED <- results_EXPOSED[results_EXPOSED$bh_q <= bhq_thres & results_EXPOSED$rel_amp >= rel_amp_thres,]$gene_name
cycling_gene_names_PROTECTED <- results_PROTECTED[results_PROTECTED$bh_q <= bhq_thres & results_PROTECTED$rel_amp >= rel_amp_thres,]$gene_name



## core clock genes in EXPOSED ONLY
coreclock_genes[coreclock_genes %in% cycling_gene_names_EXPOSED & !(coreclock_genes %in% cycling_gene_names_PROTECTED)]

## core clock genes in PROTECTED ONLY
coreclock_genes[coreclock_genes %in% cycling_gene_names_PROTECTED & !(coreclock_genes %in% cycling_gene_names_EXPOSED)]

## core clock genes in BOTH
coreclock_genes[coreclock_genes %in% cycling_gene_names_PROTECTED & coreclock_genes %in% cycling_gene_names_EXPOSED]



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 1E - BH-Q THRESHOLD PLOT


## clear environment
rm(list=ls()) 


rel_amp_thres <- 0.2


data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)



data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)


bh_q_thresholds <- seq(0,0.2,0.001)


sig_counts_EXPOSED <- c()
sig_counts_PROTECTED <- c()

for (bh_q_threshold in bh_q_thresholds){
  
  
  sig_results_EXPOSED <- results_EXPOSED[results_EXPOSED$bh_q <= bh_q_threshold & results_EXPOSED$rel_amp >= rel_amp_thres,]
  sig_results_PROTECTED <- results_PROTECTED[results_PROTECTED$bh_q <= bh_q_threshold  & results_PROTECTED$rel_amp >= rel_amp_thres,]
  
  count_EXPOSED <- dim(sig_results_EXPOSED)[1]
  count_PROTECTED <- dim(sig_results_PROTECTED)[1]
  
  sig_counts_EXPOSED <- append(sig_counts_EXPOSED, count_EXPOSED)
  sig_counts_PROTECTED <- append(sig_counts_PROTECTED, count_PROTECTED)
  
}



df <- data.frame(
  bh_q_thresholds = bh_q_thresholds,
  EXPOSED = sig_counts_EXPOSED,
  PROTECTED = sig_counts_PROTECTED
)


df_long <- data.frame(
  bh_q_thresholds = rep(df$bh_q_thresholds, 2),
  Count = c(df$EXPOSED, df$PROTECTED),
  Condition = factor(rep(c("EXPOSED", "PROTECTED"), each = length(bh_q_thresholds)))
)






  
  
p <- ggplot(df_long, aes(x = bh_q_thresholds, y = Count, color = Condition)) +
  geom_line(size = 2) +
  scale_color_manual(values = c("EXPOSED" = "#FF0000", "PROTECTED" = "#0000FF")) +
  labs(x = "BH-q Threshold", 
       y = "Number of Cycling Transcripts",
       title = "Varying BH-q Threshold") +
  theme_minimal() +
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

filename <- "plots/fig2/threshold_plot_BHQ.png"
ggsave(filename, plot = p, width = 8, height = 6, units = "in", dpi = 300)

  
  
  

########################################################################################
# WITHOUT LABELS  
p <- ggplot(df_long, aes(x = bh_q_thresholds, y = Count, color = Condition)) +
  geom_line(size = 3) +
  scale_color_manual(values = c("EXPOSED" = "#FF0000", "PROTECTED" = "#0000FF")) +
  labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
  )

filename <- "plots/fig2_NOLABELS/threshold_plot_BHQ.png"
ggsave(filename, plot = p, width = 4, height = 4, units = "in", dpi = 300)
########################################################################################





##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 1E - RELATIVE AMPLITUDE THRESHOLD PLOT


## clear environment
rm(list=ls()) 





data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)



data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)



sig_counts_EXPOSED <- c()
sig_counts_PROTECTED <- c()

bh_q_threshold <- 0.05

rel_amp_thresholds <- seq(0,0.5,0.001)

for (rel_amp_thres in rel_amp_thresholds){
  
  
  sig_results_EXPOSED <- results_EXPOSED[results_EXPOSED$bh_q <= bh_q_threshold & results_EXPOSED$rel_amp >= rel_amp_thres,]
  sig_results_PROTECTED <- results_PROTECTED[results_PROTECTED$bh_q <= bh_q_threshold  & results_PROTECTED$rel_amp >= rel_amp_thres,]
  
  count_EXPOSED <- dim(sig_results_EXPOSED)[1]
  count_PROTECTED <- dim(sig_results_PROTECTED)[1]
  
  sig_counts_EXPOSED <- append(sig_counts_EXPOSED, count_EXPOSED)
  sig_counts_PROTECTED <- append(sig_counts_PROTECTED, count_PROTECTED)
  
}




df <- data.frame(
  rel_amp_thresholds = rel_amp_thresholds,
  EXPOSED = sig_counts_EXPOSED,
  PROTECTED = sig_counts_PROTECTED
)

df_long <- data.frame( 
  rel_amp_thresholds = rep(df$rel_amp_thresholds, 2),
  Count = c(df$EXPOSED, df$PROTECTED),
  Condition = factor(rep(c("EXPOSED", "PROTECTED"), each = length(rel_amp_thresholds)))
)







p <- ggplot(df_long, aes(x = rel_amp_thresholds, y = Count, color = Condition)) +
  geom_line(size=2) +
  scale_color_manual(values = c("EXPOSED" = "red", "PROTECTED" = "blue")) +
  labs(x = "Relative Amplitude Threshold", 
       y = "Number of Cycling Genes",
       title = "Varying Relative Amplitude Threshold") +
  theme_minimal() +
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


filename <- "plots/fig2/threshold_plot_REL_AMP.png"

ggsave(filename, plot = p, width = 8, height = 6, units = "in", dpi = 300)
  
  
  

  
  
########################################################################################
# WITHOUT LABELS    
p <- ggplot(df_long, aes(x = rel_amp_thresholds, y = Count, color = Condition)) +
  geom_line(size=3) +
  scale_color_manual(values = c("EXPOSED" = "red", "PROTECTED" = "blue")) +
  labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
  )

filename <- "plots/fig2_NOLABELS/threshold_plot_REL_AMP.png"
ggsave(filename, plot = p, width = 4, height = 4, units = "in", dpi = 300)
########################################################################################


