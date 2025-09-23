
## clear environment
rm(list=ls()) 


# install.packages("ggrepel")
library(ggrepel)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("ggnewscale")
library(ggnewscale)
# install.packages("dplyr")
library(dplyr)
# install.packages("readr")
library(readr)
# install.packages("purrr")
library(purrr)

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 6A - PATHWAY PLOT


## clear environment
rm(list=ls()) 



data_label <- "EXPOSED"


###################
# GSEA  

results <- read.table(paste0('results/pathway_analysis/gsea_results/gsea_results_AMP_RATIOS/gsea_report_for_na_pos_1747673097610.tsv'), header = TRUE, sep ='\t', row.names=1)


rownames(results) <- gsub("HALLMARK_", "", rownames(results))

gsea_ES <- abs(results$ES)
names(gsea_ES) <- rownames(results)

gene_set_sizes <- results$SIZE








###################
# Enrichr  



results <- read.table(paste0('results/pathway_analysis/enrichr_results/enrichr_results_EXPOSED.txt'), header = TRUE, sep ='\t', row.names=1)


enrichr_labels <- rownames(results)
enrichr_labels <- toupper(gsub(" ", "_", enrichr_labels))





names(gsea_ES)[!(names(gsea_ES) %in% enrichr_labels)]





enrichr_labels[enrichr_labels=="IL-6/JAK/STAT3_SIGNALING"] <- "IL6_JAK_STAT3_SIGNALING"
enrichr_labels[enrichr_labels=="TNF-ALPHA_SIGNALING_VIA_NF-KB"] <- "TNFA_SIGNALING_VIA_NFKB"
enrichr_labels[enrichr_labels=="PI3K/AKT/MTOR__SIGNALING"] <- "PI3K_AKT_MTOR_SIGNALING"
enrichr_labels[enrichr_labels=="IL-2/STAT5_SIGNALING"] <- "IL2_STAT5_SIGNALING"
enrichr_labels[enrichr_labels=="G2-M_CHECKPOINT"] <- "G2M_CHECKPOINT"
enrichr_labels[enrichr_labels=="TGF-BETA_SIGNALING"] <- "TGF_BETA_SIGNALING"





enrichr_scores <- results$Combined.Score
names(enrichr_scores) <- enrichr_labels


enrichr_q_values <- results$Adjusted.P.value
names(enrichr_q_values) <- enrichr_labels

enrichr_scores <- enrichr_scores[names(gsea_ES)]
enrichr_q_values <- enrichr_q_values[names(gsea_ES)]


neg_log_enrichr_q_values <- -log10(enrichr_q_values)

labels <- names(gsea_ES)
names(labels) <- names(gsea_ES)



############################################################################################################
names(gene_set_sizes) <- names(gsea_ES)




labels <- gsub("_", " ", labels)



plot_data <- data.frame(
  x = enrichr_scores,
  y = gsea_ES,
  size_variable = gene_set_sizes,
  label = labels,
  neg_log_enrichr_q_values = neg_log_enrichr_q_values
)



# only keep terms where we have a score for both gsea and enrichr
plot_data <- plot_data[!(is.na(plot_data$x)) & !(is.na(plot_data$x)),]



plot_data$scaled_x <- (plot_data$x - min(plot_data$x)) / (max(plot_data$x) - min(plot_data$x))
plot_data$scaled_y <- (plot_data$y - min(plot_data$y)) / (max(plot_data$y) - min(plot_data$y))

plot_data$label_score <- plot_data$scaled_x + plot_data$scaled_y

plot_data$rank <- rank(-plot_data$label_score, ties.method = "first")



xlabel <- "Combined Score - Enrichr"
ylabel <- "Enrichment Score - GSEA"



plot_data$label[plot_data$label=="EPITHELIAL MESENCHYMAL TRANSITION"] <- "EPITHELIAL\nMESENCHYMAL\nTRANSITION"
plot_data$label[plot_data$label=="ESTROGEN RESPONSE EARLY"] <- "ESTROGEN\nRESPONSE\nEARLY"




  
  
  
  
  
  
scatter_plot <- ggplot(plot_data, aes(x = x, y = y, color = neg_log_enrichr_q_values)) +
  geom_point(aes(size = size_variable), alpha = 0.6) +
  scale_size_continuous(name = "Gene Set Size", range = c(1, 7)) +
  new_scale("size") +
  geom_text_repel(
    data = subset(plot_data, x > 5 & y > 0.2),
    aes(label = label),
    size = 5,
    segment.color = 'grey50',
    segment.size = 0.1,
    arrow = arrow(length = unit(0.005, 'npc')),
    force = 1,
    max.overlaps = Inf,
    box.padding = 0.75,
    point.padding = 0.5,
    show.legend = FALSE
  ) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  theme_minimal(base_size = 14) +
  xlab("Combined Score - Enrichr") +
  ylab("Enrichment Score - GSEA") +
  ggtitle(paste0("GSEA vs Enrichr -- Higher Amplitude in ", data_label)) +
  theme(legend.position = "bottom") +
  # Make the axis labels/title larger
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  scale_color_gradient(low = "black", high = "#FF0000", name = "-log10(Enrichr q-values)")


filename <- paste0("plots/fig7/gsea_amp_ratio_vs_enrichr_higher_amp_DR_", data_label, ".png")
ggsave(filename, plot = scatter_plot, width = 12, height = 8, bg = "white", dpi=300)





########################################################################################
# WITHOUT LABELS  
scatter_plot <- ggplot(plot_data, aes(x = x, y = y, color = neg_log_enrichr_q_values)) +
  geom_point(aes(size = size_variable), alpha = 0.6, show.legend = TRUE) +
  geom_text_repel(
    data = subset(plot_data, x > 10 & y > 0.2),
    aes(label = label),
    size = 9,
    segment.color = 'grey50',
    segment.size = 0.1,
    arrow = arrow(length = unit(0.005, 'npc')),
    force = 1,
    max.overlaps = Inf,
    box.padding = 0.75,
    point.padding = 0.5,
    show.legend = FALSE
  ) +
  scale_size_continuous(name = NULL, range = c(3, 10)) +
  scale_color_gradient(low = "black", high = "#FF0000", name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(order = 2)
  ) +
  theme_minimal(base_size = 14) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 14),
    plot.title = element_blank(),
    legend.title = element_blank()
  )



filename <- paste0("plots/fig7_NOLABELS/gsea_amp_ratio_vs_enrichr_higher_amp_DR_", data_label, ".png")
ggsave(filename, plot = scatter_plot, width = 14, height = 8, bg = "white", dpi=300)
########################################################################################



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 6B,C - EMT GENES, RELATIVE AMPLITUDE


## clear environment
rm(list=ls()) 




term_name <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"


term_label <- gsub("HALLMARK_", "", term_name)


read_gmt <- function(file_path) {
  gmt_data <- read_lines(file_path)
  
  gene_sets <- map(gmt_data, function(line) {
    tokens <- strsplit(line, "\t")[[1]]
    list(name = tokens[1], description = tokens[2], genes = tokens[-(1:2)])
  })
  
  return(gene_sets)
}

gene_sets <- read_gmt("pathway_analysis_inputs/h.all.v2024.1.Hs.symbols.gmt")


term <- gene_sets %>% keep(~ .x$name == term_name)
term_genes <- term[[1]]$genes






data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)

data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)






#############################################

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

#############################################

results_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_id %in% rhythmic_in_either_ids,]
results_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_id %in% rhythmic_in_either_ids,]


# next, keep MYC target genes
results_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_name %in% term_genes,]
results_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_name %in% term_genes,]

# check matching order - should be 0
# sum(results_EXPOSED$gene_id != results_PROTECTED$gene_id)



df1 <- data.frame(value = results_EXPOSED$rel_amp, group = "EXPOSED")
df2 <- data.frame(value = results_PROTECTED$rel_amp, group = "PROTECTED")

df_combined <- bind_rows(df1, df2)




  
  
p <- ggplot(df_combined, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.3, bins = 20, color = "black") +
  labs(title = paste0("Relative Amplitude Distributions - EMT-Related Genes (",length(results_EXPOSED$rel_amp)," Transcripts)" ),
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


filename <- "plots/fig7/amp_dists_EMT.png"
ggsave(filename, plot = p, width = 8, height = 6, units = "in", dpi = 300)



  
  

  
########################################################################################
# WITHOUT LABELS    
p <- ggplot(df_combined, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20, color = "black") +
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


filename <- "plots/fig7_NOLABELS/amp_dists_EMT.png"
ggsave(filename, plot = p, width = 4, height = 4, units = "in", dpi = 300)
########################################################################################


  


amp_ratios <- results_EXPOSED$rel_amp / results_PROTECTED$rel_amp


df <- data.frame(value = amp_ratios)



  

p <- ggplot(df, aes(x = log2(value))) +
  geom_histogram(bins = 20, fill = "#D9A6F9", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", , size = 1.2) +
  labs(title = paste0( "Relative Amplitude Ratio - EMT-Related Genes (",length(results_EXPOSED$rel_amp)," Transcripts)"),
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
    plot.title = element_text(size = 18)
    
  )



filename <- "plots/fig7/amp_ratios_EMT.png"
ggsave(filename, plot = p, width = 8, height = 6, units = "in", dpi = 300)

  
  
########################################################################################
# WITHOUT LABELS      
p <- ggplot(df, aes(x = log2(value))) +
  geom_histogram(bins = 20, fill = "#D9A6F9", color = "black") +
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


filename <- "plots/fig7_NOLABELS/amp_ratios_EMT.png"
ggsave(filename, plot = p, width = 4, height = 4, units = "in", dpi = 300)
########################################################################################


log2_test <- t.test(log2(amp_ratios))
print(log2_test)

# log2 dist mean
round(mean(log2(amp_ratios)),3)


