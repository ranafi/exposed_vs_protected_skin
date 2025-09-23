
## clear environment
rm(list=ls()) 

# install.packages("ggplot2")
library(ggplot2)
# install.packages("ggsignif")
library(ggsignif)
# install.packages("dplyr")
library(dplyr)
# install.packages("readr")
library(readr)
# install.packages("purrr")
library(purrr)




##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

##### CHECK - ALL RHYTHMIC

## clear environment
rm(list=ls()) 


data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)


data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)




filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)



cycling_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_id %in% rhythmic_in_either_ids,]
cycling_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_id %in% rhythmic_in_either_ids,]


# check same order - should be 0
# sum(cycling_EXPOSED$gene_id != cycling_PROTECTED$gene_id)


# 2035 genes cycling in either
length(rhythmic_in_either_ids)


# 1409 had higher rel amp in PROTECTED
sum(cycling_PROTECTED$rel_amp > cycling_EXPOSED$rel_amp)

# 1383 had higher abs amp in PROTECTED
sum(cycling_PROTECTED$abs_amp > cycling_EXPOSED$abs_amp)


# 480/2035 were differentially rhythmic
filename <- paste0("results/RDS_files/DR_ids.RDS")
DR_ids <- readRDS(filename)
sum(cycling_EXPOSED$gene_id %in% DR_ids)
sum(cycling_EXPOSED$gene_id %in% DR_ids)



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

##### CHECK - CORECLOCK

## clear environment
rm(list=ls()) 


coreclock_genes <- c( "NR1D2", "PER1", "PER2", "CRY2", "RORC", "BHLHE40", "PER3", "CRY1", "NPAS2", "NFIL3", "CIART", "TEF", "BMAL1", "CLOCK", "NR1D1", "DBP", "HLF", "RORA", "BHLHE41", "CSNK1D")
# 20 genes total (not all rhythmic)


data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)


data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)


# should both be 20
sum(results_EXPOSED$gene_name %in% coreclock_genes)
sum(unique(results_EXPOSED$gene_name) %in% coreclock_genes)



filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

rhythmic_in_either <- id_to_name$gene_name[id_to_name$gene_id %in% rhythmic_in_either_ids & id_to_name$gene_name != "na"]






cycling_core_clock_genes <- coreclock_genes[coreclock_genes %in% rhythmic_in_either]



coreclock_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_name %in% cycling_core_clock_genes,]
coreclock_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_name %in% cycling_core_clock_genes,]


# check same order - should be 0
sum(coreclock_EXPOSED$gene_id != coreclock_PROTECTED$gene_id)


# 15 core clock genes cycling in either
length(cycling_core_clock_genes)


# 13 had higher rel amp in PROTECTED
sum(coreclock_PROTECTED$rel_amp > coreclock_EXPOSED$rel_amp)
binom.test(13, 15, p = 0.5, alternative = "two.sided") # compared to 50/50
binom.test(13, 15, p = 1409/2035, alternative = "two.sided") # compared to background context



# 12 had higher abs amp in PROTECTED
sum(coreclock_PROTECTED$abs_amp > coreclock_EXPOSED$abs_amp)


# 1/15 was differentially rhythmic
filename <- paste0("results/RDS_files/DR_ids.RDS")
DR_ids <- readRDS(filename)
sum(coreclock_EXPOSED$gene_id %in% DR_ids)
sum(coreclock_PROTECTED$gene_id %in% DR_ids)
coreclock_EXPOSED$gene_name[coreclock_EXPOSED$gene_id %in% DR_ids]




##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

##### CHECK - DNA REPAIR

## clear environment
rm(list=ls()) 



term_name <- "HALLMARK_DNA_REPAIR"

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


# should both be the same, to make sure no duplicates
sum(results_EXPOSED$gene_name %in% term_genes)
sum(unique(results_EXPOSED$gene_name) %in% term_genes)





filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

rhythmic_in_either <- id_to_name$gene_name[id_to_name$gene_id %in% rhythmic_in_either_ids & id_to_name$gene_name != "na"]





cycling_term_genes <- term_genes[term_genes %in% rhythmic_in_either]



term_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_name %in% cycling_term_genes,]
term_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_name %in% cycling_term_genes,]


# check same order - should be 0
sum(term_EXPOSED$gene_id != term_PROTECTED$gene_id)


# 17 DNA repair genes cycling in either
length(cycling_term_genes)


# 14 had higher rel amp in PROTECTED
sum(term_PROTECTED$rel_amp > term_EXPOSED$rel_amp)
binom.test(14, 17, p = 0.5, alternative = "two.sided") # compared to 50/50
binom.test(14, 17, p = 1409/2035, alternative = "two.sided") # compared to background context

# 14 had higher abs amp in PROTECTED
sum(term_PROTECTED$abs_amp > term_EXPOSED$abs_amp)




# 6/17 were differentially rhythmic
filename <- paste0("results/RDS_files/DR_ids.RDS")
DR_ids <- readRDS(filename)
sum(term_EXPOSED$gene_id %in% DR_ids)
sum(term_PROTECTED$gene_id %in% DR_ids)
term_EXPOSED$gene_name[term_EXPOSED$gene_id %in% DR_ids]







##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

##### CHECK - MYC TARGETS

## clear environment
rm(list=ls()) 


# not technically a gmt file but close enough format that i can reuse this function and get correct parsing
read_gmt <- function(file_path) {
  gmt_data <- read_lines(file_path)
  
  gene_sets <- map(gmt_data, function(line) {
    tokens <- strsplit(line, "\t")[[1]]
    list(name = tokens[1], description = tokens[2], genes = tokens[-(1:2)])
  })
  
  return(gene_sets)
}

gene_sets <- read_gmt("pathway_analysis_inputs/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")

term_name <- "MYC CHEA"
term1 <- gene_sets %>% keep(~ .x$name == term_name)
term_genes1 <- term1[[1]]$genes


term_name <- "MYC ENCODE"
term2 <- gene_sets %>% keep(~ .x$name == term_name)
term_genes2 <- term2[[1]]$genes


term_genes <- unique(c(term_genes1, term_genes2))






data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)


data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)


# should both be the same, to make sure no duplicates
sum(results_EXPOSED$gene_name %in% term_genes)
sum(unique(results_EXPOSED$gene_name) %in% term_genes)


filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

rhythmic_in_either <- id_to_name$gene_name[id_to_name$gene_id %in% rhythmic_in_either_ids & id_to_name$gene_name != "na"]



cycling_term_genes <- term_genes[term_genes %in% rhythmic_in_either]



term_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_name %in% cycling_term_genes,]
term_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_name %in% cycling_term_genes,]


# check same order - should be 0
sum(term_EXPOSED$gene_id != term_PROTECTED$gene_id)


# 291 DNA repair genes cycling in either
length(cycling_term_genes)

# 259 had higher rel amp in PROTECTED
sum(term_PROTECTED$rel_amp > term_EXPOSED$rel_amp)

binom.test(259, 291, p = 0.5, alternative = "two.sided") # compared to 50/50
binom.test(259, 291, p = 1409/2035, alternative = "two.sided") # compared to background context

# 260 had higher abs amp in PROTECTED
sum(term_PROTECTED$abs_amp > term_EXPOSED$abs_amp)


# 83/291 were differentially rhythmic
filename <- paste0("results/RDS_files/DR_ids.RDS")
DR_ids <- readRDS(filename)
sum(term_EXPOSED$gene_id %in% DR_ids)
sum(term_PROTECTED$gene_id %in% DR_ids)
term_EXPOSED$gene_name[term_EXPOSED$gene_id %in% DR_ids]


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## FIG 4A - PATHWAY BAR PLOT


## clear environment
rm(list = ls())



data <- data.frame(
  Condition = c("All Rhythmic\n(2,035 Genes)", "Core Clock\n(15 Genes)",
                "DNA Repair\n(17 Genes)", "MYC Targets\n(291 Genes)"),
  Successes = c(1409, 13, 14, 259),
  Trials    = c(2035, 15, 17, 291)
)


data$Success <- data$Successes / data$Trials

get_ci <- function(x, n) binom.test(x, n)$conf.int

ci_mat <- t(mapply(get_ci, data$Successes, data$Trials))
data$Lower <- ci_mat[, 1]
data$Upper <- ci_mat[, 2]


baseline_success <- data$Successes[1]
baseline_trials  <- data$Trials[1]

pvals <- c(NA, sapply(2:nrow(data),
                      function(i) binom.test(data$Successes[i], data$Trials[i],
                                             p = baseline_success / baseline_trials)$p.value))


plot_data <- do.call(rbind, lapply(1:nrow(data), function(i) {
  data.frame(
    Condition  = rep(data$Condition[i], 2),
    Outcome    = c("Higher in PROTECTED", "Higher in EXPOSED"),
    Probability= c(data$Success[i], 1 - data$Success[i]),
    Lower      = c(data$Lower[i], 1 - data$Upper[i]),
    Upper      = c(data$Upper[i], 1 - data$Lower[i]),
    pval       = pvals[i]
  )
}))

plot_data$Condition <- factor(plot_data$Condition, levels = data$Condition)


p <- ggplot(plot_data,
            aes(x = Condition, y = Probability, fill = Outcome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           width = 0.8, alpha = 0.5) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Higher in PROTECTED" = "blue",
                               "Higher in EXPOSED"   = "red")) +
  labs(title = "Relative Amplitude",
       y = "Fraction", x = NULL) +
  theme_classic(base_size = 30)


prot_bars  <- subset(plot_data, Outcome == "Higher in PROTECTED")

ggsave("plots/fig5/binomial_bar_plot_CI.png",
       plot = p, width = 20, height = 8, units = "in", dpi = 300)




##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################


## FIG 4C,D - MYC TARGETS RELATIVE AMPLITUDE AND ACROPHASE


## clear environment
rm(list=ls()) 



# not technically a gmt file but close enough format that i can reuse this function and get correct parsing
read_gmt <- function(file_path) {
  gmt_data <- read_lines(file_path)
  
  gene_sets <- map(gmt_data, function(line) {
    tokens <- strsplit(line, "\t")[[1]]
    list(name = tokens[1], description = tokens[2], genes = tokens[-(1:2)])
  })
  
  return(gene_sets)
}

gene_sets <- read_gmt("pathway_analysis_inputs/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")

term_name <- "MYC CHEA"
term1 <- gene_sets %>% keep(~ .x$name == term_name)
term_genes1 <- term1[[1]]$genes


term_name <- "MYC ENCODE"
term2 <- gene_sets %>% keep(~ .x$name == term_name)
term_genes2 <- term2[[1]]$genes


term_genes <- unique(c(term_genes1, term_genes2))




data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)

data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)



filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

# rhythmic_in_either <- id_to_name$gene_name[id_to_name$gene_id %in% rhythmic_in_either_ids & id_to_name$gene_name != "na"]





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
  geom_histogram(position = "identity", alpha = 0.5, bins = 50, color = "black") +
  labs(title = paste0("Relative Amplitude Distributions - MYC Targets (",length(results_EXPOSED$rel_amp)," Transcripts)" ),
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


filename <- "plots/fig5/amp_dists_MYC_TARGETS.png"
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


filename <- "plots/fig5_NOLABELS/amp_dists_MYC_TARGETS.png"
ggsave(filename, plot = p, width = 4, height = 4, units = "in", dpi = 300)
########################################################################################

  
  
  







## peaktime distribution (MYC targets from ENCODE and CHEA)

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
    limits = c(0, 2*pi + 1e-5), # slight buffer here to solve out of bounds warning due to floating point rounding error
    breaks = break_angles,
    labels = paste0(hour_labels, "h")
  ) +
  scale_fill_manual(values = fill_colors) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Dataset",
    title = paste0("Peak Time Comparison - MYC Targets (",length(peaktimes1)," Transcripts)")
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



# # checks - should be false
# table(is.na(bins1))
# table(is.na(bins2))
# range(peaktimes1, na.rm = TRUE)
# range(peaktimes2, na.rm = TRUE)
# summary(data_counts1)
# summary(data_counts2)
# subset(data_counts_combined, bin_rad < 0 | bin_rad > 2 * pi | is.na(bin_rad))




ggsave(filename = "plots/fig5/rose_plot_MYC_TARGETS.png", plot = rose_plot, width = 8, height = 6, dpi = 300)



  
  
  
  




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
    limits = c(0, 2*pi + 1e-5), # slight buffer here to solve out of bounds warning due to floating point rounding error
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



# # checks - should be false
# table(is.na(bins1))
# table(is.na(bins2))
# range(peaktimes1, na.rm = TRUE)
# range(peaktimes2, na.rm = TRUE)
# summary(data_counts1)
# summary(data_counts2)
# subset(data_counts_combined, bin_rad < 0 | bin_rad > 2 * pi | is.na(bin_rad))




ggsave(filename = "plots/fig5_NOLABELS/rose_plot_MYC_TARGETS.png", plot = rose_plot, width = 8, height = 6, dpi = 300)
########################################################################################


  

