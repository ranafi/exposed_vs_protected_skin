

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






read_gmt <- function(file_path) {
  gmt_data <- read_lines(file_path)
  
  gene_sets <- map(gmt_data, function(line) {
    tokens <- strsplit(line, "\t")[[1]]
    list(name = tokens[1], description = tokens[2], genes = tokens[-(1:2)])
  })
  
  return(gene_sets)
}


make_timeseries_plots <- function(gene_ids_to_plot, text_labels=TRUE, output_path="plots/") {
  
  
  data_label_PROTECTED <- "PROTECTED"
  times_PROTECTED <- readRDS(paste("processed_data/times_", data_label_PROTECTED, ".RDS", sep=""))
  ids_PROTECTED <- readRDS(paste("processed_data/ids_", data_label_PROTECTED, ".RDS", sep=""))
  id_factors_PROTECTED <- readRDS(paste("processed_data/id_factors_", data_label_PROTECTED, ".RDS", sep=""))
  data_PROTECTED <- readRDS(paste("processed_data/gene_data_", data_label_PROTECTED, ".RDS", sep=""))
  
  data_label_EXPOSED <- "EXPOSED"
  times_EXPOSED <- readRDS(paste("processed_data/times_", data_label_EXPOSED, ".RDS", sep=""))
  ids_EXPOSED <- readRDS(paste("processed_data/ids_", data_label_EXPOSED, ".RDS", sep=""))
  id_factors_EXPOSED <- readRDS(paste("processed_data/id_factors_", data_label_EXPOSED, ".RDS", sep=""))
  data_EXPOSED <- readRDS(paste("processed_data/gene_data_", data_label_EXPOSED, ".RDS", sep=""))
  
  filename <- paste0("results/RDS_files/DR_ids.RDS")
  DR_ids <- readRDS(filename)

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
    
    y_min_PROTECTED <- min(adjusted_expression_PROTECTED)
    y_max_PROTECTED <- max(adjusted_expression_PROTECTED)
    
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
    # EXPOSED
    
    expression_EXPOSED <- unlist(as.vector(data_EXPOSED[gene_to_find,]))
    
    trig_model_EXPOSED <- lm(expression_EXPOSED ~ cos(times_EXPOSED * pi / 12) + sin(times_EXPOSED * pi / 12)  + id_factors_EXPOSED) 
    
    cos_fit_EXPOSED <- unname(trig_model_EXPOSED$coefficients[2])
    sin_fit_EXPOSED <- unname(trig_model_EXPOSED$coefficients[3])
    
    amp_EXPOSED <- sqrt(cos_fit_EXPOSED**2 + sin_fit_EXPOSED**2)
    phase_rad_fit_EXPOSED <- atan2(sin_fit_EXPOSED, cos_fit_EXPOSED)
    phase_hour_fit_EXPOSED <- unname(phase_rad_fit_EXPOSED * (12/pi))
    
    if (phase_hour_fit_EXPOSED < 0) {
      phase_hour_fit_EXPOSED <- 24 + phase_hour_fit_EXPOSED
    }
    
    c_EXPOSED <- trig_model_EXPOSED$coefficients[1]
    subject_factors_EXPOSED <- trig_model_EXPOSED$coefficients[4:length(trig_model_EXPOSED$coefficients)]
    subject_mesors_EXPOSED <- subject_factors_EXPOSED + c_EXPOSED
    subject_mesors_EXPOSED <- append(c_EXPOSED, subject_mesors_EXPOSED)
    average_mesor_EXPOSED <- mean(subject_mesors_EXPOSED)
    
    adjusted_expression_EXPOSED <- expression_EXPOSED
    
    for (sub in 1:20){
      adjusted_expression_EXPOSED[ids_EXPOSED == sub] <- adjusted_expression_EXPOSED[ids_EXPOSED == sub] - subject_mesors_EXPOSED[sub]
    }
    
    y_min_EXPOSED <- min(adjusted_expression_EXPOSED)
    y_max_EXPOSED <- max(adjusted_expression_EXPOSED)
    
    to_plot_EXPOSED <- data.frame(times = times_EXPOSED, adjusted_expression = adjusted_expression_EXPOSED)
    
    dup_to_plot_EXPOSED <- to_plot_EXPOSED
    dup_to_plot_EXPOSED$times <- dup_to_plot_EXPOSED$times+24
    
    fit_x <- seq(0, 48, by = 1)
    fit_line_EXPOSED <- data.frame(
      fit_x = fit_x,
      fit_y = cos_fit_EXPOSED * cos(fit_x * (2 * pi / 24)) + sin_fit_EXPOSED * sin(fit_x * (2 * pi / 24)) 
    )
    
    fit_line_EXPOSED_solid <- subset(fit_line_EXPOSED, fit_x <= 18)
    fit_line_EXPOSED_dotted <- subset(fit_line_EXPOSED, fit_x > 18)
    
    
    ########################################################################################################################################################################
    # plot
    
    
    y_min <- min(y_min_PROTECTED, y_min_EXPOSED)
    y_max <- max(y_max_PROTECTED, y_max_EXPOSED)
    
    
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
      
      plot_EXPOSED <- ggplot() +
        geom_point(data = to_plot_EXPOSED, aes(x = times, y = adjusted_expression), color = "black" , size   = 2.5) +
        geom_point(data = dup_to_plot_EXPOSED, aes(x = times, y = adjusted_expression), color = "gray", alpha = 0.9 , size   = 2.5) +
        geom_line(data = fit_line_EXPOSED_solid, aes(x = fit_x, y = fit_y), color = "blue", linetype = "solid", size   = 1.2) +
        geom_line(data = fit_line_EXPOSED_dotted, aes(x = fit_x, y = fit_y), color = "blue", linetype = "dotted", size   = 1.2) +
        scale_x_continuous(limits = c(0, 48), breaks = seq(0, 48, by = 12)) +
        scale_y_continuous(limits = c(y_min - 0.5, y_max + 0.5)) +
        labs(x = NULL, y = NULL, title = "EXPOSED")+
        theme_minimal() +
        theme(
          axis.text.y  = element_text(size = 16),
          axis.ticks.y = element_blank(),
          axis.text.x  = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18)
        ) +
        geom_vline(xintercept = phase_hour_fit_EXPOSED, linetype = "dashed", color = "red") +
        annotate("text", 
                 x = phase_hour_fit_EXPOSED, 
                 y = max(to_plot_EXPOSED$adjusted_expression, na.rm = TRUE), 
                 label = round(phase_hour_fit_EXPOSED, 2), 
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
      
      
      
      combined_plot <- plot_PROTECTED + plot_EXPOSED + 
        plot_layout(ncol = 2) + 
        plot_annotation(
          title = gene_name_title,
          theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
        )      
      
      
      
      if(gene_to_find %in% DR_ids){
        
        ggsave(paste0(output_path,"***timeseries_",gene_name,"_", gene_to_find, "***.png"), combined_plot, width = 10, height = 6, units = "in", dpi=300)
        
      }else{
        
        ggsave(paste0(output_path,"timeseries_",gene_name,"_", gene_to_find, ".png"), combined_plot, width = 10, height = 6, units = "in", dpi=300)
        
      }   
      
      
      
      
      
      
      
      
      
      
      
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
      
      plot_EXPOSED <- ggplot() +
        geom_point(data = to_plot_EXPOSED, aes(x = times, y = adjusted_expression), color = "black" , size   = 2.5) +
        geom_point(data = dup_to_plot_EXPOSED, aes(x = times, y = adjusted_expression), color = "gray", alpha = 0.9 , size   = 2.5) +
        geom_line(data = fit_line_EXPOSED_solid, aes(x = fit_x, y = fit_y), color = "blue", linetype = "solid", size   = 1.2) +
        geom_line(data = fit_line_EXPOSED_dotted, aes(x = fit_x, y = fit_y), color = "blue", linetype = "dotted", size   = 1.2) +
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
      
      
      
      combined_plot <- plot_PROTECTED + plot_EXPOSED + 
        plot_layout(ncol = 2) + 
        plot_annotation(
          title = gene_name_title,
          theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
        )      
      
      
      
      if(gene_to_find %in% DR_ids){
        
        ggsave(paste0(output_path,"***timeseries_",gene_name,"_", gene_to_find, "***.png"), combined_plot, width = 10, height = 6, units = "in", dpi=300)
        
      }else{
        
        ggsave(paste0(output_path,"timeseries_",gene_name,"_", gene_to_find, ".png"), combined_plot, width = 10, height = 6, units = "in", dpi=300)
        
      }
  
      
    }
    
  }
  
  
}




##########################################################################################################################################
## 15 core clock genes

coreclock_genes <- c( "NR1D2", "PER1", "PER2", "CRY2", "RORC", "BHLHE40", "PER3", "CRY1", "NPAS2", "NFIL3", "CIART", "TEF", "BMAL1", "CLOCK", "NR1D1", "DBP", "HLF", "RORA", "BHLHE41", "CSNK1D")
# 20 genes to check (not all are cycling)

filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

rhythmic_in_either <- id_to_name[id_to_name$gene_id %in% rhythmic_in_either_ids,]
core_clock_ids <- rhythmic_in_either$gene_id[rhythmic_in_either$gene_name %in% coreclock_genes]

make_timeseries_plots(core_clock_ids, output_path="plots/timeseries_expression_plots/core_clock/")


##########################################################################################################################################
## 480 differentially rhythmic genes

filename <- paste0("results/RDS_files/DR_ids.RDS")
DR_ids_to_plot <- readRDS(filename)

make_timeseries_plots(DR_ids_to_plot, output_path="plots/timeseries_expression_plots/differentially_rhythmic/")


##########################################################################################################################################
## 17 DNA repair genes

gene_sets <- read_gmt("pathway_analysis_inputs/h.all.v2024.1.Hs.symbols.gmt")

term_name <- "HALLMARK_DNA_REPAIR"
term <- gene_sets %>% keep(~ .x$name == term_name)
term_genes <- term[[1]]$genes

filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

rhythmic_in_either <- id_to_name[id_to_name$gene_id %in% rhythmic_in_either_ids,]
dna_repair_ids <- rhythmic_in_either$gene_id[rhythmic_in_either$gene_name %in% term_genes]

make_timeseries_plots(dna_repair_ids, output_path="plots/timeseries_expression_plots/DNA_repair/")


##########################################################################################################################################
## 291 MYC targets

gene_sets <- read_gmt("pathway_analysis_inputs/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")

term_name <- "MYC CHEA"
term1 <- gene_sets %>% keep(~ .x$name == term_name)
term_genes1 <- term1[[1]]$genes

term_name <- "MYC ENCODE"
term2 <- gene_sets %>% keep(~ .x$name == term_name)
term_genes2 <- term2[[1]]$genes

term_genes <- unique(c(term_genes1, term_genes2))

filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

rhythmic_in_either <- id_to_name[id_to_name$gene_id %in% rhythmic_in_either_ids,]
myc_target_ids <- rhythmic_in_either$gene_id[rhythmic_in_either$gene_name %in% term_genes]

make_timeseries_plots(myc_target_ids, output_path="plots/timeseries_expression_plots/MYC_targets/")


##########################################################################################################################################
## 462 E2F targets

gene_sets <- read_gmt("pathway_analysis_inputs/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")

term_name <- "E2F1 CHEA"
term1 <- gene_sets %>% keep(~ .x$name == term_name)
term_genes1 <- term1[[1]]$genes

term_name <- "E2F4 ENCODE"
term2 <- gene_sets %>% keep(~ .x$name == term_name)
term_genes2 <- term2[[1]]$genes

term_name <- "E2F6 ENCODE"
term3 <- gene_sets %>% keep(~ .x$name == term_name)
term_genes3 <- term3[[1]]$genes

term_genes <- unique(c(term_genes1, term_genes2, term_genes3))

filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

rhythmic_in_either <- id_to_name[id_to_name$gene_id %in% rhythmic_in_either_ids,]
e2f_target_ids <- rhythmic_in_either$gene_id[rhythmic_in_either$gene_name %in% term_genes]

make_timeseries_plots(e2f_target_ids, output_path="plots/timeseries_expression_plots/E2F_targets/")


##########################################################################################################################################
## 43 G2M checkpoint genes

gene_sets <- read_gmt("pathway_analysis_inputs/h.all.v2024.1.Hs.symbols.gmt")

term_name <- "HALLMARK_G2M_CHECKPOINT"
term <- gene_sets %>% keep(~ .x$name == term_name)
term_genes <- term[[1]]$genes

filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

rhythmic_in_either <- id_to_name[id_to_name$gene_id %in% rhythmic_in_either_ids,]
g2m_checkpoint_ids <- rhythmic_in_either$gene_id[rhythmic_in_either$gene_name %in% term_genes]

make_timeseries_plots(g2m_checkpoint_ids, output_path="plots/timeseries_expression_plots/G2M_checkpoint/")


##########################################################################################################################################
## 67 EMT genes

gene_sets <- read_gmt("pathway_analysis_inputs/h.all.v2024.1.Hs.symbols.gmt")

term_name <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
term <- gene_sets %>% keep(~ .x$name == term_name)
term_genes <- term[[1]]$genes

filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

rhythmic_in_either <- id_to_name[id_to_name$gene_id %in% rhythmic_in_either_ids,]
emt_ids <- rhythmic_in_either$gene_id[rhythmic_in_either$gene_name %in% term_genes]

make_timeseries_plots(emt_ids, output_path="plots/timeseries_expression_plots/EMT/")


##########################################################################################################################################
## 33 apical junction genes

gene_sets <- read_gmt("pathway_analysis_inputs/h.all.v2024.1.Hs.symbols.gmt")

term_name <- "HALLMARK_APICAL_JUNCTION"
term <- gene_sets %>% keep(~ .x$name == term_name)
term_genes <- term[[1]]$genes

filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

rhythmic_in_either <- id_to_name[id_to_name$gene_id %in% rhythmic_in_either_ids,]
apical_junction_ids <- rhythmic_in_either$gene_id[rhythmic_in_either$gene_name %in% term_genes]

make_timeseries_plots(apical_junction_ids, output_path="plots/timeseries_expression_plots/apical_junction/")


##########################################################################################################################################
## plots for paper figs

genes_to_plot <- c("PER3", "BMAL1", "TRPC6","PSPH","NREP","STON1","TP53","RPA2", "KLF4","SDC4",
                   "CDK4", "H2AX","PTX3","IGFBP4","VCAN","THY1")

filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

rhythmic_in_either <- id_to_name[id_to_name$gene_id %in% rhythmic_in_either_ids,]
genes_to_plot_ids <- rhythmic_in_either$gene_id[rhythmic_in_either$gene_name %in% genes_to_plot]

make_timeseries_plots(genes_to_plot_ids, output_path="plots/timeseries_expression_plots/for_paper_figs/", text_labels=FALSE)







