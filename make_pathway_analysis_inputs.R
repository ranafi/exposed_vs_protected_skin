


## clear environment
rm(list=ls()) 





##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

# MAKE PSEA INPUTS

## clear environment
rm(list=ls()) 



rel_amp_thres <- 0.2
bhq_thres <- 0.05



data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results <- readRDS(filename)

sig_results <- results[results$bh_q <= bhq_thres & results$gene_name != "na" & results$rel_amp >= rel_amp_thres,]


# check for duplicates. if > 1, that means there are duplicates
# max(table(sig_results$gene_name))


to_write <- data.frame(gene = sig_results$gene_name, phase = sig_results$phase_hour)
output_file <- paste0("pathway_analysis_inputs/psea_input_",data_label,".txt")
write.table(to_write, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)





##################################################################################################################
# shifted for plots
to_write$phase <- to_write$phase - 6
to_write$phase[to_write$phase<0] <- to_write$phase[to_write$phase<0]+24
output_file <- paste0("pathway_analysis_inputs/psea_input_shifted_for_plot_",data_label,".txt")
write.table(to_write, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
##################################################################################################################






data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results <- readRDS(filename)

sig_results <- results[results$bh_q <= bhq_thres & results$gene_name != "na" & results$rel_amp >= rel_amp_thres,]


# check for duplicates. if > 1, that means there are duplicates
# max(table(sig_results$gene_name))


to_write <- data.frame(gene = sig_results$gene_name, phase = sig_results$phase_hour)
output_file <- paste0("pathway_analysis_inputs/psea_input_",data_label,".txt")
write.table(to_write, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)





##################################################################################################################
# shifted for plots
to_write$phase <- to_write$phase - 6
to_write$phase[to_write$phase<0] <- to_write$phase[to_write$phase<0]+24
output_file <- paste0("pathway_analysis_inputs/psea_input_shifted_for_plot_",data_label,".txt")
write.table(to_write, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
##################################################################################################################



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

# MAKE GSEA INPUTS
# cycling significant scores (but writing all genes, not just cycling)


## clear environment
rm(list=ls()) 



data_label <- "EXPOSED"

filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results <- readRDS(filename)
results <- results[results$gene_name != "na",]


gene_names <- results$gene_name
sig_scores <- -log(results$p)


sig_order <- order(sig_scores, decreasing=TRUE)
sig_scores <- sig_scores[sig_order]
gene_names <- gene_names[sig_order]

to_write <- data.frame(gene=gene_names, score = sig_scores )

# removes duplicates by keeping first occurrence (which is the one with the higher score)
to_write <- to_write[!duplicated(to_write$gene), ]

output_file <- paste0("pathway_analysis_inputs/gsea_input_",data_label,".rnk")
write.table(to_write, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




data_label <- "PROTECTED"

filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results <- readRDS(filename)
results <- results[results$gene_name != "na",]


gene_names <- results$gene_name
sig_scores <- -log(results$p)


sig_order <- order(sig_scores, decreasing=TRUE)
sig_scores <- sig_scores[sig_order]
gene_names <- gene_names[sig_order]

to_write <- data.frame(gene=gene_names, score = sig_scores )

# removes duplicates by keeping first occurrence (which is the one with the higher score)
to_write <- to_write[!duplicated(to_write$gene), ]

output_file <- paste0("pathway_analysis_inputs/gsea_input_",data_label,".rnk")
write.table(to_write, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)






##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

# MAKE GSEA AMPLITUDE RATIO INPUTS
# amplitude ratios - EXPOSED / PROTECTED


## clear environment
rm(list=ls()) 


filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)



data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)
results_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_id %in% rhythmic_in_either_ids & results_PROTECTED$gene_name != "na",]



data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)
results_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_id %in% rhythmic_in_either_ids & results_EXPOSED$gene_name != "na",]


# check order - should be 0
# sum(results_PROTECTED$gene_name != results_EXPOSED$gene_name)


gene_names <- results_PROTECTED$gene_name

# check - if > 1, that means there's duplicates that need to be handled
# max(table(gene_names))


# write GSEA input
amp_ratios <- log2(results_EXPOSED$rel_amp / results_PROTECTED$rel_amp)
amp_order <- order(amp_ratios, decreasing=TRUE)
amp_ratios <- amp_ratios[amp_order]
gene_names <- gene_names[amp_order]
to_write <- data.frame(gene=gene_names, score = amp_ratios )
write.table(to_write, file = "pathway_analysis_inputs/gsea_amp_ratios.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)







##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## ENRICHR

## clear environment
rm(list=ls()) 


filename <- paste0("results/RDS_files/DR_ids.RDS")
DR_ids <- readRDS(filename)




data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)
results_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_id %in% DR_ids & results_PROTECTED$gene_name != "na",]



data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)
results_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_id %in% DR_ids & results_EXPOSED$gene_name != "na",]




to_write_EXPOSED <- results_EXPOSED$gene_name[results_EXPOSED$rel_amp > results_PROTECTED$rel_amp]
to_write_PROTECTED <- results_EXPOSED$gene_name[results_PROTECTED$rel_amp > results_EXPOSED$rel_amp]


# check - if > 1, that means there's duplicates that need to be handled
# max(table(to_write_EXPOSED))
# max(table(to_write_PROTECTED))

write.table(to_write_EXPOSED, file = "pathway_analysis_inputs/enrichr_higher_amp_EXPOSED.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(to_write_PROTECTED, file = "pathway_analysis_inputs/enrichr_higher_amp_PROTECTED.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# background - all genes tested for differential rhythmicity (meaning all genes cycling in at least one condition)


filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

filename <- paste0("results/RDS_files/id_to_name.RDS")
id_to_name <- readRDS(filename)

to_write_background <- id_to_name$gene_name[id_to_name$gene_id %in% rhythmic_in_either_ids & id_to_name$gene_name != "na"]

# check - if > 1, that means there's duplicates that need to be handled
# max(table(to_write_background))

write.table(to_write_background, file = "pathway_analysis_inputs/enrichr_background.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

