
## clear environment
rm(list=ls()) 


# install.packages("readr")
library(readr)
# install.packages("purrr")
library(purrr)


##########################################################################################################################################
## Tables S1-3

## clear environment
rm(list=ls()) 

filename <- paste0("results/RDS_files/first_pass_results_PROTECTED.RDS")
results_PROTECTED <- readRDS(filename)

rel_amp_thres <- 0.2
bhq_thres <- 0.05

sig_results_PROTECTED <- results_PROTECTED[results_PROTECTED$bh_q <= bhq_thres & results_PROTECTED$rel_amp >= rel_amp_thres,]

results_to_save <- data.frame(gene_id = sig_results_PROTECTED$gene_id, gene_name = sig_results_PROTECTED$gene_name,
                              rel_amp = sig_results_PROTECTED$rel_amp,
                              abs_amp = sig_results_PROTECTED$abs_amp,
                              phase_hour = sig_results_PROTECTED$phase_hour,
                              average_mesor = sig_results_PROTECTED$average_mesor,
                              p = sig_results_PROTECTED$p,
                              bh_q = sig_results_PROTECTED$bh_q)



write.table(results_to_save, file = "results/supplementary_tables/S1_rhythmic_transcripts_in_PROTECTED.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
## Table S1 - rhythmic transcripts in PROTECTED



filename <- paste0("results/RDS_files/first_pass_results_EXPOSED.RDS")
results_EXPOSED <- readRDS(filename)


sig_results_EXPOSED <- results_EXPOSED[results_EXPOSED$bh_q <= bhq_thres & results_EXPOSED$rel_amp >= rel_amp_thres,]


results_to_save <- data.frame(gene_id = sig_results_EXPOSED$gene_id, gene_name = sig_results_EXPOSED$gene_name,
                              rel_amp = sig_results_EXPOSED$rel_amp,
                              abs_amp = sig_results_EXPOSED$abs_amp,
                              phase_hour = sig_results_EXPOSED$phase_hour,
                              average_mesor = sig_results_EXPOSED$average_mesor,
                              p = sig_results_EXPOSED$p,
                              bh_q = sig_results_EXPOSED$bh_q)




write.table(results_to_save, file = "results/supplementary_tables/S2_rhythmic_transcripts_in_EXPOSED.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
## Table S2 - rhythmic transcripts in EXPOSED


overlap_ids <- sig_results_EXPOSED$gene_id[sig_results_EXPOSED$gene_id %in% sig_results_PROTECTED$gene_id ]

sig_results_EXPOSED_OVERLAP <- sig_results_EXPOSED[sig_results_EXPOSED$gene_id %in% overlap_ids, ]
sig_results_PROTECTED_OVERLAP <- sig_results_PROTECTED[sig_results_PROTECTED$gene_id %in% overlap_ids, ]

# ## check that they're in the same order
# any(sig_results_EXPOSED_OVERLAP$gene_id != sig_results_PROTECTED_OVERLAP$gene_id) # should be false
# sig_results_EXPOSED_OVERLAP$gene_id[400:405]
# sig_results_PROTECTED_OVERLAP$gene_id[400:405]

results_to_save <- data.frame(gene_id = sig_results_PROTECTED_OVERLAP$gene_id, gene_name = sig_results_PROTECTED_OVERLAP$gene_name,
                              rel_amp_PROTECTED = sig_results_PROTECTED_OVERLAP$rel_amp, rel_amp_EXPOSED = sig_results_EXPOSED_OVERLAP$rel_amp,
                              abs_amp_PROTECTED = sig_results_PROTECTED_OVERLAP$abs_amp, abs_amp_EXPOSED = sig_results_EXPOSED_OVERLAP$abs_amp,
                              phase_hour_PROTECTED = sig_results_PROTECTED_OVERLAP$phase_hour, phase_hour_EXPOSED = sig_results_EXPOSED_OVERLAP$phase_hour,
                              average_mesor_PROTECTED = sig_results_PROTECTED_OVERLAP$average_mesor, average_mesor_EXPOSED = sig_results_EXPOSED_OVERLAP$average_mesor,
                              p_PROTECTED = sig_results_PROTECTED_OVERLAP$p, p_EXPOSED = sig_results_EXPOSED_OVERLAP$p,
                              bh_q_PROTECTED = sig_results_PROTECTED_OVERLAP$bh_q, bh_q_EXPOSED = sig_results_EXPOSED_OVERLAP$bh_q)


write.table(results_to_save, file = "results/supplementary_tables/S3_rhythmic_transcripts_OVERLAP.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
## Table S3 - rhythmic transcripts OVERLAP


##########################################################################################################################################
## Table S4

## clear environment
rm(list=ls())

DR_ids <- readRDS("results/RDS_files/DR_ids.RDS")
results_PROTECTED <- readRDS("results/RDS_files/first_pass_results_PROTECTED.RDS")
results_EXPOSED <- readRDS("results/RDS_files/first_pass_results_EXPOSED.RDS")

results_EXPOSED_DR <- results_EXPOSED[results_EXPOSED$gene_id %in% DR_ids, ]
results_PROTECTED_DR <- results_PROTECTED[results_PROTECTED$gene_id %in% DR_ids, ]

results_to_save <- data.frame(gene_id = results_PROTECTED_DR$gene_id, gene_name = results_PROTECTED_DR$gene_name,
                              rel_amp_PROTECTED = results_PROTECTED_DR$rel_amp, rel_amp_EXPOSED = results_EXPOSED_DR$rel_amp,
                              abs_amp_PROTECTED = results_PROTECTED_DR$abs_amp, abs_amp_EXPOSED = results_EXPOSED_DR$abs_amp,
                              phase_hour_PROTECTED = results_PROTECTED_DR$phase_hour, phase_hour_EXPOSED = results_EXPOSED_DR$phase_hour,
                              average_mesor_PROTECTED = results_PROTECTED_DR$average_mesor, average_mesor_EXPOSED = results_EXPOSED_DR$average_mesor,
                              p_PROTECTED = results_PROTECTED_DR$p, p_EXPOSED = results_EXPOSED_DR$p,
                              bh_q_PROTECTED = results_PROTECTED_DR$bh_q, bh_q_EXPOSED = results_EXPOSED_DR$bh_q)

write.table(results_to_save, file = "results/supplementary_tables/S4_differentially_rhythmic_transcripts.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
## Table S4 - differentially rhythmic transcripts




##########################################################################################################################################
## Table S5-11


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



DR_ids <- readRDS("results/RDS_files/DR_ids.RDS")
id_to_name <- readRDS("results/RDS_files/id_to_name.RDS")
rhythmic_in_either_ids <- readRDS("results/RDS_files/rhythmic_in_either_ids.RDS")
id_to_name_RHYTHMIC <- id_to_name[id_to_name$gene_id %in% rhythmic_in_either_ids,]

results_PROTECTED <- readRDS("results/RDS_files/first_pass_results_PROTECTED.RDS")
results_EXPOSED <- readRDS("results/RDS_files/first_pass_results_EXPOSED.RDS")

results_PROTECTED <- results_PROTECTED[results_PROTECTED$gene_id %in% rhythmic_in_either_ids, ]
results_EXPOSED <- results_EXPOSED[results_EXPOSED$gene_id %in% rhythmic_in_either_ids, ]

DR_check <- results_PROTECTED$gene_id %in% DR_ids


full_results <- data.frame(gene_id = results_PROTECTED$gene_id, gene_name = results_PROTECTED$gene_name,
                           rel_amp_PROTECTED = results_PROTECTED$rel_amp, rel_amp_EXPOSED = results_EXPOSED$rel_amp, 
                           abs_amp_PROTECTED = results_PROTECTED$abs_amp, abs_amp_EXPOSED = results_EXPOSED$abs_amp, 
                           phase_hour_PROTECTED = results_PROTECTED$phase_hour, phase_hour_EXPOSED = results_EXPOSED$phase_hour,
                           average_mesor_PROTECTED = results_PROTECTED$average_mesor, average_mesor_EXPOSED = results_EXPOSED$average_mesor, 
                           p_PROTECTED = results_PROTECTED$p, p_EXPOSED = results_EXPOSED$p, 
                           bh_q_PROTECTED = results_PROTECTED$bh_q, bh_q_EXPOSED = results_EXPOSED$bh_q,
                           differentially_rhythmic = DR_check)





coreclock_genes <- c( "NR1D2", "PER1", "PER2", "CRY2", "RORC", "BHLHE40", "PER3", "CRY1", "NPAS2", "NFIL3", "CIART", "TEF", "BMAL1", "CLOCK", "NR1D1", "DBP", "HLF", "RORA", "BHLHE41", "CSNK1D")
coreclock_ids <- id_to_name_RHYTHMIC$gene_id[id_to_name_RHYTHMIC$gene_name %in% coreclock_genes]
coreclock_results <- full_results[full_results$gene_id %in% coreclock_ids,]

write.table(coreclock_results, file = "results/supplementary_tables/S5_core_clock.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
## Table S5 - core clock genes



gene_sets_HALLMARK <- read_gmt("pathway_analysis_inputs/h.all.v2024.1.Hs.symbols.gmt")

term_name <- "HALLMARK_DNA_REPAIR"
term <- gene_sets_HALLMARK %>% keep(~ .x$name == term_name)
term_genes <- term[[1]]$genes
dna_repair_ids <- id_to_name_RHYTHMIC$gene_id[id_to_name_RHYTHMIC$gene_name %in% term_genes]

dna_repair_results <- full_results[full_results$gene_id %in% dna_repair_ids,]
write.table(dna_repair_results, file = "results/supplementary_tables/S6_DNA_repair.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
## Table S6 - DNA repair genes




gene_sets_TF <- read_gmt("pathway_analysis_inputs/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")

term_name <- "MYC CHEA"
term1 <- gene_sets_TF %>% keep(~ .x$name == term_name)
term_genes1 <- term1[[1]]$genes
term_name <- "MYC ENCODE"
term2 <- gene_sets_TF %>% keep(~ .x$name == term_name)
term_genes2 <- term2[[1]]$genes
term_genes <- unique(c(term_genes1, term_genes2))

myc_target_ids <- id_to_name_RHYTHMIC$gene_id[id_to_name_RHYTHMIC$gene_name %in% term_genes]

myc_targets_results <- full_results[full_results$gene_id %in% myc_target_ids,]
write.table(myc_targets_results, file = "results/supplementary_tables/S7_MYC_targets.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
## Table S7 - MYC targets



term_name <- "E2F1 CHEA"
term1 <- gene_sets_TF %>% keep(~ .x$name == term_name)
term_genes1 <- term1[[1]]$genes

term_name <- "E2F4 ENCODE"
term2 <- gene_sets_TF %>% keep(~ .x$name == term_name)
term_genes2 <- term2[[1]]$genes

term_name <- "E2F6 ENCODE"
term3 <- gene_sets_TF %>% keep(~ .x$name == term_name)
term_genes3 <- term3[[1]]$genes

term_genes <- unique(c(term_genes1, term_genes2, term_genes3))
e2f_target_ids <- id_to_name_RHYTHMIC$gene_id[id_to_name_RHYTHMIC$gene_name %in% term_genes]

e2f_targets_results <- full_results[full_results$gene_id %in% e2f_target_ids,]
write.table(e2f_targets_results, file = "results/supplementary_tables/S8_E2F_targets.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
## Table S8 - E2F targets


term_name <- "HALLMARK_G2M_CHECKPOINT"
term <- gene_sets_HALLMARK %>% keep(~ .x$name == term_name)
term_genes <- term[[1]]$genes
g2m_checkpoint_ids <- id_to_name_RHYTHMIC$gene_id[id_to_name_RHYTHMIC$gene_name %in% term_genes]

g2m_checkpoint_results <- full_results[full_results$gene_id %in% g2m_checkpoint_ids,]
write.table(g2m_checkpoint_results, file = "results/supplementary_tables/S9_G2M_checkpoint.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
## Table S9 - G2M checkpoint


term_name <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
term <- gene_sets_HALLMARK %>% keep(~ .x$name == term_name)
term_genes <- term[[1]]$genes
emt_ids <- id_to_name_RHYTHMIC$gene_id[id_to_name_RHYTHMIC$gene_name %in% term_genes]

emt_results <- full_results[full_results$gene_id %in% emt_ids,]
write.table(emt_results, file = "results/supplementary_tables/S10_EMT.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
## Table S10 - EMT genes


term_name <- "HALLMARK_APICAL_JUNCTION"
term <- gene_sets_HALLMARK %>% keep(~ .x$name == term_name)
term_genes <- term[[1]]$genes
apical_junction_ids <- id_to_name_RHYTHMIC$gene_id[id_to_name_RHYTHMIC$gene_name %in% term_genes]

apical_junction_results <- full_results[full_results$gene_id %in% apical_junction_ids,]
write.table(apical_junction_results, file = "results/supplementary_tables/S11_apical_junction.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
## Table S11 - apical junction genes










