## clear environment
rm(list=ls()) 


data_label <- "EXPOSED"
# data_label <- "PROTECTED"



##########################################################################################################################################
## read in data


filename <- paste("processed_data/times_",data_label,".RDS",sep="")
times <- readRDS(filename)

filename <- paste("processed_data/ids_",data_label,".RDS",sep="")
ids <- readRDS(filename)

filename <- paste("processed_data/id_factors_",data_label,".RDS",sep="")
id_factors <- readRDS(filename)

filename <- paste("processed_data/gene_data_",data_label,".RDS",sep="")
data <- readRDS(filename)


##########################################################################################################################################

## cosinor regression 

p_values <- c()
amps <- c()
average_mesors <- c()
rel_amps <- c()

phases_rad <- c()
phases_hour <- c()



for (gene_ind in 1:dim(data)[1]){
  
  expression <- unlist(as.vector(data[gene_ind,]))
  
  trig_model <- lm(expression ~ cos(times*pi/12) + sin(times*pi/12) +id_factors) 
  base_model <- lm(expression ~ id_factors) 
  
  anova_results <- anova(trig_model,base_model)
  p_value <- anova_results[2,6]
  p_values <- append(p_values, p_value)
  
  
  cos_fit <- unname(trig_model$coefficients[2])
  sin_fit <- unname(trig_model$coefficients[3])

  

  amp_fit <- unname(sqrt(cos_fit**2 + sin_fit**2))
  
  amps <- append(amps, amp_fit)
  
  
  c <- trig_model$coefficients[1] # fit mean of first subject
  subject_factors <- trig_model$coefficients[4:length(trig_model$coefficients)] # deviations from c for every other subject
  subject_mesors <- subject_factors + c
  subject_mesors <- append(c,subject_mesors)
  average_mesor <- mean(subject_mesors)
  
  average_mesors <- append(average_mesors, average_mesor)
  
  
  
  rel_amp <- amp_fit / average_mesor
  rel_amps <- append(rel_amps, rel_amp)
  
  # phase fit
  phase_rad_fit <- atan2(sin_fit,cos_fit)
  
  
  ## convert radians to hours on 24 hour clock
  phase_hour_fit <- unname(phase_rad_fit * (24 / (2 * pi)) )
  
  if(phase_hour_fit < 0){
    phase_hour_fit <- 24 + phase_hour_fit
  }
  
  phases_rad <- append(phases_rad, phase_rad_fit)
  phases_hour <- append(phases_hour, phase_hour_fit)
  
  
  
}



bh_q <- p.adjust(p_values,method="BH")


##########################################################################################################################################

# get all gene names


# install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")


library(org.Hs.eg.db)
library(dplyr)


gene_ids <- rownames(data)

gene_names <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL",keytype = "ENSEMBL")
gene_names <- unname(gene_names)

gene_names[is.na(gene_names)] <- "na"


##########################################################################################################################################

results_to_save <- data.frame(gene_id = gene_ids, gene_name = gene_names,
                              abs_amp = amps, rel_amp = rel_amps, average_mesor = average_mesors, phase_rad = phases_rad,
                              phase_hour = phases_hour, p = p_values, bh_q = bh_q)


filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
saveRDS(results_to_save, file = filename)  


id_to_name <- results_to_save[,c(1,2)]
filename <- paste0("results/RDS_files/id_to_name.RDS")
saveRDS(id_to_name, file = filename)  






stop("done")




##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## write to results summary files


## clear environment
rm(list=ls()) 


data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)
write.table(results_EXPOSED, file = "results/first_pass_cosinor_results/all_EXPOSED.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


id_to_name <- results_EXPOSED[,c(1,2)]
filename <- paste0("results/RDS_files/id_to_name.RDS")
saveRDS(id_to_name, file = filename)  


data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)
write.table(results_PROTECTED, file = "results/first_pass_cosinor_results/all_PROTECTED.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




rel_amp_thres <- 0.2
bhq_thres <- 0.05


sig_results_EXPOSED <- results_EXPOSED[results_EXPOSED$bh_q <= bhq_thres & results_EXPOSED$rel_amp >= rel_amp_thres,]
write.table(sig_results_EXPOSED, file = "results/first_pass_cosinor_results/rhythmic_EXPOSED.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# checks
dim(sig_results_EXPOSED)[1]
dim(sig_results_EXPOSED)[1] / dim(results_EXPOSED)[1]




sig_results_PROTECTED <- results_PROTECTED[results_PROTECTED$bh_q <= bhq_thres & results_PROTECTED$rel_amp >= rel_amp_thres,]
write.table(sig_results_PROTECTED, file = "results/first_pass_cosinor_results/rhythmic_PROTECTED.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# checks
dim(sig_results_PROTECTED)[1]
dim(sig_results_PROTECTED)[1] / dim(results_PROTECTED)[1]



overlap_ids <- sig_results_EXPOSED$gene_id[sig_results_EXPOSED$gene_id %in% sig_results_PROTECTED$gene_id ]

filename <- paste0("results/RDS_files/overlap_ids.RDS")
saveRDS(overlap_ids, file = filename)  



rhythmic_in_either_ids <- unique(c(sig_results_EXPOSED$gene_id, sig_results_PROTECTED$gene_id))
filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
saveRDS(rhythmic_in_either_ids, file = filename)  


#######################
# OVERLAP


sig_results_EXPOSED_OVERLAP <- sig_results_EXPOSED[sig_results_EXPOSED$gene_id %in% overlap_ids, ]
sig_results_PROTECTED_OVERLAP <- sig_results_PROTECTED[sig_results_PROTECTED$gene_id %in% overlap_ids, ]


# ## check that they're in the same order
# any(sig_results_EXPOSED_OVERLAP$gene_id != sig_results_PROTECTED_OVERLAP$gene_id) # should be false
# sig_results_EXPOSED_OVERLAP$gene_id[400:405]
# sig_results_PROTECTED_OVERLAP$gene_id[400:405]



results_to_save <- data.frame(gene_id = sig_results_EXPOSED_OVERLAP$gene_id, gene_name = sig_results_EXPOSED_OVERLAP$gene_name,
                              abs_amp_EXPOSED = sig_results_EXPOSED_OVERLAP$abs_amp, abs_amp_PROTECTED = sig_results_PROTECTED_OVERLAP$abs_amp,
                              rel_amp_EXPOSED = sig_results_EXPOSED_OVERLAP$rel_amp, rel_amp_PROTECTED = sig_results_PROTECTED_OVERLAP$rel_amp,
                              phase_rad_EXPOSED = sig_results_EXPOSED_OVERLAP$phase_rad, phase_rad_PROTECTED = sig_results_PROTECTED_OVERLAP$phase_rad,
                              phase_hour_EXPOSED = sig_results_EXPOSED_OVERLAP$phase_hour, phase_hour_PROTECTED = sig_results_PROTECTED_OVERLAP$phase_hour,
                              average_mesor_EXPOSED = sig_results_EXPOSED_OVERLAP$average_mesor, average_mesor_PROTECTED = sig_results_PROTECTED_OVERLAP$average_mesor,
                              p_EXPOSED = sig_results_EXPOSED_OVERLAP$p, p_PROTECTED = sig_results_PROTECTED_OVERLAP$p,
                              bh_q_EXPOSED = sig_results_EXPOSED_OVERLAP$bh_q, bh_q_PROTECTED = sig_results_PROTECTED_OVERLAP$bh_q)


write.table(results_to_save, file = "results/first_pass_cosinor_results/rhythmic_OVERLAP.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

