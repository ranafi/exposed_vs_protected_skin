

## clear environment
rm(list=ls()) 


##########################################################################################################################################
# read data

data_label <- "EXPOSED"

filename <- paste("processed_data/times_",data_label,".RDS",sep="")
times_EXPOSED <- readRDS(filename)

filename <- paste("processed_data/ids_",data_label,".RDS",sep="")
ids_EXPOSED <- readRDS(filename)

filename <- paste("processed_data/gene_data_",data_label,".RDS",sep="")
data_EXPOSED <- readRDS(filename)




data_label <- "PROTECTED"

filename <- paste("processed_data/times_",data_label,".RDS",sep="")
times_PROTECTED <- readRDS(filename)

filename <- paste("processed_data/ids_",data_label,".RDS",sep="")
ids_PROTECTED <- readRDS(filename)

filename <- paste("processed_data/gene_data_",data_label,".RDS",sep="")
data_PROTECTED <- readRDS(filename)



##########################################################################################################################################

# combine

data <- merge(data_EXPOSED, data_PROTECTED, by = 'row.names', all = FALSE)
rownames(data) <- data$Row.names
data$Row.names <- NULL


conditions <- c( rep("EXPOSED",80), rep("PROTECTED",80))

ids <- c(ids_EXPOSED, ids_PROTECTED)

id_factors <- factor(ids)

times <- c(times_EXPOSED, times_PROTECTED)




##########################################################################################################################################

# only keep genes that came up cycling in one or the other



filename <- paste0("results/RDS_files/rhythmic_in_either_ids.RDS")
rhythmic_in_either_ids <- readRDS(filename)

data <- data[rownames(data) %in% rhythmic_in_either_ids,]

# gene_names <- id_to_name$gene_name[id_to_name$gene_id %in% rhythmic_in_either_ids]


##########################################################################################################################################

# cosinor analysis with indicator variable



indicator <- as.numeric(conditions == "PROTECTED")


p_values <- c()



for (gene_ind in 1:dim(data)[1]){
  
  expression <- unlist(as.vector(data[gene_ind,]))
  
  sin_term <- sin(times * pi/12)
  cos_term <- cos(times * pi/12)
  
  full_model <- lm(expression ~ 
                     sin_term + cos_term +
                     indicator * sin_term + 
                     indicator * cos_term + 
                     indicator + id_factors)
  
  reduced_model <- lm(expression ~ 
                        sin_term + cos_term +
                        indicator + id_factors)
  
  
  
  
  anova_result <- anova(reduced_model, full_model)
  
  p_value <- anova_result[2, "Pr(>F)"]
  p_values <- append(p_values, p_value)

  
}


bh_q <- p.adjust(p_values,method="BH")


##########################################################################################################################################

DR_ids <- rownames(data)[bh_q <= 0.1]

filename <- paste0("results/RDS_files/DR_ids.RDS")
saveRDS(DR_ids, file = filename)  


data_label <- "EXPOSED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_EXPOSED <- readRDS(filename)


data_label <- "PROTECTED"
filename <- paste0("results/RDS_files/first_pass_results_",data_label,".RDS")
results_PROTECTED <- readRDS(filename)


results_EXPOSED_DR <- results_EXPOSED[results_EXPOSED$gene_id %in% DR_ids, ]
results_PROTECTED_DR <- results_PROTECTED[results_PROTECTED$gene_id %in% DR_ids, ]








results_to_save <- data.frame(gene_id = results_EXPOSED_DR$gene_id, gene_name = results_PROTECTED_DR$gene_name,
                              abs_amp_EXPOSED = results_EXPOSED_DR$abs_amp, abs_amp_PROTECTED = results_PROTECTED_DR$abs_amp,
                              rel_amp_EXPOSED = results_EXPOSED_DR$rel_amp, rel_amp_PROTECTED = results_PROTECTED_DR$rel_amp,
                              phase_rad_EXPOSED = results_EXPOSED_DR$phase_rad, phase_rad_PROTECTED = results_PROTECTED_DR$phase_rad,
                              phase_hour_EXPOSED = results_EXPOSED_DR$phase_hour, phase_hour_PROTECTED = results_PROTECTED_DR$phase_hour,
                              average_mesor_EXPOSED = results_EXPOSED_DR$average_mesor, average_mesor_PROTECTED = results_PROTECTED_DR$average_mesor,
                              p_EXPOSED = results_EXPOSED_DR$p, p_PROTECTED = results_PROTECTED_DR$p,
                              bh_q_EXPOSED = results_EXPOSED_DR$bh_q, bh_q_PROTECTED = results_PROTECTED_DR$bh_q)


write.table(results_to_save, file = "results/differential_cosinor_results/differentially_rhythmic_results.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





