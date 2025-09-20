## clear environment
rm(list=ls()) 




# install.packages("BiocManager")
# BiocManager::install("edgeR")
library(edgeR)

# needed for ComBat-seq
# install.packages("BiocManager")
# BiocManager::install("sva")
library(sva)


#######################
# sample label info

# photoprotected
# A, C, E, G
# 12, 18, 0, 6

# photoexposed
# B, D, F, H
# 12, 18, 0, 6
#######################




##########################################################################################################################################
## read in raw counts data

data_0 <- read.table('counts_data/Time 00 UV protected vs UV exposed/raw_counts T00.csv', header = TRUE, sep =',',row.names=1)
colnames(data_0) <- sub("X","",colnames(data_0) )
PROTECTED_0 <- data_0[, seq(1, ncol(data_0), by = 2)]
EXPOSED_0 <- data_0[, seq(2, ncol(data_0), by = 2)]


data_6 <- read.table('counts_data/Time 06 UV protected vs UV exposed/raw_counts T06.csv', header = TRUE, sep =',',row.names=1)
colnames(data_6) <- sub("X","",colnames(data_6) )
PROTECTED_6 <- data_6[, seq(1, ncol(data_6), by = 2)]
EXPOSED_6 <- data_6[, seq(2, ncol(data_6), by = 2)]


data_12 <- read.table('counts_data/Time 12 UV protected vs UV exposed/raw_counts T12.csv', header = TRUE, sep =',',row.names=1)
colnames(data_12) <- sub("X","",colnames(data_12) )
PROTECTED_12 <- data_12[, seq(1, ncol(data_12), by = 2)]
EXPOSED_12 <- data_12[, seq(2, ncol(data_12), by = 2)]


data_18 <- read.table('counts_data/Time 18 UV protected vs UV exposed/raw_counts T18.csv', header = TRUE, sep =',',row.names=1)
colnames(data_18) <- sub("X","",colnames(data_18) )
PROTECTED_18 <- data_18[, seq(1, ncol(data_18), by = 2)]
EXPOSED_18 <- data_18[, seq(2, ncol(data_18), by = 2)]



##########################################################################################################################################

# merge PROTECTED timepoints

data_PROTECTED <- merge(PROTECTED_0, PROTECTED_6, by = 'row.names', all = FALSE)
rownames(data_PROTECTED) <- data_PROTECTED$Row.names
data_PROTECTED$Row.names <- NULL

data_PROTECTED <- merge(data_PROTECTED, PROTECTED_12, by = 'row.names', all = FALSE)
rownames(data_PROTECTED) <- data_PROTECTED$Row.names
data_PROTECTED$Row.names <- NULL

data_PROTECTED <- merge(data_PROTECTED, PROTECTED_18, by = 'row.names', all = FALSE)
rownames(data_PROTECTED) <- data_PROTECTED$Row.names
data_PROTECTED$Row.names <- NULL

times_PROTECTED <- c( rep(0,20), rep(6,20), rep(12,20), rep(18,20))
subs_PROTECTED <- rep( seq(1,20),4 )
sub_factors_PROTECTED <- factor(subs_PROTECTED)


# merge EXPOSED timepoints

data_EXPOSED <- merge(EXPOSED_0, EXPOSED_6, by = 'row.names', all = FALSE)
rownames(data_EXPOSED) <- data_EXPOSED$Row.names
data_EXPOSED$Row.names <- NULL

data_EXPOSED <- merge(data_EXPOSED, EXPOSED_12, by = 'row.names', all = FALSE)
rownames(data_EXPOSED) <- data_EXPOSED$Row.names
data_EXPOSED$Row.names <- NULL

data_EXPOSED <- merge(data_EXPOSED, EXPOSED_18, by = 'row.names', all = FALSE)
rownames(data_EXPOSED) <- data_EXPOSED$Row.names
data_EXPOSED$Row.names <- NULL

times_EXPOSED <- c( rep(0,20), rep(6,20), rep(12,20), rep(18,20))
subs_EXPOSED <- rep( seq(1,20),4 )
sub_factors_EXPOSED <- factor(subs_EXPOSED)



##########################################################################################################################################
## combine exposed and protected to normalize together

# check and make sure rownames match -- should be false
# any(rownames(data_exposed) != rownames(data_protected))

data <- merge(data_EXPOSED, data_PROTECTED, by = 'row.names', all = FALSE)
rownames(data) <- data$Row.names
data$Row.names <- NULL


skin_types <- c( rep("EXPOSED",80), rep("PROTECTED",80))





##########################################################################################################################################
# read in batch info, and reorder according to order of samples in data


batch_info <- read.table('counts_data/batch_info.csv', header = TRUE, sep =',',row.names=1)
sample_ids <- colnames(data)
batch_info_ordered <- batch_info[sample_ids,]



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## data processing steps

# start with raw counts
# filter low count genes (using raw counts)
# Combat-seq batch correction
# EdgeR normalization with batch-corrected counts (TMM to calculate effective library sizes, then CPM)
# result - unlogged CPM data
# mean expression filter >= 5


#########################################################################################################
## sources:

# ComBat-seq
# https://github.com/zhangyuqing/ComBat-seq
# biological variable (skin condition, exposed or protected) should be included as "group", and signal will be preserved
# "In ComBat-seq, user may specify biological covariates, whose signals will be preserved in the adjusted data. 
# If the user would like to specify one biological variable, they may use the group parameter"
# confirmed with real-data example from repo: https://github.com/zhangyuqing/ComBat-seq/blob/master/real_data_application/gfrn_DE.R


# https://github.com/zhangyuqing/ComBat-seq/issues/3
# Combat-seq is appropriate to use before EdgeR normalizaiton.,
# filtering out low expression genes should be done before Combat-seq batch correction




#########################################################################################################



# filter low expression genes
keep <- filterByExpr(data)
filtered_counts <- data[keep,]
# before: 57500 genes
# after: 17102 genes



batch_corrected_counts <- ComBat_seq(
  counts = as.matrix(filtered_counts),
  batch = batch_info_ordered$Library_preparation_batch,
  group = factor(skin_types)
)




dge_batch_corrected <- DGEList(counts = batch_corrected_counts)
dge_batch_corrected <- calcNormFactors(dge_batch_corrected)
# TMM done at this point
# result is dge_batch_corrected$samples$norm.factors, which is then used to calculate effective library sizes for CPM normalization



CPM_batch_corrected <- cpm(dge_batch_corrected, log=FALSE)
# log=FALSE is actually the default here, but I wrote it explicitly anyway as a reminder

# colSums(CPM) # columns sum to exactly a million if we don't do the calcNormFactors(), close to a million if we do calcNormFactors() because of TMM adjustment

# sum(CPM_batch_corrected<0) # negative number check - should be 0


# mean expression filter
CPM_batch_corrected <- CPM_batch_corrected[rowMeans(CPM_batch_corrected) >=5, ]
# before: 17102 genes
# after: 12471 genes


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################


# separate into EXPOSED and PROTECTED again, and save


cpm_EXPOSED <- CPM_batch_corrected[,skin_types=="EXPOSED"]
cpm_PROTECTED <- CPM_batch_corrected[,skin_types=="PROTECTED"]


saveRDS(times_PROTECTED, file = "processed_data/times_PROTECTED.RDS")  
saveRDS(subs_PROTECTED, file = "processed_data/ids_PROTECTED.RDS")  
saveRDS(sub_factors_PROTECTED, file = "processed_data/id_factors_PROTECTED.RDS")  
saveRDS(cpm_PROTECTED, file = "processed_data/gene_data_PROTECTED.RDS")  


saveRDS(times_EXPOSED, file = "processed_data/times_EXPOSED.RDS")  
saveRDS(subs_EXPOSED, file = "processed_data/ids_EXPOSED.RDS")  
saveRDS(sub_factors_EXPOSED, file = "processed_data/id_factors_EXPOSED.RDS")  
saveRDS(cpm_EXPOSED, file = "processed_data/gene_data_EXPOSED.RDS")  
