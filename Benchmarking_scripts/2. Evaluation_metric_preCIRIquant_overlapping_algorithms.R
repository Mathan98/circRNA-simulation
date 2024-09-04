
### Calculate the evaluation metrics from the overlapping output of the standalone tools (Overlapping algorithm)


library(dplyr)
library(rtracklayer)

## The folders below are in the ( Overlap_CIRI2_No_Stringent folder ) in the main directory
# For the overlapping, we use the results from CIRI2 No Stringent analysis


                                                ### Folders
# CIRI2-CIRCexplorer2 Individual                                              -> Overlap_CIRI2_CE2
# CIRI2-Find_circ Individual                                                  -> Overlap_CIRI2_FC
# Find_circ-CIRCexplorer2 Individual                                          -> Overlap_FC_CE2
# CIRI2-CIRCexplorer2-find_circ_2 Individual (Detected by more than 2 tools)  -> Overlap_CIRI2_CE2_FC_2
# CIRI2-CIRCexplorer2-find_circ_3 Individual (Detected by all 3 tools)        -> Overlap_CIRI2_CE2_FC_3



# Superdepth_Overlap_CIRI2_CIRCexplorer2              -> Superdepth_Overlap_CIRI2_CE2
# Superdepth_Overlap_CIRI2_find_circ                  -> Superdepth_Overlap_CIRI2_FC
# Superdepth_Overlap_find_circ_CIRCexplorer2          -> Superdepth_Overlap_FC_CE2
# Superdepth_Overlap_CIRI2_CIRCexplorer2_find_circ_2  -> Superdepth_Overlap_CIRI2_CE2_FC_2
# Superdepth_Overlap_CIRI2_CIRCexplorer2_find_circ_3  -> Superdepth_Overlap_CIRI2_CE2_FC_3


## Path to directories containing the folders above and the true positive (tsv files)
# E.g. coverage circ 2x lin 100x

circRNA_sim_dir <- 'C:/Users/matha/Documents/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage2_100'
circRNA_predict_dir <- 'C:/Users/matha/Documents/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage2_100/'


#############################################################################################################
##### Evaluation metrics

# A function defined to calculate precision, recall, F1 score and other metrics
calculate_metrics <- function(conf_matrix) {
  
  TP = confusion_table[2,2]
  TN = confusion_table[1,1]
  FN = confusion_table[2,1]
  FP = confusion_table[1,2]
  precision = TP / (TP + FP)
  recall = TP / (TP + FN)
  FDR = FP/(FP + TP)
  F1_score = (2*(precision * recall)/(precision + recall))
  data.frame(TP = TP,
             TN = TN,
             FP = FP,
             FN = FN,
             precision = precision, 
             recall_sensitivity = recall,
             FDR = FDR,
             F1_score = F1_score,
             TP_Percent = TP / sum(TP, FP),
             TN_Percent = 0,
             #TN_Percent = TN / nrow(known_negative_gene),
             FP_Percent = FP / sum(TP, FP),
             FN_Percent = FN / nrow(circRNAs)
  )
  
}

#####################################################################################################################################
#################

## CIRI2-CIRCexplorer2 Individual

metrics_CIRI2_CE2 <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/Overlap_CIRI2_No_Stringent/Overlap_CIRI2_CE2'), pattern = "bed")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, sep = '\t')
  
  if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F,sep = ',')}
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    ## Predicted
    predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/Overlap_CIRI2_No_Stringent/Overlap_CIRI2_CE2/Sim_', i, '.bed'), header = F, sep = '\t')
    predicted_positive$V4 <- trimws(predicted_positive$V4)
    predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$V4)
    
    
    
    # Full set of possible IDs for gene_ids and circRNA_ids
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circRNA_ID)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circRNA_ID)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circRNA_ID, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circRNA_ID)
    
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_CE2[[paste0('CIRI2|CE2_Ind_', i)]] <- calculate_metrics(confusion_table)
}


# Combine the results into a single data frame
metrics_df_CIRI2_CE2_Ind <- do.call(rbind, metrics_CIRI2_CE2)


#####################################################################################################################################

## CIRI2-Find_circ Individual

metrics_CIRI2_FC <- list()

for (i in sprintf("%03d", 001:016)) {
  
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, sep = '\t')
  
  if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F,sep = ',')}
  
  ## Predicted
  predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/Overlap_CIRI2_No_Stringent/Overlap_CIRI2_FC/Sim_', i, '.bed'), header = F, sep = '\t')
  predicted_positive$V4 <- trimws(predicted_positive$V4)
  predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$V4)
  
  
  
  # Full set of possible IDs for gene_ids and circRNA_ids
  all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circRNA_ID)
  
  
  true_positive <- intersect(circRNAs$V1, predicted_positive$circRNA_ID)
  true_negative <- integer()
  false_positive <- setdiff(predicted_positive$circRNA_ID, circRNAs$V1)
  false_negative <- setdiff(circRNAs$V1, predicted_positive$circRNA_ID)
  
  
  confusion_table = matrix(nrow = 2, ncol = 2)
  if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
  if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
  if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
  if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
  
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_FC[[paste0('CIRI2|FC_Ind_', i)]] <- calculate_metrics(confusion_table)
}


# Combine the results into a single data frame
metrics_df_CIRI2_FC_Ind <- do.call(rbind, metrics_CIRI2_FC)


#####################################################################################################################################

## Find_circ-CIRCexplorer2 Individual 

metrics_FC_CE2 <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/Overlap_CIRI2_No_Stringent/Overlap_FC_CE2'), pattern = "bed")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, sep = '\t')
  
  if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F,sep = ',')}
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    ## Predicted
    predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/Overlap_CIRI2_No_Stringent/Overlap_FC_CE2/Sim_', i, '.bed'), header = F, sep = '\t')
    predicted_positive$V4 <- trimws(predicted_positive$V4)
    predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$V4)
    
    
    
    # Full set of possible IDs for gene_ids and circRNA_ids
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circRNA_ID)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circRNA_ID)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circRNA_ID, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circRNA_ID)
    
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_FC_CE2[[paste0('FC|CE2_Ind_', i)]] <- calculate_metrics(confusion_table)
}


# Combine the results into a single data frame
metrics_df_FC_CE2_Ind <- do.call(rbind, metrics_FC_CE2)


#####################################################################################################################################

## CIRI2-CIRCexplorer2-find_circ_2 Individual

metrics_CIRI2_CE2_FC_2 <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/Overlap_CIRI2_No_Stringent/Overlap_CIRI2_CE2_FC_2'), pattern = "bed")
lf <- gsub("[^0-9]", "", lf)


for (i in sprintf("%03d", 001:016)) {
  
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, sep = '\t')
  
  if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F,sep = ',')}
  
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    ## Predicted
    predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/Overlap_CIRI2_No_Stringent/Overlap_CIRI2_CE2_FC_2/Sim_', i, '.bed'), header = F, sep = '\t')
    predicted_positive$V4 <- trimws(predicted_positive$V4)
    predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$V4)
    
    
    
    # Full set of possible IDs for gene_ids and circRNA_ids
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circRNA_ID)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circRNA_ID)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circRNA_ID, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circRNA_ID)
    
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_CE2_FC_2[[paste0('CIRI2|CE2|FC|2_Ind_', i)]] <- calculate_metrics(confusion_table)
}


# Combine the results into a single data frame
metrics_df_CIRI2_CE2_FC_2_Ind <- do.call(rbind, metrics_CIRI2_CE2_FC_2)


#####################################################################################################################################

## CIRI2-CIRCexplorer2-find_circ_3 Individual

metrics_CIRI2_CE2_FC_3 <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/Overlap_CIRI2_No_Stringent/Overlap_CIRI2_CE2_FC_3'), pattern = "bed")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, sep = '\t')
  
  if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F,sep = ',')}
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    ## Predicted
    predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/Overlap_CIRI2_No_Stringent/Overlap_CIRI2_CE2_FC_3/Sim_', i, '.bed'), header = F, sep = '\t')
    predicted_positive$V4 <- trimws(predicted_positive$V4)
    predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$V4)
    
    
    
    # Full set of possible IDs for gene_ids and circRNA_ids
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circRNA_ID)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circRNA_ID)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circRNA_ID, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circRNA_ID)
    
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_CE2_FC_3[[paste0('CIRI2|CE2|FC|3_Ind_', i)]] <- calculate_metrics(confusion_table)
}


# Combine the results into a single data frame
metrics_df_CIRI2_CE2_FC_3_Ind <- do.call(rbind, metrics_CIRI2_CE2_FC_3)


#####################################################################################################################################
#####################################################################################################################################

## Superdepth_Overlap_CIRI2_CIRCexplorer2

circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F, sep = '\t')
if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F,sep = ',')}

# Predicted positive ids
predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/Overlap_CIRI2_No_Stringent/Superdepth_Overlap_CIRI2_CE2/Sim_Superdepth.bed'), header = F, sep = '\t')
predicted_positive$V4 <- trimws(predicted_positive$V4)
predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$V4)


# Full set of possible IDs
all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circRNA_ID)


true_positive <- intersect(circRNAs$V1, predicted_positive$circRNA_ID)
true_negative <- integer()
false_positive <- setdiff(predicted_positive$circRNA_ID, circRNAs$V1)
false_negative <- setdiff(circRNAs$V1, predicted_positive$circRNA_ID)


confusion_table = matrix(nrow = 2, ncol = 2)
if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}


# Loop through the list and calculate metrics for each confusion matrix
metrics_df_CIRI2_CE2_SD <- calculate_metrics(confusion_table)
rownames(metrics_df_CIRI2_CE2_SD) <- "CIRI2|CE2_SD"

#####################################################################################################################################

## Superdepth_Overlap_CIRI2_find_circ

circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F, sep = '\t')
if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F,sep = ',')}

# Predicted positive ids
predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/Overlap_CIRI2_No_Stringent/Superdepth_Overlap_CIRI2_FC/Sim_Superdepth.bed'), header = F, sep = '\t')
predicted_positive$V4 <- trimws(predicted_positive$V4)
predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$V4)


# Full set of possible IDs
all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circRNA_ID)


true_positive <- intersect(circRNAs$V1, predicted_positive$circRNA_ID)
true_negative <- integer()
false_positive <- setdiff(predicted_positive$circRNA_ID, circRNAs$V1)
false_negative <- setdiff(circRNAs$V1, predicted_positive$circRNA_ID)


confusion_table = matrix(nrow = 2, ncol = 2)
if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}


# Loop through the list and calculate metrics for each confusion matrix
metrics_df_CIRI2_FC_SD <- calculate_metrics(confusion_table)
rownames(metrics_df_CIRI2_FC_SD) <- "CIRI2|FC_SD"

#####################################################################################################################################

## Superdepth_Overlap_find_circ_CIRCexplorer2

circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = T, sep = '\t')
if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = T,sep = ',')}


# Predicted positive ids
if (file.size(paste0(circRNA_predict_dir, 'Overlap_CIRI2_No_Stringent/Superdepth_Overlap_FC_CE2/Sim_Superdepth.bed')) == 0){
  
  confusion_table = matrix(nrow = 2, ncol = 2)
  confusion_table[1,2] <- 0
  confusion_table[2,1] <- nrow(circRNAs)
  confusion_table[2,2] <- 0
  confusion_table[1,1] <- 0
  
} else {
  
  predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/Overlap_CIRI2_No_Stringent/Superdepth_Overlap_FC_CE2/Sim_Superdepth.bed'), header = F, sep = '\t')
  predicted_positive$V4 <- trimws(predicted_positive$V4)
  predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$V4)
  
  
  # Full set of possible IDs
  all_ids_circRNA <- c(circRNAs$Circ_ID, predicted_positive$circRNA_ID)
  
  
  true_positive <- intersect(circRNAs$Circ_ID, predicted_positive$circRNA_ID)
  true_negative <- integer()
  false_positive <- setdiff(predicted_positive$circRNA_ID, circRNAs$Circ_ID)
  false_negative <- setdiff(circRNAs$Circ_ID, predicted_positive$circRNA_ID)
  
  
  confusion_table = matrix(nrow = 2, ncol = 2)
  if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
  if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
  if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
  if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
  
}

# Loop through the list and calculate metrics for each confusion matrix
metrics_df_FC_CE2_SD <- calculate_metrics(confusion_table)
rownames(metrics_df_FC_CE2_SD) <- "FC|CE2_SD"

#####################################################################################################################################

## Superdepth_Overlap_CIRI2_CIRCexplorer2_find_circ_2

circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F, sep = '\t')
if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F,sep = '\t')}


# Predicted positive ids
predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/Overlap_CIRI2_No_Stringent/Superdepth_Overlap_CIRI2_CE2_FC_2/Sim_Superdepth.bed'), header = F, sep = '\t')
predicted_positive$V4 <- trimws(predicted_positive$V4)
predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$V4)


# Full set of possible IDs
all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circRNA_ID)


true_positive <- intersect(circRNAs$V1, predicted_positive$circRNA_ID)
true_negative <- integer()
false_positive <- setdiff(predicted_positive$circRNA_ID, circRNAs$V1)
false_negative <- setdiff(circRNAs$V1, predicted_positive$circRNA_ID)


confusion_table = matrix(nrow = 2, ncol = 2)
if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}


# Loop through the list and calculate metrics for each confusion matrix
metrics_df_CIRI2_CE2_FC_2_SD <- calculate_metrics(confusion_table)
rownames(metrics_df_CIRI2_CE2_FC_2_SD) <- "CIRI2|CE2|FC|2_SD"

#####################################################################################################################################

## Superdepth_Overlap_CIRI2_CIRCexplorer2_find_circ_3

circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F, sep = '\t')
if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F,sep = ',')}


# Predicted positive ids
if (file.size(paste0(circRNA_predict_dir, 'Overlap_CIRI2_No_Stringent/Superdepth_Overlap_CIRI2_CE2_FC_3/Sim_Superdepth.bed')) == 0){
  
  confusion_table = matrix(nrow = 2, ncol = 2)
  confusion_table[1,2] <- 0
  confusion_table[2,1] <- nrow(circRNAs)
  confusion_table[2,2] <- 0
  confusion_table[1,1] <- 0
  
} else {
  # Predicted positive ids
  predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/Overlap_CIRI2_No_Stringent/Superdepth_Overlap_CIRI2_CE2_FC_3/Sim_Superdepth.bed'), header = F, sep = '\t')
  predicted_positive$V4 <- trimws(predicted_positive$V4)
  predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$V4)
  
  
  # Full set of possible IDs
  all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circRNA_ID)
  
  
  true_positive <- intersect(circRNAs$V1, predicted_positive$circRNA_ID)
  true_negative <- integer()
  false_positive <- setdiff(predicted_positive$circRNA_ID, circRNAs$V1)
  false_negative <- setdiff(circRNAs$V1, predicted_positive$circRNA_ID)
  
  
  confusion_table = matrix(nrow = 2, ncol = 2)
  if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
  if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
  if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
  if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
}

# Loop through the list and calculate metrics for each confusion matrix
metrics_df_CIRI2_CE2_FC_3_SD <- calculate_metrics(confusion_table)
rownames(metrics_df_CIRI2_CE2_FC_3_SD) <- "CIRI2|CE2|FC|3_SD"

#####################################################################################################################################

### Summarize all evaluation metrics (CIRI2 Superdepth, CIRI2 Individuals,
###                                   CE2 Superdepth, CE2 Individuals, 
###                                   FC Superdepth, FC Individual)

Evaluation_metrics_df <- rbind(metrics_df_CIRI2_CE2_Ind, metrics_df_CIRI2_CE2_SD, 
                               metrics_df_CIRI2_FC_Ind, metrics_df_CIRI2_FC_SD,
                               metrics_df_FC_CE2_Ind, metrics_df_FC_CE2_SD, 
                               metrics_df_CIRI2_CE2_FC_2_Ind, metrics_df_CIRI2_CE2_FC_2_SD, 
                               metrics_df_CIRI2_CE2_FC_3_Ind, metrics_df_CIRI2_CE2_FC_3_SD)

xlsx::write.xlsx2(Evaluation_metrics_df, paste0(circRNA_predict_dir, "Evaluation_metrics_pre_CIRIquant_Overlap.xlsx"))

#######################