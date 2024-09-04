### Calculate the evaluation metrics from output of the standalone tools


## Load libraries
library(dplyr)
library(rtracklayer)

#############################################################################################################

## Output folders in the main directory
### Folders
# CIRIquant CIRI2-CIRCexplorer2 Individual              -> CIRIquant_Overlap_CIRI2_CE2_Individual
# CIRIquant CIRI2-find_circ Individual                  -> CIRIquant_Overlap_CIRI2_FC_Individual
# CIRIquant find_circ-CIRCexplorer2 Individual          -> CIRIquant_Overlap_FC_CE2_Individual
# CIRIquant CIRI2-CIRCexplorer2-find_circ_2 Individual  -> CIRIquant_Overlap_CIRI2_CE2_FC_2_Individual
# CIRIquant CIRI2-CIRCexplorer2-find_circ_3 Individual  -> CIRIquant_Overlap_CIRI2_CE2_FC_3_Individual

# CIRIquant CIRI2-CIRCexplorer2 Superdepth              -> CIRIquant_Overlap_CIRI2_CE2_Superdepth
# CIRIquant CIRI2-find_circ Superdepth                  -> CIRIquant_Overlap_CIRI2_FC_Superdepth
# CIRIquant find_circ-CIRCexplorer2 Superdepth          -> CIRIquant_Overlap_FC_CE2_Superdepth
# CIRIquant CIRI2-CIRCexplorer2-find_circ_2 Superdepth  -> CIRIquant_Overlap_CIRI2_CE2_FC_2_Superdepth
# CIRIquant CIRI2-CIRCexplorer2-find_circ_3 Superdepth  -> CIRIquant_Overlap_CIRI2_CE2_FC_3_Superdepth

#############################################################################################################

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

##########################################################################################################

## CIRIquant CIRI2-CIRCexplorer2 Individual

metrics_CIRI2_CE2_Individual <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_Overlap_CIRI2_CE2_Individual'), pattern = "gtf")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  # Known true circRNA positive ids
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                        sep = '\t') 
  
  if (length(colnames(circRNAs)) == 1) {
    circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                          sep = ',')
  }
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    # Predicted positive ids
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_Overlap_CIRI2_CE2_Individual/Sim_', i, '.gtf')))
    predicted_positive$circ_id <- gsub("\\|", "-", predicted_positive$circ_id)
    
    
    # Full set of possible IDs
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)
    
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- nrow(circRNAs)} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_CE2_Individual[[paste0('CIRI2|CE2_Ind_', i)]] <- calculate_metrics(confusion_table)
  
  
}


# Combine the results into a single data frame
metrics_df_CIRI2_CE2_Individual <- do.call(rbind, metrics_CIRI2_CE2_Individual)
metrics_df_CIRI2_CE2_Individual[is.na(metrics_df_CIRI2_CE2_Individual)] <- 0



##########################

## CIRIquant CIRI2-CIRCexplorer2 Superdepth

metrics_CIRI2_CE2_Superdepth <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_Overlap_CIRI2_CE2_Superdepth'), pattern = "gtf")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  # Known true circRNA positive ids
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                        sep = '\t') 
  
  if (length(colnames(circRNAs)) == 1) {
    circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                          sep = ',')
  }
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    # Predicted positive ids
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_Overlap_CIRI2_CE2_Superdepth/Sim_', i, '.gtf')))
    predicted_positive$circ_id <- gsub("\\|", "-", predicted_positive$circ_id)
    
    
    # Full set of possible IDs
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)
    
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- nrow(circRNAs)} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_CE2_Superdepth[[paste0('CIRI2|CE2_SD_', i)]] <- calculate_metrics(confusion_table)
  
  
}


# Combine the results into a single data frame
metrics_df_CIRI2_CE2_Superdepth <- do.call(rbind, metrics_CIRI2_CE2_Superdepth)
metrics_df_CIRI2_CE2_Superdepth[is.na(metrics_df_CIRI2_CE2_Superdepth)] <- 0

##########################

## CIRIquant CIRI2-find_circ Individual

metrics_CIRI2_FC_Individual <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_Overlap_CIRI2_FC_Individual'), pattern = "gtf")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  # Known true circRNA positive ids
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                        sep = '\t') 
  
  if (length(colnames(circRNAs)) == 1) {
    circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                          sep = ',')
  }
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    # Predicted positive ids
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_Overlap_CIRI2_FC_Individual/Sim_', i, '.gtf')))
    predicted_positive$circ_id <- gsub("\\|", "-", predicted_positive$circ_id)
    
    
    # Full set of possible IDs
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)
    
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- nrow(circRNAs)} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_FC_Individual[[paste0('CIRI2|FC_Ind_', i)]] <- calculate_metrics(confusion_table)
  
  
}


# Combine the results into a single data frame
metrics_df_CIRI2_FC_Individual <- do.call(rbind, metrics_CIRI2_FC_Individual)
metrics_df_CIRI2_FC_Individual[is.na(metrics_df_CIRI2_FC_Individual)] <- 0

##########################

## CIRIquant CIRI2-find_circ Superdepth

metrics_CIRI2_FC_Superdepth <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_Overlap_CIRI2_FC_Superdepth'), pattern = "gtf")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  # Known true circRNA positive ids
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                        sep = '\t') 
  
  if (length(colnames(circRNAs)) == 1) {
    circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                          sep = ',')
  }
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    # Predicted positive ids
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_Overlap_CIRI2_FC_Superdepth/Sim_', i, '.gtf')))
    predicted_positive$circ_id <- gsub("\\|", "-", predicted_positive$circ_id)
    
    
    # Full set of possible IDs
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)
    
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- nrow(circRNAs)} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_FC_Superdepth[[paste0('CIRI2|FC_SD_', i)]] <- calculate_metrics(confusion_table)
  
  
}


# Combine the results into a single data frame
metrics_df_CIRI2_FC_Superdepth <- do.call(rbind, metrics_CIRI2_FC_Superdepth)
metrics_df_CIRI2_FC_Superdepth[is.na(metrics_df_CIRI2_FC_Superdepth)] <- 0

##########################

## CIRIquant find_circ-CIRCexplorer2 Individual

metrics_FC_CE2_Individual <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_Overlap_FC_CE2_Individual'), pattern = "gtf")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  # Known true circRNA positive ids
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                        sep = '\t') 
  
  if (length(colnames(circRNAs)) == 1) {
    circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                          sep = ',')
  }
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    # Predicted positive ids
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_Overlap_FC_CE2_Individual/Sim_', i, '.gtf')))
    predicted_positive$circ_id <- gsub("\\|", "-", predicted_positive$circ_id)
    
    
    # Full set of possible IDs
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- nrow(circRNAs)} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_FC_CE2_Individual[[paste0('FC|CE2_Ind_', i)]] <- calculate_metrics(confusion_table)
  
  
}


# Combine the results into a single data frame
metrics_df_FC_CE2_Individual <- do.call(rbind, metrics_FC_CE2_Individual)
metrics_df_FC_CE2_Individual[is.na(metrics_df_FC_CE2_Individual)] <- 0

##########################

## CIRIquant find_circ-CIRCexplorer2 Superdepth

metrics_FC_CE2_Superdepth <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_Overlap_FC_CE2_Superdepth'), pattern = "gtf")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  # Known true circRNA positive ids
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                        sep = '\t') 
  
  if (length(colnames(circRNAs)) == 1) {
    circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                          sep = ',')
  }
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    # Predicted positive ids
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_Overlap_FC_CE2_Superdepth/Sim_', i, '.gtf')))
    predicted_positive$circ_id <- gsub("\\|", "-", predicted_positive$circ_id)
    
    
    # Full set of possible IDs
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- nrow(circRNAs)} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_FC_CE2_Superdepth[[paste0('FC|CE2_SD_', i)]] <- calculate_metrics(confusion_table)
  
  
}


# Combine the results into a single data frame
metrics_df_FC_CE2_Superdepth <- do.call(rbind, metrics_FC_CE2_Superdepth)
metrics_df_FC_CE2_Superdepth[is.na(metrics_df_FC_CE2_Superdepth)] <- 0




##########################

## CIRIquant CIRI2-CIRCexplorer2-find_circ_2 Individual

metrics_CIRI2_CE2_FC_2_Individual <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_Overlap_CIRI2_CE2_FC_2_Individual'), pattern = "gtf")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  # Known true circRNA positive ids
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                        sep = '\t') 
  
  if (length(colnames(circRNAs)) == 1) {
    circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                          sep = ',')
  }
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    # Predicted positive ids
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_Overlap_CIRI2_CE2_FC_2_Individual/Sim_', i, '.gtf')))
    predicted_positive$circ_id <- gsub("\\|", "-", predicted_positive$circ_id)
    
    
    # Full set of possible IDs
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- nrow(circRNAs)} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_CE2_FC_2_Individual[[paste0('CIRI2|CE2|FC|2_Ind_', i)]] <- calculate_metrics(confusion_table)
  
  
}


# Combine the results into a single data frame
metrics_df_CIRI2_CE2_FC_2_Individual <- do.call(rbind, metrics_CIRI2_CE2_FC_2_Individual)
metrics_df_CIRI2_CE2_FC_2_Individual[is.na(metrics_df_CIRI2_CE2_FC_2_Individual)] <- 0



##########################

## CIRIquant CIRI2-CIRCexplorer2-find_circ_2 Superdepth

metrics_CIRI2_CE2_FC_2_Superdepth <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_Overlap_CIRI2_CE2_FC_2_Superdepth'), pattern = "gtf")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  # Known true circRNA positive ids
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                        sep = '\t') 
  
  if (length(colnames(circRNAs)) == 1) {
    circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                          sep = ',')
  }
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    # Predicted positive ids
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_Overlap_CIRI2_CE2_FC_2_Superdepth/Sim_', i, '.gtf')))
    predicted_positive$circ_id <- gsub("\\|", "-", predicted_positive$circ_id)
    
    
    # Full set of possible IDs
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- nrow(circRNAs)} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_CE2_FC_2_Superdepth[[paste0('CIRI2|CE2|FC|2_SD_', i)]] <- calculate_metrics(confusion_table)
  
  
}


# Combine the results into a single data frame
metrics_df_CIRI2_CE2_FC_2_Superdepth <- do.call(rbind, metrics_CIRI2_CE2_FC_2_Superdepth)
metrics_df_CIRI2_CE2_FC_2_Superdepth[is.na(metrics_df_CIRI2_CE2_FC_2_Superdepth)] <- 0




##########################

## CIRIquant_CIRI2_CE2_FC_3_Individual

metrics_CIRI2_CE2_FC_3_Individual <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_Overlap_CIRI2_CE2_FC_3_Individual'), pattern = "gtf")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  # Known true circRNA positive ids
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                        sep = '\t') 
  
  if (length(colnames(circRNAs)) == 1) {
    circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                          sep = ',')
  }
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    # Predicted positive ids
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_Overlap_CIRI2_CE2_FC_3_Individual/Sim_', i, '.gtf')))
    predicted_positive$circ_id <- gsub("\\|", "-", predicted_positive$circ_id)
    
    
    # Full set of possible IDs
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- nrow(circRNAs)} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_CE2_FC_3_Individual[[paste0('CIRI2|CE2|FC|3_Ind_', i)]] <- calculate_metrics(confusion_table)
  
  
}


# Combine the results into a single data frame
metrics_df_CIRI2_CE2_FC_3_Individual <- do.call(rbind, metrics_CIRI2_CE2_FC_3_Individual)
metrics_df_CIRI2_CE2_FC_3_Individual[is.na(metrics_df_CIRI2_CE2_FC_3_Individual)] <- 0






##########################

## CIRIquant CIRI2-CIRCexplorer2-find_circ_3 Superdepth

metrics_CIRI2_CE2_FC_3_Superdepth <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_Overlap_CIRI2_CE2_FC_3_Superdepth'), pattern = "gtf")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  
  # Known true circRNA positive ids
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                        sep = '\t') 
  
  if (length(colnames(circRNAs)) == 1) {
    circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, 
                          sep = ',')
  }
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    # Predicted positive ids
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_Overlap_CIRI2_CE2_FC_3_Superdepth/Sim_', i, '.gtf')))
    predicted_positive$circ_id <- gsub("\\|", "-", predicted_positive$circ_id)
    
    
    # Full set of possible IDs
    all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- nrow(circRNAs)} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CIRI2_CE2_FC_3_Superdepth[[paste0('CIRI2|CE2|FC|3_SD_', i)]] <- calculate_metrics(confusion_table)
  
  
}


# Combine the results into a single data frame
metrics_df_CIRI2_CE2_FC_3_Superdepth <- do.call(rbind, metrics_CIRI2_CE2_FC_3_Superdepth)
metrics_df_CIRI2_CE2_FC_3_Superdepth[is.na(metrics_df_CIRI2_CE2_FC_3_Superdepth)] <- 0



###########################################################################################

### Summarize all evaluation metrics (All Overlaps )

Evaluation_metrics_df <- rbind(metrics_df_CIRI2_CE2_Individual, metrics_df_CIRI2_CE2_Superdepth,
                               metrics_df_CIRI2_FC_Individual, metrics_df_CIRI2_FC_Superdepth,
                               metrics_df_FC_CE2_Individual, metrics_df_FC_CE2_Superdepth,
                               metrics_df_CIRI2_CE2_FC_2_Individual, metrics_df_CIRI2_CE2_FC_2_Superdepth, 
                               metrics_df_CIRI2_CE2_FC_3_Individual, metrics_df_CIRI2_CE2_FC_3_Superdepth)

xlsx::write.xlsx2(Evaluation_metrics_df, paste0(circRNA_predict_dir,"/Evaluation_metrics_CIRIquant_Overlaps.xlsx"))

# #######################
