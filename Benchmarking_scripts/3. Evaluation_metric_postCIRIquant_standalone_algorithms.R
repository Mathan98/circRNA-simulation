### Calculate the evaluation metrics from output of the standalone tools


## Load libraries
library(dplyr)
library(rtracklayer)

#############################################################################################################

## Output folders in the main directory
### Folders
# CIRIquant_CIRI2_Individual_NoStringent  -> CIRIquant_CIRI2_Individual_No_Stringent (Folder)
# CIRIquant_CIRI2_Individual_Stringent    -> CIRIquant_CIRI2_Individual_Stringent
# CIRIquant_CIRCexplorer2_Individual      -> CIRIquant_CIRCexplorer2_Individual
# CIRIquant_find_circ_Ind                 -> CIRIquant_find_circ_Individual

# CIRIquant_CIRI2_Superdepth_noStringent  -> CIRIquant_CIRI2_Superdepth_No_Stringent
# CIRIquant_CIRI2_Superdepth_Stringent    -> CIRIquant_CIRI2_Superdepth_Stringent
# CIRIquant_CIRCexplorer2_Superdepth      -> CIRIquant_CIRCexplorer2_Superdepth
# CIRIquant_Find_circ_Superdepth          -> CIRIquant_find_circ_Superdepth


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

##################################################################################################################333


###############
## CIRIquant_CIRI2_Individual_NoStringent

metrics_CIRI2_Ind_No_Stringent <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_CIRI2_Individual_No_Stringent'), pattern = "gtf")
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
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_CIRI2_Individual_No_Stringent/Sim_', i, '.gtf')))
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
  metrics_CIRI2_Ind_No_Stringent[[paste0('CIRI2_Ind No Stringent_', i)]] <- calculate_metrics(confusion_table)
  
  
}

# Combine the results into a single data frame
metrics_df_CIRI2_Ind_No_Stringent <- do.call(rbind, metrics_CIRI2_Ind_No_Stringent)
metrics_df_CIRI2_Ind_No_Stringent[is.na(metrics_df_CIRI2_Ind_No_Stringent)] <- 0


################

## CIRIquant_CIRI2_Individual_Stringent

metrics_CIRI2_Ind_Stringent <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_CIRI2_Individual_Stringent'), pattern = "gtf")
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
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_CIRI2_Individual_Stringent/Sim_', i, '.gtf')))
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
  metrics_CIRI2_Ind_Stringent[[paste0('CIRI2_Ind Stringent_', i)]] <- calculate_metrics(confusion_table)
  
  
}

# Combine the results into a single data frame
metrics_df_CIRI2_Ind_Stringent <- do.call(rbind, metrics_CIRI2_Ind_Stringent)
metrics_df_CIRI2_Ind_Stringent[is.na(metrics_df_CIRI2_Ind_Stringent)] <- 0


##########################

## CIRIquant_CIRI2_Superdepth_noStringent

metrics_CIRI2_SD_No_Stringent <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_CIRI2_Superdepth_No_Stringent'), pattern = "gtf")
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
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_CIRI2_Superdepth_No_Stringent/Sim_', i, '.gtf')))
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
  metrics_CIRI2_SD_No_Stringent[[paste0('CIRI2_SD No Stringent_', i)]] <- calculate_metrics(confusion_table)
  
  
}

# Combine the results into a single data frame
metrics_df_CIRI2_SD_No_Stringent <- do.call(rbind, metrics_CIRI2_SD_No_Stringent)
metrics_df_CIRI2_SD_No_Stringent[is.na(metrics_df_CIRI2_SD_No_Stringent)] <- 0


##########################

## CIRIquant_CIRI2_Superdepth_Stringent

metrics_CIRI2_SD_Stringent <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_CIRI2_Superdepth_Stringent'), pattern = "gtf")
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
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_CIRI2_Superdepth_Stringent/Sim_', i, '.gtf')))
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
  metrics_CIRI2_SD_Stringent[[paste0('CIRI2_SD Stringent_', i)]] <- calculate_metrics(confusion_table)
  
}

# Combine the results into a single data frame
metrics_df_CIRI2_SD_Stringent <- do.call(rbind, metrics_CIRI2_SD_Stringent)
metrics_df_CIRI2_SD_Stringent[is.na(metrics_df_CIRI2_SD_Stringent)] <- 0


##########################

## CIRIquant_CIRCexplorer2_Superdepth

metrics_CE2_SD <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_CIRCexplorer2_Superdepth'), pattern = "gtf")
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
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_CIRCexplorer2_Superdepth/Sim_', i, '.gtf')))
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
  metrics_CE2_SD[[paste0('CE2_SD_', i)]] <- calculate_metrics(confusion_table)
  
  
}

# Combine the results into a single data frame
metrics_df_CE2_SD <- do.call(rbind, metrics_CE2_SD)
metrics_df_CE2_SD[is.na(metrics_df_CE2_SD)] <- 0

##########################

## CIRIquant_CIRCexplorer2_Individual

metrics_CE2_Individual <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_CIRCexplorer2_Individual'), pattern = "gtf")
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
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_CIRCexplorer2_Individual/Sim_', i, '.gtf')))
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
  metrics_CE2_Individual[[paste0('CE2_Ind_', i)]] <- calculate_metrics(confusion_table)
  
  
}

# Combine the results into a single data frame
metrics_df_CE2_Individual <- do.call(rbind, metrics_CE2_Individual)
metrics_df_CE2_Individual[is.na(metrics_df_CE2_Individual)] <- 0

##########################

## CIRIquant_Find_circ_Superdepth

metrics_FC_Superdepth <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_find_circ_Superdepth'), pattern = "gtf")
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
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_find_circ_Superdepth/Sim_', i, '.gtf')))
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
  metrics_FC_Superdepth[[paste0('FC_SD_', i)]] <- calculate_metrics(confusion_table)
  
  
}

# Combine the results into a single data frame
metrics_df_FC_Superdepth <- do.call(rbind, metrics_FC_Superdepth)
metrics_df_FC_Superdepth[is.na(metrics_df_FC_Superdepth)] <- 0

##########################

## CIRIquant_find_circ_Individual

metrics_FC_Ind <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRIquant_find_circ_Individual'), pattern = "gtf")
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
    predicted_positive <- as.data.frame(rtracklayer::import(paste0(circRNA_predict_dir, 'CIRIquant_find_circ_Individual/Sim_', i, '.gtf')))
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
  metrics_FC_Ind[[paste0('FC_Ind_', i)]] <- calculate_metrics(confusion_table)
  
  
}

# Combine the results into a single data frame
metrics_df_FC_Ind <- do.call(rbind, metrics_FC_Ind)
metrics_df_FC_Ind[is.na(metrics_df_FC_Ind)] <- 0




###########################################################################################

### Summarize all evaluation metrics

Evaluation_metrics_df <- rbind(metrics_df_CIRI2_Ind_No_Stringent, metrics_df_CIRI2_Ind_Stringent,
                               metrics_df_CIRI2_SD_No_Stringent, metrics_df_CIRI2_SD_Stringent, 
                               metrics_df_CE2_Individual, metrics_df_CE2_SD, 
                               metrics_df_FC_Ind,metrics_df_FC_Superdepth)

xlsx::write.xlsx2(Evaluation_metrics_df, paste0(circRNA_predict_dir, "Evaluation_metrics_CIRIquant.xlsx"))

###########################################################################################
