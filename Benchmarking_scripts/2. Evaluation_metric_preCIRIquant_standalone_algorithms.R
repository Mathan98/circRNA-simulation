
### Calculate the evaluation metrics from output of the standalone tools


library(dplyr)
library(rtracklayer)

## Output folders in the main directory
                                   ### Folders
# CIRI2 Individual No Stringent -> CIRI2_Ind_No_Stringent (Folder)
# CIRI2 Individual Stringent    -> CIRI2_Ind_Stringent
# CIRCexplorer2 Individual      -> CIRCexplorer2_Individual
# find_circ Individual          -> find_circ_Individual

# CIRI2 Superdepth No Stringent -> CIRI2_superdepth_noStringent
# CIRI2 Superdepth No Stringent -> CIRI2_superdepth_Stringent
# CIRCexplorer2 Superdepth      -> CIRCexplorer2_Superdepth
# find_circ Superdepth          -> find_circ_Superdepth


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
###############

##CIRI2 Individual No Stringent

metrics_CIRI2_Ind_No_Stringent <- list()

for (i in sprintf("%03d", 001:016)) {
  
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, sep = '\t')
  
  if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F,sep = ',')}
  
  ## Predicted
  predicted_positive <- read.csv2(paste0(circRNA_predict_dir, '/CIRI2_Ind_No_Stringent/Sim_', i, '.ciri'), header = T, sep = '\t')
  predicted_positive$circRNA_ID <- trimws(predicted_positive$circRNA_ID)
  predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$circRNA_ID)
  
  
  
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
  metrics_CIRI2_Ind_No_Stringent[[paste0('CIRI2_Ind No Stringent_', i)]] <- calculate_metrics(confusion_table)
}


# Combine the results into a single data frame
metrics_df_CIRI2_Ind_No_Stringent <- do.call(rbind, metrics_CIRI2_Ind_No_Stringent)

#####################################################################################################################################
###############

##CIRI2 Individual Stringent

# Predicted positive ids
metrics_CIRI2_Ind_Stringent <- list()

for (i in sprintf("%03d", 001:016)) {
  
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, sep = '\t')
  
  if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F,sep = ',')}
  
  ## Predicted
  predicted_positive <- read.csv2(paste0(circRNA_predict_dir, '/CIRI2_Ind_Stringent/Sim_', i, '.ciri'), header = T, sep = '\t')
  predicted_positive$circRNA_ID <- trimws(predicted_positive$circRNA_ID)
  predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$circRNA_ID)
  
  
  
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
  metrics_CIRI2_Ind_Stringent[[paste0('CIRI2_Ind Stringent_', i)]] <- calculate_metrics(confusion_table)
}


# Combine the results into a single data frame
metrics_df_CIRI2_Ind_Stringent <- do.call(rbind, metrics_CIRI2_Ind_Stringent)

#####################################################################################################################################
###############

## CIRCexplorer2 Individual

metrics_CE2_Individual <- list()

lf <- list.files(path=paste0(circRNA_sim_dir, '/CIRCexplorer2_Individual'), pattern = "txt")
lf <- gsub("[^0-9]", "", lf)

for (i in sprintf("%03d", 001:016)) {
  ## Sim circRNAs
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, sep = '\t')
  
  
  if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F,sep = ',')}
  
  if (!i %in% lf ){
    confusion_table = matrix(nrow = 2, ncol = 2)
    confusion_table[1,2] <- 0
    confusion_table[2,1] <- nrow(circRNAs)
    confusion_table[2,2] <- 0
    confusion_table[1,1] <- 0
    
  } else {
    
    ## predicted circRNAs
    predicted_positive <- read.delim2(paste0(circRNA_predict_dir, "CIRCexplorer2_Individual/circularRNA_known_", i, ".txt"), header = F)
    
    predicted_positive$circ_id <- paste0(predicted_positive$V1, ":", predicted_positive$V2+1, "-", predicted_positive$V3)
    
    
    # Full set of possible IDs for gene_ids and circRNA_ids
    all_ids_circRNA <- c(circRNAs$Circ_ID, predicted_positive$circ_id)
    
    
    true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
    true_negative <- integer()
    false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
    false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)
    
    
    confusion_table = matrix(nrow = 2, ncol = 2)
    if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
    if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
    if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
    if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}
    
  }
  
  # Loop through the list and calculate metrics for each confusion matrix
  metrics_CE2_Individual[[paste0('CE2_Ind_', i)]] <- calculate_metrics(confusion_table)
  
  
}


# Combine the results into a single data frame
metrics_df_CE2_Individual <- do.call(rbind, metrics_CE2_Individual)

#####################################################################################################################################
###############

## find_circ Individual


metrics_FC_Ind <- list()

for (i in sprintf("%03d", 001:016)) {
  
  circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F, sep = '\t')
  
  if (length(colnames(circRNAs)) == 1) {circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/Sim_', i, '_tp.csv'), header = F,sep = ',')}
  
  # Predicted positive ids
  predicted_positive <- read.delim2(paste0(circRNA_predict_dir, "find_circ_Individual/Sim_", i, "_circ_candidates_bsj1.bed"), header = F)
  predicted_positive$circRNA_ID <- paste0(predicted_positive$V1, ":", predicted_positive$V2+1, "-", predicted_positive$V3)
  
  
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
  metrics_FC_Ind[[paste0('FC_Ind_', i)]] <- calculate_metrics(confusion_table)
}

# Combine the results into a single data frame
metrics_df_FC_Ind <- do.call(rbind, metrics_FC_Ind)



#####################################################################################################################################
###############

### CIRI2 Superdepth No Stringent

circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F, sep = '\t')

# Predicted positive ids
predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/CIRI2_superdepth_noStringent/Sim_All.ciri'))
predicted_positive$circRNA_ID <- trimws(predicted_positive$circRNA_ID)
predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$circRNA_ID)


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
metrics_df_CIRI2_SD_No_Stringent <- calculate_metrics(confusion_table)
rownames(metrics_df_CIRI2_SD_No_Stringent) <- "CIRI2_SD No Stringent"



##############################################################################################################
###############

## CIRI2 Superdepth Stringent

circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F, sep = '\t')

# Predicted positive ids
predicted_positive <- read.delim2(paste0(circRNA_predict_dir, '/CIRI2_superdepth_Stringent/Sim_All.ciri'))
predicted_positive$circRNA_ID <- trimws(predicted_positive$circRNA_ID)
predicted_positive$circRNA_ID <- gsub("\\|", "-", predicted_positive$circRNA_ID)


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
metrics_df_CIRI2_SD_Stringent <- calculate_metrics(confusion_table)
rownames(metrics_df_CIRI2_SD_Stringent) <- "CIRI2_SD Stringent"

####################################################################################################################################
###############

## CIRCexplorer2 Superdepth

circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F, sep = '\t')

# Predicted positive ids
predicted_positive <- read.delim2(paste0(circRNA_predict_dir, "CIRCexplorer2_Superdepth/circularRNA_known_All.txt"), header = F)
predicted_positive$circ_id <- paste0(predicted_positive$V1, ":", predicted_positive$V2+1, "-", predicted_positive$V3)


# Full set of possible IDs for gene_ids and circRNA_ids
all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)


true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
true_negative <- integer()
false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)


confusion_table = matrix(nrow = 2, ncol = 2)
if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}


# Loop through the list and calculate metrics for each confusion matrix
metrics_df_CE2_SD <- calculate_metrics(confusion_table)
rownames(metrics_df_CE2_SD) <- "CE2_SD"

####################################################################################################################################
###############

## find_circ Superdepth

circRNAs <- read.csv2(paste0(circRNA_sim_dir, '/SimALL_tp.csv'), header = F, sep = '\t')

# Predicted positive ids
predicted_positive <- read.delim2(paste0(circRNA_predict_dir, "find_circ_Superdepth/Sim_All_circ_candidates_bsj1.bed"), header = F)
predicted_positive$circ_id <- paste0(predicted_positive$V1, ":", predicted_positive$V2+1, "-", predicted_positive$V3)


# Full set of possible IDs for gene_ids and circRNA_ids
all_ids_circRNA <- c(circRNAs$V1, predicted_positive$circ_id)


true_positive <- intersect(circRNAs$V1, predicted_positive$circ_id)
true_negative <- integer()
false_positive <- setdiff(predicted_positive$circ_id, circRNAs$V1)
false_negative <- setdiff(circRNAs$V1, predicted_positive$circ_id)


confusion_table = matrix(nrow = 2, ncol = 2)
if (length(false_positive) == 0) {confusion_table[1,2] <- 0} else {confusion_table[1,2] <- length(false_positive)}
if (length(false_negative) == 0) {confusion_table[2,1] <- 0} else {confusion_table[2,1] <- length(false_negative)}
if (length(true_positive) == 0) {confusion_table[2,2] <- 0} else {confusion_table[2,2] <- length(true_positive)}
if (length(true_negative) == 0) {confusion_table[1,1] <- 0} else {confusion_table[1,1] <- length(true_negative)}

# Loop through the list and calculate metrics for each confusion matrix
metrics_df_FC_SD <- calculate_metrics(confusion_table)
rownames(metrics_df_FC_SD) <- "FC_SD"


###########################################################################################

### Summarize all evaluation metrics 

Evaluation_metrics_df <- rbind(metrics_df_CIRI2_Ind_No_Stringent, metrics_df_CIRI2_Ind_Stringent, metrics_df_CE2_Individual, metrics_df_FC_Ind,
                               metrics_df_CIRI2_SD_No_Stringent, metrics_df_CIRI2_SD_Stringent, metrics_df_CE2_SD, metrics_df_FC_SD)

xlsx::write.xlsx2(Evaluation_metrics_df, paste0(circRNA_predict_dir, "Evaluation_metrics_pre_CIRIquant_allSamples.xlsx"))

#######################


