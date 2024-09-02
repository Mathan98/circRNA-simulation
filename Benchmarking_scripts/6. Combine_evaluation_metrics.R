
library(xlsx)
library(readxl)


## Path to directories containing the folders above and the true positive (tsv files)
# E.g. coverage circ 2x lin 100x

circRNA_sim_dir <- 'C:/Users/matha/Documents/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage2_100'
circRNA_predict_dir <- 'C:/Users/matha/Documents/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage2_100/'

################################################################################################################################

## Pre-CIRIquant
# Combine standalone and overlapping evaluation metrics excel file

## Evaluation metrics combined
Evaluation_metrics_Standalone_df <- read_excel(paste0(circRNA_predict_dir, "/Evaluation_metrics_pre_CIRIquant_allSamples.xlsx"))
Evaluation_metrics_Overlaps_df <- read_excel(paste0(circRNA_predict_dir, "Evaluation_metrics_pre_CIRIquant_Overlap.xlsx"))

Evaluation_metrics_df <- rbind(Evaluation_metrics_Standalone_df, Evaluation_metrics_Overlaps_df)
colnames(Evaluation_metrics_df)[1] <- "Methods"

xlsx::write.xlsx2(Evaluation_metrics_df %>% as.data.frame(), paste0(circRNA_predict_dir,"/Evaluation_metrics_preCIRIquant_combine.xlsx"), row.names = F)



################################################################################################################################

## Post-CIRIquant
# Combine standalone and overlapping evaluation metrics excel file

Evaluation_metrics_Standalone_df <- read_excel(paste0(circRNA_predict_dir, "/Evaluation_metrics_CIRIquant.xlsx"))
Evaluation_metrics_Overlaps_df <- read_excel(paste0(circRNA_predict_dir, "Evaluation_metrics_CIRIquant_Overlaps.xlsx"))

Evaluation_metrics_df <- rbind(Evaluation_metrics_Standalone_df, Evaluation_metrics_Overlaps_df)
colnames(Evaluation_metrics_df)[1] <- "Methods"

xlsx::write.xlsx2(Evaluation_metrics_df%>% as.data.frame(), paste0(circRNA_predict_dir,"/Evaluation_metrics_CIRIquant_combine.xlsx"), row.names = F)