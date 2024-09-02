### ANOVA for simulated Dataset

library(rtracklayer)
library(dplyr)
library(tidyr)
library(rstatix)
library(readxl)
library(reshape2) 
library(ggpubr)
library(xlsx)
library(multcompView)
library(patchwork)

## Path to directories containing the folders above and the true positive (tsv files)
# E.g. coverage circ 2x lin 100x

circRNA_sim_dir <- 'C:/Users/matha/Documents/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage2_100'
circRNA_predict_dir <- 'C:/Users/matha/Documents/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage2_100/'

##################################################################################################
## Create df to store data

Anova_table_Recall <- as.data.frame(matrix(NA, 0, 3))
colnames(Anova_table_Recall) <- c("Method", "Dataset", "Recall")


##################################################################################################
## Pre CIRIquant
Evaluation_metrics <- read_excel(paste0(circRNA_predict_dir, '/Evaluation_metrics_preCIRIquant_combine.xlsx'))


## Post CIRIquant
# Evaluation_metrics <- read.xlsx(paste0(circRNA_predict_dir, 'Evaluation_metrics_CIRIquant_combine.xlsx'), sheetIndex = 1)

##################################################################################################

colnames(Evaluation_metrics)[1] <- "Methods_ID"

# Split the ID column based on "_"
split_column <- strsplit(Evaluation_metrics$Methods_ID, "_")

# Create new columns based on the split results
Evaluation_metrics$Methods <- sapply(split_column, function(x) x[1])
Evaluation_metrics$Type <- sapply(split_column, function(x) x[2])
Evaluation_metrics$ID <- sapply(split_column, function(x) x[3])


## Create dataframes for Recall, Recall and F1 score

# Recall
Anova_table_Recall[Evaluation_metrics$Methods_ID,1] <- paste0(Evaluation_metrics$Methods, "_", Evaluation_metrics$Type)
Anova_table_Recall[Evaluation_metrics$Methods_ID,2] <- Evaluation_metrics$ID
Anova_table_Recall[Evaluation_metrics$Methods_ID,3] <- Evaluation_metrics$Recall


## Rename columns and standardise strategy acronyms
colnames(Anova_table_Recall)[1] <- "Strategies"

Anova_table_Recall$Strategies <- gsub("\\_s*\\d+$", "", Anova_table_Recall$Strategies)
Anova_table_Recall$Strategies <- gsub("_", " ", Anova_table_Recall$Strategies)
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2 Ind Stringent"] <- "C2_Str_IND"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2 Ind No Stringent"]<- "C2_NStr_IND"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2 SD No Stringent"]<- "C2_NStr_SD"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2 SD Stringent"]<- "C2_Str_SD"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CE2 Ind"]<- "CE2_IND"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CE2 SD"]<- "CE2_SD"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "FC Ind"]<- "FC_IND"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "FC SD"]<- "FC_SD"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2|CE2 Ind"]<- "C2CE2_IND"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2|CE2 SD"]<- "C2CE2_SD"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2|FC Ind"]<- "C2FC_IND"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2|FC SD"]<- "C2FC_SD"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "FC|CE2 Ind"]<- "CE2FC_IND"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "FC|CE2 SD"]<- "CE2FC_SD"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2|CE2|FC|2 Ind"]<- "C2CE2FC_2_IND"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2|CE2|FC|2 SD"]<- "C2CE2FC_2_SD"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2|CE2|FC|3 Ind"]<- "C2CE2FC_3_IND"
Anova_table_Recall$Strategies[Anova_table_Recall$Strategies == "CIRI2|CE2|FC|3 SD"]<- "C2CE2FC_3_SD"

## Factorise the strategies
Strategies_list <-c("C2_NStr_IND", "C2_NStr_SD", "C2_Str_IND", "C2_Str_SD", "CE2_IND", "CE2_SD",
                    "FC_IND", "FC_SD", "C2CE2_IND", "C2CE2_SD", "C2FC_IND", "C2FC_SD", "CE2FC_IND", "CE2FC_SD",
                    "C2CE2FC_2_IND", "C2CE2FC_2_SD", "C2CE2FC_3_IND", "C2CE2FC_3_SD")


##################################################################################################
#### Statistical tests
Anova_table_Recall$Strategies <- factor(Anova_table_Recall$Strategies, levels = Strategies_list)
Anova_table_Recall$Recall <- as.numeric(Anova_table_Recall$Recall)
res.aov <- Anova_table_Recall %>% anova_test(Recall ~ Strategies )
model_1 <- aov(Recall ~ Strategies, data=Anova_table_Recall)
Tukey= TukeyHSD(model_1)
Cld=multcompLetters4(model_1, Tukey)

#add letter to table
Cld= as.data.frame.list(Cld$Strategies)


data_summ = group_by(Anova_table_Recall, Strategies) %>% 
  summarise(mean=mean(Recall), sd=sd(Recall)) %>% 
  arrange(desc(mean))
data_summ$Cld=toupper(Cld$Letters)

pwc <- Anova_table_Recall %>% tukey_hsd(Recall ~ Strategies)


##################################################################################################
## Plot the ANOVA boxplot

ggboxplot(Anova_table_Recall, x = "Strategies", y = "Recall", color = "Strategies", bxp.errorbar =T, 
          size=0.4, width = 0.5,
          label.rectangle = T, fill = "Strategies" ,  colour='black') +
  #stat_pvalue_manual(pwc, hide.ns = T, y.position = 560, step.increase = 0.1, label = "{p.adj}",
  #fontface = "bold", color = "dark red") +
  labs(#subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc),
    #    title = "Sequencing Error Rate: 5 %"
  ) +
  xlab("Analysis") +
  ylab("Recall Score") +
  geom_text(data = data_summ, aes(label=Cld,x=Strategies, y=mean +sd), 
            position=position_dodge2(0.75), vjust = -2, size = 8, color = "red",
            fontface = 'bold') +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 40),
        plot.subtitle = element_text(size = 12, face = "bold.italic", colour="black", hjust=0.5),
        plot.caption = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 17, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 30, colour = "red", hjust = 0.5, vjust = 1.5),
        legend.title = element_text(size = 20, face = "bold", hjust=0.5, colour = "red"),
        legend.text = element_text(size = 15, face = "bold"),
        legend.position = "right",
        legend.box.background = element_rect(colour = "black"),
        legend.box.margin = margin(t = 1, l = 1),
        legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm')) + 
  ylim(0,1)




################################################################################
