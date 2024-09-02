### Simulation Overlaps analysis
library(dplyr)
library(ggVennDiagram)


setwd('C:/Users/matha/Documents/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage2_100/')


#############################################################################################################
# Create folders

dir.create("Overlap_CIRI2_CE2")
dir.create("Overlap_CIRI2_FC")
dir.create("Overlap_FC_CE2")
dir.create("Overlap_CIRI2_CE2_FC_2")
dir.create("Overlap_CIRI2_CE2_FC_3")
dir.create("Overlap_CIRI2_No_Stringent")

#############################################################################################################

### Individual analysis

## Create a dictionary for the circRNAs detected in CIRI2, CIRCexplorer2, find_circ

for (j in sprintf("%03d", 001:016)) {
  
  ## CIRI2_No_Stringent
  CIRI2_read <- read.delim2(paste0("./CIRI2_Ind_No_Stringent/Sim_", j, ".ciri"), header = T)
  CIRI2_read_1 <- CIRI2_read %>% 
    mutate(CIRI2_Stringent = 
             paste0(CIRI2_read$chr, ":", CIRI2_read$circRNA_start-1, "-", CIRI2_read$circRNA_end, " ", CIRI2_read$strand))
  
  
  ## CE2
  if (is.na(file.size(paste0('./CIRCexplorer2_Individual/circularRNA_known_', j, '.txt')))) {
    CE2_read <- data.frame(CircExplorer2 = character())
  } else {
    CE2_read <- read.delim2(paste0('./CIRCexplorer2_Individual/circularRNA_known_', j, '.txt'), header = F) %>% 
      mutate(CircExplorer2 = paste0(V1, ":", V2, "-", V3, " ", V6))
  }
  
  
  ## FC
  FC_read <- read.delim2(paste0('./find_circ_Individual/Sim_', j, '_circ_candidates_bsj1.bed'), header = FALSE) %>% 
    mutate(Find_circ = paste0(V1, ":", V2, "-", V3, " ", V6))
  
  
  ## Create a data frame by merging all discovered circRNAs to do a table with 1s and 0s
  circRNA_list <- data.frame(circRNA = unique(c(CIRI2_read_1$CIRI2_Stringent, CE2_read$CircExplorer2, FC_read$Find_circ)))
  
  for (i in 1:nrow(circRNA_list)) {
    
    ## CIRI2
    if (circRNA_list$circRNA[i] %in% CIRI2_read_1$CIRI2_Stringent) {
      circRNA_list$CIRI2[i] <- 1
    } 
    else (circRNA_list$CIRI2[i] <- 0)
    
    ## CE2
    if (circRNA_list$circRNA[i] %in% CE2_read$CircExplorer2) {
      circRNA_list$CE2[i] <- 1
    } 
    else (circRNA_list$CE2[i] <- 0)
    
    ## FC
    if (circRNA_list$circRNA[i] %in% FC_read$Find_circ) {
      circRNA_list$FC[i] <- 1
    } 
    else (circRNA_list$FC[i] <- 0)
  }
  
  
#############################################################################################################  
  ## Overlap of Individual results between two tools only
  
  # CIRI2-CIRCexplorer2
  
  CIRI2_CE2 <- circRNA_list[(circRNA_list$CIRI2 + circRNA_list$CE2 == 2),] %>% select(1)
  CIRI2_CE2 <-  CIRI2_read_1[(CIRI2_read_1$CIRI2_Stringent %in% CIRI2_CE2$circRNA),]
  if (!nrow(CIRI2_CE2) == 0) {
    CIRI2_CE2 <- data.frame(chr = CIRI2_CE2$chr, start = CIRI2_CE2$circRNA_start, end = CIRI2_CE2$circRNA_end, 
                            circRNA = CIRI2_CE2$circRNA_ID, 
                            gene = ".", strand = CIRI2_CE2$strand)
  }
  
  # CIRI2-find_circ
  
  CIRI2_FC <- circRNA_list[(circRNA_list$CIRI2 + circRNA_list$FC == 2),] %>% select(1)
  CIRI2_FC <-  CIRI2_read_1[(CIRI2_read_1$CIRI2_Stringent %in% CIRI2_FC$circRNA),]
  if (!nrow(CIRI2_FC) == 0) {
    CIRI2_FC <- data.frame(chr = CIRI2_FC$chr, start = CIRI2_FC$circRNA_start, end = CIRI2_FC$circRNA_end, 
                           circRNA = CIRI2_FC$circRNA_ID, 
                           gene = ".", strand = CIRI2_FC$strand)
  }
  
  # find_circ-CIRCexplorer2
  
  FC_CE2 <- circRNA_list[(circRNA_list$FC + circRNA_list$CE2 == 2),] %>% select(1)
  FC_CE2 <-  FC_read[(FC_read$Find_circ %in% FC_CE2$circRNA),]
  if (!nrow(FC_CE2) == 0) {
    FC_CE2 <- data.frame(chr = FC_CE2$V1, start = as.numeric(FC_CE2$V2) + 1, end = FC_CE2$V3, 
                         circRNA = paste0(paste0(FC_CE2$V1, ":", as.numeric(FC_CE2$V2) + 1, "|", FC_CE2$V3)), 
                         gene = ".", strand = FC_CE2$V6)
  }
  
  
  
  ## Overlap of results using 3 tools
  
  # CIRI2-CIRCexplorer2-find_circ, where circRNAs are predicted by at least 2 tools
  
  CIRI2_CE2_FC_2 <- circRNA_list[(circRNA_list$CIRI2 + circRNA_list$CE2 + circRNA_list$FC >= 2),] %>% select(1)
  CIRI2_CE2_FC_2 <-  CIRI2_read_1[(CIRI2_read_1$CIRI2_Stringent %in% CIRI2_CE2_FC_2$circRNA),]
  if (!nrow(CIRI2_CE2_FC_2) == 0) {
    CIRI2_CE2_FC_2 <- data.frame(chr = CIRI2_CE2_FC_2$chr, start = CIRI2_CE2_FC_2$circRNA_start, end = CIRI2_CE2_FC_2$circRNA_end, 
                                 circRNA = CIRI2_CE2_FC_2$circRNA_ID, 
                                 gene = ".", strand = CIRI2_CE2_FC_2$strand)
  }
  
  
  # CIRI2-CIRCexplorer2-find_circ, where circRNAs are predicted by all 3 tools
  
  CIRI2_CE2_FC_3 <- circRNA_list[(circRNA_list$CIRI2 + circRNA_list$CE2 + circRNA_list$FC == 3),] %>% select(1)
  CIRI2_CE2_FC_3 <-  CIRI2_read_1[(CIRI2_read_1$CIRI2_Stringent %in% CIRI2_CE2_FC_3$circRNA),]
  if (!nrow(CIRI2_CE2_FC_3) == 0) {
    CIRI2_CE2_FC_3 <- data.frame(chr = CIRI2_CE2_FC_3$chr, start = CIRI2_CE2_FC_3$circRNA_start, end = CIRI2_CE2_FC_3$circRNA_end, 
                                 circRNA = CIRI2_CE2_FC_3$circRNA_ID, 
                                 gene = ".", strand = CIRI2_CE2_FC_3$strand)
  }
  
  #############################################################################################################  
  
  ### Write those dataframes into bed files
  
  ## Overlap2 BED file
  write.table(CIRI2_CE2, file = paste0("./Overlap_CIRI2_CE2/Sim_", j, ".bed"), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(CIRI2_FC, file = paste0("./Overlap_CIRI2_FC/Sim_", j, ".bed"), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(FC_CE2, file = paste0("./Overlap_FC_CE2/Sim_", j, ".bed"), col.names = F, row.names = F, quote = F, sep = "\t")
  
  ## Overlap3 BED file
  
  write.table(CIRI2_CE2_FC_2, file = paste0("./Overlap_CIRI2_CE2_FC_2/Sim_", j, ".bed"), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(CIRI2_CE2_FC_3, file = paste0("./Overlap_CIRI2_CE2_FC_3/Sim_", j, ".bed"), col.names = F, row.names = F, quote = F, sep = "\t")
}


## Move those folders to ( Overlap_CIRI2_No_Stringent )
system('bash -c "mv Overlap_CIRI2_CE2 Overlap_CIRI2_No_Stringent"')
system('bash -c "mv Overlap_CIRI2_FC Overlap_CIRI2_No_Stringent"')
system('bash -c "mv Overlap_FC_CE2 Overlap_CIRI2_No_Stringent"')
system('bash -c "mv Overlap_CIRI2_CE2_FC_2 Overlap_CIRI2_No_Stringent"')
system('bash -c "mv Overlap_CIRI2_CE2_FC_3 Overlap_CIRI2_No_Stringent"')


###########################################################################################################################################################################

## Make new folders

dir.create("Superdepth_Overlap_CIRI2_CE2")
dir.create("Superdepth_Overlap_CIRI2_FC")
dir.create("Superdepth_Overlap_FC_CE2")
dir.create("Superdepth_Overlap_CIRI2_CE2_FC_2")
dir.create("Superdepth_Overlap_CIRI2_CE2_FC_3")


#############################################################################################################

### Superdepth analysis

## CIRI2_No_Stringent
CIRI2_read <- read.delim2(paste0("./CIRI2_superdepth_noStringency/Sim_All.ciri"), header = T)
CIRI2_read_1 <- CIRI2_read %>% 
  mutate(CIRI2_Stringent = 
           paste0(CIRI2_read$chr, ":", CIRI2_read$circRNA_start-1, "-", CIRI2_read$circRNA_end, " ", CIRI2_read$strand))


## CE2
if (file.size(paste0('./CIRCexplorer2_Superdepth/circularRNA_known_All.txt')) == 0) {
  CE2_read <- data.frame(CircExplorer2 = character())
} else {
  CE2_read <- read.delim2(paste0('./CIRCexplorer2_Superdepth/circularRNA_known_All.txt'), header = F) %>% 
    mutate(CircExplorer2 = paste0(V1, ":", V2, "-", V3, " ", V6))
}


## FC
FC_read <- read.delim2(paste0('./find_circ_Superdepth/Sim_All_circ_candidates_bsj1.bed'), header = FALSE) %>% 
  mutate(Find_circ = paste0(V1, ":", V2, "-", V3, " ", V6))



## Create a data frame by merging all discovered circRNAs to do a table with 1s and 0s
circRNA_list <- data.frame(circRNA = unique(c(CIRI2_read_1$CIRI2_Stringent, CE2_read$CircExplorer2, FC_read$Find_circ)))

for (i in 1:nrow(circRNA_list)) {
  
  ## CIRI2
  if (circRNA_list$circRNA[i] %in% CIRI2_read_1$CIRI2_Stringent) {
    circRNA_list$CIRI2[i] <- 1
  } 
  else (circRNA_list$CIRI2[i] <- 0)
  
  ## CE2
  if (circRNA_list$circRNA[i] %in% CE2_read$CircExplorer2) {
    circRNA_list$CE2[i] <- 1
  } 
  else (circRNA_list$CE2[i] <- 0)
  
  ## FC
  if (circRNA_list$circRNA[i] %in% FC_read$Find_circ) {
    circRNA_list$FC[i] <- 1
  } 
  else (circRNA_list$FC[i] <- 0)
}

#############################################################################################################  
## Overlap of Superdepth results between two tools only

# CIRI2-CIRCexplorer2

CIRI2_CE2 <- circRNA_list[(circRNA_list$CIRI2 + circRNA_list$CE2 == 2),] %>% select(1)
CIRI2_CE2 <-  CIRI2_read_1[(CIRI2_read_1$CIRI2_Stringent %in% CIRI2_CE2$circRNA),]
if (!nrow(CIRI2_CE2) == 0) {
  CIRI2_CE2 <- data.frame(chr = CIRI2_CE2$chr, start = CIRI2_CE2$circRNA_start, end = CIRI2_CE2$circRNA_end, 
                          circRNA = CIRI2_CE2$circRNA_ID, 
                          gene = ".", strand = CIRI2_CE2$strand)
}


# CIRI2-find_circ

CIRI2_FC <- circRNA_list[(circRNA_list$CIRI2 + circRNA_list$FC == 2),] %>% select(1)
CIRI2_FC <-  CIRI2_read_1[(CIRI2_read_1$CIRI2_Stringent %in% CIRI2_FC$circRNA),]
if (!nrow(CIRI2_FC) == 0) {
  CIRI2_FC <- data.frame(chr = CIRI2_FC$chr, start = CIRI2_FC$circRNA_start, end = CIRI2_FC$circRNA_end, 
                         circRNA = CIRI2_FC$circRNA_ID, 
                         gene = ".", strand = CIRI2_FC$strand)
}


# find_circ-CIRCexplorer2

FC_CE2 <- circRNA_list[(circRNA_list$FC + circRNA_list$CE2 == 2),] %>% select(1)
FC_CE2 <-  FC_read[(FC_read$Find_circ %in% FC_CE2$circRNA),]
if (!nrow(FC_CE2) == 0) {
  FC_CE2 <- data.frame(chr = FC_CE2$V1, start = as.numeric(FC_CE2$V2) + 1, end = FC_CE2$V3, 
                       circRNA = paste0(paste0(FC_CE2$V1, ":", as.numeric(FC_CE2$V2) + 1, "|", FC_CE2$V3)), 
                       gene = ".", strand = FC_CE2$V6)
}


#######################
## Overlap of results using 3 tools

# CIRI2-CIRCexplorer2-find_circ, where circRNAs are predicted by at least 2 tools

CIRI2_CE2_FC_2 <- circRNA_list[(circRNA_list$CIRI2 + circRNA_list$CE2 + circRNA_list$FC >= 2),] %>% select(1)
CIRI2_CE2_FC_2 <-  CIRI2_read_1[(CIRI2_read_1$CIRI2_Stringent %in% CIRI2_CE2_FC_2$circRNA),]
if (!nrow(CIRI2_CE2_FC_2) == 0) {
  CIRI2_CE2_FC_2 <- data.frame(chr = CIRI2_CE2_FC_2$chr, start = CIRI2_CE2_FC_2$circRNA_start, end = CIRI2_CE2_FC_2$circRNA_end, 
                               circRNA = CIRI2_CE2_FC_2$circRNA_ID, 
                               gene = ".", strand = CIRI2_CE2_FC_2$strand)
}


# CIRI2-CIRCexplorer2-find_circ, where circRNAs are predicted by all 3 tools

CIRI2_CE2_FC_3 <- circRNA_list[(circRNA_list$CIRI2 + circRNA_list$CE2 + circRNA_list$FC == 3),] %>% select(1)
CIRI2_CE2_FC_3 <-  CIRI2_read_1[(CIRI2_read_1$CIRI2_Stringent %in% CIRI2_CE2_FC_3$circRNA),]
if (!nrow(CIRI2_CE2_FC_3) == 0) {
  CIRI2_CE2_FC_3 <- data.frame(chr = CIRI2_CE2_FC_3$chr, start = CIRI2_CE2_FC_3$circRNA_start, end = CIRI2_CE2_FC_3$circRNA_end, 
                               circRNA = CIRI2_CE2_FC_3$circRNA_ID, 
                               gene = ".", strand = CIRI2_CE2_FC_3$strand)
}


#############################################################################################################

## Overlap2 BED file
write.table(CIRI2_CE2, file = paste0("./Superdepth_Overlap_CIRI2_CE2/Sim_Superdepth.bed"), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(CIRI2_FC, file = paste0("./Superdepth_Overlap_CIRI2_FC/Sim_Superdepth.bed"), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(FC_CE2, file = paste0("./Superdepth_Overlap_FC_CE2/Sim_Superdepth.bed"), col.names = F, row.names = F, quote = F, sep = "\t")


## Overlap3 BED file

write.table(CIRI2_CE2_FC_2, file = paste0("./Superdepth_Overlap_CIRI2_CE2_FC_2/Sim_Superdepth.bed"), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(CIRI2_CE2_FC_3, file = paste0("./Superdepth_Overlap_CIRI2_CE2_FC_3/Sim_Superdepth.bed"), col.names = F, row.names = F, quote = F, sep = "\t")


###########################################################################################################################################################################

## Make new folders
system('bash -c "mv Superdepth_Overlap_CIRI2_CE2 Overlap_CIRI2_No_Stringent"')
system('bash -c "mv Superdepth_Overlap_CIRI2_FC Overlap_CIRI2_No_Stringent"')
system('bash -c "mv Superdepth_Overlap_FC_CE2 Overlap_CIRI2_No_Stringent"')
system('bash -c "mv Superdepth_Overlap_CIRI2_CE2_FC_2 Overlap_CIRI2_No_Stringent"')
system('bash -c "mv Superdepth_Overlap_CIRI2_CE2_FC_3 Overlap_CIRI2_No_Stringent"')

#############################################################################################################






