
## Install packages and dependencies

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("polyester", quietly = TRUE))
  BiocManager::install("polyester")

if (!require("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")

if (!require("GenomicFeatures", quietly = TRUE))
  BiocManager::install("GenomicFeatures")

if (!require("dplyr", quietly = TRUE))
  BiocManager::install("dplyr")

if (!require("ballgown", quietly = TRUE))
  BiocManager::install("ballgown")

if (!require("stringr", quietly = TRUE))
  BiocManager::install("stringr")

if (!require("glue", quietly = TRUE))
  BiocManager::install("glue")

if (!require("rtracklayer", quietly = TRUE))
  BiocManager::install("rtracklayer")


## Load packages

library(polyester)
library(Biostrings)
library(GenomicFeatures)
library(dplyr)
library(ballgown)
library(stringr)
library(glue)
library(rtracklayer)


# Read in the GTF file; in this case human gencode v25
human_gtf <- rtracklayer::import('human_gencode_vch38.gtf') 


# Retain only chr1 to chr 22, chrX and chrY
human_gtf <- human_gtf[seqnames(human_gtf) %in% c(paste0("chr", 1:22), "chrX", "chrY")]

# Get unique gene_id
human_gene_id <- unique(human_gtf$gene_id) %>% as.data.frame()


###############################################################################################################

dir.create('linearRNA_simulated_reads')

## if simulating linear RNA coverage of 100X with sequencing error 1; otherwise, 40X, 20X, 10X or 1,2,5
dir.create('linearRNA_simulated_reads/lin_reads_chrAll_coverage100_error1') ## if simulating linear RNA coverage of 100X


# Randomly select transcripts from gene_id
for (i in sprintf("%03d", 001:016)) {
  
  ## Set the range of numbers to sample from; in this case arbitrary
  start_num <- 600
  end_num <- 800
  
  num_transcripts <- sample(start_num:end_num, 1)
  
  selected_gene_ID <- list(human_gene_id[sample(nrow(human_gene_id), num_transcripts),])
  
  ## Subset the GTF to keep only the sampled transcripts
  regex <- glue_collapse(selected_gene_ID[[1]], "|")
  human_gtf_1 <- human_gtf[str_detect(human_gtf$gene_id, regex),]
  
  
  ## Write the subset gtf files
  export(human_gtf_1, paste0("linearRNA_simulated_reads/lin_reads_chrAll_coverage100_error1/Sim_", i, '.gtf'))
  
}


###############################################################################################################


# Simulate reads from polyester using the subsetted GTF

for (i in sprintf("%03d", 001:016)) {
  
  ## Scan for transcripts in GTF file
  n.trx <- sum("transcript" == scan(file = paste0('Sim_', i, '.gtf'),
                                    what = "character"),
               na.rm = T)
  
  ## LogFC table of 1 (no change in expression in this case)
  fold_changes = matrix(1,
                        ncol = 1,
                        nrow = n.trx)
  
  # ## create directory
  results.path <- paste0('C:\\linearRNA_simulated_reads\\lin_reads_chrAll_coverage100_error1\\SIM_', i)
  dir.create(path = results.path, recursive = T, showWarnings = T)
  
  
  ## simulate reads
  simulate_experiment(seqpath = 'ChrAll_fa',
                      gtf = paste0('linearRNA_simulated_reads/lin_reads_chrAll_coverage100_error1/Sim_', i, '.gtf'),
                      # reads_per_transcript = readspertx,
                      reads_per_transcript = 100, ## linear RNA coverage {100X => 100; 40X => 40; 20X => 20; 10X => 10}
                      num_reps = c(1, 1), #c(1),
                      readlen = 120,
                      fraglen = 350,
                      fragsd = 70,
                      fold_changes = fold_changes,
                      error_rate = 0.01, ## sequencing error {1% => 0.01; 2% => 0.02; 5% => 0.05}
                      # meanmodel = TRUE,
                      outdir = results.path,
                      strand_specific = F,
                      paired = T)
  
}


list_files <- list.files(path = "linearRNA_simulated_reads/lin_reads_chrAll_coverage100_error1", 
                         pattern = "^sample_02", recursive = T)
file.remove(list_files)
