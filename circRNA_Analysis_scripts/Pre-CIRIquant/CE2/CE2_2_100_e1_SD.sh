#! /bin/bash

#SBATCH --job-name=CE2_5_40_SD
#SBATCH --output=CE2_5_40_SD_output.log
#SBATCH --error=CE2_5_40_SD_error.log
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --ntasks=12
#SBATCH --hint=multithread
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=17170725@siswa.um.edu.my

module load miniconda
source activate CIRCexplorer2_Latest

#bowtie-build /scr/user/mk_98/Simulated_reads/Homo_sapiens.fa /scr/user/mk_98/Simulated_reads/bowtie_index/Homo_sapiens
#bowtie2-build /scr/user/mk_98/Simulated_reads/Homo_sapiens.fa /scr/user/mk_98/Simulated_reads/bowtie2_index/Homo_sapiens

mkdir /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRCexplorer2_Superdepth


#tophat_align
tophat2 -o /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRCexplorer2_Superdepth/tophat_fusion_All -p 10 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search /scr/user/mk_98/Simulated_reads/bowtie_index/Homo_sapiens /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/Sim_All_1.fq /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/Sim_All_2.fq

CIRCexplorer2 parse --pe -t TopHat-Fusion /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRCexplorer2_Superdepth/tophat_fusion_All/accepted_hits.bam -b /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRCexplorer2_Superdepth/back_spliced_junction_All.bed > /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRCexplorer2_Superdepth/CIRCexplorer_All_parse.log

CIRCexplorer2 annotate -r /scr/user/mk_98/Simulated_reads/hg38_ref_all.txt -g /scr/user/mk_98/Simulated_reads/Homo_sapiens.fa -b /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRCexplorer2_Superdepth/back_spliced_junction_All.bed -o /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRCexplorer2_Superdepth/circularRNA_known_All.txt > /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRCexplorer2_Superdepth/CIRCexplorer2_All_annotate.log

