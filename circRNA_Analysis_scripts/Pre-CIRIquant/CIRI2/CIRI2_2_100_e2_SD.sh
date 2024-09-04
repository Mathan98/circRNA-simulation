#! /bin/bash

#SBATCH --job-name=CIRI2_2_100_SD
#SBATCH --output=CIRI2_2_100_SD_output.log
#SBATCH --error=CIRI2_2_100_SD_error.log
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --ntasks=12
#SBATCH --hint=multithread
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=17170725@siswa.um.edu.my

module load miniconda
source activate CIRI2_New

mkdir /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_superdepth_noStringency
mkdir /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_superdepth_Stringent

/home/user/mk_98/.conda/envs/CIRI/bin/bwa mem -t 10 -T 19 /scr/user/mk_98/Simulated_reads/Homo_sapiens.fa /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/Sim_All_1.fq /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/Sim_All_2.fq 1> /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_superdepth_noStringency/Sim_All_unmapped.sam 2> /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_superdepth_noStringency/Sim_All_unmapped.log

## CIRIquant_noStringency
perl /home/user/mk_98/CIRIquant/libs/CIRI2.pl -t 10 -I /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_superdepth_noStringency/Sim_All_unmapped.sam -O /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_superdepth_noStringency/Sim_All.ciri -F /scr/user/mk_98/Simulated_reads/Homo_sapiens.fa -A /scr/user/mk_98/Simulated_reads/human_gencode_vch38.gtf -0 -G /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_superdepth_noStringency/CIRI2_All.log

## CIRIquant_Stringent
perl /home/user/mk_98/CIRIquant/libs/CIRI2.pl -t 10 -I /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_superdepth_noStringency/Sim_All_unmapped.sam -O /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_superdepth_Stringent/Sim_All.ciri -F /scr/user/mk_98/Simulated_reads/Homo_sapiens.fa -A /scr/user/mk_98/Simulated_reads/human_gencode_vch38.gtf -G /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_superdepth_Stringent/CIRI2_All.log
