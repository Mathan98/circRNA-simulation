#! /bin/bash

#SBATCH --job-name=FC_5_40_SD
#SBATCH --output=FC_5_40_SD_output.log
#SBATCH --error=FC_5_40_SD_error.log
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --ntasks=12
#SBATCH --hint=multithread
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=17170725@siswa.um.edu.my


module load miniconda
source activate CIRI2_New

## Individual
mkdir /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRIquant_find_circ_Superdepth

for i in {001..016}
do
CIRIquant --config '/scr/user/mk_98/Simulated_reads/GRCh38.yml' -1 /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/Sim_"$i"_1.fq -2 /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/Sim_"$i"_2.fq -t 10 -p Sim_"$i" -o /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRIquant_find_circ_Superdepth/Sim_"$i" --circ /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/find_circ_Superdepth/Sim_All_circ_candidates_bsj1.bed --tool find_circ --bsj-file /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRIquant_find_circ_Superdepth/Sim_BSJ_"$i".txt
done;
