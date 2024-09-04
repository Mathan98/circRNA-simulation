#! /bin/bash

#SBATCH --job-name=GPU1__CIRI2_error1_SD_NS
#SBATCH --output=GPU1__CIRI2_error1_SD_NS_output.log
#SBATCH --error=GPU1__CIRI2_error1_SD_NS_error.log
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --hint=multithread
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=17170725@siswa.um.edu.my


module load miniconda
source activate CIRI2_New

## Superdepth_No Stringent
mkdir /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRIquant_CIRI2_Superdepth_No_Stringent

for i in {006..016}
do
CIRIquant --config '/scr/user/mk_98/Simulated_reads/GRCh38.yml' -1 /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/Sim_"$i"_1.fq -2 /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/Sim_"$i"_2.fq -v -t 10 -p Sim_"$i" -o /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRIquant_CIRI2_Superdepth_No_Stringent/Sim_"$i" --circ /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRI2_superdepth_noStringency/Sim_All.ciri --tool CIRI2 --bsj-file /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRIquant_CIRI2_Superdepth_No_Stringent/Sim_BSJ_"$i".txt
done;

