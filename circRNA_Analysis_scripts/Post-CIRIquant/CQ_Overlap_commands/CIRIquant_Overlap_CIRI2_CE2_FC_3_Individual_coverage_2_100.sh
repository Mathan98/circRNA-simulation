#! /bin/bash

#SBATCH --job-name=Overlap_CIRI2_CE2_FC_3_coverage2_coverage100
#SBATCH --output=Overlap_CIRI2_CE2_FC_3_coverage2_coverage100_Ind_output.log
#SBATCH --error=Overlap_CIRI2_CE2_FC_3_coverage2_coverage100_Ind_error.log
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --ntasks=12
#SBATCH --hint=multithread
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=17170725@siswa.um.edu.my


module load miniconda
source activate CIRI2_New

## Individual
mkdir /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRIquant_Overlap_CIRI2_CE2_FC_3_Individual

for i in {001..016}
do
CIRIquant --config '/scr/user/mk_98/Simulated_reads/GRCh38.yml' -1 /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/Sim_"$i"_1.fq -2 /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/Sim_"$i"_2.fq -t 10 -p Sim_"$i" -o /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRIquant_Overlap_CIRI2_CE2_FC_3_Individual/Sim_"$i" --bed /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/Overlap_CIRI2_No_Stringent/Overlap_CIRI2_CE2_FC_3/Sim_"$i".bed --bsj-file /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/CIRIquant_Overlap_CIRI2_CE2_FC_3_Individual/Sim_BSJ_"$i".txt
done;