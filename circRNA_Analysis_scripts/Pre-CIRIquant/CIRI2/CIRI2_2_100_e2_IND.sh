#! /bin/bash

#SBATCH --job-name=CIRI2_2_100_IND
#SBATCH --output=CIRI2_2_100_IND_output.log
#SBATCH --error=CIRI2_2_100_IND_error.log
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

mkdir /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_Ind_No_Stringent
mkdir /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_Ind_Stringent

for i in {001..016}
do
/home/user/mk_98/.conda/envs/CIRI/bin/bwa mem -t 10 -T 19 /scr/user/mk_98/Simulated_reads/Homo_sapiens.fa /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/Sim_"$i"_1.fq /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/Sim_"$i"_2.fq 1> /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_Ind_No_Stringent/Sim_"$i"_unmapped.sam 2> /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_Ind_No_Stringent/Sim_"$i"_unmapped.log

## CIRIquant_noStringency
perl /home/user/mk_98/CIRIquant/libs/CIRI2.pl -t 10 -I /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_Ind_No_Stringent/Sim_"$i"_unmapped.sam -O /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_Ind_No_Stringent/Sim_"$i".ciri -F /scr/user/mk_98/Simulated_reads/Homo_sapiens.fa -A /scr/user/mk_98/Simulated_reads/human_gencode_vch38.gtf -0 -G /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_Ind_No_Stringent/CIRI2_"$i".log

## CIRIquant_Stringent
perl /home/user/mk_98/CIRIquant/libs/CIRI2.pl -t 10 -I /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_Ind_No_Stringent/Sim_"$i"_unmapped.sam -O /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_Ind_Stringent/Sim_"$i".ciri -F /scr/user/mk_98/Simulated_reads/Homo_sapiens.fa -A /scr/user/mk_98/Simulated_reads/human_gencode_vch38.gtf -G /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRI2_Ind_Stringent/CIRI2_"$i".log
done