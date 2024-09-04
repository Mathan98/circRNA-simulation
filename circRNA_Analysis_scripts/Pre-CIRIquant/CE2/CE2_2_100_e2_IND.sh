#! /bin/bash

#SBATCH --job-name=CE2_2_100_e2_IND
#SBATCH --output=CE2_2_100_e2_IND_output.log
#SBATCH --error=CE2_2_100_e2_IND_error.log
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

mkdir /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRCexplorer2_Individual

for i in {001..016}
do
#tophat_align
tophat2 -o /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRCexplorer2_Individual/tophat_fusion_"$i" -p 10 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search /scr/user/mk_98/Simulated_reads/bowtie_index/Homo_sapiens /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/Sim_"$i"_1.fq /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/Sim_"$i"_2.fq

CIRCexplorer2 parse --pe -t TopHat-Fusion /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRCexplorer2_Individual/tophat_fusion_"$i"/accepted_hits.bam -b /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRCexplorer2_Individual/back_spliced_junction_"$i".bed > /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRCexplorer2_Individual/CIRCexplorer_"$i"_parse.log

CIRCexplorer2 annotate -r /scr/user/mk_98/Simulated_reads/hg38_ref_all.txt -g /scr/user/mk_98/Simulated_reads/Homo_sapiens.fa -b /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRCexplorer2_Individual/back_spliced_junction_"$i".bed -o /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRCexplorer2_Individual/circularRNA_known_"$i".txt > /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage10_20/CIRCexplorer2_Individual/CIRCexplorer2_"$i"_annotate.log
done
