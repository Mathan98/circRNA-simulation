#! /bin/bash

#SBATCH --job-name=FC_2_100_IND
#SBATCH --output=FC_2_100_IND_output.log
#SBATCH --error=FC_2_100_IND_error.log
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --ntasks=12
#SBATCH --hint=multithread
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=17170725@siswa.um.edu.my

module load miniconda
source activate Find_circ

#bowtie2-build /scr/user/mk_98/Simulated_reads/Homo_sapiens.fa /scr/user/mk_98/Simulated_reads/bowtie2_index/Homo_sapiens

## find_circ Individual
mkdir /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/find_circ_Individual

for i in {001..016}
do

## Find_circ Superdepth
bowtie2 -p 10 --very-sensitive --score-min=C,-15,0 --mm -x /scr/user/mk_98/Simulated_reads/bowtie2_index/Homo_sapiens -q -U <(cat /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/Sim_"$i"_1.fq /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/Sim_"$i"_2.fq) 2> /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/find_circ_Individual/Sim_"$i"_bowtie2.log | samtools view -hbuS - | samtools sort -T sorted -o /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/find_circ_Individual/Sim_"$i".bam

##get the unmapped and pipe through unmapped2anchors.py
samtools view -hf 4 /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/find_circ_Individual/Sim_"$i".bam | samtools view -Sb - | python2.7 /scr/user/mk_98/Simulated_reads/find_circ-master/unmapped2anchors.py -  | bowtie2 -p 10 --reorder --score-min=C,-15,0 -q -x /scr/user/mk_98/Simulated_reads/bowtie2_index/Homo_sapiens -U - 2> /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/find_circ_Individual/Sim_"$i"_bowtie2_unmapped.log | python2.7 /scr/user/mk_98/Simulated_reads/find_circ-master/find_circ.py -G "/scr/user/mk_98/Simulated_reads/Homo_sapiens.fa" --prefix=Sim_"$i"_find_circ --name=Sim_"$i" --stats=/scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/find_circ_Individual/Sim_"$i"_stats.txt --reads=/scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/find_circ_Individual/Sim_"$i"_spliced_reads.fa > /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/find_circ_Individual/Sim_"$i"_splice_sites.bed
 
grep CIRCULAR /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/find_circ_Individual/Sim_"$i"_splice_sites.bed | grep -v chrM | awk '$5>=1' | grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | awk '$9 == 40 && $10 == 40'| python2.7 /scr/user/mk_98/Simulated_reads/find_circ-master/maxlength.py 100000 > /scr/user/mk_98/Simulated_reads/Latest_constant_coverage_2024/circlin_reads_chrAll_error1_coverage5_40/find_circ_Individual/Sim_"$i"_circ_candidates_bsj1.bed
done
