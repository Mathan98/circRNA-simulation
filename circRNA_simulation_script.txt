Command:
perl CIRI_simulator.pl -O {Dataset_name_ID} -G human_gencode_vch38.gtf -C {circRNA_coverage} -LC 0 -R 1 -LR 1 -L 120 -E {sequencing_error} -D ChrAll_fa/ -CHR1 0 -M 320 -M2 550 -PM 0 -S 70 -S2 70 -SE 0 -PSI 0 -RN {Dataset_ID}


## CircRNA Coverage ##
======================
2X => 2
5X => 5
10X => 10
20X => 20

## Sequencing error ##
======================
1% => 1
2% => 2
5% => 5




###########################
For example:

\\\
if circRNA coverage 2X, sequencing error 1% and want to simulate 16 replicates
==============================================================================
for i in {001..016}
do
perl CIRI_simulator.pl -O Sim_"$i"_circ -G human_gencode_vch38.gtf -C 2 -LC 0 -R 1 -LR 1 -L 120 -E 1 -D ChrAll_fa/ -CHR1 0 -M 320 -M2 550 -PM 0 -S 70 -S2 70 -SE 0 -PSI 0 -RN "$i"
done


\\\
if circRNA coverage 2X, sequencing error 5% and want to simulate 16 replicates
==============================================================================
for i in {001..016}
do
perl CIRI_simulator.pl -O Sim_"$i"_circ -G human_gencode_vch38.gtf -C 2 -LC 0 -R 1 -LR 1 -L 120 -E 5 -D ChrAll_fa/ -CHR1 0 -M 320 -M2 550 -PM 0 -S 70 -S2 70 -SE 0 -PSI 0 -RN "$i"
done