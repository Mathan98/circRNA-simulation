# circRNA-simulation
The "circRNA-simulation" repository provides a comprehensive framework for simulating circular RNA (circRNA) datasets.

This perl script is originally made by CIRI-full authors for simulating circRNAs using a reference fasta file and a gtf annotation file.

Originally, there is a lack of manual to use this simulation tool. Hence, we describe in detail here regarding the output file.

Moreover, the perl script deposited here is **MODIFIED** to include additional functions. Modifications include:
i) manual read id designation (maybe useful when simulation involves pooling of data) (-RN parameter)
ii) new tab-delimited file that describes BSJ count number for each circRNA simulated.
iii) circRNA reads are defined as BSJs **ONLY** if they overlap the left and right sequence of a BSJ by **5 nucleotides** (as oppose to 20 in the original). This is according to the default parameter of the CIRIquant tool. Can be manually changed at line 472-479

Output file: 
**.out file**
1st row = chromosome id, gene id, transcript id, circRNA id
2nd row = isoform num (if there is more than one isoform other than "isoform1") separated by "_" with the length of simulated circRNA (in number format), exon information separated by commas ( exon1:exon2!strand ) 
3rd row until the next 1st row = ">", read number for that circRNA, read number in fastq file
                 = if "**" indicates BSJ read id, with the number beside as either read 1 or read 2


**bsj_counts.out**

1st column = circRNA id
2nd column = BSJ read count simulated



**References**

Gao, Y., Wang, J., & Zhao, F. (2015). CIRI: an efficient and unbiased algorithm for de novo circular RNA identification. Genome Biology, 16(1), 4. https://doi.org/10.1186/s13059-014-0571-3 

Zheng, Y., Ji, P., Chen, S., Hou, L., & Zhao, F. (2019). Reconstruction of full-length circular RNAs enables isoform-level quantification. Genome Medicine, 11(1), 2. https://doi.org/10.1186/s13073-019-0614-1 

Zhang, J., Chen, S., Yang, J., & Zhao, F. (2020). Accurate quantification of circular RNAs identifies extensive circular isoform switching events. Nature Communications, 11(1), 90. https://doi.org/10.1038/s41467-019-13840-9 





