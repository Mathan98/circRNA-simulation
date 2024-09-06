# circRNA-simulation
The "circRNA-simulation" repository provides a comprehensive framework for simulating circular RNA (circRNA) and linear RNA reads in a RNA-seq format (fastq files). We also included **circRNA analysis runs (CIRI2, CIRCexplorer2 and find_circ)**, **benchmarking** and **visualization** scripts.

The ciri-simulator perl script is originally made by CIRI authors for simulating circRNAs by using reference fasta file and gtf annotation file. Originally, there is a lack of manual to use this simulation tool and intepret its output. Hence, we describe in detail here regarding the use and output file.


Moreover, the perl script deposited here is **MODIFIED** to include additional functions. Modifications include:
```
i) manual read id designation (maybe useful when simulation involves pooling of data) (-RN parameter)

ii) new tab-delimited file that describes BSJ count number for each circRNA simulated.

iii) circRNA reads are defined as BSJs **ONLY** if they overlap the left and right sequence of a BSJ by **5 nucleotides** (as oppose to 20 in the original perl script). The overlapping 5 nucleotides parameter is according to the default parameter of CIRIquant. This value can be manually changed at line 472-479 of the perl script.
```

## A total of 4 output files: circRNA fastq, linearRNA fastq, out file and BSJ reads file

### Output file (**.out file**) has 3 rows for each circRNA read:
```

1st row = chromosome id, gene id, transcript id, circRNA id

2nd row = isoform num (if there is more than one isoform other than "isoform1") separated by "_" with the length of simulated circRNA (in number format), exon information separated by commas ( exon1:exon2!strand ) 

3rd row until the next 1st row = ">", read number for that circRNA, read number in fastq file
                 = if "**" indicates BSJ read, with the number next to it describing either read 1 or read 2

```

### bsj_counts.out
```
1st column = circRNA id

2nd column = BSJ read count simulated

```

### Example run
```

```


## References

Gao, Y., Wang, J., & Zhao, F. (2015). CIRI: an efficient and unbiased algorithm for de novo circular RNA identification. Genome Biology, 16(1), 4. https://doi.org/10.1186/s13059-014-0571-3 

Zheng, Y., Ji, P., Chen, S., Hou, L., & Zhao, F. (2019). Reconstruction of full-length circular RNAs enables isoform-level quantification. Genome Medicine, 11(1), 2. https://doi.org/10.1186/s13073-019-0614-1 

Zhang, J., Chen, S., Yang, J., & Zhao, F. (2020). Accurate quantification of circular RNAs identifies extensive circular isoform switching events. Nature Communications, 11(1), 90. https://doi.org/10.1038/s41467-019-13840-9 





