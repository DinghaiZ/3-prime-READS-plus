# 3-prime-READS-plus
Pipeline for analyzing 3' end RNA-seq (3'READS+) data

3'READS+ is a 3' end RNA-seq protocol that produces sequencing reads with variable number of T-streches in the beginning (5') of the reads. Steming from partial digestion of the poly(A) tail and absent from the genome, these T-streches can be used to identify cleavage and polyadenylation sites in the genome. However, long T-streches can also mislead alignment of the reads to the genome. Moreover, due to micro-heterogeneity, the mapped reads that are close to each other in the genome need to be clustered into cleavage and polyadenylation (pA) sites.

To solve the above issues, this pipeline will trim 5' T-streches while recording T-strech length, map cleaned reads to the genome, and use the recorded T-strech length and genomic alignment result to identify PASS (PolyA Site Supporting) reads, which are reads containing >= 2 extra Ts originated from the poly(A) tail but not from the genome. Each PASS read comes from one cleavage and polyadenylation (CS) site in the genome. PASS reads within a 24-nt window are then clustered to define pA sites. Number of PASS reads mapped to each CS and pA sites are counted for each input sample. In addition, this pipeline also generate a fastq QC report, a PASS read QC report, and create genome browser tracks for visualizing PASS reads in the genome (optional).

To use this pipeline, please follow the following steps. 
