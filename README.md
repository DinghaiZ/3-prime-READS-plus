# 3-prime-READS-plus
## Pipeline for analyzing 3' end RNA-seq (3'READS+) data 

[3'READS+](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5029459/) is a 3' end RNA-seq protocol that produces sequencing reads with variable number of T-streches in the 5' end of the reads. Steming from partial digestion of the poly(A) tail and absent from the genome, these T-streches can be used to identify cleavage and polyadenylation sites in the genome. However, long T-streches can also mislead alignment of the reads to the genome. Moreover, due to micro-heterogeneity, the mapped reads that are close to each other in the genome need to be clustered into cleavage and polyadenylation (pA) sites.

To solve the above issues, this pipeline (see 3-prime-READS-plus.ipynb) will trim 5' T-streches while recording T-strech length, map (using [STAR](https://github.com/alexdobin/STAR)) cleaned reads to the genome, and use the recorded T-strech length and genomic alignment result to identify PASS (PolyA Site Supporting) reads, which are defined as reads containing >= 2 extra Ts originated from the poly(A) tail but not from the genome. Each PASS read comes from a cleavage and polyadenylation site (CPS) in the genome. PASS reads within a 24-nt window are then clustered to define a pA site (PAS). Besides creating a fastq QC report, a PASS read QC report, and genome browser tracks for visualizing PASS reads (if http server is available), this pipeline counts number of PASS reads mapped to each CPS and PAS in each input sample. 


## Hardware/Software Requirements


## Installation


## Quick Start
To use this pipeline, please follow the following steps. 

https://github.com/bcbio/bcbio-nextgen
