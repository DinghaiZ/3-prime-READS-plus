# 3-prime-READS-plus
## Pipeline for analyzing 3' end RNA-seq (3'READS+) data 

[3'READS+](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5029459/) is a 3' end RNA-seq protocol that produces sequencing reads with variable number of T-streches in the 5' end of the reads:![3'READS+ Flow Chart](images/flowchart.jpg) 

Steming from partial digestion of the poly(A) tail and absent from the genome, these T-streches can be used to identify cleavage and polyadenylation sites in the genome:![3'READS+ Read Structure](images/read.jpg) 

However, long T-streches can also mislead alignment of the reads to the genome. Moreover, due to micro-heterogeneity, the mapped reads that are close to each other in the genome need to be clustered into cleavage and polyadenylation (pA) sites.

To solve the above issues, the step 1 [note book](https://github.com/DinghaiZ/3-prime-READS-plus/blob/master/3%60READS%2B%20Step-1.ipynb) of this pipeline will trim 5' T-streches while recording T-strech length, map (using [STAR](https://github.com/alexdobin/STAR)) cleaned reads to the genome, and use the recorded T-strech length and genomic alignment result to identify PASS (PolyA Site Supporting) reads, which are defined as reads containing >= 2 extra Ts originated from the poly(A) tail but not from the genome. Each PASS read comes from a cleavage and polyadenylation site (CPS) in the genome. PASS reads within a 24-nt window are then clustered to define a pA site (PAS). The numbers of PASS reads mapped to genome-wide CPS and PAS in different input sample are then used for further analysis. 


## Quick Start
Step 1. Run the notebook [note book](https://github.com/DinghaiZ/3-prime-READS-plus/blob/master/3%60READS%2B%20Step-1.ipynb). 
Step 2. pA site annotation (Coming soon)
Step 3. APA and DE analysis (Coming soon)

