# 3-prime-READS-plus
## Pipeline for analyzing 3' end RNA-seq (3'READS+) data 

About 70% of mRNA genes in eukaryotes contain multiple cleavage and polyadenylation sites (PAS), resulting in alternative cleavage and polyadenylation (APA) isoforms with different coding sequences and/or variable 3′ untranslated regions (3′UTRs). [3'READS+](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5029459/) is a 3' end RNA-seq protocol for accurate identification of PAS in the genome, even in A-rich regions:

<p align="center">
  <img src="images/flowchart.jpg" width="500" height="500">
</p>

Steming from incomplete digestion of the poly(A) tail and absent from the genome, the T-streches at the 5' end of the reads can be used to identify cleavage and polyadenylation sites in the genome:

<p align="center">
  <img src="images/read.jpg" width="350" height="210" class="center">
</p>

However, long T-streches can also mislead alignment of the reads to the genome. In addition, due to micro-heterogeneity, the mapped reads that are close to each other in the genome need to be clustered into cleavage and polyadenylation (pA) sites.

To solve the above issues, this pipeline will trim 5' T-streches while recording T-strech length, map (using [STAR](https://github.com/alexdobin/STAR)) cleaned reads to the genome, and use the recorded T-strech length and genomic alignment result to identify PASS (PolyA Site Supporting) reads, which are defined as reads containing >= 2 extra Ts originated from the poly(A) tail but not from the genome. Each PASS read comes from a cleavage and polyadenylation site (CPS) in the genome. Due to microheterogeneity during cleavage and polyadenylation, the CPSs tend to form clusters in a small window in the genome. Therefore PASS reads within a 24-nt window are clustered to define a pA site (PAS). The numbers of PASS reads mapped to genome-wide CPS and PAS in different input sample are then used for further analysis: 

<p align="center">
  <img src="images/pipeline.png" width="500" height="250" class="center">
</p>

**The pipeline has the following two parts:**

**[Part 1](https://github.com/DinghaiZ/3-prime-READS-plus/blob/master/projects/project_1/experiment_1/notebooks/Part-1.ipynb). Read QC, generate PASS read count matrix, and visualize PASS and nonPASS reads in UCSC genome browser.**

**[Part 2](https://github.com/DinghaiZ/3-prime-READS-plus/blob/master/projects/project_1/experiment_1/notebooks/Part-2.ipynb). pA site annotation and feature extraction.**  

**What are in each folder?**
*modules*: Definitions of functions and classes used in the pipeline.
*notebooks*: Template Jupyter notebooks that should be copied into each /projects/project_name/experiment_name folder, edited, and run for each experiment under different projects.
*projects*: A tree-like directory containing one *project_name/experiment_name* folder for each experiment under different projects. Within each *project_name/experiment_name* folder, there are three subfolders: *notebooks* (code for project and experiment-specific analysis), *data* (data for project and experiment-specific analysis), and *results* (analysis results). 
*tests*: Code for testing my pipeline during development. Not needed for running the notebook.
*images*: Images displayed in the README and notebooks on github. Not needed for running the notebook.

## Quick Start

**Before each analysis, please do the following:**
1. git clone the repo to your computer.
2. Copy the notebooks folder into the 'projects/project_name/experiment_name' folder. You can name the projects and experiments the way you want.
3. Edit the analysis configuration sections (at the beginning) of the notebooks for project and experiment-specific analysis.
4. Make sure that 3rd party softwares have been installed on your system and are in your path.
5. Run the notebooks and have fun!



