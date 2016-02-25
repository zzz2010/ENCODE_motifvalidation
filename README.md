ENCODE_motifvalidation pipeline
====
Three types of motif enrichment scores are computed by overlapping the motif instances with the given ChIP-seq peak locations, which includes global enrichment z-score (compare the actual motif with the shuffled-version motif), positional bias z-score (compare peak center to the peak flanking region, +/-100bp), and peak rank bias z-score (compare the high signal value regions to the low signal value regions). Then, a combined enrichment score can be derived from taking average of three enrichment scores listed above.  Forth, the known motifs are grouped as 282 clusters by their PWM similarities and each motif cluster is ranked by the highest combined enrichment score of that motif cluster. Therefore, each known motif is assigned to the ranking of the motif cluster it belong to.

We developed a Bayesian approach assuming the enrichment ranking distribution of the known motif following a mixture of two negative-binomial distributions, which corresponding to two cases: 1. antibody targeting corresponding TF of that motif, and  2. Antibody is targeting other TFs. If these two negative-binomial distributions for the tested motif is known, we can derived the accept/reject probability of the tested antibody given the enrichment ranking of that motif. 

The parameters of two negative-binomial distribution is estimated from the previous ChIP-seq data of that TF and all other TFs (stored in pipeline_script/cistrome_model.pickle). This accepted probability calculation also assumes the prior probabilities of passing validation is 50%, and takes into account that different TFs may share the same motif and one TF may use multiple motifs. 

Software Requirment
----
- add "bin" folder to your enviroment PATH , and these are several utility scripts created by Pouya Kheradpour
- Python 2.7, Numpy, Scipy, sklearn
- Bedtool
- R-3.1

Input
-----
- BED format peak file with at least 5 column, and assume the 5th column is the peak score which usually corresponding to the 7th column in the ENCODE NarrowPeak format.
- reference genome, default is hg19
- motif instances are defined and can be downloaded in http://compbio.mit.edu/encode-motifs/  

Run Example
---
'cd pipeline_script; make test'


This pipeline has only been tested in RedHat 6.7


Author
---
Zhizhuo Zhang (zhizhuo@mit.edu)
MIT Kellis Lab 
