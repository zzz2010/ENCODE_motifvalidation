ENCODE_motifvalidation pipeline
====
Three types of motif enrichment scores are computed by overlapping the motif instances with the given ChIP-seq peak locations, which includes global enrichment z-score (compare the actual motif with the shuffled-version motif), positional bias z-score (compare peak center to the peak flanking region, +/-100bp), and peak rank bias z-score (compare the high signal value regions to the low signal value regions). Then, a combined enrichment score can be derived from taking average of three enrichment scores listed above.  Forth, the known motifs are grouped by their PWM similarities and each motif group is ranked by the highest combined enrichment scores within group. Therefore, each known motif is assigned to the ranking of the PWM group it belong to.
The accepted probability of the tested TF is calculated by a Bayesian approach comparing its known motifs' ranking in this ChIP-seq library and those known motifs in previous published TF ChIP-seq data.   This Bayesian approach assumes the known motif ranking distribution following a mixture of two negative-binomial distributions, which corresponding to cases antibody targeting corresponding TF or not. The parameters of two negative-binomial distribution is estimated from the previous ChIP-seq data of that TF and all other TFs. This accepted probability calculation assumes the prior probabilities of passing validation is 50%, and takes into account that different TFs may share the same motif and one TF may use multiple motif. 

Software Requirment
----
- add "bin" folder to your enviroment PATH , these are several scripts by Pouya Kheradpour
- Python 2.7, Numpy, Scipy, sklearn
- Bedtool
- R-3.1

Input
-----
- BED format peak file with at least 5 column, and assume the 5th column is the peak score which usually corresponding to the 7th column in the ENCODE NarrowPeak format.
- reference genome, default is hg19
- motif instance file which can be download in http://compbio.mit.edu/encode-motifs/matches.txt.gz

Run Example
---
'cd pipeline_script; make test'


This pipeline has only been tested in RedHat 6.7


Author
---
Zhizhuo Zhang (zhizhuo@mit.edu)
MIT Kellis Lab 
