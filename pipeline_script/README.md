
Run Demo
===========
- type "make test" or type "python2.7 pipeline.py test/K562_CTCF_UW.bed"
- output will be "/tmp/motifpipeline/K562_CTCF_UW"


Commandline Workflow
==============

Reformat Peak file
-----------------
* center 100bp window for each peak
* sort the peak based on the score and compute the peak rank


Motif Overlap
--------------
* given reference genome, overlap pouya motif instance with the peak file


Compute Global Enrichment score
---------------
Enricher.sh (developped by Pouya )


Compute Positional and Peak-Rank Enrichment score
-----------------
computeCentdistScore.py 


Compute Avg-combine Score
---------------
combineAvgScore.py


Compute Accept and Reject Probability
------------------------
AcceptReject.py


Generate Html output
--------------------------
makehtml_png.py
