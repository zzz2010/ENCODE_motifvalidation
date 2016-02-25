data/name-mapping.txt:~/compbio/projects/centdist2/result/2015-02-12/CentdistScoreForPouyaPipelineOutput/name-mapping.txt
	cut -f 2,3 $<|sed 's/_/\t/1'|cut -f 1,3|sed 's/_/\t/1'|cut -f 1,2|sort -u > $@

test:pipeline.py
	rm -rf  /tmp/motifpipeline/K562_CTCF_UW 
	python2.7 pipeline.py test/K562_CTCF_UW.bed

enricher:
	zsh /local/scratch/pipeline_script/Enricher.zzz.sh -d /tmp/motifpipeline/K562_CTCF_UW/K562_CTCF_UW.bed.global_enrichment_score -o hg19 -i /tmp/motifpipeline/K562_CTCF_UW/K562_CTCF_UW.bed.origwidth -C 0 -p  /local/scratch/pipeline_script/../insts/hg19/ -P 0 -Z 1.96 -N 0 -n 0.0 -k 268435455  -K 8
