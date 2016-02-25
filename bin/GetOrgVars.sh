#!/bin/sh

# use with: eval $(GetOrgVars.sh $ORG)
# or with GetOrgVars.sh $ORG $VAR

VARS="NUMSP BLTREE ALIGNFILE SPECIES NUMOTHER OTHERREGIONS OTHERBLTREE[1] OTHERREGIONS[1] ORGREGIONDIR DEFMATCHMASK"

if [ -z "$1" ]; then 
	echo ""
	exit
fi

case "$1" in
	dm3)
		NUMSP=12
		SPECIES="fly"
		BLTREE="/broad/compbio/pouyak/data/flies/dm3/ucsc12way/ptree.txt"
		ALIGNFILE="/broad/compbio/pouyak/data/flies/dm3/ucsc12way/alignment-files.txt"
		ORGREGIONDIR="/broad/compbio/pouyak/motifs/flies/data/regions/dm3"
		;;
	ce6)
		NUMSP=6
		SPECIES="worm"
		BLTREE="/broad/compbio/pouyak/data/worms/ce6/ucsc6way/ptree.txt"
		ALIGNFILE="/broad/compbio/pouyak/data/worms/ce6/ucsc6way/alignment-files.txt"
		ORGREGIONDIR="/broad/compbio/pouyak/motifs/worms/data/regions/ce6"
		;;
	ce10)
		NUMSP=7
		SPECIES="worm"
		BLTREE="/broad/compbio/pouyak/data/worms/ce10/ucsc7way/ptree.txt"
		ALIGNFILE="/broad/compbio/pouyak/data/worms/ce10/ucsc7way/alignment-files.txt"
		ORGREGIONDIR="/broad/compbio/pouyak/motifs/worms/data/regions/ce10"
		;;
	hg18old)
		NUMSP=28
		SPECIES="human"
		BLTREE="/home/zzhang/motifpipeline/insts/hg18/ucsc28way/ptree.txt"
		ALIGNFILE="/home/zzhang/motifpipeline/insts/hg18/ucsc28way/alignment-files.txt"
		ORGREGIONDIR="/home/zzhang/motifpipeline/insts/hg18"
		;;
	hg18)
		NUMSP=32
		SPECIES="human"
		BLTREE="/home/zzhang/motifpipeline/insts/hg18/ucsc32way/ptree.txt"
		ALIGNFILE="/home/zzhang/motifpipeline/insts/hg18/ucsc32way/alignment-files.txt"
		NUMOTHER=1
		OTHERREGIONS="/home/zzhang/motifpipeline/insts/hg18/ucsc32way/region-chrX.txt"
		OTHERBLTREE[1]="/home/zzhang/motifpipeline/insts/hg18/ucsc32way/ptree-chrX.txt"
		OTHERREGIONS[1]="/home/zzhang/motifpipeline/insts/hg18/ucsc32way/region-chrX.txt"
		ORGREGIONDIR="/home/zzhang/motifpipeline/insts/hg18"
		DEFMATCHMASK=268435455
		;;
	hg19)
		NUMSP=29
		SPECIES="human"
		BLTREE="/home/zzhang/motifpipeline/insts/hg19/ucsc29way/ptree.txt"
		ALIGNFILE="/home/zzhang/motifpipeline/insts/hg19/ucsc29way/alignment-files.txt"
		NUMOTHER=1
		OTHERREGIONS="/home/zzhang/motifpipeline/insts/hg19/ucsc29way/region-chrX.txt"
		OTHERBLTREE[1]="/home/zzhang/motifpipeline/insts/hg19/ucsc29way/ptree-chrX.txt"
		OTHERREGIONS[1]="/home/zzhang/motifpipeline/insts/hg19/ucsc29way/region-chrX.txt"
		ORGREGIONDIR="/home/zzhang/motifpipeline/insts/hg19"
		;;
	hg19all)
		NUMSP=46
		SPECIES="human"
		BLTREE="/home/zzhang/motifpipeline/insts/hg19/ucsc46way/ptree.txt"
		ALIGNFILE="/home/zzhang/motifpipeline/insts/hg19/ucsc46way/alignment-files.txt"
		NUMOTHER=1
		OTHERREGIONS="/home/zzhang/motifpipeline/insts/hg19/ucsc46way/region-chrX.txt"
		OTHERBLTREE[1]="/home/zzhang/motifpipeline/insts/hg19/ucsc46way/ptree-chrX.txt"
		OTHERREGIONS[1]="/home/zzhang/motifpipeline/insts/hg19/ucsc46way/region-chrX.txt"
		ORGREGIONDIR="/home/zzhang/motifpipeline/insts/hg19"
		;;
	hg38)
		NUMSP=20
		SPECIES="human"
		BLTREE="/home/zzhang/motifpipeline/insts/hg38/ucsc20way/ptree.txt"
		ALIGNFILE="/home/zzhang/motifpipeline/insts/hg38/ucsc20way/alignment-files.txt"
		ORGREGIONDIR="/home/zzhang/motifpipeline/insts/hg38"
		;;
	mm9)
		NUMSP=21
		SPECIES="mouse"
		ALIGNFILE="/home/zzhang/motifpipeline/insts/mm9/ucsc21way/alignment-files.txt"
		ORGREGIONDIR="/home/zzhang/motifpipeline/insts/mm9"
		BLTREE="/home/zzhang/motifpipeline/insts/mm9/ucsc21way/ptree.txt"
		DEFMATCHMASK=1048575
		;;
	mm10)
		NUMSP=60
		SPECIES="mouse"
		ALIGNFILE="/home/zzhang/motifpipeline/insts/mm10/ucsc60way/alignment-files.txt"
		ORGREGIONDIR="/home/zzhang/motifpipeline/insts/mm10"
		BLTREE="/home/zzhang/motifpipeline/insts/mm10/ucsc60way/ptree.txt"
		DEFMATCHMASK=1099511627775
		;;
	AgamP3)
		NUMSP=21
		SPECIES="mosquito"
		ALIGNFILE="/broad/compbio/pouyak/data/mosqs/AgamP3/rmwater/alignment-files.txt"
		ORGREGIONDIR="/broad/compbio/pouyak/motifs/mosqs/data/regions/AgamP3"
		BLTREE="/broad/compbio/pouyak/data/mosqs/AgamP3/rmwater/ptree.txt"
		;;
	taeGut2)
		NUMSP=1
		SPECIES="zebrafinch"
		ALIGNFILE="/home/zzhang/motifpipeline/insts/taeGut2/noalign/alignment-files.txt"
		ORGREGIONDIR="/home/zzhang/motifpipeline/insts/taeGut2"
		BLTREE="/home/zzhang/motifpipeline/insts/taeGut2/noalign/ptree.txt"
		;;

	*)echo "exit 1"; exit ;; 
esac

if [ -z "$DEFMATCHMASK" -a ! -z "$NUMSP" ]; then
	DEFMATCHMASK=$(bc <<< "2 ^ ($NUMSP - 1) - 1")
fi

if [ ! -z "$2" ]; then
	eval echo \${$2}
else
	for v in $VARS; do
		eval test -n \"\${$v+1}\" && eval echo $v=\${$v}
	done
fi

