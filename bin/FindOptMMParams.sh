#!/bin/sh -x

WINDOWS=$(echo 0 5 `seq 10 10 100` `seq 120 20 500` | tr '\n' ' ')
FLIPS="0 1"
CUTS=`seq -w 0.0 0.1 1.01 | tr '\n' ' '`
CONFIDZ=0
NUMMOTIFSKEEP=10
QUEUE=""
HSTR=0
MATCHMASKS=-1
CLUSTCORR=0.8
REGWINDOW=1000
CREATEMATCHES=1
USEFLATDIR=0
NUMSHUFBINS=1
NUMOTHER=0
OTHERREGIONS="/dev/null"
TOLDIFF=0.2
NUMNODES=500

########################################## PARSE ARGUMENTS
if [ $# -eq 0 ]; then
		echo >&2 "USAGE: $0 [OPTIONS] 
 -m     File listing motifs to consider [required]
 -C     Changes from previous run for motifs. NewName OldName [optional]
 -r     File containing regions to use [required]
 -o     Organism. See GetOrgVars.sh for details. (support for -a, -n, -B, -e is deprecated)
 -E     Directory to use for storing region files. [default: -d]
 -S     Do all analysis only on the first strand (0/1) [default: $HSTR]
 -x     Create matches file (0/1) [default: 1]
 -k     Space separated list of species subsets to use (-1 indicates all species) [default $MATCHMASKS]
 -z     Z value used to correct confidence measurement [default: $CONFIDZ]
 -F     Flips to use (always just 0 for -S 1; empty to exit after producing controls) [default: $FLIPS]
 -w     Windows to use [default: $WINDOWS]
 -W     Surrounding region to extract to allow for window scanning [default: $REGWINDOW]
 -D     Use flat output directories (i.e. do not include matchmask; MATCHMASKS must contain only one species) [default: $USEFLATDIR]
 -f     Cutoffs to use on confidence [default: $CUTS]
 -q     Queue to use [default: $QUEUE]
 -d     Output directory [required]
 -N     Number of nodes to use in submissions [default: $NUMNODES]

Parameters applying to control motifs:
 -M     File listing motifs in addition to -m to avoid in controls (- disables avoiding known; -- avoids only self) [default: -m motifs]
 -U     Number of control motifs to keep [default: $NUMMOTIFSKEEP]
 -h     Number of bins to use when shuffling motifs [default: $NUMSHUFBINS]
 -c     Correlation cutoff used by motif-match for clustering of control motifs [default: $CLUSTCORR]
 -t     Toldiff paramter (-1 to disable) [default: $TOLDIFF]

Deprecated (use -o):
 -a     Alignment file to extract regions from [required]
 -n     Number of species in alignment file [required]
 -B     Branch length file in parent tree format [required]
"
	exit 1
fi

while getopts m:r:a:n:d:z:M:r:S:U:k:c:w:W:x:B:D:F:q:E:h:o:C:t:f:N: o
do      case "$o" in
		m)		MOTIFFILE="$OPTARG";;
		C)		MOTIFCHANGESFILE="$OPTARG";;
		M)		CMOTIFFILE="$OPTARG";;
		r)		REGIONSFILE="$OPTARG";;
		E)		RDIR="$OPTARG";;
		c)		CLUSTCORR="$OPTARG";;
		t)		TOLDIFF="$OPTARG";;
		h)		NUMSHUFBINS="$OPTARG";;
		k)		MATCHMASKS="$OPTARG";;
		U)		NUMMOTIFSKEEP="$OPTARG";;
		z)		CONFIDZ="$OPTARG";;
		S)		HSTR="$OPTARG";;
		w)		WINDOWS="$OPTARG";;
		F)		FLIPS="$OPTARG";;
		W)		REGWINDOW="$OPTARG";;
		x)		CREATEMATCHES="$OPTARG";;
		D)		USEFLATDIR="$OPTARG";;
		d)		DIR="$OPTARG";;
		q)		QUEUE="$OPTARG";;
		o)		ORG="$OPTARG";;
		a)		ALIGNFILE="$OPTARG";;
		B)		BLTREE="$OPTARG";;
		n)		NUMSP="$OPTARG";;
		f)		CUTS="$OPTARG";;
		N)		NUMNODES="$OPTARG";;
		[?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
		esac
done

eval $(GetOrgVars.sh $ORG)

if [ -z "$REGIONSFILE" -o -z "$NUMSP" -o -z "$BLTREE" -o -z "$ALIGNFILE" -o -z "$RDIR$DIR" ]; then
	echo >&2 "ERROR: -r -n -d/-E -B -a are required."; exit 1
fi

###################### SETUP REGIONS
if [ -z "$RDIR" ]; then
	RDIR=$DIR
fi

mkdir -p $RDIR

# clean up regions file, this isn't strictly required but can make a small difference when
# HSTR is on (it does not permit the use of regions that are overlapping in opposite
# directions because grep-overlap -t does not match "." with anything)
if [ ! -f $RDIR/regions.txt ]; then 
	rtool namerepl -r Region $REGIONSFILE | rtool prec-overlap -s $HSTR -r Region > $RDIR/regions.txt
fi

if [ ! -f $RDIR/regions-long.txt -a ! -z "$FLIPS" ]; then
	rtool window -r $REGWINDOW $RDIR/regions.txt | prec-overlap Region > $RDIR/regions-long.txt
fi

########################
if [ -z "$DIR" ]; then 
	# this will happen when only regions are being setup
	exit 0
fi

if [ -z "$MOTIFFILE" -o  -z "$DIR" -o -z "$REGIONSFILE" -o -z "$NUMSP" -o -z "$BLTREE" -o -z "$ALIGNFILE" ]; then
	echo >&2 "ERROR: -m -r -a -n -B -d -a are required."; exit 1
fi

# its important that WINDOWS be exactly a list of space separated values without a trailing space for use with branch-sn
WINDOWS=$(echo $WINDOWS | awk '{for (i=1; i<=NF; i++) L = L " " $i}; END{print substr(L,2)}')

mkdir -p $DIR

###### CREATE CONTROLS FILE
if [ -f "$MOTIFCHANGESFILE" ]; then
	# when there are changes, save the old files so we can pass them off and delete
	# everything else
	mv $DIR/motifs-toscan.txt{,.old}
	mv $DIR/final-controls.txt{,.old}
	mv $DIR/controls{,.old}
	rm -rf $DIR/wmatches-flip{0,1}.txt.gz
	for MATCHMASK in $MATCHMASKS; do
		rm -rf $DIR/$MATCHMASK
	done
	mkdir -p $DIR/controls
	mv $DIR/controls.old $DIR/controls/old
fi

if [ ! -f $DIR/final-controls.txt ]; then 
	# also deals with injection of old motifs (only if $MOTIFCHANGESFILE is found)
	CreateControlMotifs.sh \
		-d $DIR/controls -r $RDIR/regions.txt \
		-o $ORG -S $HSTR -t $TOLDIFF \
		-e $([ "$CMOTIFFILE" = "-" -o "$CMOTIFFILE" = "--" ] && echo "/dev/null" || echo $CMOTIFFILE) \
		-E $([ "$CMOTIFFILE" = "-" ] && echo 0 || ([ "$CMOTIFFILE" = "--" ] && echo 2 || echo 1)) \
		-c $CLUSTCORR -h $NUMSHUFBINS -N $NUMMOTIFSKEEP \
		-C "$MOTIFCHANGESFILE" -m $DIR/final-controls.txt.old -M $DIR/motifs-toscan.txt.old \
		-O $DIR/motifs-toscan.txt -n $NUMNODES \
		< $MOTIFFILE > $DIR/final-controls.txt
fi

if [ -z "$FLIPS" ]; then
	# we just want control motifs... exit
	exit 0
fi

##### WINDOW ANALYSIS
# do one run for each of flips/not flips with -v -2 flag which outputs verbose 
# info that is sufficient to produce values for all cutoffs
MAXWINDOW=$(echo $WINDOWS | awk -vRS="[\n\t ]" '($1>=X){X=$1} END{print X}')

if [ $HSTR -eq 1 ]; then
	FLIPS=0
fi

for f in $FLIPS; do
	if [ ! -f $DIR/wmatches-flip$f.txt.gz ] ; then
		rmc -N wmatches-flip$f -v $DIR/wmatches-flip$f.rmc -c "extract-mfa $NUMSP $ALIGNFILE < %s | motif-match -n $NUMSP -w $MAXWINDOW -s $f -m $DIR/motifs-toscan.txt -v -2 | grep-overlap -t$HSTR -r1 -o - $RDIR/regions.txt | gzip -c" -i $RDIR/regions-long.txt -o $DIR/wmatches-flip$f.txt.gz -e - -n $NUMNODES &
		PIDS="$PIDS $!"
	fi
done

#### WAIT ON SEQUENCE FILE AND WINDOW ANALYSIS
if [ ! -z "$PIDS" ]; then
	wait $PIDS
fi

## SPLIT OFF DIFFERENT MATCH MASKS HERE
LPIDS=""
for MATCHMASK in $MATCHMASKS; do
	# find actual match mask for MATCHMASK=-1 so that the directory name is nice
	if [ $MATCHMASK -eq -1 ]; then
		MATCHMASK=$(echo "2 ^ ($NUMSP - 1 ) - 1" | bc)
	fi

	if [ $USEFLATDIR -eq 1 ]; then 
		MDIR=$DIR
	else
		MDIR=$DIR/$MATCHMASK
	fi

	mkdir -p $MDIR/{bsn,mmc,matches,bsub/matches}

	unset PIDS

	{
		## CREATE THE INDIVIDUAL WINDOW FILES (NOW DONE ALL AT ONCE)
		for f in $FLIPS; do
			if [ ! -f $MDIR/bsn/$MAXWINDOW-$f ]; then
				psub -N bsn-$f -b -v $MDIR/bsub/bsn-$f "grep-overlap -v -r1 $DIR/wmatches-flip$f.txt.gz $OTHERREGIONS | cut -d' ' -f1,7 | branch-sn -mm 2 -w '$WINDOWS' -k $MATCHMASK -f $DIR/final-controls.txt -bp $BLTREE -n $NUMSP -z $CONFIDZ -o $MDIR/bsn/%d-$f -" &
				PIDS="$PIDS $!"
			fi
		done

		if [ ! -z "$PIDS" ]; then
			wait $PIDS
		fi

		# create a table with columns: 1) motif, 2) window-flip, 3) cutoff, 4) confidence (min above cutoff), 5) branch length, 6) number of instances
		if [ ! -f $MDIR/w-table.txt ]; then
			for w in $WINDOWS; do
				for f in $FLIPS; do
					for c in $CUTS; do
						awk -vC=$c -vW=$w-$f '{
							for (i=2; i<=NF; i++)
							{
								split($i, A, ";")
								if (A[1]>= C)
								{
									print $1, W, C, A[1], A[3], A[2]
									break
								}
							}
						}' $MDIR/bsn/$w-$f
					done
				done
			done > $MDIR/w-table.txt
		fi

		# make sure the w-table is ok (basically a  check of the wmatches files), otherwise exit
		if ! awk '$3 != "0.0"{next}; $2=="0-0"{X[$1]=$6}; $2=="0-1"{Y[$1]=$6}; END{for (i in X) if (i in Y) if (X[i] != Y[i]) exit 1}' $MDIR/w-table.txt
		then
			echo "Inconsistency in $MDIR/w-table.txt!" >&2
			exit 1
		fi

		PIDS=""
		for c in $CUTS; do
			# use the table above to create motif-match files at different confidences using the best possible window
			if [ ! -f $MDIR/mmc/$c ]; then
				awk -vC=$c -vCF=$DIR/final-controls.txt '
					($3 == C) && $6 > M[$1]{
						split($2, WF, "-")
						M[$1]=$6
						W[$1]=WF[1]
						F[$1]=WF[2]
						B[$1]=$5
						N[$1]=$6
					}
					
					END{
						while (getline < CF)
							for (i=1; i<=NF; i++) 
								CN[$1,i]=$i
							
						for (m in M) 
							for (i=1; (m,i) in CN; i++)
								print CN[m,i], m (i>1 ? ("_C" (i-1)) : ""), W[m], F[m], B[m], N[m]
					}' $MDIR/w-table.txt > $MDIR/mmc/$c
			fi

			# create matches files by converting from the wmatches files
			if [ $CREATEMATCHES -eq 1 ]; then
				if [ ! -f $MDIR/matches/$c.gz ]; then
					#NOTE: motif-match-conv knows to only open flip1 if necessary
					#NOTE: there is no longer an option to "USEORIGREG"! 
					# ->see FindOptMMParams.sh-pre-fast-final-matches for what it used to do
					{
						psub -N matches-$c -b -v $MDIR/bsub/matches/$c -q "$QUEUE" "
							motif-match-conv -i0 $DIR/wmatches-flip0.txt.gz \
											 -i1 $DIR/wmatches-flip1.txt.gz \
											 -m $MDIR/mmc/$c -bp $BLTREE    \
											 -k $MATCHMASK -n $NUMSP        \
								| grep-overlap -v -r1 - $OTHERREGIONS       \
								| gzip -c" > $MDIR/matches/$c.gz &
						PIDS2="$!"

						# also create match files for other regions (each with own tree)
						for i in $(seq 1 $NUMOTHER); do
							psub -N matches-$c-o$i -b -v $MDIR/bsub/matches/$c-$i -q "$QUEUE" "
								motif-match-conv -i0 $DIR/wmatches-flip0.txt.gz  \
												 -i1 $DIR/wmatches-flip1.txt.gz  \
												 -m $MDIR/mmc/$c -bp ${OTHERBLTREE[$i]} \
												 -k $MATCHMASK -n $NUMSP         \
									| grep-overlap -o -r1 - ${OTHERREGIONS[$i]}  \
									| gzip -c" > $MDIR/matches/$c-$i.gz &
							PIDS2="$PIDS2 $!"
						done
						wait $PIDS2

						# merge all the files together
						if [ $NUMOTHER -gt 0 ]; then
							rtool merge $MDIR/matches/$c{,-*}.gz | gzip -c > $MDIR/matches/$c.gz.new
							mv $MDIR/matches/$c.gz.new $MDIR/matches/$c.gz
							rm $MDIR/matches/$c-*.gz
						fi
					} &
					PIDS="$PIDS $!"
				fi
			fi
		done

		if [ ! -z "$PIDS" ]; then
			wait $PIDS
		fi
	} & 
	LPIDS="$LPIDS $!"

	# only allow one MATCHMASK if USEFLATDIR=1
	if [ $USEFLATDIR -eq 1 ] ; then 
		break
	fi
done

if [ ! -z "$LPIDS" ]; then 
	wait $LPIDS
fi

