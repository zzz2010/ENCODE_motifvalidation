#!/bin/sh

EXCLUDESELF=0
EXCLUDEMOTIFS=/dev/null
NUMSHUFBINS=1
NUMMOTIFSSCAN=100
NUMMOTIFSKEEP=10
HSTR=0
TOLDIFF=0.2
CLUSTCORR=0.8
NUMNODES=500

# NOTE: does not exclude EXTRAREGIONS as previous code did; if this is desired filter
#       input regions appropriately

if [ $# -eq 0 ]; then
		echo >&2 "USAGE: $0 [OPTIONS] < INPUTMOTIFS > CONTROLFILE
  -d     Work directory [required]
  -r     Background regions to scan in [required]
  -O     File to put motifs-toscan.txt output in [required]
  -e     Exclude motifs (sort against this list) [optional]
  -E     Exclude self motifs (0/1/2; 2 indicates only exclude the original motif; 1 excludes all input motifs for all other input motifs)  [default: 0]
  -h     Number of bins to use when shuffling motifs [default: $NUMSHUFBINS]
  -N     Number of control motifs to keep [default: $NUMMOTIFSKEEP]
  -o     Organism. See GetOrgVars.sh for details [required]
  -S     Do all analysis only on the first strand (0/1) [default: $HSTR]
  -c     Correlation cutoff used by motif-match for clustering of control motifs [default: $CLUSTCORR]
  -t     Difference tolerated in number of counts (-1 to disable) [default: $TOLDIFF]
  -C     Changes from previous run for motifs. NewName OldName [optional]
  -m     Output from previous run [optional]
  -M     -O from previous run [optional]
  -n     Number of nodes to use for submission [default: $NUMNODES]
"
	exit 1
fi

while getopts d:r:e:E:h:N:o:S:c:C:m:M:t:O:x:n: o
do      case "$o" in
		d)		DIR="$OPTARG";;
		r)		REGIONS="$OPTARG";;
		e)		EXCLUDEMOTIFS="$OPTARG";;
		E)		EXCLUDESELF="$OPTARG";;
		h)		NUMSHUFBINS="$OPTARG";;
		N)		NUMMOTIFSKEEP="$OPTARG";;
		o)		ORG="$OPTARG";;
		S)		HSTR="$OPTARG";;
		c)		CLUSTCORR="$OPTARG";;
		C)		MOTIFCHANGESFILE="$OPTARG";;
		m)		OLDOUTPUT="$OPTARG";;
		M)		OLDMOTIFSTOSCAN="$OPTARG";;
		n)		NUMNODES="$OPTARG";;
		t)		TOLDIFF="$OPTARG";;
		O)		MOTIFSTOSCAN="$OPTARG";;
		x)		LISTFILE="$OPTARG";;
        [?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
        esac
done

if [ -z "$LISTFILE" ]; then 
	mkdir -p $DIR

	eval $(GetOrgVars.sh $ORG)

	cat > $DIR/input.txt

	###### SETUP EXCLUDE MOTIFS
	(
		if [ $EXCLUDESELF -eq 1 ]; then 
			cat $DIR/input.txt
		fi
		cat $EXCLUDEMOTIFS
	) > $DIR/exclude.txt


	###### INDICATES MOTIFS ARE IN PWM FORMAT
	PWMINPUT=$(awk 'NR==1{print substr($1,1,1) == ">"}' $DIR/input.txt)


	###### CLEAN UP NON-PWM MOTIFS (ONLY TAKE FIRST COLUMN; UNIQUE)
	if [ $PWMINPUT -eq 0 ]; then
		for i in input exclude; do 
			awk '{print $1}' $DIR/$i.txt | sort -u > $DIR/$i.new
			mv $DIR/$i.{new,txt}
		done
	fi


	###### PRODUCE POTENTIAL CONTROL MOTIFS
	# also includes input motifs, but excludes all motifs that we are keeping from an old run
	# TODO: test to make sure this works with new MotifShuffle.awk
	if [ ! -f $DIR/shuffles.txt ]; then 
		if [ ! -f "$MOTIFCHANGESFILE" ]; then 
			cat	
		elif [ $PWMINPUT -eq 1 ]; then 
			awk -vF=$MOTIFCHANGESFILE 'BEGIN{while (getline < F) X[$1]}; /^>/{On = !(gensub(/_[0-9]*mer$/, "", 1, substr($1,2)) in X)}; On'
		else
			awk -vF=$MOTIFCHANGESFILE 'BEGIN{while (getline < F) X[$1]}; !($1 in X)'
		fi < $DIR/input.txt | MotifShuffle.awk -vOrigName='%s' -vScrambleName='%s_shuffled_%s' -vNumBin=$NUMSHUFBINS -vNumRand=$NUMMOTIFSSCAN -vReplace=0 > $DIR/shuffles.txt
	fi


	##### FOR EACH MOTIF INDICATE LIST OF SHUFFLE MOTIFS
	if [ ! -f $DIR/shuffles-list.txt ]; then 
		cat $DIR/shuffles.txt | if [ $PWMINPUT -eq 1 ]; then
			awk '/^>/{print substr($1,2), substr($1,2)}'
		else
			cat
		fi | awk '
				L != "" && substr($2,1,length(L)) == L{
					X = X OFS $1
					next
				}
				
				{
					if (X != "") print X
					X = $1
					L = $2 "_shuffled_"
				}
				
				END{if (X != "") print X}
				' > $DIR/shuffles-list.txt
	fi


	###### CONVERT MOTIFS FOR SCANNING AGAINST GENOME
	if [ ! -f  $DIR/shuffles-toscan.txt ]; then 
		if [ $PWMINPUT -eq 1 ]; then
			motif-util -om l -im $DIR/shuffles.txt
		else
			awk '{print $1}' $DIR/shuffles.txt | sort -u
		fi > $DIR/shuffles-toscan.txt
	fi

	###### GET COUNTS FOR EACH OF THE MOTIFS
	# note use of -h 6 for much lower memory usage
	if [ ! -f $DIR/counts.txt ]; then 
		rmc -N ccm-counts -c "extract-mfa -k 1 $NUMSP $ALIGNFILE < %s | motif-match -h 6 -k 0 -n 1 -m $DIR/shuffles-toscan.txt -t ${HSTR/0/3} | awk '{X[\$1]++}; END{for (x in X) print x, X[x]}'" -i $REGIONS -n $NUMNODES -e - -v $DIR/counts.rmc | awk '{X[$1]+=$2}; END{for (x in X) print x, X[x]}' > $DIR/counts.txt
	fi


	###### GET LIST OF MOTIFS THAT PASS THE FREQUENCY CRITERION 
	if [ ! -f $DIR/freqmatch-list.txt ]; then 
		awk -vCF=$DIR/counts.txt -vTD=$TOLDIFF '
			BEGIN{while (getline < CF) C[$1]=$2}

			{
				L=$1
				for (i=2; i<=NF; i++)
					if (TD < -0.5 || ((C[$i] >= C[$1] * (1 - TD)) && (C[$i] <= C[$1] * (1 + TD))))
						L = L OFS $i
				print L
			}' $DIR/shuffles-list.txt > $DIR/freqmatch-list.txt
	fi


	###### FROM freqmatch-list CHOOSE THE FINAL LIST OF MOTIFS
	if [ ! -f $DIR/final-list.txt ]; then 
		rmc -N ccm-final -c "CreateControlMotifs.sh -d $DIR -c $CLUSTCORR -S $HSTR -E $EXCLUDESELF -N $NUMMOTIFSKEEP -x %s" -i $DIR/freqmatch-list.txt -n $NUMNODES -o $DIR/final-list.txt -e - -v $DIR/final-list.rmc
	fi

	##### CREATE FINAL MOTIF LIST TO SCAN FOR
	if [ $PWMINPUT -eq 1 ]; then
		awk -vCF=$DIR/final-list.txt '
			BEGIN{
				while (getline < CF)
					for (i=1; i<=NF; i++) O[$i]
			}
			
			{
				L=length($1)
				On = $2 in O
				if (On) 
					print
				for (i=0; i<=L; i++)
				{
					getline
					if (On) print
				}
			}
		' $DIR/shuffles-toscan.txt

		if [ -f "$MOTIFCHANGESFILE" ] ; then 
			awk -vF=$MOTIFCHANGESFILE -vOFS="\t" '
				BEGIN{
					while (getline < F)
						Old[$2] = $1
					FS="\t"
				}

				$1 ~ /^XXX*$/{
					c = gensub(/_[0-9]*mer(_shuffled_[0-9]*)?$/, "", 1, $2)
					On = c in Old
					if (On)
						print $1, Old[c] substr($2,length(c)+1)
					next
				}
				On
				' $OLDMOTIFSTOSCAN
		fi
	else
		tr '\t' '\n' < $DIR/final-list.txt | sort -u | awk 'NF!=0'

		if [ -f "$MOTIFCHANGESFILE" ] ; then 
			awk -vF=$MOTIFCHANGESFILE -vN=$DIR/final-list.txt -vOFS="\t" '
				BEGIN{
					while (getline < F)
						Old[$2] = $1

					while (getline < N)
						for (i=1; i<=NF; i++)
							Seen[$i]
					FS="\t"
				}
				
				$1 in Old{
					for (i=1; i<=NF; i++)
						if (!($i in Seen))
						{
							Seen[$i]
							print $i
						}
				}
				' $OLDOUTPUT
		fi
	fi > $MOTIFSTOSCAN

	######## OUTPUT THE FINAL CONTROL MOTIFS LIST
	cat $DIR/final-list.txt

	if [ -f "$MOTIFCHANGESFILE" ] ; then 
		if [ $PWMINPUT -eq 1 ]; then
			awk -vF=$MOTIFCHANGESFILE -vOFS="\t" '
				BEGIN{
					while (getline < F)
						Old[$2] = $1
					FS="\t"
				}

				{c = gensub(/_[0-9]*mer$/, "", 1, $1)}
				c in Old{
					L = ""
					for (i=1; i<=NF; i++)
						L = L OFS Old[c] substr($i,length(c)+1)
					print substr(L,2)
				}
				' $OLDOUTPUT
		else
			awk -vF=$MOTIFCHANGESFILE -vOFS="\t" '
				BEGIN{
					while (getline < F)
						Old[$2] = $1
				}
				
				$1 in Old' $OLDOUTPUT
		fi
	fi

else 
	ALLOWREV=$(expr 1 - $HSTR)
	TMPDIR=$(mktemp -d /local/scratch/pk.ccm.XXXXXX)
	PWMINPUT=$(awk 'NR==1{print substr($1,1,1) == ">"}' $DIR/input.txt)
	cat $LISTFILE | while read i; do
		m=${i/ */}
		cm=$(echo $i | awk '{print substr($0, length($1) + 2)}')

		if [ $PWMINPUT -eq 1 ]; then 
			awk -vML="$cm" 'BEGIN{split(ML,A,/ /); for (i in A) X[A[i]]}; /^>/{On=substr($1,2) in X}; On' $DIR/shuffles.txt > $TMPDIR/shuffled-motifs.txt
		fi

		# immediately exclude all shuffles within $CLUSTCORR of the original motif
		if [ $EXCLUDESELF -eq 2 ]; then 
			if [ $PWMINPUT -eq 1 ]; then
				awk -vM="$m" '/^>/{On=substr($1,2) == M}; On' $DIR/shuffles.txt | motif-util -am $TMPDIR/shuffled-motifs.txt -im - -c $CLUSTCORR -r $ALLOWREV -oo 1 | cut -f4 | tr ';' '\n' | awk 'ARGIND==1{X[$1]}; ARGIND==2 && /^>/{On=!(substr($1,2) in X)}; On' - $TMPDIR/shuffled-motifs.txt > $TMPDIR/shuffled-motifs-new.txt
				mv $TMPDIR/shuffled-motifs-new.txt $TMPDIR/shuffled-motifs.txt
				cm=$(awk '/^>/{L = L " " substr($1,2)}; END{print substr(L,2)}' $TMPDIR/shuffled-motifs.txt)
			else
				cm=$(echo $cm | motif-util -a - -ic $m -c $CLUSTCORR -r $ALLOWREV | cut -f4 | tr ';' '\n' | awk -vCM="$cm" '
					{X[$1]}
					
					END{
						n=split(CM,A,/ /)
						for (i=1; i<=n; i++)
							if (!(A[i] in X))
								L = L " " A[i]
							print substr(L,2)
					}')
			fi
		fi

		# (1) sort all the potential controls for motif m by their similarity to known motifs
		# note: motifs that match /exactly/ are NO LONGER excluded
		if [ -s $DIR/exclude.txt ]; then 
			if [ $PWMINPUT -eq 1 ]; then
				motif-util -im $TMPDIR/shuffled-motifs.txt -am $DIR/exclude.txt -c -1 -r $ALLOWREV -oo 1
			else
				echo $cm | motif-util -i - -a $DIR/exclude.txt -c -1 -r $ALLOWREV 
			fi | sort -k3,3g | cut -f1
		else
			echo $cm | tr ' ' '\n' | rand-lines - 
		fi > $TMPDIR/top-order.txt

		# (2) cluster the potential controls
		# (3) choose up to $NUMMOTIFSKEEP motifs in ascending order of similarity to known motifs, but only choose one from each motif cluster
		if [ $PWMINPUT -eq 1 ]; then
			motif-util -im $TMPDIR/shuffled-motifs.txt -c $CLUSTCORR -r $ALLOWREV -oo 1 -ot 0
			rm $TMPDIR/shuffled-motifs.txt
		else
			echo $cm | motif-util -i - -c $CLUSTCORR -r $ALLOWREV -ot 0
		fi | cut -f2 | awk -F';' -vOF=$TMPDIR/top-order.txt -vMt=$m -vNMK=$NUMMOTIFSKEEP -vOFS="\t" '
			{for (i=1; i<=NF; i++) G[$i]=NR}
			
			END{
				L=Mt
				while (getline < OF)
					if (!(G[$1] in M) && X < NMK)
					{
						M[G[$1]]
						++X
						L=L OFS $1
					}
				print L
			}
		'
		rm $TMPDIR/top-order.txt
	done
	rm -rf $TMPDIR
fi

