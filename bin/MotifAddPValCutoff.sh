#!/bin/sh

SUFFIX=""
KEEPWITHMAX=0

if [ $# -eq 0 ]; then
		echo >&2 "USAGE: $0 [OPTIONS] <PWM FILE>
 -k     k-mer equivalent length (may be list) [required]
 -s     Suffix to add to motif names (run through printf with kmer) [default: $SUFFIX]
 -m     Keep motifs with information content less than -k with the maximum score (1 to keep only first >= MAX; 2 to keep all) [default: 0]
"
	exit 1
fi

while getopts k:s:m: o
do		case "$o" in
	k)		KMERS="$OPTARG";;
	s)		SUFFIX="$OPTARG";;
	m)		KEEPWITHMAX="$OPTARG";;
	[?])    echo >&2 "ERROR: command line parameter not recognized."; exit 1;;
	esac
done
shift $(expr $OPTIND - 1)
MOTIFFILE=$1

if [ -z "$KMERS" ]; then
	echo >&2 "ERROR: -k is required."; exit 1
fi	

TMPFILE1=$(mktemp /tmp/pk.mapvc.XXXXXX)
TMPFILE2=$(mktemp /tmp/pk.mapvc.XXXXXX)

for m in $(seq 1 $(awk '/^>/{++X}; END{print X}' $MOTIFFILE)); do 
	awk -vM=$m '/>/{MX++}; MX==M' $MOTIFFILE > $TMPFILE1
	
	motif-util -im $TMPFILE1 -om l | awk 'NR>2{for (i=1; i<=4; i++) X[i] = X[i] " " $(i+1)}; END{for (i=1; i<=4; i++) print substr(X[i],2)}' > $TMPFILE2

	if [ $KEEPWITHMAX != 0 ]; then 
		MAX=$(awk '{for (i=1; i<=NF; i++) if (!(i in X) || $i > X[i]) X[i]=$i}; END{for (i=1; i in X; i++) M+=X[i]; print M}' $TMPFILE2)
	fi
	
	for K in $(echo $KMERS); do 
		TFMpvalue-pv2sc -a 0.25 -t 0.25 -c 0.25 -g 0.25 -w  -m $TMPFILE2 -p $(awk -vK=$K 'BEGIN{print 4**(-K)}') | awk -vSuf=$(printf $SUFFIX $K) -vX="$MAX" '
			ARGIND==1 && FNR==1 && $3 != "0"{X = $2}
			X == ""{exit}
			ARGIND==2 && FNR==1{$0 = $1 Suf "\t" X}
			ARGIND==2
		' - $TMPFILE1
	done | if [ $KEEPWITHMAX = 1 ]; then 
		awk -F"\t" -vX="$MAX" '/^>/{On = ($2 != X) || !($2 in L); L[$2]}; On'
	else
		cat
	fi
done

rm -f $TMPFILE1 $TMPFILE2

