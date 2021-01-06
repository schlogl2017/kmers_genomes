#!/bin/bash
#Counting k-mers for different k on phiX
echo "k,unique,distinct,total"
for k in {1..15}; do
	kat hist -o phiX_$k.hist -m $k phiX.fasta >/dev/null 2>&1
	egrep -v '^#' phiX_$k.hist|awk '{if ($1==1) u=$2; d+=$2; t+=$2*$1;}\
	END{print "'$k',"u","d","t}'
done
