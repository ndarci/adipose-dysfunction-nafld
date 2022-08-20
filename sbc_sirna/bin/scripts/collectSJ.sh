#!/bin/bash
#
# collect only unannotated
# Remove MT splice reads
# Remove non-canonical junctions
# Remove junctions supported by less than 3 reads

mkdir -p ../../data/genome_for_pass2

awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} 
{if($6 == 0 && $1 != "MT" && $5>0 && $7 >= 3) {print $1,$2,$3,strChar[$4]}}' \
	../../data/aligned/pass1/*/*_SJ.out.tab \
	> ../../data/genome_for_pass2/all.SJ.out
