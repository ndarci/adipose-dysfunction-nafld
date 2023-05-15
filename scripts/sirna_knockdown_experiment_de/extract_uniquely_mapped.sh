#!/bin/bash

outfile=../../data/de_cov.txt
echo "sample unique_map_pct" > $outfile

{
	read
	while IFS='	' read -r sample fastq
        do
		bamdir=../../data/aligned/pass2/${sample}/
		ump=$( grep "Uniquely mapped reads %" ${bamdir}/${sample}_Log.final.out | cut -f 2 | cut -d '%' -f 1)
		echo $sample $ump >> $outfile
	done
} < ../../data/sample_fastq.txt


