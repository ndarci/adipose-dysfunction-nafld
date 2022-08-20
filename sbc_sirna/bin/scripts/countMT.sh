#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=8G,h_rt=2:00:00
#$ -m n
#$ -r n
#$ -o countMT.sh.log

. /u/local/Modules/default/init/modules.sh
module load samtools

outfile=../../data/aligned/pass2/mtreadproportions.txt
echo "sample mtreads totreads" > $outfile

{
	read
	while IFS='	' read -r sample fastq
	do
		bamdir=../../data/aligned/pass2/${sample}/
		mtreads=$( samtools view ${bamdir}/${sample}_Aligned.out.bam | grep chrM | wc -l )
		totreads=$( grep "Uniquely mapped reads number" ${bamdir}/${sample}_Log.final.out | cut -f 2 )
		echo ${sample} ${mtreads} ${totreads} >> $outfile
	done
} < ../../data/sample_fastq.txt



