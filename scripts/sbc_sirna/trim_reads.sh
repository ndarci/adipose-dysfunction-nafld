#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=1G,h_vmem=8G,h_rt=04:00:00,highp
#$ -m n
#$ -r n
#$ -o trim_reads.sh.log

. /u/local/Modules/default/init/modules.sh
module load python

cutadapt="/u/home/n/nikodm/.local/bin/cutadapt"

polya='AAAAAAAAAA'
r1adapt='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
contam=${polya}${r1adapt}

{
	read
	while IFS='	' read -r sample fastq
	do
		fastqout=$( echo $fastq | sed 's/.fastq.gz/_trimmed.fastq.gz/g' )
		$cutadapt \
			--cores 8 \
			--adapter "${contam};min_overlap=10" \
			--minimum-length 20 \
			--output $fastqout \
			$fastq
	done
} < ../../data/sample_fastq.txt

