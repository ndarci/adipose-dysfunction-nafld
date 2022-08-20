#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=12G,h_rt=04:00:00
#$ -m n
#$ -r n
#$ -o qc_rawfastq.sh.log

fastqc="/u/home/n/nikodm/bin/fastqc"
multiqc="/u/home/n/nikodm/.local/bin/multiqc"

cd ../../data/fastq/

$fastqc S*/*.fastq.gz

multiqc */*001_fastqc* --filename multiqc_rawfastq

multiqc */*trimmed_fastqc* --filename multiqc_trimmedfastq
