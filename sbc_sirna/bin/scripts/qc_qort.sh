#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=16G,h_rt=02:00:00
#$ -m n
#$ -r n
#$ -o qc_qort.sh.log

. /u/local/Modules/default/init/modules.sh
module load java/jre-1.8.0_281

qort="../QoRTs-STABLE.jar"

fastq="../../data/fastq/S10/Baseline_CCDC80_1_rep3_S10_R1_001.fastq.gz"
bam="../../data/aligned/pass2/S10/S10_Aligned.out.bam"

java -jar $qort QC \
	--generatePlots \
	--rawfastq $fastq \
	--singleEnded \
	$bam \
	../../data/gencode.v19.annotation.gtf \
	../../data/QC_QoRT/
