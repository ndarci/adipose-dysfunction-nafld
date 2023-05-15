#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 12
#$ -l h_data=6000M,h_rt=10:00:00
#$ -m n
#$ -r n
#$ -o run_featureCounts.sh.log

featurecounts="/u/project/pajukant/nikodm/bulk_practice_METSIM/bin/subread-1.6.2-Linux-x86_64/bin/featureCounts"

cd ../../

touch data/merged_sample_paths.txt
ls data/aligned/pass2/*/*.readNameSorted.bam > data/merged_sample_paths.txt

samplefile="data/merged_sample_paths.txt"
gtf="/u/project/pajukant/nikodm/sbc_sirna/data/gencode.v19.annotation.gtf"
outdir="data/featureCounts/"
out="${outdir}/gencode19.featureCounts.txt"
mkdir -p $outdir

$featurecounts \
	-p \
	-T 12 \
	-B \
	-a $gtf \
	--extraAttributes gene_name,gene_type \
	-o $out \
	-s 0 \
	$( cat $samplefile )

