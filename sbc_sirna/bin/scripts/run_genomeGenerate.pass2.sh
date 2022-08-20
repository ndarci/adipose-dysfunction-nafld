#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 12
#$ -l h_data=4000M,h_rt=6:00:00
#$ -v QQAPP=openmp
#$ -m n
#$ -r n
#$ -o run_genomeGenerate.pass2.sh.log

star="/u/project/pajukant/nikodm/bulk_practice_METSIM/bin/STAR-2.5.2b/bin/Linux_x86_64/STAR"

refFasta="../../data/GRCh37.primary_assembly.genome.fa"
gencodeAnno="../../data/gencode.v19.annotation.gtf"
rlMinusOne=85
genomedir="../../data/genome_for_pass2"
mkdir -p $genomedir

$star --runMode genomeGenerate \
  --runThreadN 12 \
  --genomeDir $genomedir \
  --genomeFastaFiles $refFasta \
  --sjdbGTFfile $gencodeAnno \
  --sjdbFileChrStartEnd ${genomedir}/all.SJ.out \
  --sjdbOverhang $rlMinusOne
