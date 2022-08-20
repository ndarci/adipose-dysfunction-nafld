#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=16G,h_rt=12:00:00
#$ -m n
#$ -r n
#$ -o permuteBestSubsets.sh.log

. /u/local/Modules/default/init/modules.sh
module load R/4.0.2
Rscript permuteBestSubsets.R
~
