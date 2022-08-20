#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=16G,h_rt=1:00:00
#$ -m n
#$ -r n
#$ -o countPolyA.sh.log

. /u/local/Modules/default/init/modules.sh
module load R/4.0.2
Rscript countPolyA.R
