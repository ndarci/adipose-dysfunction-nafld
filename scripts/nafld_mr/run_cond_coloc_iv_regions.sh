#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=32G,h_rt=8:00:00
#$ -m n
#$ -r n
#$ -o run_cond_coloc_iv_regions.sh.log

. /u/local/Modules/default/init/modules.sh
module load R/4.0.2
Rscript run_cond_coloc_iv_regions.R
