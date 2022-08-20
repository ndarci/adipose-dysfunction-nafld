#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=36G,h_rt=00:15:00,highp
#$ -m n
#$ -r n
#$ -o construct_network.sh.log

. /u/local/Modules/default/init/modules.sh
module load R/4.0.2

tissue=$1

Rscript construct_network.R $tissue
