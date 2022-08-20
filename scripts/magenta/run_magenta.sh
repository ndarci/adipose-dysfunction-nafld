#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=8G,h_rt=12:00:00,highp
#$ -m n
#$ -r n
#$ -o run_magenta.sh.log

. /u/local/Modules/default/init/modules.sh
module load matlab

# test for gwas hit enrichment in every combination of interesting gwases and gene sets
for gwas in TG NAFLD_imp ALT GGT
do
	for geneset in SFRP2_cor all_SBC_cor adi_aware_DE
	do
		matlab -nodisplay -nodesktop -nosplash -nojvm -r \
			"\
			GWAS_SNP_score_file_name='${gwas}_gwas_UKB_3cols';\
			GeneSet_db_file_name='${geneset}_geneset';\
			exp_label='${gwas}_gwas_${geneset}_geneset';\
			Run_MAGENTA_vs2_July_2011\
			"
	done
done


