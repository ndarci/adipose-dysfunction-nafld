#!/bin/bash

gwasfile="/u/project/pajukant/zmiao/UKBiobank/GWAS/NAFLD_imp/NAFLD_imp_bgen.stats.gz"
outfile="../../data/nafld_gwas_ivonly.txt"

zcat $gwasfile | \
head -1 > $outfile

while IFS= read -r snp
do
	zcat $gwasfile | \
	grep $snp >> $outfile
done < ../../data/iv_snps_cond_colocalized_ldpruned_just_rsID.txt

