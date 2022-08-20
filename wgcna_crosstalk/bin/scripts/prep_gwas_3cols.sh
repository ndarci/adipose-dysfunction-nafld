#!/bin/bash

# cut each relevant GWAS down to 3-column version

# iterate over each file path
for gwas in /u/project/pajukant/zmiao/UKBiobank/GWAS/TG/TG_bgen.stats.gz \
			/u/project/pajukant/zmiao/UKBiobank/GWAS/NAFLD_imp/NAFLD_imp_bgen.stats.gz \
			/u/project/pajukant/zmiao/UKBiobank/GWAS/ALT/ALT_bgen.stats.gz \
			/u/project/pajukant/zmiao/UKBiobank/GWAS/GGT/GGT_bgen.stats.gz
do
	# extract phenotype name from file path
	pheno=$(echo ${gwas} | awk -F "/" '{print $8}')

	# take the relevant columns and chop off header
	zcat $gwas | \
	cut -f 2,3,14 | \
	sed 1d > \
	../MAGENTA_software_package_vs2_July2011/${pheno}_gwas_UKB_3cols
done
