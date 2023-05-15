#!/bin/bash

plink="/u/home/n/nikodm/bin/PLINK_v1.9/plink"

$plink \
	--bfile /u/project/pajukant/nikodm/kobs_prs/data/kobs_geno_rsID \
	--snp rs2845885 \
	--recode A \
	--out ../../data/vegfb_mr_snp_genotypes_text

$plink \
	--bfile /u/project/pajukant/nikodm/kobs_prs/data/kobs_geno_rsID \
	--snps rs2845885, rs467812, rs6442038, rs6465120, rs6502821, rs6714157 \
	--recode A \
	--out ../../data/all_iv_snp_genotypes_text