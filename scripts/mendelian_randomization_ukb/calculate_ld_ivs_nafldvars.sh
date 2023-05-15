#!/bin/bash

# point to PLINK tool
plink="/u/home/n/nikodm/bin/PLINK_v1.9/plink"
# point to imputed METSIM genotype data
infile="/u/project/pajukant/pajukant/datashare/METSIM_10K/imputed/hrc1.1_phased_rsq0.3/metsim.10K_imputed.hrc1.1.rsq0.3.hwe1e-6.maf0.005.grm0.0625"

snpfile="../../data/iv_candidate_snps_plus_nafld_gwashits.txt"
outfile="../../data/ld_matrix_ivs_nafldhits"

# generate an LD matrix (r2 values) with plink
$plink \
--bfile $infile \
--extract $snpfile \
--write-snplist \
--r2 square \
--out $outfile

