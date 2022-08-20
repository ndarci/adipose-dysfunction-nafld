#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=8G,h_rt=01:00:00
#$ -m n
#$ -r n
#$ -o runPLINKldmatrix.sh.log

direc=$1

# forward
if [[ "$direc" == "f" ]]; then
	genelist="../../data/iv_candidate_genes_ensembl.txt"
	ovtag="TG_noNAFLD_"
# backward
elif [[ "$direc" == "b" ]]; then
	genelist="../../data/iv_candidate_genes_ensembl_backward.txt"
	ovtag="NAFLD_noTG_"
else
    echo "Invalid input arguments"
fi

# point to PLINK tool
plink="/u/home/n/nikodm/bin/PLINK_v1.9/plink"
# point to imputed METSIM genotype data
infile="/u/project/pajukant/pajukant/datashare/METSIM_10K/imputed/hrc1.1_phased_rsq0.3/metsim.10K_imputed.hrc1.1.rsq0.3.hwe1e-6.maf0.005.grm0.0625"

# loop thru candidate IV regions
while IFS= read -r gene
do
	snpfile="../../data/ld_matrix/${ovtag}${gene}_SNPoverlap_SNPlist.txt"
	outfile="../../data/ld_matrix/ld_matrix_${gene}"
	# generate an LD matrix (r, not r^2) with plink
	$plink \
	--bfile $infile \
	--extract $snpfile \
	--write-snplist \
	--r square \
	--out $outfile
done < $genelist
