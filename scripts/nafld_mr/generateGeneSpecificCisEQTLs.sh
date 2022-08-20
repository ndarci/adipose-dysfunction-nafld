#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=4G,h_rt=06:00:00
#$ -m n
#$ -r n
#$ -o generateGeneSpecificCisEQTLs.sh.log

direc=$1

# forward direction
if [[ "$direc" == "f" ]]; then
	genelist="/u/project/pajukant/nikodm/kobs_limmaVoom/data/DE_genes_filter_except_secreted_all.txt"
	cis="/u/project/pajukant/nikodm/kobs_coloc/data/KOBS_eQTL_BaselineAll_Adipose_cis.txt"
	outprefix="../../data/ciseqtl_gene/KOBS_eQTL_BaselineAll_Adipose_cis_"
# backward direction
elif [[ "$direc" == "b" ]]; then
	genelist="/u/project/pajukant/nikodm/kobs_limmaVoom/data_liver/DE_genes_liver_not_adipose_all.txt"
	cis="/u/project/pajukant/seunglee/KOBS/liver/bulkRNA/matrix.eqtl/KOBS_liver_hg19_covar_peer10.cis.txt"
	outprefix="../../data/ciseqtl_gene/KOBS_eQTL_BaselineAll_Liver_cis_"
else
	echo "Invalid input arguments"
fi

while IFS= read -r gene
do
	# get ciseQTL header line
	head -1 $cis > ${outprefix}${gene}.txt
	# get every entry for this gene
	grep $gene $cis >> ${outprefix}${gene}.txt
done < $genelist


