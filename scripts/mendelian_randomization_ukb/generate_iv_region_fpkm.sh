#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=4G,h_rt=00:30:00
#$ -m n
#$ -r n
#$ -o generate_iv_region_fpkm.sh.log

direc=$1

# forward
if [[ "$direc" == "f" ]]; then
        ivgenes="../../data/iv_candidate_genes_ensembl.txt"
        fpkm="/u/project/pajukant/dzpan29/Juno/KOBS/KOBS_Invrs_FPKM_Baseline_All_peer25_invnorm.csv"
        outfile="../../data/KOBS_Invrs_FPKM_Baseline_All_peer25_invnorm_IVregions.csv"
# backward
elif [[ "$direc" == "b" ]]; then
        ivgenes="../../data/iv_candidate_genes_ensembl_backward.txt"
        fpkm="/u/project/pajukant/seunglee/KOBS/liver/bulkRNA/featureCounts/KOBS_liver_hg19_FPKM_nonzero_INT_correctedID.csv"
        outfile="../../data/KOBS_liver_hg19_FPKM_nonzero_INT_correctedID_IVregions_backward.csv"
else
        echo "Invalid input arguments"
fi

# get header of fpkm file
head -1 $fpkm > $outfile
while IFS= read -r gene
do
        # extract every gene entry into subset file
        grep $gene $fpkm >> $outfile
done < $ivgenes
