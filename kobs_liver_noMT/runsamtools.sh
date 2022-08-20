#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=8G,h_rt=8:00:00
#$ -t 1-288:15
#$ -m n
#$ -r n
#$ -o runsamtools.sh.log.$TASK_ID

# load samtools module
. /u/local/Modules/default/init/modules.sh
module load samtools

# find the range of samples that this job will process
start=$SGE_TASK_ID
step=$SGE_TASK_STEPSIZE
stop=$(( $start + $step - 1))
last=$SGE_TASK_LAST

# account for running over the last line of the sample ID file
if [ "$stop" -gt "$last" ]
then
    stop=$last
fi

# link the file with sample IDs
dsfile="sampleIDs.txt"

# link the source of the .bam files
bamdirec="/u/project/pajukant/malvarez/KOBS_liver_align/data/processed/STAR_output/gencode26.sample"

# link the output directory
outdir="./kobs_liver_bam_noMT/"

# run samtools on each sample to produce new bams with MT reads removed
for i in $( seq $start $stop ); do
	# get sample ID
    sample=$( awk -v I=$i 'NR == I {print $1}' $dsfile )
    # get corresponding input bam
    infile="${bamdirec}/${sample}/Aligned.sortedByCoord.out.bam"
    # create an output directory and file name for this sample
    mkdir -p "$outdir"
    mkdir -p "${outdir}/${sample}"
    outfile="${outdir}/${sample}/Aligned.sortedByCoord.out.noMT.bam"
    # call samtools
	samtools idxstats $infile \
	| cut -f 1 \
	| grep -v MT \
	| xargs samtools view -b $infile \
	> $outfile
done

