#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 4
#$ -l h_data=8G,h_rt=16:00:00
#$ -t 1-288:30
#$ -m n
#$ -r n
#$ -o runpicard.sh.log.$TASK_ID

# import java module (used for picard)
. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77

# define the samples this job will process
start=$SGE_TASK_ID
step=$SGE_TASK_STEPSIZE
stop=$(( $start + $step - 1))
last=$SGE_TASK_LAST

if [ "$stop" -gt "$last" ]
then
    stop=$last
fi

# point to relevant programs and files
picard="/u/project/pajukant/malvarez/KOBS_liver_align/programs/picard.jar"
dsfile="sampleIDs.txt"
refflat="refFlat_GRCh37.txt"
strand="FIRST_READ_TRANSCRIPTION_STRAND"

for i in $( seq $start $stop ); do
	# get sample ID
    sample=$( awk -v I=$i 'NR == I {print $1}' $dsfile )
    # point to ins and outs
    infile="kobs_liver_bam_noMT/${sample}/Aligned.sortedByCoord.out.noMT.bam"
    outfile="kobs_liver_bam_noMT/${sample}/picard_RNAmetrics_noMT.txt"
    # call picard
    java -Xmx8g -jar $picard CollectRnaSeqMetrics \
        I=$infile \
        REF_FLAT=$refflat \
        STRAND=$strand \
        O=$outfile
done
