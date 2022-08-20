#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 12
#$ -l h_data=2000M,h_rt=4:00:00
#$ -m n
#$ -r n
#$ -t 1-1
#$ -o readName_sort.sh.log.$TASK_ID

. /u/local/Modules/default/init/modules.sh
module load java/jre-1.8.0_281

cd ../../

samplefile="data/sample_fastq.txt"

picard="/u/project/pajukant/nikodm/bulk_practice_METSIM/bin/picard.jar"

n_datasets=$( wc -l $samplefile |  awk '{print $1}' )
ntasks=1
#ntasks=24
step=$(( n_datasets / ntasks ))
task=$(( $SGE_TASK_ID ))
start=$(( ($task - 1) * $step ))
stop=$(( $start + $step - 1 ))
if [[ $SGE_TASK_ID == $ntasks ]]; then
        stop=$(( $n_datasets - 1 ))
fi

for i in $( seq $start $stop ); do
        # Since header, add 1
        i=$((i + 1))
        sample=$( awk -v I=$i 'NR == I {print $1}' $samplefile )
	bamdir="data/aligned/pass2/${sample}"
        
	ibam=${bamdir}/${sample}_Aligned.out.bam
        obam=${bamdir}/${sample}_Aligned.out.readNameSorted.bam

        echo $sample

        java -Xmx16g -jar $picard SortSam \
                I=$ibam \
                O=$obam \
                SORT_ORDER=queryname
done
