#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 1
#$ -l h_data=24G,h_rt=4:00:00
#$ -m n
#$ -r n
#$ -t 1-1
#$ -o coordinate_sort.sh.log.$TASK_ID

. /u/local/Modules/default/init/modules.sh
module load java/jre-1.8.0_281

# cd ../../

pass2_dir="../../data/aligned/pass2"

dsfile="../../data/sample_fastq.txt"

picard="/u/project/pajukant/nikodm/bulk_practice_METSIM/bin/picard.jar"

n_datasets=$( wc -l $dsfile |  awk '{print $1}' )
ntasks=1
# ntasks=24
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
        sample=$( awk -v I=$i 'NR == I {print $1}' $dsfile )
        # instr=$( awk -v I=$i 'NR == I {print $2}' $dsfile )
        # runn=$( awk -v I=$i 'NR == I {print $3}' $dsfile )
        # flowcell=$( awk -v I=$i 'NR == I {print $4}' $dsfile )
        # lane=$( awk -v I=$i 'NR == I {print $5}' $dsfile )

        ibam="${pass2_dir}/${sample}/${sample}_Aligned.out.bam"
        obam="${pass2_dir}/${sample}/${sample}_Aligned.out.coordinateSorted.bam"

        echo $sample

        java -Xmx16g -jar $picard SortSam \
                I=$ibam \
                O=$obam \
                SORT_ORDER=coordinate
done
