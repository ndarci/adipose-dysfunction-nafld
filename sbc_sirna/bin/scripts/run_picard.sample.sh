#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=24000M,h_rt=2:00:00
#$ -t 1-1
#$ -m n
#$ -r n
#$ -o run_picard.sample.sh.log.$TASK_ID

. /u/local/Modules/default/init/modules.sh
module load java/jre-1.8.0_281

cd ../../

picard="/u/project/pajukant/nikodm/bulk_practice_METSIM/bin/picard.jar"

samplefile="data/sample_fastq.txt"

pass2dir="data/aligned/pass2"

n_datasets=$( wc -l $samplefile | awk '{print $1}')
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
        i=$(( i + 1 ))
        sample=$( awk -v I=$i 'NR == I {print $1}' $samplefile )

        outdir="data/picard/${sample}"
		mkdir -p $outdir

        echo $sample

        ibam="${pass2dir}/${sample}/${sample}_Aligned.out.coordinateSorted.bam"
        out="${outdir}/${sample}_picard.RNA_Metrics"

        java -Xmx16g -jar $picard CollectRnaSeqMetrics \
        I=$ibam \
        REF_FLAT=data/refFlat.txt.gz \
        STRAND=NONE \
        O=$out
done
