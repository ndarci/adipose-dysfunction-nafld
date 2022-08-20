README for kobs_liver_noMT, part of adipose dysfunction project that generates Picard RNA-seq metrics for the KOBS liver data (used downstream as DE covariates)

# get the list of sample IDs
```
./gen_liver_sampleIDs.sh
```

# run samtools to remove MT reads from the original bam files
```
qsub runsamtools.sh
```

# run picard to generate the RNA-seq metrics
```
qsub runpicard.sh
```

# merge sample-level picard results into one file
```
Rscript merge_picard.R
```

After all of this, copy picardRNAmetrics_merged_kobs_liver_noMT.txt to ~/project-pajukant/kobs_limmaVoom/data_liver
