# adiDys siRNA knockdown experiment data pipeline
Niko Darci-Maher 2021

All R analysis done in R version 4.0.2

# Transfer data from Kuopio, Finland
raw fastq files from bulk cell line RNA-seq downloaded manually from OneDrive and transferred to hoffman2

# Trim the raw data
```bash
qsub trim_reads.sh
```

## QC the raw and trimmed data
QC script will generate a MultiQC HTML report that can be viewed in browser
```bash
cd bin/scripts/
qsub qc_rawfastq.sh
```

# Download reference genome and annotation
Using main assembly files ("CHR") for GENCODE release 19/GRCh37
```bash
./download_refs.sh
```

# Map to the human genome (PASS 1)
Pass 1 maps to existing transcripts

## Generate genome with STAR
```bash
qsub run_genomeGenerate.pass1.sh
```

## Create key file
The script needs to be pointed to the raw data directory to create a key file used throughout the downstream analysis
```bash
python3 get_seq_key.py
```

## Align reads with STAR
```bash
qsub map.py 1
```

# Map to the human genome (PASS 2)
Pass 2 adds new splice junctions present in the data to the reference and re-maps to that updated reference

## Generate new genome with STAR (including new SJs)
```bash
./collectSJ.sh
qsub run_genomeGenerate.pass2.sh
```

## Align reads to new genome with STAR (including new SJs)
Only output uniquely mapped reads
```bash
qsub map.py 2
qsub countMT.sh
Rscript visualizeMT.R
```

## Generate sorted versions of the alignments for downstream
Can run these simultaneously
```bash
qsub coordinate_sort.sh
qsub readName_sort.sh
```

# Get some QC stats with Picard CollectRNASeqMetrics and QoRT
```bash
qsub qc_qort.sh
qsub run_picard.sample.sh
```

# Quantify expression with subread's featureCounts
```bash
qsub run_featureCounts.sh
```

# QC the entire mapping and quantification process at once
```bash
./multiqc_veryend.sh
```

# Check that the expression looks as expected
```bash
./extract_uniquely_mapped.sh
Rscript expression_sanity_checks.R
```

# Download adipogenesis marker genes from wikipathways
https://www.wikipathways.org/index.php/Pathway:WP236
also download srebf1 pathway genes: https://www.wikipathways.org/index.php/Pathway:WP2706
Use biomaRt to convert Entrez IDs to ensembl IDs
```
Rscript convert_markergenes_to_ensembl.R
Rscript convert_srebf1genes_to_ensembl.R
```

# Run DE analysis between knockdown and controli at each timepoint separately
```bash
Rscript run_limmavoom_pertimepoint.R
Rscript explore_de_pertimepoint_results.R
```

# write out DE results in a nice table for supplement
```bash
Rscript write_supplement_table_knockdownde.R
```








