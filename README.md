# Cross-tissue omics analysis discovers adipose genes encoding secreted proteins in obesity-related non-alcoholic fatty liver disease

This repo contains the code used for computational data analysis in Darci-Maher et al. 2022.

The scripts used for each section of the analysis are inside the directory `scripts`, each with their own sub-directory. 

This README details the order in which to run each script, assuming that the user has access to a high-performance computing cluster and the data used in the project (not provided here).

## Data
Unfortunately, due to privacy laws concerning the sharing of personally identifiable human data, we are unable to provide the KOBS, METSIM, and UKB cohort data on this repo. The GTEx, HPA, and WikiPathways data are freely available online as described below.

## Search for adipose-liver crosstalk using WGCNA

```bash
cd scripts/wgcna_crosstalk
```

### Get a list of valid samples for the analysis
```bash
Rscript get_adipose_liver_overlap_samples.R
```

### Perform data normalization and correction per developer suggestions
```bash
Rscript prep_expression_data.R adipose
Rscript prep_expression_data.R liver
Rscript check_expr_genes_overlap.R
```

### Use builtin WGCNA functions to run QC on the cleaned data
```bash
Rscript initial_wgcna_qc_checks.R adipose
Rscript initial_wgcna_qc_checks.R liver
```

### Generate independent correlation networks in each tissue
```bash
Rscript construct_network.R adipose
Rscript construct_network.R liver
```

### Correlate modules in each network with relevant phenotypes
```bash
Rscript correlate_modules_phenotypes.R adipose
Rscript correlate_modules_phenotypes.R liver
```

### Correlate modules across adipose-liver
```bash
Rscript correlate_modules_crosstissue.R
```

### Generate heatmap correlation plots
```bash
Rscript arrange_sup_fig_heatmaps.R
```

### Investigate biological function of correlated modules

* Plug module gene lists into webgestalt with all expressed genes as background
* Get all expressed genes from limmaVoom folder

```bash
Rscript explore_xtissue_modules.R
```

### Write out a table for the supplement with important modules and gene names
```bash
Rscript write_supplement_table_wgcna.R
cd ../../
```

## Run DE analysis in KOBS adipose and liver data

### Prep the liver data

```
cd scripts/kobs_liver_noMT
```

### Get the list of sample IDs
```
./gen_liver_sampleIDs.sh
```

### Run samtools to remove MT reads from the original bam files
```
qsub runsamtools.sh
```

### Run picard to generate the RNA-seq metrics
```
qsub runpicard.sh
```

### Merge sample-level picard results into one file
```
Rscript merge_picard.R
cd ../../
```

After all of this, copy picardRNAmetrics_merged_kobs_liver_noMT.txt to `scripts/kobs_limmaVoom/data_liver`

### Set up the data for DE analysis

```
cd scripts/kobs_limmaVoom
```

Read in the GENCODE v26 .gtf file, convert it to a dataframe, and filter for our DE genes
```{bash}
Rscript extractAnnotations.R
```

### Consolidate covariate data from various sources and define groups to use for DE (based on liver histology measurements)
```{bash}
Rscript mergeCovs_defGrps.R
```

### Run gene-level DE analysis
Use LIMMA to run several DE experiments for different group definitions
Run identical analyses for adipose and liver
```{bash}
Rscript runLimmaVoom.R adipose
Rscript runLimmaVoom.R liver
```

## Select serum biomarker candidates (SBCs) from DE results
Note: when dowloading PANTHER results, added header manually and copy-pasted to each file
```{bash}
Rscript collectMetadata.R
Rscript test_de_ctm_enrich.R
```

### Generate DE plots for the paper
```
Rscript plot_kobs_de.R
```

### Run pairwise correlation of expression values for SBCs
```{bash}
Rscript runSBCcorrelation.R
```

## Run best subsets analysis
After generating the models, validate their significance with a permutation test
```{bash}
Rscript runBestSubsets.R
Rscript permuteBestSubsets.R
```

### Write supplemental tables with relevant gene stats from DE and best subsets
```bash
Rscript write_supplement_table_kobsde.R
Rscript write_supplement_table_bestsubsets.R
cd ../../
```

## Align and quantify gene expression data from siRNA knockdown experiment

```
cd scripts/sbc_sirna
```

### Trim the raw data
```bash
qsub trim_reads.sh
```

### QC the raw and trimmed data
QC script will generate a MultiQC HTML report that can be viewed in browser
```bash
qsub qc_rawfastq.sh
```

### Download reference genome and annotation
Using main assembly files ("CHR") for GENCODE release 19/GRCh37
```bash
./download_refs.sh
```

### Map to the human genome (PASS 1)
Pass 1 maps to existing transcripts

#### Generate genome with STAR
```bash
qsub run_genomeGenerate.pass1.sh
```

#### Create key file
The script needs to be pointed to the raw data directory to create a key file used throughout the downstream analysis
```bash
python3 get_seq_key.py
```

#### Align reads with STAR
```bash
qsub map.py 1
```

### Map to the human genome (PASS 2)
Pass 2 adds new splice junctions present in the data to the reference and re-maps to that updated reference

#### Generate new genome with STAR (including new SJs)
```bash
./collectSJ.sh
qsub run_genomeGenerate.pass2.sh
```

#### Align reads to new genome with STAR (including new SJs)
Only output uniquely mapped reads
```bash
qsub map.py 2
qsub countMT.sh
Rscript visualizeMT.R
```

#### Generate sorted versions of the alignments for downstream
Can run these simultaneously
```bash
qsub coordinate_sort.sh
qsub readName_sort.sh
```

### Get some QC stats with Picard CollectRNASeqMetrics and QoRT
```bash
qsub qc_qort.sh
qsub run_picard.sample.sh
```

### Quantify expression with subread's featureCounts
```bash
qsub run_featureCounts.sh
```

### QC the entire mapping and quantification process at once
```bash
./multiqc_veryend.sh
```

### Check that the expression looks as expected
```bash
./extract_uniquely_mapped.sh
Rscript expression_sanity_checks.R
```

## Run DE analysis on knockdown expression data

Download adipogenesis marker genes from wikipathways: 
https://www.wikipathways.org/index.php/Pathway:WP236


Use biomaRt to convert Entrez IDs to ensembl IDs

```
Rscript convert_markergenes_to_ensembl.R
Rscript convert_srebf1genes_to_ensembl.R
```

### Run DE analysis between knockdown and controli at each timepoint separately
```bash
Rscript run_limmavoom_pertimepoint.R
Rscript explore_de_pertimepoint_results.R
```

### Write out DE results in a nice table for supplement
```bash
Rscript write_supplement_table_knockdownde.R
cd ../../
```

### Generate knockdown plots for the paper
```bash
Rscript generate_knockdown_plots.R
```

## Correlate gene expression across adipose and liver to look for ligand-receptor interactions with the SBCs
```bash
cd scripts/wgcna_crosstalk
Rscript correlate_sbcs_ligands.R 
Rscript explore_ligand_corrs.R 
Rscript get_ligand_cor_summarystats.R
cd ../../
```

## Run Mendelian Randomization to test for a causal effect of cis-regulatory SNPs for adipose aware DE genes on NAFLD

### Run MAGENTA tool to check for enrichment of GWAS hits in adipose aware DE genes
To run this section:
* Download MAGENTA from the Broad institute website
* Move the MAGENTA folder inside `scripts`
* Move run_magenta.sh into the MAGENTA folder with the MATLAB scripts

```
cd scripts/magenta
```

```bash 
./prep_gwas_3cols.sh 
Rscript format_magenta_geneset_entrez.R
```

Here, `cd` into the MAGENTA folder (will have a unique name based on version)
```
qsub run_magenta.sh
```

Now, return to `scripts/magenta`
```
Rscript collect_magenta_results.R
```

### Forward MR direction: adipose IVs -> TG -> NAFLD

```
cd scripts/nafld_mr
```

### Break huge eqtl result file into small pieces for necessary genes
```bash
qsub generateGeneSpecificCisEQTLs.sh f
```

### Select candidate regions for IV SNPs
```bash
Rscript select_iv_snp_set.R f
```

### Get FPKMs for candidate region genes
```bash
qsub generate_iv_region_fpkm.sh f
```

### Set up LD matrix for conditional coloc
```bash
qsub runPLINKldmatrix.sh f
```

### Run conditional coloc on every candidate IV region
```bash
qsub run_cond_coloc_iv_regions.sh f
```

### Prune IV candidates based on LD and output the final list for MR
```bash
Rscript generate_iv_nafldhit_snplist.R
./calculate_ld_ivs_nafldvars.sh
Rscript ld_prune_ivs.R
```

### Get ready for the actual MR analysis by collecting and formatting data
```bash
./extract_ivs_nafldgwas.sh
Rscript prep_mr_inputs.R
```

### Run MR method 1: MR-PRESSO
```
Rscript run_mrpresso.R
```

### Run MR methods 2-4 with MendelianRandomization
```
Rscript run_MendelianRandomization.R
```

### Make a nice plot of colocalized VEGFB region above MR effect size plot
```
Rscript plot_vegfb_coloc_and_MR_effectsizes.R
```

### Backward direction: liver IVs -> NAFLD -> TG

### Break huge eqtl result file into small pieces for necessary genes
```bash
qsub generateGeneSpecificCisEQTLs.sh b
```

### Select candidate regions for IV SNPs
```bash
Rscript select_iv_snp_set.R b
```

### Get FPKMs for candidate region genes
```bash
qsub generate_iv_region_fpkm.sh b
```

### Set up LD matrix for conditional coloc
```bash
qsub runPLINKldmatrix.sh b
```

### Run conditional coloc on every candidate IV region
```bash
qsub run_cond_coloc_iv_regions.sh b
```

















