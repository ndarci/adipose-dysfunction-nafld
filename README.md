# Cross-tissue omics analysis discovers adipose genes encoding secreted proteins in obesity-related non-alcoholic fatty liver disease

This repo contains the code used for computational data analysis in Darci-Maher et al. 2022.

The scripts used for each section of the analysis are inside the directory `scripts`, each with their own sub-directory. 

This README details the order in which to run each script, assuming that the user has access to a high-performance computing cluster and the data used in the project (not provided here).

## Data
Unfortunately, due to privacy laws concerning the sharing of personally identifiable human data, we are unable to provide the KOBS cohort data on this repo.

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

### Overlap DE results with gene metadata to identify serum biomarker candidates (SBCs)
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

### Run best subsets analysis to find genes that predict steatosis and NASH
After generating the models, validate their significance with a permutation test
```{bash}
Rscript runBestSubsets.R
Rscript permuteBestSubsets.R
```

### Write a table for the supplement with relevant gene stats
```bash
Rscript write_supplement_table_kobsde.R
Rscript write_supplement_table_bestsubsets.R
```



















## Correlate gene expression across adipose and liver to look for ligand-receptor interactions with the SBCs

```bash
Rscript correlate_sbcs_ligands.R 
Rscript explore_ligand_corrs.R 
Rscript get_ligand_cor_summarystats.R
```

### Run MAGENTA tool to check for enrichment of gwas hits in correlated liver genes
```bash 
./prep_gwas_3cols.sh 
Rscript format_magenta_geneset_entrez.R
qsub run_magenta.sh
Rscript collect_magenta_results.R
```
