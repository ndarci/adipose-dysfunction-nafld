# Cross-tissue omics analysis discovers adipose genes encoding secreted proteins in obesity-related non-alcoholic fatty liver disease

This repo contains the code used for computational data analysis in Darci-Maher et al. 2023, published in *eBioMedicine* in June 2023.

The scripts used for each section of the analysis are inside the directory `scripts`, each with their own sub-directory. 

This README details the order in which to run each script, assuming that the user has access to a high-performance computing cluster and the data used in the project (not provided here).




# Data
Unfortunately, due to privacy laws, we are unable to provide the KOBS, METSIM, or UK Biobank cohort data on this repo. The GTEx, HPA, and WikiPathways data are freely available online as described in the manuscript data sharing statement.


# Figures 
The code used to generate each figure from the manuscript is available in the scripts listed here. These scripts are also integrated into the full analysis pipeline below.

## Figure 1: Study design flowchart

This figure was drawn using Adobe Illustrator and does not include data.

## Figure 2: KOBS adipose DE results for NASH

```bash
cd scripts/
Rscript kobs_de/plot_kobs_de.R
```

## Figure 3: SBC filtering results

This figure was drawn using Adobe Illustrator and does not include data.

## Figure 4: SBC gene-gene expression correlation and best subsets analysis

This script produces the heatmap in Figure 4a. The diagram in Figure 4b was drawn using Adobe Illustrator. 

```bash
Rscript kobs_de/runSBCcorrelation.R
```

## Figure 5: Preadipoctye siRNA knockdown experiment DE results
```bash
Rscript sirna_knockdown_experiment_de/generate_knockdown_plots.R 
```

## Figure 6: HepG2 recombinant protein experiment DE results
```bash
Rscript hepg2_experiment_de/plot_hepg2_de.R
```

## Figure 7: Colocalization and mendelian randomization results
```bash
Rscript mendelian_randomization_ukb/plot_vegfb_coloc_and_MR_effectsizes.R
```

## Supplementary Figure S1: WGCNA module correlation heatmap
```bash
Rscript wgcna_crosstalk/correlate_modules_phenotypes.R 
Rscript wgcna_crosstalk/arrange_sup_fig_heatmaps.R 
```

## Supplementary Figure S2: Key WGCNA module KEGG pathway enrichment results
```bash
Rscript wgcna_crosstalk/plot_pathway_enrichment_coolmodules.R 
```

## Supplementary Figure S3: KOBS adipose DE results for steatosis and fibrosis
```bash
Rscript kobs_de/plot_kobs_de.R 
```

## Supplementary Figure S4: Preadipocyte cell ORO lipid staining images and quantification
```bash
Rscript mendelian_randomization_ukb/plot_vegfb_coloc_and_MR_effectsizes.R
cd ../
```








# Full analysis pipeline

Follow the instructions below, in order, to reproduce all of the data analysis conducted in this project.

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
Plug module gene lists into WebGestalt with all expressed genes as background

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
```bash
cd scripts/kobs_liver_data_prep
```

### Get the list of sample IDs
```bash
./gen_liver_sampleIDs.sh
```

### Run samtools to remove MT reads from the original bam files
```bash
qsub runsamtools.sh
```

### Run picard to generate the RNA-seq metrics
```bash
qsub runpicard.sh
```

### Merge sample-level picard results into one file
```bash
Rscript merge_picard.R
cd ../../
```

After all of this, copy picardRNAmetrics_merged_kobs_liver_noMT.txt to `scripts/kobs_de/data_liver`

### Set up the data for DE analysis

```bash
cd scripts/kobs_de
```

Read in the GENCODE v26 .gtf file, convert it to a dataframe, and filter for our DE genes
```bash
Rscript extractAnnotations.R
```

### Consolidate covariate data from various sources and define groups to use for DE (based on liver histology measurements)
```bash
Rscript mergeCovs_defGrps.R
```

### Run gene-level DE analysis
Use LIMMA to run several DE experiments for different group definitions

Run identical analyses for adipose and liver
```bash
Rscript runLimmaVoom.R adipose
Rscript runLimmaVoom.R liver
```

## Select serum biomarker candidates (SBCs) from DE results
Note: when dowloading PANTHER results, added header manually and copy-pasted to each file
```bash
Rscript collectMetadata.R
Rscript test_de_ctm_enrich.R
```

### Generate DE plots for the paper
```bash
Rscript plot_kobs_de.R
```

### Run pairwise correlation of expression values for SBCs
```bash
Rscript runSBCcorrelation.R
```













## Run best subsets analysis
After generating the models, validate their significance with a permutation test
```bash
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
cd scripts/sirna_knockdown_experiment_de
```

### Create key file for the raw data
```bash
python3 get_seq_key.py /u/project/pajukant/nikodm/sbc_sirna/data/fastq/ raw
```

### Trim the raw data
```bash
qsub trim_reads.sh
```

### Create key file for the trimmed data
```bash
python3 get_seq_key.py /u/project/pajukant/nikodm/sbc_sirna/data/fastq/ trimmed
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

### Run DE analysis on knockdown expression data

Download adipogenesis marker genes from wikipathways: 
https://www.wikipathways.org/index.php/Pathway:WP236


Use biomaRt to convert Entrez IDs to ensembl IDs

```bash
Rscript convert_markergenes_to_ensembl.R
Rscript convert_srebf1genes_to_ensembl.R
```

### Run DE analysis between knockdown and scrambled control at each timepoint separately
```bash
Rscript run_limmavoom_pertimepoint.R
Rscript explore_de_pertimepoint_results.R
```

### Write out DE results in a table for the supplement
```bash
Rscript write_supplement_table_knockdownde.R
```

### Generate knockdown plots for the paper
```bash
Rscript generate_knockdown_plots.R
```

### Generate plot of the ORO lipid staining absorbance readings
```bash
Rscript plot_oro_absorbance.R
cd ../../
```


















## Align and quantify gene expression from the SBC recombinant protein HepG2 cell experiment

This pipeline is almost identical to the one in the sirna_knockdown_experiment_de folder, with a few small changes documented below.

### Find some likely downstream targets for CCDC80 and SOD3 for qPCR calibration testing
```bash
cd scripts/hepg2_experiment_de
Rscript collect_qPCR_targets.R
```

### Create key file for the raw data
```bash
python ../sirna_knockdown_experiment_de/get_seq_key.py /u/project/pajukant/nikodm/hepatocyte_rnaseq/data/fastq/ raw
```

### Trim the raw data
```bash
../sirna_knockdown_experiment_de/trim_reads.sh
```

### Create key file for the trimmed data
```bash
python ../sirna_knockdown_experiment_de/get_seq_key.py /u/project/pajukant/nikodm/hepatocyte_rnaseq/data/fastq/ trimmed
```

### QC the raw and trimmed data
QC script will generate a MultiQC HTML report that can be viewed in browser
```bash
qsub ../sirna_knockdown_experiment_de/qc_rawfastq.sh
```

### Map to the human genome (PASS 1)
Pass 1 maps to existing transcripts

#### Generate genome with STAR
```bash
qsub ../sirna_knockdown_experiment_de/run_genomeGenerate.pass1.sh
```

#### Align reads with STAR
```bash
qsub ../sirna_knockdown_experiment_de/map.py 1
```

### Map to the human genome (PASS 2)
Pass 2 adds new splice junctions present in the data to the reference and re-maps to that updated reference

#### Generate new genome with STAR (including new SJs)
```bash
/../sirna_knockdown_experiment_de/collectSJ.sh
qsub ../sirna_knockdown_experiment_de/run_genomeGenerate.pass2.sh
```

#### Align reads to new genome with STAR (including new SJs)
Only output uniquely mapped reads, and count MT reads
```bash
qsub ../sirna_knockdown_experiment_de/map.py 2
qsub ../sirna_knockdown_experiment_de/countMT.sh
Rscript visualizeMT_hepg2.R
```

### Generate sorted versions of the alignments for downstream
Can run these simultaneously
```bash
qsub ../sirna_knockdown_experiment_de/coordinate_sort.sh
qsub ../sirna_knockdown_experiment_de/readName_sort.sh
```

### Get some QC stats with Picard CollectRNASeqMetrics
```bash
qsub ../sirna_knockdown_experiment_de/run_picard.sample.sh
```

### Quantify expression with subread's featureCounts
```bash
qsub ../sirna_knockdown_experiment_de/run_featureCounts.sh
```

### QC the entire mapping and quantification process at once
```bash
../sirna_knockdown_experiment_de/multiqc_veryend.sh
```

### Check that the expression looks as expected
```bash
../sirna_knockdown_experiment_de/extract_uniquely_mapped.sh
Rscript expression_sanity_checks_hepg2.R
```

### Find the genes we will test for DE in the rna-seq data
```bash
Rscript collect_DE_targets.R
```

### Run the DE analyis comparing SBC protein treatment with control
```bash
Rscript run_limmavoom_proteinamt.R
Rscript explore_de_results_hepg2.R
```

### Plot the DE results
```bash
Rscript plot_hepg2_de.R
```

### Make a pretty DE table for the supplement
```bash
Rscript write_supplement_table_hepg2de.R
cd ../../
```










## Correlate SBC adipose expression with liver WGCNA modules that correlate with NAFLD traits in liver

### Run the correlation between SBC adipose expression and liver modules
```bash
cd scripts/wgcna_crosstalk
Rscript correlate_sbcs_ligands.R
```

### Investigate the implications of these correlations
```bash
Rscript explore_ligand_corrs.R
Rscript get_ligand_cor_summarystats.R
Rscript explore_SBC_livermodule_corrs.R
```

### Generate a supplementary table and plot highlighting the correlation results
```bash
Rscript write_supplement_table_sbc_livermodule_corr.R
Rscript plot_pathway_enrichment_coolmodules.R
cd ../../
```
















## Run Mendelian Randomization to test for a causal effect of cis-regulatory SNPs for adipose aware DE genes on NAFLD

### Run MAGENTA tool to check for enrichment of GWAS hits in adipose aware DE genes
To run this section:
* Download MAGENTA from the Broad institute website
* Move the MAGENTA folder inside `scripts`
* Move run_magenta.sh into the MAGENTA folder with the MATLAB scripts

```bash
cd scripts/magenta
```

```bash 
./prep_gwas_3cols.sh 
Rscript format_magenta_geneset_entrez.R
```

Here, `cd` into the MAGENTA folder (will have a unique name based on version)
```bash
qsub run_magenta.sh
```

Now, return to `scripts/magenta`
```bash
Rscript collect_magenta_results.R
cd ../../
```




### Forward MR direction: adipose IVs -> TG -> NAFLD

```bash
cd scripts/mendelian_randomization_ukb
```

### Break huge eQTL result file into small pieces for necessary genes
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
```bash
Rscript run_mrpresso.R
```

### Run MR methods 2-4 with MendelianRandomization
```bash
Rscript run_MendelianRandomization.R
```

### Make a plot of colocalized VEGFB region above MR effect size plot
```bash
Rscript plot_vegfb_coloc_and_MR_effectsizes.R
```

### Backward direction: liver IVs -> NAFLD -> TG

### Break huge eQTL result file into small pieces for necessary genes
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

No need to go any further because we observed 0 valid liver IVs at this point















## Run a regression analysis to quantify the additional variance in NAFLD explained by VEGFB compared to TG alone

### Get the genotypes for important SNPs in R compatible format
```bash
./get_vegfb_snp_textformat.sh
```

### Fit VEGFB + TG to NAFLD and compare to the fit with TG alone
```bash
Rscript build_vegfb_tg_models.sh
```

### Write supplementary tables for the VEGFB and TG model results
```bash
Rscript write_supplement_table_vegfb_tg_models.R
cd ../../
```







