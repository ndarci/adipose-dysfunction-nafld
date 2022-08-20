README for gene-level DE analysis with limma-voom

# Set up the data for DE analysis
Read in the GENCODE v26 .gtf file, convert it to a dataframe, and filter for our DE genes
```{bash}
Rscript extractAnnotations.R
```
Consolidate covariate data from various sources and define groups to use for DE (based on liver histology measurements)
```{bash}
Rscript mergeCovs_defGrps.R
```
# Run gene-level DE analysis
Use LIMMA to run several DE experiments for different group definitions
run identical analyses for adipose and liver
```{bash}
Rscript runLimmaVoom.R adipose
Rscript runLimmaVoom.R liver
```
# Overlap DE results with gene metadata to identify serum biomarker candidates (SBCs)
Note: when dowloading PANTHER results, added header manually and copy-pasted to each file
```{bash}
Rscript collectMetadata.R
Rscript test_de_ctm_enrich.R
```

# generate DE plots for the paper
```
Rscript plot_kobs_de.R
```

# Run pairwise correlation of expression values for SBCs
```{bash}
Rscript runSBCcorrelation.R
```

# run best subsets analysis to find genes that predict steatosis and NASH
After generating the models, validate their significance with a permutation test
```{bash}
Rscript runBestSubsets.R
Rscript permuteBestSubsets.R
```

# write a table for the supplement with relevant gene stats
```bash
Rscript write_supplement_table_kobsde.R
Rscript write_supplement_table_bestsubsets.R
```


--below here not used in the paper--


# look for isoQTL genes in DE results
```{bash}
Rscript compare_isoQTL_DE.R
```

# Compare sleuth results to LIMMA results
```{bash}
Rscript analyzeSleuthResults.R
```

# Check for SBCs in the T2D DE results
```{bash}
Rscript analyze_T2D_DE.R
```

# Check for associations between SBC gene expression in adipose tissue and T2D/TG
```{bash}
Rscript associateSBCs_T2D_TG.R
```

# Find enrichment of tissue types in DE results
```{bash}
Rscript runTissueEnrich.R
```
