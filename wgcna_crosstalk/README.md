# Searching for adipose-liver crosstalk using WGCNA

The WGCNA arm of the adipose dysfunction project

# Get a list of valid samples for the analysis
```bash
cd bin/scripts/
Rscript get_adipose_liver_overlap_samples.R
```

# Perform data normalization and correction per developer suggestions
```bash
Rscript prep_expression_data.R adipose
Rscript prep_expression_data.R liver
Rscript check_expr_genes_overlap.R
```

# Use builtin WGCNA functions to run QC on the cleaned data
```bash
Rscript initial_wgcna_qc_checks.R adipose
Rscript initial_wgcna_qc_checks.R liver
```

# Generate independent correlation networks in each tissue
```bash
Rscript construct_network.R adipose
Rscript construct_network.R liver
```

# Correlate modules in each network with relevant phenotypes
```bash
Rscript correlate_modules_phenotypes.R adipose
Rscript correlate_modules_phenotypes.R liver
```

# Correlate modules across adipose-liver
```bash
Rscript correlate_modules_crosstissue.R
```

# Investigate biological function of correlated modules

* Plug module gene lists into webgestalt with all expressed genes as background
* Get all expressed genes from limmaVoom folder
* plug important module gene lists into panther.db

```bash
Rscript overlap_panther_xtissuemods.R A_lightyellow
Rscript overlap_panther_xtissuemods.R L_saddlebrown
```

```bash
Rscript explore_xtissue_modules.R
```


```bash
Rscript convert_ensembl_to_entrez.R
Rscript run_functional_enrichment.R
```

# write out a table for the supplement with important modules and gene names
```bash
Rscript write_supplement_table_wgcna.R
```


# Side thing: directly correlate gene expression across adipose and liver to look for ligand-receptor interactions with the SBCs

```bash
Rscript correlate_sbcs_ligands.R 
Rscript explore_ligand_corrs.R 
Rscript get_ligand_cor_summarystats.R
```

# Run MAGENTA tool to check for enrichment of gwas hits in correlated liver genes
```bash 
./prep_gwas_3cols.sh 
Rscript format_magenta_geneset_entrez.R
qsub run_magenta.sh
Rscript collect_magenta_results.R
```







