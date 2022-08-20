# Run Mendelian Randomization to test for a causal effect of cis-regulatory SNPs for adipose aware DE genes on NAFLD

# Forward direction: adipose IVs -> TG -> NAFLD

## Break huge eqtl result file into small pieces for necessary genes
```bash
qsub generateGeneSpecificCisEQTLs.sh f
```

## Select candidate regions for IV SNPs
```bash
Rscript select_iv_snp_set.R f
```

## Get FPKMs for candidate region genes
```bash
qsub generate_iv_region_fpkm.sh f
```

## Set up LD matrix for conditional coloc
```bash
qsub runPLINKldmatrix.sh f
```

## Run conditional coloc on every candidate IV region
```bash
qsub run_cond_coloc_iv_regions.sh f
```

## Prune IV candidates based on LD and output the final list for MR
```bash
Rscript generate_iv_nafldhit_snplist.R
./calculate_ld_ivs_nafldvars.sh
Rscript ld_prune_ivs.R
```

## Get ready for the actual MR analysis by collecting and formatting data
```bash
./extract_ivs_nafldgwas.sh
Rscript prep_mr_inputs.R
```

## Run MR method 1: MR-PRESSO
```
Rscript run_mrpresso.R
```

## run MR methods 2-4 with MendelianRandomization
```
Rscript run_MendelianRandomization.R
```

# make a nice plot of colocalized VEGFB region above MR effect size plot
```
Rscript plot_vegfb_coloc_and_MR_effectsizes.R
```

# alternative thing: choose GWAS snps instead of cis-eQTLs as IVs
Substitute the following scripts for their analogues in the above pipeline, otherwise run the normal scripts
```bash
generate_iv_nafldhit_snplist_choosegwas.R
ld_prune_ivs_choosegwas.R
prep_mr_inputs_choosegwas.R
run_mrpresso_choosegwas.R
run_MendelianRandomization_choosegwas.R
```


# Backward direction: liver IVs -> NAFLD -> TG

## Break huge eqtl result file into small pieces for necessary genes
```bash
qsub generateGeneSpecificCisEQTLs.sh b
```

## Select candidate regions for IV SNPs
```bash
Rscript select_iv_snp_set.R b
```

## Get FPKMs for candidate region genes
```bash
qsub generate_iv_region_fpkm.sh b
```

## Set up LD matrix for conditional coloc
```bash
qsub runPLINKldmatrix.sh b
```



## Run conditional coloc on every candidate IV region
```bash
qsub run_cond_coloc_iv_regions.sh b
```











