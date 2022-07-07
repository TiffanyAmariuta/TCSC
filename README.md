# TCSC (Tissue co-regulation score regression)

TCSC is a statistical genetics method to identify causal tissues in diseases and complex traits. We leverage TWAS and GWAS summary statistics while explicitly modeling the genetic co-regulation of genes across tissues.

## Using TCSC

### Dependencies

R v.3.6.1

### Input Data 

1. Signed GWAS summary statistics (formatted for S-LDSC; Finucane 2015 Nat Genet)
```
sumstats/
```

2. Gene expression prediction models
GTEx models for tissues downsampled or meta-analyzed such that N = 320 for all tissues
```
weights/
```

3. TWAS summary statistics computed using FUSION (Gusev 2016 Nat Genet; http://gusevlab.org/projects/fusion/) 
```
twas_statistics/
```

### Running TCSC

1. Using gene expression prediction models for cis-heritable genes (GCTA p < 0.01, GCTA h2 estimation > 0, cross validation R2 > 0), construct gene and tissue co-regulation scores. 

```
analysis/MakeCoregulationScores.R
```

2. Prepare input files for the analysis of all traits. 
```
analysis/TCSC_setup.R
```

3. Run TCSC regression for each trait. 
```
analysis/Run_TCSC.R $trait
```


