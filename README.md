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

For GTEx v8 European samples, these can be downloaded from http://gusevlab.org/projects/fusion/. In our primary analysis, for tissues with sample size greater than 320, we subsampled to 320 individuals per tissue; for groups of tissues with sample size less than 320 and highly correlation of marginal eQTL effect sizes, we meta-analyzed gene expression models such that each meta-tissue has an effective sample size of 320 individuals. 

While we cannot provide gene expression prediction models on Github due to file quantity, we provide a list of cis heritable genes per tissue in each of the above analyses. 
```
weights/
```

3. TWAS summary statistics computed using FUSION (Gusev 2016 Nat Genet; http://gusevlab.org/projects/fusion/)

While we cannot provide TWAS summary statistics for all of our analyses on Github due to file quantity, we provide one representative set of TWAS summary statistics for BMI.
```
twas_statistics/
```

### Running TCSC

1. Using gene expression prediction models for cis-heritable genes (GCTA p < 0.01, GCTA h2 estimation > 0, cross validation R2 > 0), impute gene expression into reference panel, e.g. 489 European individuals in 1KG. 

We use FUSION (http://gusevlab.org/projects/fusion/) to first compute score files for each gene expression model. Then we use plink, as done in FUSION, to impute gene expression into 1KG individuals.
```
predicted_expression/
```

2. Construct gene and tissue co-regulation scores using squared correlations of predicted expression and applying bias correction (using gene model accuracy) when tissue t = t'. 

```
analysis/MakeCoregulationScores.R
```

3. Prepare input files for the analysis of all traits. 
```
analysis/TCSC_setup.R
```

4. Run TCSC regression for each trait. 
```
analysis/Run_TCSC.R $trait
```


