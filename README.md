# TCSC (Tissue co-regulation score regression)

TCSC is a statistical genetics method to identify causal tissues in diseases and complex traits. We leverage TWAS and GWAS summary statistics while explicitly modeling the genetic co-regulation of genes across tissues.

## Using TCSC

### Dependencies

R v.3.6.1

TCSC uses TWAS statistics computed from the FUSION software (http://gusevlab.org/projects/fusion/). We have decided that TCSC as an R package would not be compatible with future maintenance and bug fixes made to the TWAS software, which exists as a Github repo. 

### Quick Start: Apply TCSC to new GWAS summary statistics and default gene expression data from GTEx
1. Preprocess your GWAS summary statistics such that they are in LDSC format

See directions at https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format
Refer to your genome-wide summary statistics using the path and name of your trait as suggested below
```
your_genomewide_sumstats=your_path/${trait}.sumstats
```

2. Clone the TCSC github repository
```
git clone https://github.com/TiffanyAmariuta/TCSC
```

3. Clone FUSION github repository

```
git clone https://github.com/gusevlab/fusion_twas
```

4. Download and unzip LD reference files for 1000 Genomes European population
```
wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
tar xjvf LDREF.tar.bz2
```

5. Download and unzip TCSC gene expression weight files
```
mkdir TCSC_weight_files
cd TCSC_weight_files
for tissue in `cat ../TCSC/analysis/TissuesA.txt`
do
wget https://storage.googleapis.com/broad-alkesgroup-public/TCSC/GeneExpressionModels/eQTL_samplesize_320/${tissue}.tar.gz
tar zxvf ${tissue}.tar.gz #creates directory TCSC_weight_files/weights/v8_320EUR/META_${tissue}/
done

cd ../
mkdir -p TCSC_weight_files/v8_allEUR/
cd TCSC_weight_files/v8_allEUR/
for tissue in `cat ../TCSC/analysis/TissuesB.txt`
do
wget https://storage.googleapis.com/broad-alkesgroup-public/TCSC/GeneExpressionModels/eQTL_samplesize_all/${tissue}.tar.gz
tar zxvf ${tissue}.tar.gz #creates directory TCSC_weight_files/weights/v8_allEUR_${tissue}_blup/
done

cd ../v8_320EUR/
for tissue in `cat ../TCSC/analysis/TissuesC.txt`
do
wget https://storage.googleapis.com/broad-alkesgroup-public/TCSC/GeneExpressionModels/eQTL_samplesize_320/${tissue}.post.meta.tar.gz
tar zxvf ${tissue}.tar.gz #creates directory TCSC_weight_files/weights/v8_320EUR/META_${tissue}/
done

cd ../
```

6. Run FUSION to perform TWAS analysis on provided gene expression data

```
#A. GTEx tissues that we downsampled to N = 320 individuals 
for tissue in `cat TCSC/analysis/TissuesA.txt`
do
for chr in {1..22}
do
Rscript fusion_twas-master/FUSION.assoc_test.R --sumstats $your_genomewide_sumstats --weights TCSC/weights/320EUR_metatissues/${tissue}.pos --weights_dir TCSC_weight_files/weights/v8_320EUR --ref_ld_chr LDREF/1000G.EUR. --chr $chr --out results_320/v8_320EUR.${trait}/v8_320EUR.${trait}.${tissue}.${chr}.dat
done
done

#B. GTEx tissues with small sample size that we could not downsample
for tissue in `cat TCSC/analysis/TissuesB.txt`
do
for chr in {1..22}
do
Rscript fusion_twas-master/FUSION.assoc_test.R --sumstats $your_genomewide_sumstats --weights TCSC/weights/allEUR_tissues/v8_allEUR_${tissue}_blup.pos --weights_dir TCSC_weight_files/weights --ref_ld_chr LDREF/1000G.EUR. --chr $chr --out results_all/v8_allEUR.${trait}/v8_allEUR.${trait}.${tissue}.${chr}.dat
done
done

#C. Meta-analyzed GTEx tissues which had small sample size and high genetic correlation to another GTEx tissue 
for tissue in `cat TCSC/analysis/TissuesC.txt`
do
for chr in {1..22}
do
Rscript TCSC/analysis/FUSION.assoc_test.meta.R --sumstats $your_genomewide_sumstats --weights TCSC/weights/320EUR_metatissues/${tissue}.pos --weights_dir TCSC_weight_files/weights/v8_320EUR --ref_ld_chr LDREF/1000G.EUR. --chr $chr --out results_320/v8_320EUR.${trait}/v8_320EUR.${trait}.${tissue}.${chr}.dat
done
done
```

7. Aggregate TWAS results across chromosomes

```
type=large
for tissue in `cat TCSC/analysis/TissuesA.txt`
do
Rscript concat_chrs.R ${trait},${tissue},${type}
done

type=small
for tissue in `cat TCSC/analysis/TissuesB.txt`
do
Rscript concat_chrs.R ${trait},${tissue},${type}
done

type=large
for tissue in `cat TCSC/analysis/TissuesC.txt`
do
Rscript concat_chrs.R ${trait},${tissue},${type}
done

```

8. Run TCSC using your GWAS sumstats and default GTEx gene expression data 

Obtain common SNP heritability from LDSC and set $h2g to that value. Set $N to the GWAS sample size. 

```
Rscript TCSC/analysis/Run_TCSC.R ${trait},${h2g},${N}
```

### Description of Input Data Provided in Repo

1. Signed GWAS summary statistics (formatted for S-LDSC; Finucane 2015 Nat Genet)
```
TCSC/sumstats/
```

2. Gene expression prediction models

In our primary analysis, for tissues with sample size greater than 320, we subsampled to 320 individuals per tissue; for groups of tissues with sample size less than 320 and high correlation of marginal eQTL effect sizes, we meta-analyzed gene expression models such that each meta-tissue has an effective sample size of 320 individuals. These gene expression prediction models can be found here: https://alkesgroup.broadinstitute.org/TCSC/GeneExpressionModels/eQTL_samplesize_320/

For tissues with small sample size (< 320) that could not be reasonably meta-analyzed with another tissue, we used the full set of European samples for that tissue. These gene expression prediction models can be found here: https://alkesgroup.broadinstitute.org/TCSC/GeneExpressionModels/eQTL_samplesize_all/

You can find a list of significantly cis-heritable genes (GCTA p < 0.01, GCTA h2 estimation > 0, cross validation R2 > 0) per tissue in each of the above analyses here: 
```
TCSC/weights/
```

3. TWAS summary statistics computed using FUSION (Gusev 2016 Nat Genet; http://gusevlab.org/projects/fusion/)

While we cannot provide TWAS summary statistics for all of our analyses on Github due to file quantity, we provide one representative set of TWAS summary statistics for BMI. Other TWAS summary statistics generated by our study can be found here: https://alkesgroup.broadinstitute.org/TCSC/TWAS_sumstats/
```
TCSC/twas_statistics/
```

4. Co-regulation scores

You can find the co-regulation score matrix (where columns correspond to tissues and rows correspond to significantly cis-heritable gene-tissue pairs) for the primary TCSC analysis here: 
```
TCSC/coregulation_scores/CoregulationMatrix_320orlessGTEx_062122.txt.gz
```

5. Predicted expression in 489 European individuals from 1000 Genomes
```
TCSC/predicted_expression/
```

6. Data for brain-specific analysis (TCSC power to analyze brain-related traits can be improved by restricting analysis to brain tissues):
```
TCSC/analysis/TCSC_setup_BrainSpecific.R
TCSC/analysis/InputCoreg_BrainSpecific.RData
TCSC/analysis/Run_TCSC_BrainSpecific.R
```

### Code for running the RTC Coloc method with simulated data
```
TCSC/simulation_analysis/
```

### Our steps for preparing inputs to TCSC and running TCSC

1. Using gene expression prediction models for significantly cis-heritable genes (GCTA p < 0.01, GCTA h2 estimation > 0, cross validation R2 > 0), impute gene expression into reference panel, e.g. 489 European individuals in 1000 Genomes. 

We use FUSION (http://gusevlab.org/projects/fusion/) to first compute score files for each gene expression model. Then we use plink, as done in FUSION, to impute gene expression into 1KG individuals.
```
TCSC/predicted_expression/
```

2. Construct gene and tissue co-regulation scores using squared correlations of predicted expression and applying bias correction (using gene model accuracy) when tissue t = t'. 

```
TCSC/analysis/MakeCoregulationScores_N320.R
```

3. Prepare input files for the analysis of all traits. 
```
TCSC/analysis/TCSC_setup.R
```

4. Run TCSC regression for each trait. 
```
TCSC/analysis/Run_TCSC.R ${trait},${h2g},${N}
```


