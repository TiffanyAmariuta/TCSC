#!/usr/bin/env python
import argparse as ap
import sys
import numpy as np
import pandas as pd
import scipy.linalg as linalg
from numpy.linalg import multi_dot as mdot
from pandas_plink import read_plink
from scipy import stats
from sklearn import linear_model as lm
import random
import warnings
import scipy
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
#from sklearn.utils.testing import ignore_warnings #not loading 
from sklearn.exceptions import ConvergenceWarning
import time
import subprocess, os
import gzip
import os.path  #new
from os import path  #new
from statsmodels.regression.linear_model import OLS

mvn = stats.multivariate_normal

def sim_geno(L, n): 
    """
    removed first argument with is L (works with regular LD matrix too), added p to indicate number of snps.
    Sample genotypes from an MVN approximation.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param n: int the number of genotypes to sample

    :return: numpy.ndarray n x p centered/scaled genotype matrix
    """
    p, p = L.shape
    Z = L.dot(np.random.normal(size=(n, p)).T).T
    #Z = np.random.normal(size=(n, p)).T #snps x people
    #Z = Z.T #I added this. people x snps
    Z -= np.mean(Z, axis=0)
    Z /= np.std(Z, axis=0) #173 x 4 in original, ine too after transpose.
    return Z

def sim_gwas(ngwas, Z_gwas, b_qtls, ve1, ve2, alpha1, alpha2):
    """
    Simulate a GWAS using `ngwas` individuals such that genetics explain `var_explained` of phenotype.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param ngwas: int the number of GWAS genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene
    :param var_explained: float the amount of phenotypic variance explained by genetic component of gene expression

    :return: (pandas.DataFrame, float) estimated GWAS beta and standard error, causal GE effect
    """
    nGenes = b_qtls.shape[1]
    complextrait1 = []
    complextrait2 = []
    for g in range(nGenes):
        ve1_g = ve1[g]
        ve2_g = ve2[g]
        if ve1_g == 0:
            y1 = np.random.normal(0, 1, ngwas)
        if ve2_g == 0:
            y2 = np.random.normal(0, 1, ngwas)
        if ve1_g > 0:
            gwas_expr = np.dot(Z_gwas, b_qtls[:,g])
            xalpha1 = gwas_expr #* alpha1[g]
            y1 = sim_trait(xalpha1, ve1_g)[0]
        if ve2_g > 0:
            gwas_expr = np.dot(Z_gwas, b_qtls[:,g])
            xalpha2 = gwas_expr #* alpha2[g]
            y2 = sim_trait(xalpha2, ve2_g)[0]
        complextrait1.append(y1)
        complextrait2.append(y2)
    complextrait1 = np.array(complextrait1)
    complextrait2 = np.array(complextrait2)
    return (complextrait1, complextrait2) 

def sim_trait(g, h2g): #use this twice, first simulate GE as a trait in the GE cohort, then simulate complex trait based on GE in the GWAS cohort 
    """
    Simulate a complex trait as a function of latent genetic value and env noise.

    :param g: numpy.ndarray of latent genetic values
    :param h2g: float the heritability of the trait in the population

    :return: (numpy.ndarray, float) simulated phenotype, sd of Y
    """
    n = len(g)
    if h2g > 0:
        s2g = np.var(g, ddof=1)
        s2e = s2g * ( (1.0 / h2g ) - 1 )
        e = np.random.normal(0, np.sqrt(s2e), n)
        y = g + e #appears to add noise to the genetic values, adding the same about of noise if the seed is the same 
    else:
        e = np.random.normal(0, 1, n)
        y = e
    # standardize
    y -= np.mean(y)
    y_std = np.std(y)
    y /= y_std
    y_std = y_std.item()
    return y, y_std #should we add these together to get a complex trait? 

def fit_lasso_katie(Z, y):
    """
    Infer eqtl coefficients using LASSO regression for simeQTL function above. Uses the PLINK-style coordinate descent algorithm
    that is bootstrapped by the current h2g estimate.
    :param Z:  numpy.ndarray n x p genotype matrix
    :param y: numpy.ndarray gene expression for n individuals
    :return: (numpy.ndarray, float, float) tuple of the LASSO coefficients, the r-squared score, and log-likelihood
    """
    n, p = Z.shape
    alphas=np.logspace(-2,1.2,50,base=10)
    kf=KFold(n_splits=5)
    CV_r2_allpenalty = []
    CV_r2_allpenalty_std = []
    CV_r2_allpenalty_sem = []
    for penalty in alphas:
        CV_r2=[]
        for train_index, test_index in kf.split(Z):
            X_train, X_test = Z[train_index], Z[test_index]
            y_train, y_test = y[train_index], y[test_index]
            y_train_std=(y_train-np.mean(y_train))/np.std(y_train)
            #this seems to be not quite standardizing correctly
            y_test_std=(y_test-np.mean(y_test))/np.std(y_test)
            lasso = lm.Lasso(fit_intercept=True)
            lasso.set_params(alpha=penalty,tol=1e-4,max_iter=1000)
            lasso.fit(X_train,y_train_std)
            r2 = lasso.score(X_test, y_test_std)
            CV_r2.append(r2)
        CV_r2_allpenalty.append(np.mean(CV_r2))
        CV_r2_allpenalty_std.append(np.std(CV_r2))
        CV_r2_allpenalty_sem.append(scipy.stats.sem(CV_r2))
    #print(CV_r2_allpenalty)
    besti=np.argmax(CV_r2_allpenalty)
    bestalpha=alphas[besti]
    lasso = lm.Lasso(fit_intercept=True)
    lasso.set_params(alpha=bestalpha,tol=1e-4,max_iter=1000)
    lasso.fit(Z, y)
    r2_train = lasso.score(Z, y)
    coef = lasso.coef_
    return coef, CV_r2_allpenalty[besti], r2_train

def sim_eqtl(g, ti, sim, Z_qtl, nqtl, b_qtls, eqtl_h2):
    """
    removed first arg which was LD.
    Simulate an eQTL study using `nqtl` individuals, can do multiple genes.
    :param LD: LD matrix (from real or fake genotypes? katie doesn't simulate LD matrix here, so supply a real LD matrix of chosen size, how many variants ? MxM)
    :param nqtl: int the number of eQTL-panel genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene: SNPs x genes
    :param eqtl_h2: list of expression variance explained by linear model of SNPs for each gene (not used) rather min(1,b@LD_qtl[np.ix_(nearSNPsi,nearSNPsi)]@b)
    :param geneLocs: a numpy array of the start position of each gene
    :param SNPloc: a list of the location of each SNP in basepair units
    :param eqtlDist: a  integer that's how far away from the start of the gene eQTLs can be ( I use 50,000 for my sims)
    :return:  (numpy.ndarray x 5) Arrays of eqtl effect sizes and heritability and cross validation accuracies
    """
    allqtlshape = b_qtls.shape
    nGenes = b_qtls.shape[1] #should be snps by genes, not just one vector
    fullNSNPs = b_qtls.shape[0] #LD.shape[0]
    n, p = [float(x) for x in Z_qtl.shape] #LD_qtl = np.dot(Z_qtl.T, Z_qtl) / n

    tempdir = "twas_sim/Sims_biascorrec_minhsqthresh/temps_ggcoreg/temp_gene"+str(g)+"/temp_gene"+str(g)+"_tissue"+str(ti)+"_sim"+str(sim)+"_nqtl"+str(nqtl)
    b = b_qtls[:,0]
    nonzeroindex = np.where(b != 0)[0]
    gexpr = sim_trait(np.dot(Z_qtl,b), eqtl_h2)[0]
    coef, r2, r2_train = fit_lasso_katie(Z_qtl, gexpr)
    if r2_train>0 and r2>0:
        h2g,hsq_se,hsq_p = get_hsq(Z_qtl,gexpr,tempdir,"fusion_twas-master/gcta_nr_robust")
    else:
        h2g=0
        hsq_se=10
        hsq_p=1
    return (coef,h2g,hsq_p,r2,r2_train,gexpr)

def get_hsq(geno, gene_expr, out, gcta):
    nindv, nsnp = geno.shape
    grm = np.dot(geno, geno.T) / nsnp
    write_pheno(gene_expr, out)
    write_grm(grm, nsnp, out)
    run_gcta(gcta, out)
    clean_up(out)
    hsq_f = '{}.hsq'.format(out)
    if os.path.exists(hsq_f):
        hsq_out = pd.read_table(hsq_f)
        hsq = hsq_out['Variance'][0]
        hsq_se = hsq_out['SE'][0]
        #Should also be able to get p-value like Sasha
        hsq_p = hsq_out['Variance'][8]
        command = 'rm {}'.format(hsq_f)
        subprocess.run(command.split())
        return hsq, hsq_se, hsq_p
    else:
        return 0,10000,1
    return None

def run_gcta(gcta, out):
    command = '{} --grm-gz {} --pheno {} --reml --reml-no-constrain --thread-num 4 --out {}'\
        .format(gcta, out, out+'.phen', out)
    subprocess.run(command.split())
    command = 'rm {}.log'.format(out)
    subprocess.run(command.split())
    command = 'rm {}.phen'.format(out)
    subprocess.run(command.split())

def write_pheno(phen, out):
    out_f = open(out+'.phen', 'w')
    nindv = phen.shape[0]
    for i in range(nindv):
        out_f.write('{}\t{}\t{}\n'.format(i, i, phen[i]))
    out_f.close()

def clean_up(out):
    command = 'rm {}.grm.gz'.format(out)
    subprocess.run(command.split())
    command = 'rm {}.grm.id'.format(out)
    subprocess.run(command.split())

def write_grm(grm, nsnp, out):
    out_id_f = open(out+'.grm.id', 'w')
    out_f = gzip.open(out+'.grm.gz', 'w')
    nindv = grm.shape[0]
    for i in range(nindv):
        for j in range(0, i+1):
            out_f.write('{}\t{}\t{}\t{}\n'.format(i+1, j+1, nsnp,
                grm[i,j]).encode('utf-8'))
        out_id_f.write('{}\t{}\n'.format(i, i))
    out_f.close()
    out_id_f.close()

# Create simulated 1Kg cohort for imputation 

bim, fam, G = read_plink("TCSC/simulation_analysis/1KG_HM3_chr1", verbose=False)
G = G.T
bim = np.array(bim)

np.random.seed(12345)    
# estimate LD for population from PLINK data
n, p = [float(x) for x in G.shape]
p_int = int(p)
mafs = np.mean(G, axis=0) / 2
G -= mafs * 2
G /= np.std(G, axis=0)
# regularize so that LD is PSD (positive semi definite)
LD = np.dot(G.T, G) / n + np.eye(p_int) * 0.1 #this is the LD matrix.
# compute cholesky decomp for faster sampling/simulation
L = linalg.cholesky(LD, lower=True)
Z_qtl = pd.DataFrame(sim_geno(L,1500))
filename = "TCSC/simulation_analysis/Simulated_eQTL_Cohort_for_ImputedExpression.txt.gz" 
Z_qtl.round(4).to_csv(filename, sep="\t", index=False,mode = "a", header = False)


      
