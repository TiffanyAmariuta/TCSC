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

# SIMULATE RUNNING AN EQTL STUDY + DERIVE WEIGHTS FOR CIS SNP GENE MODELS 

sim = int(sys.argv[1])
print(sim)
ti = int(sys.argv[2])
print(ti)

nGenesTot = 1000
nvar = 5 
samplesizes = [100,200,300,500,1000,1500]
filename = "TCSC/simulation_analysis/eQTLs_1000sims_Nov/Nov_sims_eQTLeffectsizes_0.75rg_SNPpolygen5_tissue%s_sim%s.txt.gz" % (ti+1,sim) 
b_qtls_df = pd.read_csv(filename, sep = "\t", header = None) #flag
b_qtls = np.array(b_qtls_df) #flag
nCausalGenes = 100

filename = "TCSC/simulation_analysis/Nov_cish2_pergene_toachievebettercoreg.txt"
gene_cish2 = pd.read_csv(filename, sep = "\t", header = None)
gene_cish2 = np.array(gene_cish2)

genenum = []
for i in range(1,nGenesTot+1):
    a = [i]*nvar
    genenum = np.concatenate((genenum,a))

bim, fam, G_all = read_plink("1KG_HM3_chr1", verbose=False)
G_all = G_all.T
bim = np.array(bim)

filename = "TCSC/simulation_analysis/Simulated_eQTL_Cohort_for_ImputedExpression.txt.gz"
z_eqtl = np.array(pd.read_csv(filename, sep = "\t", header = None))

#import code
#code.interact(local=locals()) #remember need to activate python 3 otherwise libraries don't work. 
count_g = 0
count_g_n = np.zeros(len(samplesizes))
#count_g_n100 = 0
#count_g_n200 = 0
#count_g_n300 = 0 
#count_g_n500 = 0 
#count_g_n1000 = 0 
#count_g_n1500 = 0

for gene in range(1,nGenesTot+1):
    np.random.seed(int(sim*gene))    
    x = np.where(genenum == gene)[0]
    min_snp = int(b_qtls[np.min(x),2]) #make 1-indexed values 0-indexed for python. is this the same per sim? 
    max_snp = int(b_qtls[np.max(x),3])
    p_int = int(max_snp-(min_snp-1))
    b_qtls_pass = np.zeros([p_int,1])
    pos = b_qtls[x,0]
    pos = pos.astype(int)
    sequence = np.array(range(min_snp,max_snp+1))
    x_var = []
    for j in pos:
        a = np.where(sequence == j)[0]
        x_var = np.concatenate((x_var,a))
    x_var = x_var.astype(int)
    b_qtls_pass[x_var,0] = b_qtls[x,1] #need new x-var where chunk is the gene.
    for nqtl in samplesizes:
        print(nqtl)
        Z_qtl = z_eqtl[range(nqtl),(min_snp-1):(max_snp)] #do this for each sample size and assign name to use later. 
        h2 = gene_cish2[gene-1,0]
        if np.sum(b_qtls_pass) != 0: 
            coef, h2g, hsq_p, r2all, r2train, gexpr = sim_eqtl(gene, ti, sim, Z_qtl, nqtl, b_qtls_pass, float(h2))
        else: #skipping gcta step  
            h2g = 0 
            hsq_p = 1 
            r2all = 0 
            r2train = 0
            gexpr = sim_trait(np.dot(Z_qtl,b_qtls_pass), 0)[0] #just random normal, why waste time saving this info? just increases file size by 2x. not terrible.  
        totalGE = np.zeros([1,nqtl])
        totalGE[0,:] = gexpr    
        xxx = [nqtl, gene, h2g, hsq_p, r2all, r2train]
        xxx = np.transpose(pd.DataFrame(xxx))
        filename = "TCSC/simulation_analysis/weights/Nov_0.75_ggcoreg/h2_r2_Nov_Group1_Tissue%s_sim%s.txt.gz" % (ti,sim) #every gene x sample size 
        if ((gene==1) & (nqtl == 100)): #only for very first time, not for each sample size 
            pd.DataFrame(xxx).to_csv(filename, sep="\t", index=False, header = False)
        else: 
            pd.DataFrame(xxx).to_csv(filename, sep="\t", index=False, mode="a", header = False)
        filename = "TCSC/simulation_analysis/totalexpression/Nov_0.75_ggcoreg/%s/TotExp_Nov_Group1_Tissue%s_sim%s.txt.gz" % (nqtl,ti,sim) #every gene
        totalGE = pd.DataFrame(totalGE)
        if gene == 1: #because totalexpression file per sample size. 
            totalGE.round(4).to_csv(filename, sep="\t", index=False, header = False)
        else:
            totalGE.round(4).to_csv(filename, sep="\t", index=False, mode = "a", header = False)
        if h2g > 0 and hsq_p < 0.01:  
            count_g = count_g + 1
            count_g_n[samplesizes.index(nqtl)] = count_g_n[samplesizes.index(nqtl)] + 1 
            xxx = np.transpose(pd.DataFrame(np.append(nqtl,coef))) #since each gene has a different number of cis SNPs, this will not be a matrix. different number of columns each time. 
            filename = "TCSC/simulation_analysis/weights/Nov_0.75_ggcoreg_coef/coef_Nov_Group1_Tissue%s_sim%s.txt.gz" % (ti,sim)
            if count_g == 1: #first cis-h2 gene found.
                xxx.round(6).to_csv(filename, sep="\t", index=False, header = False)
            else:
                xxx.round(6).to_csv(filename, sep="\t", index=False, mode = "a", header = False)
            predexp = np.transpose(pd.DataFrame(np.dot(z_eqtl[range(500),(min_snp-1):(max_snp)],coef))) #t(people on rows), want people on columns and genes on rows. 
            filename = "TCSC/simulation_analysis/predExp/Nov_0.75_ggcoreg/%s/predExp_Nov_Group1_sim%s_Tissue%s.txt.gz" % (nqtl,sim,ti)
            if count_g_n[samplesizes.index(nqtl)] == 1: #first cis-h2 gene found per nqtl
                predexp.round(6).to_csv(filename, sep="\t", index=False, header = False)
            else:
                predexp.round(6).to_csv(filename, sep="\t", index=False, mode = "a", header = False)

      
