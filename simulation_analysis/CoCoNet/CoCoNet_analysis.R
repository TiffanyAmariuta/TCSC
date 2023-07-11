sim <- as.numeric(commandArgs(trailingOnly=TRUE))
#####
#simulations:
#1. 1000 genes, 10 tissues; done
#2. for each tissue, do kmeans on genes (gene x person matrix?), get 10-15 clusters (hard clustering) 
#3. tissue-specific gene adjacency matrix: block diagonal matrix for each tissue separately, 1 if belong to the same cluster. 0 else. 
#4. pick tissue 1 as causal tissue 
#5. make outcome variable using this one tissue-specific gene adjacency matrix: MVN mean 0 cov matrix Asigma_1^2 + Isigma_0^2 where A is TSGAM
#6. p = 0.02, sigma1 = 0.13, sigma2 = 0.87
#7. apply model, pick tissue with highest log likelihood, can we use CoCoNet() function for this? 

#2. 
library(data.table)
library(MASS)
library(CoCoNet)
set.seed(1)
#sims <- 15
tissues <- 0:9
samplesizes <- c(100,200,300,500,1000,1500)
ngenes <- 1000
causalgenes <- as.matrix(fread("TCSC/simulation_analysis/Nov_Y1_alpha_nCausalGenes100.txt.gz", header = F))

    outcomevar <- causalgenes[,sim] #in pratice, gene-level gwas association statistics 
    outcomevar_binary <- ifelse(causalgenes[,sim] != 0, 1, 0)

    m <- matrix(0,length(samplesizes),length(tissues))
    for (samp in 1:length(samplesizes)){
        print(samplesizes[samp])
        for (ti in tissues){
            print(ti)
            ge <- as.matrix(fread(paste0("TCSC/simulation_analysis/SEG/totalge/",samplesizes[samp],"/TotExp_Nov_Tissue",ti,"_sim",sim,".txt.gz"), header = F))
            #ge <- fread(paste0("totalexpression/Nov_0.75_ggcoreg/",samplesizes[samp],"/TotExp_Nov_Group1_Tissue",ti,"_sim",sim,".txt.gz"), header = F)
            kobj <- kmeans(ge, centers = 10) #clusters on the rows (genes), cols are people 
            clusters <- kobj$cluster
            adjmat <- matrix(0,nrow(ge),nrow(ge))
            for (i in 1:10){
                w <- which(clusters == i)
                adjmat[w,w] <- 1
            }
            #if(ti == 0){outcomevar <- mvrnorm(1, mu = rep(0,nrow(adjmat)), Sigma = adjmat*0.13^2 + diag(nrow(adjmat))*0.87^2 )} #pergene, per disease (so just 1)
            result = CoCoNet(outcomevar_binary[1:ngenes], max_path = 1, adjmat)
            print(result$loglikelihood)
            m[samp,ti+1] <- result$loglikelihood
        }
    }
    filename <- paste0("TCSC/simulation_analysis/CoCoNet/results/res_",sim,".txt")
    write.table(m, file = filename, row.names = F, col.names = F, sep = "\t", quote = F) 
    system(paste0("gzip TCSC/simulation_analysis/CoCoNet/results/res_",sim,".txt"))

#afterward 
sims <- 1000 
res_mat <- data.frame()
for (sim in 1:sims){
  tryCatch({
  res <- as.matrix(fread(paste0("TCSC/simulation_analysis/CoCoNet/results/res_",sim,".txt.gz")))
  res_mat <- rbind(res_mat, cbind(1:6,res))
  },error=function(cond){message(paste("Didn't finish ", sim))})
}
power <- c()
fpr <- c()
for (i in 1:6){
w <- which(res_mat[,1] == i)
res_mat2 <- res_mat[w,-1]
performance <- sapply(1:nrow(res_mat2), function(x) which.max(as.numeric(res_mat2[x,]))) #1 1 1 2
power[i] <- length(which(performance == 1))/length(performance)*100
fpr[i] <- 100-power[i]
}
####

#https://lulushang.org/docs/Projects/CoCoNet/Reproduce

#load("CoCoNet/outcome_tissue_scale.RData") #double, one column per trait. rows = genes. already for one tissue only? 
#load("CoCoNet/tissue_net.RData") #very large. 

#gene by gene correlation, squared? how do they compute it? follow tutorial.

#scaled gene level effect sizes, twas res? how do they compute it? follow tutorial.

#A = tissue_net[[1]] #for one tissue
#result = CoCoNet(outcome_tissue_scale[,3], max_path = 1, A) #genes on rows 1 column for disease. 


#from https://lulushang.org/docs/Projects/CoCoNet/Reproduce
#library(data.table)
#load("CoCoNet/outcome_cell_scale.RData")
#gene_table = fread("CoCoNet/genetable.txt.gz", header=T)
#motif = fread("CoCoNet/motif.txt.gz") # same as in GTEx tissue PANDA input, downloaded from same website 
#ppi = fread("CoCoNet/ppi.txt.gz") # same as in GTEx tissue PANDA input, downloaded from same website 
#expr = read.table("CoCoNet/GTEx_droncseq_hip_pcf/GTEx_droncseq_hip_pcf.umi_counts.txt.gz",header=T) # GTEx single cell expression genes x cells
#cell_anno = read.csv("CoCoNet/GTEx_droncseq_hip_pcf/cell_annotation.csv",header=T) # cell type annotation for GTEx single cell expression cells x 6 columns 
#library(devtools)
#devtools::install_github('xzhoulab/CoCoNet')

