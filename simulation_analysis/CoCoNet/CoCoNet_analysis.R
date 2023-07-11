sim <- as.numeric(commandArgs(trailingOnly=TRUE))
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
            kobj <- kmeans(ge, centers = 10) #clusters on the rows (genes), cols are people 
            clusters <- kobj$cluster
            adjmat <- matrix(0,nrow(ge),nrow(ge))
            for (i in 1:10){
                w <- which(clusters == i)
                adjmat[w,w] <- 1
            }
            result = CoCoNet(outcomevar_binary[1:ngenes], max_path = 1, adjmat)
            print(result$loglikelihood)
            m[samp,ti+1] <- result$loglikelihood
        }
    }
    filename <- paste0("TCSC/simulation_analysis/CoCoNet/results/res_",sim,".txt")
    write.table(m, file = filename, row.names = F, col.names = F, sep = "\t", quote = F) 
    system(paste0("gzip TCSC/simulation_analysis/CoCoNet/results/res_",sim,".txt"))


    

