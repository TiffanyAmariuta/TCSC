sim <- as.numeric(commandArgs(trailingOnly=TRUE))
set.seed(sim)

library(rolypoly)
library(data.table)

ngenes <- 1000
samplesizes <- c(100,200,300,500,1000,1500)
tissues <- 0:9
ve <- 0.1

#did once: 
genes <- fread("TCSC/simulation_analysis/1000genes_Nov.txt", header = T)
annot <- cbind(genes$CHR,genes$START,genes$STOP,paste0("g",1:ngenes))
colnames(annot) <- c("chrom","start","end","label")
#write.table(annot, file = "RolyPoly/annot.txt", row.names = F, col.names = T, sep = "\t", quote = F) #was problem with levels, so wrote and reread file. 
annot <- fread("TCSC/simulation_analysis/RolyPoly/annot.txt", header = T)
annot <- as.data.frame(annot)
#tutorial: sim_block_annotation

#sims <- 1:1000
#tutorial: sim_expression_data_normalized
#res_mat <- matrix(0,length(samplesizes),3)
#want matrix where each row is a sim, samp pair. 
res_mat <- matrix(0,0,2+length(tissues))
#for (sim in sims){
    print(sim)
    sumstats <- as.data.frame(fread(paste0("TCSC/simulation_analysis/RolyPoly/sumstats/sumstats_10000_vgene",ve,"_sim",sim,".txt.gz"), header = T)) #tutorial: sim_gwas_data
    for (samp in 1:length(samplesizes)){
        ge_alltissues <- matrix(0,ngenes,length(tissues))
        for (ti in tissues){
            ge <- as.matrix(fread(paste0("TCSC/simulation_analysis/SEG/totalge/",samplesizes[samp],"/TotExp_Nov_Tissue",ti,"_sim",sim,".txt.gz"), header = F))
            ge_alltissues[,ti+1] <- rowSums(ge)/ncol(ge)
        }
        colnames(ge_alltissues) <- paste0("Tissues",tissues)
        rownames(ge_alltissues) <- paste0("g",1:ngenes)
#        if(sim == 1 & samp == 1){
            rp <- rolypoly_roll(
            gwas_data = sumstats,
            block_annotation = annot,
            block_data = data.frame(ge_alltissues),
            ld_folder = "RolyPoly/EUR_LD_FILTERED_NONAN_R"
            )
#        }else{
#            rp_update <- rolypoly_roll(
#            rolypoly = rp,
#            block_data = data.frame(ge_alltissues)
#            )
#        }  
        pvals <- rp$bootstrap_results$bp_value[-1]
        #inferred_causal <- tissues[which(pvals < 0.05)]        
        res_mat <- rbind(res_mat, c(sim, samp, pvals))
#        res_mat[samp,1] <- res_mat[samp,1] + length(which(inferred_causal == 0))
#        res_mat[samp,2] <- res_mat[samp,2] + length(which(inferred_causal %in% c(0,1,2)))
#        res_mat[samp,3] <- res_mat[samp,3] + length(which(inferred_causal %in% c(1:9)))
    }
  write.table(res_mat, file = paste0("TCSC/simulation_analysis/RolyPoly/res/rolypoly_res_10000_ve0.1_sim",sim,".txt"), row.names=F, col.names = F, sep = "\t", quote = F)
#} 
#write.table(res_mat, file = "rolypoly_res.txt", row.names=F, col.names = F, sep = "\t", quote = F, append = T) #ve = 0.05

library(data.table)
res_mat <- data.frame()
sims <- 1000 
for (sim in 1:sims){
  tryCatch({
  res <- as.matrix(fread(paste0("TCSC/simulation_analysis/RolyPoly/res/rolypoly_res_10000_ve0.1_sim",sim,".txt")))
  },error=function(cond){message(paste("Didn't finish ", sim))})
  res_mat <- rbind(res_mat, res)
}
res_mat <- as.matrix(res_mat) #each row: sim, samp, p values
samplesizes <- c(100, 200, 300, 500, 1000, 1500)
power <- c()
fpr <- c()
for (samp in 1:length(samplesizes)){
    rmat <- res_mat[which(res_mat[,2] == samp),]
    power[samp] <- length(which(rmat[,3] < 0.05))/length(rmat[,3])*100
    fpr[samp] <- length(which(rmat[,4:12] < 0.05))/length(rmat[,4:12])*100
}
print(power)
print(fpr)
#pretty good. 

#ROC curve. Pick one sample size; n = 300, because this will affect our analysis of Figure 5. 
pvals <- seq(0,1,0.001) #1000 points
for (samp in 1:6){
sensitivity <- c()
oneminus_specificity <- c()
#samp <- 3
rmat <- res_mat[which(res_mat[,2] == samp),]
for (p in 1:length(pvals)){
    sensitivity[p] <- length(which(rmat[,3] < pvals[p]))/length(rmat[,3])
    oneminus_specificity[p] <- length(which(rmat[,4:12] < pvals[p]))/length(rmat[,4:12])
}
a <- cbind(oneminus_specificity,sensitivity)
colnames(a) <- c("oneminus_specificity","sensitivity")
write.table(a, file = paste0("TCSC/simulation_analysis/RolyPoly_ROCcurve_samp",samp,".txt"), row.names = F, col.names = T, sep = "\t", quote = F)
}
#from tutorial
#ld_path <- system.file("extdata", "example_ld", package = "rolypoly")
#ld_data <- readRDS(paste(ld_path, '/1.Rds', sep = ''))
#rp <- rolypoly_roll(
#  gwas_data = sim_gwas_data,
#  block_annotation = sim_block_annotation,
#  block_data = sim_expression_data_normalized,
#  ld_folder = ld_path
#)

#add new ge (for other samps) without having to link SNP data again
#rp <- rolypoly_roll(
  # some new set of expression data
 # block_data = ge_nextsamp,
#)

simple_auc <- function(TPR, FPR){
  # inputs already sorted, best scores first 
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}
library(data.table)
for (samp in 1:6){
y <- fread(paste0("TCSC/simulation_analysis/RolyPoly_ROCcurve_samp",samp,".txt"), header=T)
print(simple_auc(y$sensitivity,y$oneminus_specificity))
}
