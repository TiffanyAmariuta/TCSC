sim <- as.numeric(commandArgs(trailingOnly=TRUE))
library(data.table)

ngwas <- 10000
ngenes <- 1000
nvar <- 5
snps <- fread("TCSC/simulation_analysis/1KG_HM3_chr1.bim", header = F)
QTL_geno_all <- as.matrix(fread(paste0("TCSC/simulation_analysis/Simulated_eQTL_Cohort_for_ImputedExpression.txt.gz"), header = F))
tissues <- 0:9
samplesizes <- c(100,200,300,500,1000,1500)
boundary <- 500000
genes <- fread("TCSC/simulation_analysis/1000genes_Nov.txt", header = T)
colnames(QTL_geno_all) <- snps$V2
print(sim)
for (samp in 1:length(samplesizes)){
    print(samp)
    QTL_geno <- QTL_geno_all[1:samplesizes[samp],]
    for (ti in tissues){
       print(ti)
       #load total expression
       tot_GE <- as.matrix(fread(paste0("TCSC/simulation_analysis/totalexpression/Nov_0.75_ggcoreg/",samplesizes[samp],"/TotExp_Nov_Group1_Tissue",ti,"_sim",sim,".txt.gz"),header = F)) #note, this is an output of the main TCSC simulations
       pval_mat <- data.frame()
       for (gene in 1:nrow(tot_GE)){
           snps_in_range <- which(snps$V4 >= genes$START[gene]-boundary & snps$V4 <= genes$STOP[gene]+boundary)
           pvals <- sapply(1:length(snps_in_range), function(x) cor.test(tot_GE[gene,], QTL_geno[,snps_in_range[x]])$p.val)
           dump <- cbind(pvals,gene,snps_in_range)
           pval_mat <- rbind(pval_mat,dump)
       }
       eqtls <- which(as.numeric(pval_mat[,1]) < 0.05/nrow(pval_mat))
       filename <- paste0("TCSC/simulation_analysis/RTC_Coloc/output/Ongen_files/eqtl_tissue",ti,"/SigeQTLs_samp",samplesizes[samp],"_sim",sim,".txt") #note, need to make this directory
       write.table(pval_mat[eqtls,c(2:3)], file = filename, row.names = F, col.names = F, sep = "\t", quote = F) #gene and snp combo. 
       system(paste0("gzip -f ",filename))
   }
}
#1.39911e-06 threshold, same regardless of tissue or samp, use this threshold for conditional analysis
#eqtls identified here should be coming up in gwas coloc step next.   
