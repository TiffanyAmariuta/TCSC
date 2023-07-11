sim_vegene_vesnp <- commandArgs(trailingOnly=TRUE)
s <- strsplit(sim_vegene_vesnp, split = "_")[[1]]
sim <- as.numeric(s[1])
ve_gene <- s[2]
ve_snp <- s[3]
ngwas <- s[4]

library(data.table)

#Make sumstats with this format: 
#SNP     A1      A2      N       CHISQ   Z #for ldsc 
#rs3766192       T       C       25604.000       2.604996        -1.614
#for rolypoly
#chrom  pos       rsid        beta    se        maf

snps <- fread("TCSC/simulation_analysis/1KG_HM3_chr1.bim", header = F)
frq <- fread("TCSC/simulation_analysis/1000G.EUR.QC.hm3_noMHC.1.frq", header = T)
frq <- frq[match(snps$V2,frq$SNP),]

sumstat_names <- c()
if(as.numeric(ve_gene) != 0){sumstat_names <- c(sumstat_names, paste0("vgene",ve_gene,"_sim",sim,".txt"))}
if(as.numeric(ve_snp) != 0){sumstat_names <- c(sumstat_names, paste0("vgene",ve_gene,"_vsnp",ve_snp,"_sim",sim,".txt"))}

for (i in 1:length(sumstat_names)){
    sumstats <- cbind(snps$V2, snps$V5, snps$V6, ngwas, 0, 0) #LDSC SEG
    sumstats_rolypoly <- cbind(snps$V1,snps$V4,snps$V2,0,0,frq$MAF) #RolyPoly
    mysumstats <- as.matrix(fread(paste0("TCSC/simulation_analysis/GWASsumstats/sumstats_",ngwas,"_",sumstat_names[i],".gz")), header = F)
    z <- as.numeric(mysumstats[,2])/as.numeric(mysumstats[,3])
    beta <- as.numeric(mysumstats[,2])
    se <- as.numeric(mysumstats[,3])
    chisq <- z^2
    sumstats_rolypoly[,5] <- se
    sumstats_rolypoly[,4] <- beta
    sumstats[,6] <- z
    sumstats[,5] <- chisq 

    colnames(sumstats) <- c("SNP","A1","A2","N","CHISQ","Z")
    write.table(sumstats, file = paste0("TCSC/simulation_analysis/SEG/sumstats/sumstats_",ngwas,"_",sumstat_names[i]), row.names = F, col.names = T, sep = "\t", quote = F)
    system(paste0("gzip -f TCSC/simulation_analysis/SEG/sumstats/sumstats_",ngwas,"_",sumstat_names[i]))

    colnames(sumstats_rolypoly) <- c("chrom","pos","rsid","beta","se","maf")
    write.table(sumstats_rolypoly, file = paste0("TCSC/simulation_analysis/RolyPoly/sumstats/sumstats_",ngwas,"_",sumstat_names[i]), row.names = F, col.names = T, sep = "\t", quote = F)
    system(paste0("gzip -f TCSC/simulation_analysis/RolyPoly/sumstats/sumstats_",ngwas,"_",sumstat_names[i]))
}
