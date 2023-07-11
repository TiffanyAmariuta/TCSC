sim <- as.numeric(commandArgs(trailingOnly=TRUE))

library(data.table)
library(MASS)
set.seed(sim)
ngenes <- 1000
samplesizes <- c(100,200,300,500,1000,1500)
causalgenes <- as.matrix(fread("TCSC/simulation_analysis/Nov_Y1_alpha_nCausalGenes100.txt.gz", header = F))
causalgenes <- which(causalgenes[,sim] != 0)
tissues <- 0:9

for (samp in 1:length(samplesizes)){
    for (ti in tissues){
        ge <- as.matrix(fread(paste0("TCSC/simulation_analysis/totalexpression/Nov_0.75_ggcoreg/",samplesizes[samp],"/TotExp_Nov_Group1_Tissue",ti,"_sim",sim,".txt.gz")))
        if(ti %in% c(0,1,2)){
            randomothergenes <- sample(c(1:ngenes)[-causalgenes], size = 200, replace = F)
            de_in_tissue <- c(causalgenes,randomothergenes)
            de_down_tissue <- sample(c(1:ngenes)[-de_in_tissue], size = 100, replace = F)
        }else{
            de_in_tissue <- sample(c(1:ngenes), size = 300, replace = F)
            de_down_tissue <- c(1:ngenes)[-de_in_tissue]
        }
        addDE <- runif(n = length(de_in_tissue)*samplesizes[samp],min=0, max = 10)
        minusDE <- runif(n = length(de_down_tissue)*samplesizes[samp],min=0, max = 10)
        ge[de_in_tissue,] <- ge[de_in_tissue,]+matrix(addDE,ncol=samplesizes[samp])
        ge[de_down_tissue,] <- ge[de_down_tissue,]-matrix(minusDE,ncol=samplesizes[samp])
        filename <- paste0("TCSC/simulation_analysis/SEG/totalge/",samplesizes[samp],"/TotExp_Nov_Tissue",ti,"_sim",sim,".txt")
        write.table(ge, file = filename, row.names = F, col.names = F, sep = "\t", quote = F)
        system(paste0("gzip -f ",filename))
        assign(paste0("ge_ti",ti,"_nqtl",samp),ge)
    }
}

#tstats
modules <- cbind(c(1,1,2,4,4,5,7,7,7,8,8,9),c(2,3,3,5,6,6,8,9,10,9,10,10))
modules <- modules - 1

for (samp in 1:length(samplesizes)){
    print(samp)
    for (ti in tissues){
        print(ti)
        immune <- c(ti,modules[which(modules[,1] == ti),2])
        ge_focal <- get(paste0("ge_ti",ti,"_nqtl",samp))
        compare <- tissues[-match(immune,tissues)]

        tstats <- c()
        for (g in 1:ngenes){
            y_other <- sapply(1:length(compare), function(x) get(paste0("ge_ti",compare[x],"_nqtl",samp))[g,])
            y <- c(ge_focal[g,],y_other)
            lab <- c(rep(1,ncol(ge_focal)), rep(-1,length(y)-ncol(ge_focal)))
            mod <- lm(y ~ lab)
            res <- summary(mod)
            tstats[g] <- res$coefficients[2,3]
            #take the top 10% as SEG for this tissue. 
        }
        s <- sort(tstats,decreasing =T)
        thresh <- s[0.1*length(s)]
        seg <- which(tstats > thresh)
        write.table(seg, file = paste0("TCSC/simulation_analysis/SEG/tstats/SEG_Ti",ti,"_nqtl",samp,"_sim",sim,".txt"), row.names = F, col.names = F, sep = "\t", quote = F)
    }
}

system(paste0("gzip -f TCSC/simulation_analysis/SEG/tstats/*_sim",sim,".txt"))

#annotate SNPs 
library(data.table)
snps <- fread("TCSC/simulation_analysis/1000G.EUR.QC.hm3_noMHC.1.frq", header = F) #98K 1000G EUR SNPs 
#write.table(snps$V2, file = "list.txt", row.names = F, col.names = F, sep = "\t", quote = F)
snpannot <- fread("TCSC/simulation_analysis/sample_annot.1.annot.gz", header=T) #obtain any sample annot from s-ldsc 
snpannot$Chromatin <- 0
w <- which(snpannot$SNP %in% snps$V2)
snpannot <- snpannot[w,]
genes <- fread("TCSC/simulation_analysis/1000genes_Nov.txt", header = T)
for (ti in tissues){
    for (samp in 1:length(samplesizes)){
        y <- as.numeric(fread(paste0("TCSC/simulation_analysis/SEG/tstats/SEG_Ti",ti,"_nqtl",samp,"_sim",sim,".txt.gz"), header = F)$V1)
        genessubset <- genes[y,]
        snpannot_fill <- snpannot
        for (i in 1:nrow(genessubset)){
            w <- which(snpannot_fill$BP > genessubset$START[i] - 100000 & snpannot_fill$BP < genessubset$STOP[i] + 100000)
            snpannot_fill$Chromatin[w] <- 1
        }
        write.table(snpannot_fill, paste0("TCSC/simulation_analysis/SEG/annots/SEG_Ti",ti,"_nqtl",samp,"_sim",sim,".annot"), row.names = F, col.names = T, sep = "\t", quote = F)
    }
}
system(paste0("gzip -f TCSC/simulation_analysis/SEG/annots/*_sim",sim,".annot"))



