library(data.table)
set.seed(1)
sims <- 1000
ve <- 0.1
ngwas <- 10000
ngenes <- 1000
nvar <- 5
snps <- fread("TCSC/simulation_analysis/1KG_HM3_chr1.bim", header = F)

#distance of every SNP to nearest TSS
#put into bins (assign bin number)
#remove gwas variants from selection but retain info on which bin they came from
#pick random variants for null analysis. 
genes <- fread("TCSC/simulation_analysis/1000genes_Nov.txt", header = T)
s <- sapply(1:nrow(snps), function(x) min(abs(snps$V4[x]-genes$START))) #distance to closest tss 
#lets do 50 bins
binorder <- sort(rep(c(1:50),20))
snps <- cbind(snps,s)
snps <- snps[order(as.numeric(snps$s),decreasing = F),]
snps <- cbind(snps,binorder)
#return to original order: 
snps <- snps[order(as.numeric(snps$V4),decreasing = F),]
#double check 
snpscheck <- fread("TCSC/simulation_analysis/1KG_HM3_chr1.bim", header = F)
all(snpscheck$V4 == snps$V4) #T

for (sim in 1:sims){
    #load summary statistics, and see which have p < 0.05/(1e6*0.08)
    gwas <- fread(paste0("SEG/sumstats/sumstats_",ngwas,"_vgene",ve,"_sim",sim,".txt.gz"), header = T) 
    p <- pnorm(q = abs(gwas$Z), lower.tail = F)
    w <- which(p < 5e-8)
    if(length(w) > 0){ 
        filename <- paste0("TCSC/simulation_analysis/output/Ongen_files/gwas_snps/Nov_sim",sim,"_ve",ve,".txt") #need to make this output directory somewhere. 
        write.table(w, file = filename, row.names = F, col.names = F, sep = "\t", quote = F)
        system(paste0("gzip -f ",filename))
   
        #allnonsigvar <- c(1:nrow(gwas))[-w]
        possiblesnps <- snps[-w,] #snps[allnonsigvar,]
        #nullvariants <- sapply(1:length(w), function(x) which.min(abs(w[x] - allnonsigvar))) #null variants were next to other variants. 
        #do random null variants, matched based on distance to TSS and MAF (doesn't matter for us) 
        s <- sapply(1:length(w), function(x) sample(possiblesnps$V2[possiblesnps$binorder == snps$binorder[w[x]]], size = 1))
        nullvariants <- match(s,snps$V2)
        filename <- paste0("TCSC/simulation_analysis/output/Ongen_files/gwas_null_snps/Nov_sim",sim,"_ve",ve,".txt")
        write.table(nullvariants, file = filename, row.names = F, col.names = F, sep = "\t", quote = F)
        system(paste0("gzip -f ",filename))
    }
}

  
