#Features
#1. neighbor genes share 3 eqtls and have same cish2 parameter (doesn't depend on effect sizes), what if calculate h2 like before, but introduce constraint. 
#2. 1/2 genes are not cish2, but it is not random; only in pairs. 
#3. cish2 drawn from exponential distribution 
#4. neighbor genes share eqtl effect sizes too. 

library(data.table)
library(MASS)


nsims <- 1000
set.seed(12345)

y <- fread("TCSC/simulation_analysis/1000genes_Nov.txt", header = T)
snps <- fread("TCSC/simulation_analysis/1KG_HM3_chr1.bim", header = T)
#FIND NEIGHBOR GENES WITH OVERLAPPING CIS WINDOW WITH > 5 SNPS
res <- matrix(0,nrow(y)-1,5)
for (i in c(1:(nrow(y)-1))){
#which snps are in the cis regions of both genes? 
cis_region_i <- c(y$START[i]-500000,y$STOP[i]+500000)
cis_region_iplusone <- c(y$START[i+1]-500000,y$STOP[i+1]+500000)
w <- intersect(which(snps$V4 > cis_region_i[1] & snps$V4 < cis_region_i[2]),which(snps$V4 > cis_region_iplusone[1] & snps$V4 < cis_region_iplusone[2]))
if(length(w) > 5){ #end of the first gene minus the start of the second gene (overlapping space)
#pick 5 snps that are in this window. 
three_eQTL <- sample(w, size = 3, replace = F)
res[i,] <- c(i,i+1,three_eQTL)
}else{res[i,] <- c(i,i+1,rep(NA,3))}
}
#don't want gene to appear more than once in res, so remove duplicate occurences, since can have gene to the left or gene to the right. can't have both. 
res2 <- res
j <- 1
while(j <= nrow(res)){
m1 <- match(j, res2[,1])
m2 <- match(j, res2[,2])
m <- c(m1,m2)
m <- m[which(!is.na(m))]
if(length(m) > 1){
res2 <- res2[-m[1],]
}
j <- j + 1
}
#remove rows with NAs, because later we will search for gene index in this matrix. 
w <- which(is.na(res2[,3]))
res2 <- res2[-w,]
write.table(res2[,c(1:2)], file = "TCSC/simulation_analysis/Nov_coreg_genepairs.txt", row.names=F, col.names = F, sep = "\t", quote = F)

#ASSIGN NEIGHBOR CO-REG GENES THE SAME H2 SO THEY HAVE A CHANCE AT BEING DETECTED AS CISH2
#if don't do this intentionally, cish2 genes in sims are never in pairs so no co-reg. 
maxh2 <- 0.5 #was 0.6 when most things looked good, dec 12 check
mihh2 <- 0.01
cish2_exp <- rexp(n = nrow(res2), rate = 1/0.15)
cish2_exp[cish2_exp < mihh2] <- mihh2
cish2_exp[cish2_exp > maxh2] <- maxh2 #makes gcta faster (slower when h2 closer to 1)
h2mat <- cbind(res2[,c(1:2)],cish2_exp)
h2_pergene <- rep(0,ngenes)
m <- sapply(1:ngenes, function(x) which(h2mat[,c(1,2)] == x, arr.ind = T)[1])
w <- which(!is.na(m))
m <- m[w]
h2_pergene[w] <- h2mat[m,3]
w <- which(h2_pergene == 0) #ones that are not in a pair
cish2_exp <- rexp(n = length(w), rate = 1/0.15)
cish2_exp[cish2_exp < mihh2] <- mihh2
cish2_exp[cish2_exp > maxh2] <- maxh2 #mean = 0.15
h2_pergene[w] <- cish2_exp
write.table(h2_pergene, file = "TCSC/simulation_analysis/Nov_cish2_pergene_toachievebettercoreg.txt", row.names=F, col.names = F, sep = "\t", quote = F)

#SET UP COVARIANCE OF EQTLS ACROSS TISSUES
ntiss <- 10
nvar <- 5
ngenes <- nrow(y)
modules <- cbind(c(1,1,2,4,4,5,7,7,7,8,8,9),c(2,3,3,5,6,6,8,9,10,9,10,10))
modules2 <- modules[,c(2,1)]
N <- 3
mu <- rep(0,ntiss)
rho_module <- 0.88 #was 0.8; was 0.85 for v2 with cish2 <- 0.025; value among 3, value among next 3, 0.8 among next four. is this the same as it was before? 
rho_outside_module <- 0.82 #was 0.74; was 0.8 for v2
cish2 <- 0.075 
sigma <- matrix(sqrt(cish2*cish2)*rho_outside_module,ntiss,ntiss)
diag(sigma) <- cish2
sigma[modules] <- sqrt(cish2*cish2)*rho_module
sigma[modules2] <- sqrt(cish2*cish2)*rho_module
not_modules <- which(sigma == sqrt(cish2*cish2)*rho_outside_module, arr.ind = T) #has duplicates in pairs 66 instead of 33 
for (i in 1:nrow(not_modules)){
mymin <- min(not_modules[i,])
mymax <- max(not_modules[i,])
not_modules[i,] <- c(mymin,mymax)
}
not_modules <- unique(not_modules) #now 33 

#for every gene pick 3 snps within 50kb of gene body to be the same for gene i across all tissues, note position
#3/5 snps shared between each of 10 tissues (no tissue-specific eQTLs) opportunity to change to make sims more realistic for RTC Coloc. Should be fine. Already some tissue-specific ones. 
#for every gene pick a unique 2 snps within 1 Mb of gene start to be differnet for gene i across all tissues, note position. 
#if gene is 2nd in a pair, just use those eqtl effect sizes for all 5 variants. (3 will have pos, 2 will not) need co-reg btwn genes to be explicit, position only doesn't work.
summary_mat <- matrix(0,nsims,5)
for (sim in 1:nsims){
  print(sim)
  rg_mod <- c()
  rg_out <- c()
 
  #CHOOSE WHICH GENES ARE AND ARE NOT CISH2 IN THIS SIM
  #select 1/2 pairs of genes to be not cish2 to better match % of detected cish2 from total in gtex. 
  immune_genes <- as.matrix(fread("TCSC/simulation_analysis/Nov_Y1_alpha_nCausalGenes100.txt.gz", header = F))
  immune_genes <- which(immune_genes[,sim]!=0) #use for causal tissue and tissues in module. Otherwise instead of a random half of genes not being cis h2, 
  gene_cish2_pertissue_matrix <- matrix(0,ngenes,ntiss)
  gene_cish2_pertissue_matrix_key <- matrix(1,ngenes,ntiss)
  s <- sapply(1:length(immune_genes), function(x) which(res2[,c(1,2)] == immune_genes[x], arr.ind = T)[1]) #protected pairs
  s <- s[which(!is.na(s))]
  s_zero <- sample(c(1:nrow(res2))[-s], size = 0.5*nrow(res2), replace = F)
  #intersect(s_zero,s) #check that 0 overlap 
  s_zero_genes <- unique(res2[s_zero,c(1,2)])
  #s <- sample(c(1:ngenes)[-immune_genes],size = 0.5*ngenes, replace = F)
  gene_cish2_pertissue_matrix_key[s_zero_genes,c(1:3)] <- 0 #causal tissue + module of related tissues
  for (ti in 4:ntiss){ #iterate over tissues that are not in module with causal tissue. 
  s_zero <- sample(c(1:nrow(res2))[-s], size = 0.5*nrow(res2), replace = F)
  s_zero_genes <- unique(res2[s_zero,c(1,2)])
  #s <- sample(c(1:ngenes),size = 0.5*ngenes, replace = F)
  gene_cish2_pertissue_matrix_key[s_zero_genes,ti] <- 0
  }
  
  betas <- matrix(0,0,ntiss)
  for (gene in 1:ngenes){ 
    #draw from a mvn with number of tissues corresponding to tissues for which the gene is actually going to be heritable? 
    #or can we just zero out the matrix w/o much of a problem? 
    w <- which(res2[,c(1:2)] == gene,arr.ind = T)[2]
    if(is.na(w)){w <- 0}
    if(w != 2){ #if gene is first in a pair or is not in the list 
    s <- 0 #well now we want average h2 to be 0.1, so change this to like 0.01 (changed March 9 from 0.1 to 0.01 below)
    while(s < 0.01){ #min h2 allowed is 0.1 in causal tissue because if gene is causal and it has extremely low h2, it will never be detectable. 
    mvn1 <- mvrnorm(N, mu = mu, Sigma = sigma )
    uvn1 <- matrix(rnorm(n = (nvar-N)*ntiss, mean = 0, sd = sqrt(cish2)),nrow = nvar-N, ncol = ntiss)
    gene_betas <- rbind(mvn1,uvn1)
    gene_cish2_pertissue_matrix[gene,] <- sapply(1:ncol(gene_betas), function(x) sum(gene_betas[,x]^2)) #want these numbers to be approx the same. 
    s <- min(gene_cish2_pertissue_matrix[gene,]) #if any cish2 for a tissue is < 0.01, rerun. surprisingly rare, so this is fast. 
    }
    } #ifelse, still want to assess if gene needs to be zeroed out.  but want same gene_betas matrix. 
    zerooutgene <- which(gene_cish2_pertissue_matrix_key[gene,] == 0)
    if(length(zerooutgene) > 0){gene_betas[,zerooutgene] <- 0} #tissues on columns 
    betas <- rbind(betas, gene_betas)
    rg_mod <- c(rg_mod, mean(sapply(1:nrow(modules), function(x) cor(mvn1[,modules[x,1]],mvn1[,modules[x,2]]))))
    rg_out <- c(rg_out, mean(sapply(1:nrow(not_modules), function(x) cor(mvn1[,not_modules[x,1]],mvn1[,not_modules[x,2]])))) #need not_modules here as is; e.g. with dup across diagonal.
  }
  
  summary_mat[sim,1] <- mean(gene_cish2_pertissue_matrix[,1]) #0.1 is target  
  summary_mat[sim,2] <- mean(gene_cish2_pertissue_matrix[,-1]) #0.1 is target 
  summary_mat[sim,3] <- mean(rg_mod) #want around 0.8
  summary_mat[sim,4] <- mean(rg_out) #want around 0.74 #this is unaffected by mistake about 66 instead of 33 pairs. because it's the mean.  
  summary_mat[sim,5] <- (mean(rg_mod)*nrow(modules) + mean(rg_out)*nrow(not_modules)/2)/(nrow(not_modules)/2 + nrow(modules)) #0.75 (this one is a bit affected), recompute this column. 
  print(summary_mat[sim,])
  
  #PICK WHICH SNPS WILL BE CAUSAL FOR EACH GENE
  #col1: chunk-specific gene number (out of total genes in chunk)
  #col2: chunk number for variant/gene pair 
  #col3: variant position for non-zero effects (out of total variants in 1000G chr1) 
  position_matrix <- matrix(0,nvar*nrow(y),ntiss) #want a different one for each simulation. 
  #pick 60% of genes to actually have their eqtl defined in a way that creates gene co-reg. tune this param. dec 12 check
  s <- sample(nrow(res2), size = ceiling(0.35*nrow(res2)), replace = F) #changed to smaller number to decrease co-reg
  res_ds <- res2[s,]
  #length(unique(res_ds[,c(1,2)]))/length(unique(res2[,c(1,2)])) #60%
  #check if j is found in s, 
  snp_border_up <- c()
  snp_border_down <- c()
  count <- 0
  for (j in 1:nrow(y)){
  m1 <- match(j, res_ds[,1])
  m2 <- match(j, res_ds[,2])
  m <- c(m1,m2)
  m <- m[which(!is.na(m))]
  w <- which(snps$V4 > y$START[j] - 500000 & snps$V4 < y$STOP[j] + 500000)
  if(length(m) > 0){ #snps decided
  w_shared <- res_ds[m,c(3,4,5)]}else{
  w_shared <- sample(w,size = 3, replace = F)
  }
  position_matrix[((count+1):(count + nvar))[1:3],1:ntiss] <- matrix(rep(w_shared,ntiss),nrow = 3)
  w_leftover <- w[-which(w %in% w_shared)]
  for (m in 1:ntiss){
  w_unique <- sample(w_leftover,size=2,replace = F)
  position_matrix[((count+1):(count + nvar))[4:5],m] <- w_unique
  w_leftover <- w_leftover[-which(w_leftover %in% w_unique)]
  }
  #w <- which(snps$V4 > y$START[j] - 500000 & snps$V4 < y$STOP[j] + 500000)
  snp_border_up <- c(snp_border_up, rep(min(w),nvar))
  snp_border_down <- c(snp_border_down, rep(max(w),nvar))  
  count <- count + nvar
  } #double check, looks good, position_matrix[5*4+1,] (end of 4th gene + 1 = start of 5th gene, matches position_matrix[5*5+1,] = start of 6th gene

  for (i in 1:ntiss){
    beta_matrix <- cbind(position_matrix[,i],betas[,i],snp_border_up,snp_border_down)
    write.table(beta_matrix, file = paste0("TCSC/simulation_analysis/eQTLs_1000sims_Nov/Nov_sims_eQTLeffectsizes_0.75rg_SNPpolygen",nvar,"_tissue",i,"_sim",sim,".txt"), quote=F, sep = "\t", row.names = F, col.names = F)
    system(paste0("gzip -f TCSC/simulation_analysis/eQTLs_1000sims_Nov/Nov_sims_eQTLeffectsizes_0.75rg_SNPpolygen",nvar,"_tissue",i,"_sim",sim,".txt"))
  }
}







