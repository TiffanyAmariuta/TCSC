
ConditionalAnalysis <- function(QTL_geno, tot_GE){
pmin <- 0
QTL_geno_iter <- QTL_geno
colnames(QTL_geno_iter) <- 1:ncol(QTL_geno)
topSNP <- NULL
condition <- NULL
thresh <- 1.39911e-06 #0.05 #changed from 0.01 for Nov 
independent_signals <- c()
#FORWARD STAGE
while(pmin < thresh & is.matrix(QTL_geno_iter)){
eqtl_res <- matrix(0,ncol(QTL_geno_iter),1)
for (i in 1:ncol(QTL_geno_iter)){ #cis snps 
if(is.null(topSNP)){
mod <- summary(lm(tot_GE ~ QTL_geno_iter[,i]))
}else{mod <- summary(lm(tot_GE ~ QTL_geno_iter[,i] + condition))}
eqtl_res[i,1] <- coef(mod)[2,4]
}
topSNP <- which.min(eqtl_res) #5
independent_signals <- c(independent_signals, colnames(QTL_geno_iter)[topSNP])
#then condition next round on this one 
pmin <- min(eqtl_res[,1])
w <- which(eqtl_res[,1] < thresh & eqtl_res[,1] > pmin) #snps that are also significant, test for 2ndary association. 
condition <- cbind(condition,QTL_geno_iter[,topSNP])
if(length(w) > 0){QTL_geno_iter <- QTL_geno_iter[,w]}else{QTL_geno_iter <- rep(0,nrow(QTL_geno))} #could be any number. 
}
#BACKWARD STAGE only if there are multiple signals: finally test each independent signal against all other independent signals. Remove if no longer independent 
independent_signals <- as.numeric(independent_signals)
return(independent_signals)
}

MostSigeQTL <- function(QTL_geno, tot_GE){
eqtl_res <- matrix(0,ncol(QTL_geno),1)
for (i in 1:ncol(QTL_geno)){
mod <- summary(lm(tot_GE ~ QTL_geno[,i]))
eqtl_res[i,1] <- coef(mod)[2,4]
}
topSNP <- which.min(eqtl_res[,1])
return(topSNP)
}

RTC_old <- function(tot_GE, independent_signals, QTL_geno, chunkspecific_SNPindex_forGWASvariant, is){
if(length(independent_signals) > 1){
starting_pheno <- resid(lm(tot_GE ~ QTL_geno[,independent_signals[-is]]))}else{
starting_pheno <- tot_GE
}
pseudo_pheno <- matrix(0,nrow(QTL_geno),ncol(QTL_geno))  #people x snps (in chunk) 300 x 238 for example
for (i in 1:ncol(QTL_geno)){ #correcting for each of the variants BUT ALSO the other independent eQTLs. 
pseudo_pheno[,i] <- resid(lm(starting_pheno ~ QTL_geno[,i]))
}
eqtl_res <- matrix(0,ncol(QTL_geno),1)
for (i in 1:ncol(QTL_geno)){
mod <- summary(lm(pseudo_pheno[,i] ~ QTL_geno[,independent_signals[is]])) #all N pseudo phenos against the same variant. 
eqtl_res[i,] <- coef(mod)[-1,4]
}
#need to find rank of each. 
if(all(eqtl_res[,1] == 1)){rank <- nrow(eqtl_res)}else{
snp_ID <- rep("qtl",ncol(QTL_geno)) #has to be all cis variants because gwas snp is not necessarily exactly the same snp as the eQTL, therefore need a p value for every cis window variant. 
snp_ID[chunkspecific_SNPindex_forGWASvariant] <- "gwas_snp"
a <- cbind(eqtl_res,snp_ID)
a <- a[order(as.numeric(a[,1]),decreasing = T),]
rank <- which(a[,2] == "gwas_snp")}
rtc <- (ncol(QTL_geno) - rank)/ncol(QTL_geno)
mod <- summary(lm(tot_GE ~ QTL_geno[,independent_signals[is]]))
betas <- as.numeric(coef(mod)[,1])
return(list = c(rtc, betas))
}

RTC <- function(tot_GE_gene, independent_signals, QTL_geno, chunkspecific_SNPindex_forGWASvariant, is){
#if more than one eqtl, regress out effects of all other eqtl and do rtc separately.
if(length(independent_signals) > 1){
    starting_pheno <- resid(lm(tot_GE_gene ~ QTL_geno[,independent_signals[-is]]))
    tot_GE_gene <- starting_pheno
}

pseudo_pheno <- sapply(1:ncol(QTL_geno), function(x) resid(lm(tot_GE_gene ~ QTL_geno[,x])))
eqtl_res <- sapply(1:ncol(pseudo_pheno), function(x) cor.test(pseudo_pheno[,x],QTL_geno[,independent_signals[is]])$p.val)
a <- cbind(eqtl_res,1:ncol(pseudo_pheno))
a <- a[order(as.numeric(a[,1]),decreasing = T),]
rank <- which(a[,2] == chunkspecific_SNPindex_forGWASvariant) #rank
rtc <- (ncol(QTL_geno) - rank)/ncol(QTL_geno)
mod <- summary(lm(tot_GE_gene ~ QTL_geno[,independent_signals[is]]))
betas <- as.numeric(coef(mod)[,1])
return(list = c(rtc, betas))
}


RTC_null <- function(QTL_geno){
rtc <- c()
for (tr in 1:200){ #maybe we should only pick causal, rather than linked, b/c our variants don't have linkage. 
print(tr)
twovar_c <- sample(ncol(QTL_geno), size = 2, replace = F) #convention will be gwas is first, eqtl is second
beta_eqtlc <- 3
int_eqtlc <- 1
totexp_c <- QTL_geno[,twovar_c[2]]*beta_eqtlc + int_eqtlc #gene expression phenotype that has nothing to do with gwas snp. 
pseudo_pheno <- sapply(1:ncol(QTL_geno), function(x) resid(lm(totexp_c ~ QTL_geno[,x]))) #people x snps (in chunk) 300 x 238 for example, correcting for each of the variants BUT ALSO the other independent eQTLs. 
eqtl_res <- sapply(1:ncol(QTL_geno), function(x) cor.test(pseudo_pheno[,x],QTL_geno[,twovar_c[2]])$p.val) #check that p=1 when x = twovar_c[2], yes. 
snp_ID <- rep("qtl",ncol(QTL_geno)) #has to be all cis variants because gwas snp is not necessarily exactly the same snp as the eQTL, therefore need a p value for every cis window variant. 
snp_ID[twovar_c[1]] <- "gwas_snp" #the gwas proxy variant 
a <- cbind(eqtl_res,snp_ID)
a <- a[order(as.numeric(a[,1]),decreasing = T),]
rank <- which(a[,2] == "gwas_snp")
rtc <- c(rtc, (ncol(QTL_geno) - rank)/ncol(QTL_geno))
}
return(rtc) #since there's no linkage, the rank might be predictable, e.g. always in the middle. whats the mean rank, what's the sd, maybe we can make a fake distribution from this? since we have no genetic architecture in our sims as far as linkage.  
}

RTC_null_nolinkage <- function(QTL_geno){
rtc <- rnorm(mean = 0.5, sd = 0.3, n = 200) 
rtc[rtc < 0] <- 0.5
return(rtc)
}

RTC_alt <- function(QTL_geno){ #rank = 1 all the time, b/c we have no linkage, therefore need to select the actual variant itself. or can say in perfect linkage. therefore rtc = 0.9964664, depends on ncol QTL_geno so actually we can leave this. 
rtc <- c()
for (tr in 1:200){
print(tr)
onevar_c <- sample(ncol(QTL_geno), size = 1) #convention will be gwas is first, eqtl is second
beta_eqtlc <- 3
int_eqtlc <- 1
totexp_c <- QTL_geno[,onevar_c]*beta_eqtlc + int_eqtlc
pseudo_pheno <- sapply(1:ncol(QTL_geno), function(x) resid(lm(totexp_c ~ QTL_geno[,x])))
eqtl_res <- sapply(1:ncol(QTL_geno), function(x) cor.test(pseudo_pheno[,x],QTL_geno[,onevar_c])$p.val)
snp_ID <- rep("qtl",ncol(QTL_geno)) #has to be all cis variants because gwas snp is not necessarily exactly the same snp as the eQTL, therefore need a p value for every cis window variant. 
snp_ID[onevar_c] <- "gwas_snp" #the gwas proxy variant 
a <- cbind(eqtl_res,snp_ID)
a <- a[order(as.numeric(a[,1]),decreasing = T),]
rank <- which(a[,2] == "gwas_snp")
rtc <- c(rtc, (ncol(QTL_geno) - rank)/ncol(QTL_geno))
}
return(rtc)
}

RTC_alt_nolinkage <- function(QTL_geno){
rank <- sample(1:10, size = 200, replace = T)   
rtc <- (ncol(QTL_geno) - rank)/ncol(QTL_geno) #didn't change anything, ok. 
return(rtc)
}


