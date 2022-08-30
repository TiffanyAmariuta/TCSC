#please set working directory with setwd() to your local path to TCSC/

trait <- commandArgs(trailingOnly=TRUE) #format: example: UKB_460K.body_BMIz as in first column in TCSC/sumstats/alltraits.txt

library(data.table)
library(Hmisc)

traits <- fread("sumstats/alltraits.txt", header=T)
m <- match(trait,traits$Trait_Identifier)
N <- traits$N[m]
h2g <- traits$h2g[m]

y <- fread("analysis/TissueGroups.txt", header = T)
tissues <- y$Tissues

get_cov_alpha1alpha1_multitissue <- function(alpha_z,CoRegMat,nGWAS,weights1,weights2,weights3,covar){ 
y <- (alpha_z^2)/nGWAS
w1 <- 1/weights1
w2 <- 1/weights2
w3 <- 1/weights3
mod <- summary(glm(y ~ CoRegMat, weights = w1*w2*w3))
cov_b1b1 <- coef(mod)[-1,1]  
cov_b1b1 <- cov_b1b1[1:ncol(CoRegMat)]
return(cov_b1b1)
}

load("analysis/InputCoreg_BrainSpecific.RData")
gtex <- fread("analysis/gene_annotation.txt.gz", header = F, sep = "\t")

#### trait specific analysis #### 
g <- grep("Brain",tissues)
tissues <- tissues[g]
alpha_z <- c()
for (i in 1:length(tissues)){
transcript_key <- fread(paste0("weights/heritablegenes/Nall/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1
z <- fread(paste0("twas_statistics/allEUR_GTEx/",trait,"/Marginal_alphas_",trait,"_",tissues[i],".txt.gz"), header = F)$V1
m <- match(transcript_key,gtex$V7)
genetype <- gtex$V8[m]
alpha_z <- c(alpha_z,z[which(genetype == "protein_coding")])
}
remove_genes1 <- unique(which(is.na(rowSums(X)))) 
remove_genes2 <- unique(which(rowSums(X) == 0))
mhc <- which(chrs == 6 & as.numeric(starts) > 29000000 & as.numeric(starts) < 33000000)
remove_genes <- which(alpha_z^2 > max(80,0.001*N) | is.na(alpha_z))
w <- unique(c(remove_genes1,remove_genes2,mhc,remove_genes))
if(length(w)>0){
alpha_z <- alpha_z[-w]
transcripts <- transcripts[-w]
starts <- starts[-w]
chrs <- chrs[-w]
tissueassign <- tissueassign[-w] #add this to save
X <- X[-w,]
}
N_tissuespecific <- as.numeric(table(tissueassign))
a_transcripts <- table(transcripts)
weights3 <- sapply(1:length(transcripts), function(x) as.numeric(a_transcripts)[match(transcripts[x], names(a_transcripts))])  #number of tissues in which each gene is cis heritable
weights2 <- sapply(1:nrow(X), function(x) sum(X[x,-tissueassign[x]])) + 1 #these weights are generated using co-regulation scores without bias correction, because we don't want to downweight genes because higher accuracy
expected_true_cish2_genes <- N_tissuespecific #the number of GCTA-detected cis heritable genes is an estimate of the true number of cis heritable genes at infinite sample size

#### set up jackknife #### 
a <- cbind(1:length(starts),as.numeric(chrs),as.numeric(starts),unlist(sapply(1:length(N_tissuespecific), function(x) rep(x,N_tissuespecific[x]))))
a <- a[order(a[,2],a[,3],decreasing = F),]
chunks <- 200
size_groups <- floor(length(starts)/chunks)
size_group5 <- length(starts) - (chunks-1)*size_groups
group_assignment <- c()
for (j in 1:chunks){
if(j == chunks){
group_assignment <- c(group_assignment, rep(j,size_group5))
}else{
group_assignment <- c(group_assignment, rep(j,size_groups))
}
}
a <- cbind(a, group_assignment)
a <- a[order(a[,1], decreasing = F),]
a <- cbind(a,transcripts)

mean_chisq <- mean(alpha_z^2)
mean_coreg <- mean(rowSums(X)) 
crude_h2_est <- (mean_chisq - 1)/(N*mean_coreg)
weights1 <- (1 + N*crude_h2_est*rowSums(X))^2
variance_mat <- matrix(0,1+length(tissues),5) #tissue, h2_ge_t, jackknife SE, P, fdr across tissues 
variance_mat[,1] <- c(tissues,"AllTissues")
variance_mat[1:length(tissues),2] <- get_cov_alpha1alpha1_multitissue(alpha_z,X,N,weights1,weights2,weights3,covar) * expected_true_cish2_genes / h2g
variance_mat[(length(tissues)+1),2] <- sum(as.numeric(variance_mat[1:length(tissues),2]))

chunks <- 200
jk <- matrix(0,nrow = chunks,ncol = length(tissues))
jk_weights <- matrix(0,nrow = chunks,ncol = 1)
jk_sum <- matrix(0,nrow = chunks,ncol = 1)
for (chunk in 1:chunks){
print(paste0("processing jackknife chunk ",chunk," out of 200"))
remove_genes <- which(a[,5] == chunk)
tab <- table(a[remove_genes,4]) #freq of tissues
subtract_genes <- rep(0,length(tissues))
m <- match(1:length(tissues), names(tab))
w <- which(!is.na(m))
subtract_genes[w] <- as.numeric(tab)[m[w]]
N_tissuespecific_jk <- sapply(1:length(N_tissuespecific), function(x) N_tissuespecific[x] - subtract_genes[x])
alpha_z_jk <- alpha_z[-remove_genes]
covar <- rep(1,length(alpha_z_jk))
mean_chisq <- mean(alpha_z_jk^2)
X_jk <- X[-remove_genes,]
mean_coreg <- mean(rowSums(X_jk)) 
crude_h2_est <- (mean_chisq - 1)/(N*mean_coreg)
weights1_jk <- (1 + N*crude_h2_est*rowSums(X_jk))^2 
weights2_jk <- weights2[-remove_genes] 
weights3_jk <- weights3[-remove_genes] #upweights tissue specific genes
jk[chunk,] <- get_cov_alpha1alpha1_multitissue(alpha_z_jk,X_jk,N,weights1_jk,weights2_jk,weights3_jk,covar) * N_tissuespecific_jk / h2g
jk_sum[chunk,1] <- sum(as.numeric(jk[chunk,]))
jk_weights[chunk,1] <- sum(1/(weights1_jk*weights2_jk*weights3_jk))
} 

variance_mat[1:ncol(jk),3] <- sapply(1:ncol(jk), function(x) sqrt(wtd.var(jk[,x], jk_weights[,1]))*sqrt(length(which(jk_weights[,1] != 0)))) 
variance_mat[nrow(variance_mat),3] <- sqrt(wtd.var(jk_sum[,1],jk_weights[,1]))*sqrt(length(which(jk_weights[,1] != 0)))

variance_mat[,4] <- pnorm(q = 0, mean = as.numeric(variance_mat[,2]), sd = as.numeric(variance_mat[,3])) 
variance_mat[1:ncol(jk),5] <- p.adjust(as.numeric(variance_mat[1:ncol(jk),4]), method = "fdr")
variance_mat[nrow(variance_mat),5] <- NA

colnames(variance_mat) <- c("Tissue","h2ge_t","JK_SE","P","FDRP")
write.table(variance_mat, file = paste0("results/TCSC_",trait,"_BrainSpecificAnalysis.txt"), row.names = F, col.names = T, sep = "\t", quote = F)

