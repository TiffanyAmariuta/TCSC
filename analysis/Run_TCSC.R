trait <- commandArgs(trailingOnly=TRUE)

library(data.table)
library(Hmisc)

traits <- fread("TCSC/sumstats/alltraits.txt", header= F)$V1
N <- traits$N[match(trait,traits$Trait_Identifier)]

y <- fread("TCSC/analysis/TissueGroups.txt", header = T)
tissues <- unique(y$MetaTissue)

small_tissues <- c(3,5,11,12,14,22,24,27:30,33:35,37:38) #tissues with gene expression sample size < 320
normal_tissues <- c(1:length(tissues))[-small_tissues] #tissues with gene expression sample size = 320

get_cov_alpha1alpha1_multitissue <- function(alpha_z,CoRegMat,nGWAS,weights1,weights2,weights3){ 
y <- (alpha_z^2)/nGWAS
w1 <- 1/weights1
w2 <- 1/weights2
w3 <- 1/weights3
mod <- summary(glm(y ~ CoRegMat, weights = w1*w2*w3)) #to constrain intercept, use summary(glm(I(y+1) ~ CoRegMat, weights = w1*w2*w3))
cov_b1b1 <- coef(mod)[-1,1]  
cov_b1b1 <- cov_b1b1[1:ncol(CoRegMat)]
return(cov_b1b1)
}

load("TCSC/analysis/InputCoreg_TCSC.RData")
gtex <- fread("TCSC/analysis/gene_annotation.txt.gz", header = F, sep = "\t")

#### trait specific analysis #### 
alpha_z <- c()
for (i in 1:length(tissues)){
if(i %in% small_tissues){
transcript_key <- fread(paste0("TCSC/weights/heritablegenes/Nall/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1
keep <- fread(paste0("TCSC/weights/heritablegenes/Nall/TranscriptsIn",tissues[i],"Model_keep.txt"), header = F)$V1
z <- fread(paste0("TCSC/twas_statistics/allEUR_GTEx/",trait,"/Marginal_alphas_",trait,"_",tissues[i],".txt.gz"), header = F)$V1
w <- which(transcript_key %in% keep) 
transcript_key <- transcript_key[w] 
z <- z[w]
m <- match(transcript_key,gtex$V7)
genetype <- gtex$V8[m]
w <- which(genetype == "protein_coding")
z <- z[w]
alpha_z <- c(alpha_z,z)
}else{
transcript_key <- fread(paste0("TCSC/weights/heritablegenes/N320/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1
keep <- fread(paste0("TCSC/weights/heritablegenes/N320/TranscriptsIn",tissues[i],"Model_keep.txt"), header = F)$V1 
z <- fread(paste0("TCSC/twas_statistics/320EUR_GTEx/",trait,"/Marginal_alphas_",trait,"_",tissues[i],".txt.gz"), header = F)$V1
w <- which(transcript_key %in% keep)
transcript_key <- transcript_key[w]
z <- z[w]
m <- match(transcript_key,gtex$V7)
genetype <- gtex$V8[m]
w <- which(genetype == "protein_coding")
z <- z[w]
alpha_z <- c(alpha_z,z)
}
}
alpha_z <- alpha_z[-remove_genes]
w <- which(alpha_z^2 > max(80,0.001*N) | is.na(alpha_z)) #trait specific qc 
if(length(w)>0){
alpha_z <- alpha_z[-w]
transcripts <- transcripts[-w]
starts <- starts[-w]
chrs <- chrs[-w]
tissueassign <- tissueassign[-w]
X <- X[-w,]
}
N_tissuespecific <- as.numeric(table(tissueassign))
a_transcripts <- table(transcripts)
weights3 <- sapply(1:length(transcripts), function(x) as.numeric(a_transcripts)[match(transcripts[x], names(a_transcripts))])
weights2 <- sapply(1:nrow(X), function(x) sum(X[x,-tissueassign[x]])) + 1
expected_true_cish2_genes <- N_tissuespecific

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

mean_chisq <- mean(alpha_z^2, na.rm = T)
totcoreg <- rowSums(X)
crude_h2_est <- (mean_chisq - 1)/(N*mean(totcoreg)) #should this be average co-reg score or average sum across tissues co-reg score? 
weights1 <- (1 + N*crude_h2_est*totcoreg)^2 #bias corrected but not scaled X 

variance_mat <- matrix(0,1+length(tissues),6) #estimate, jackknife SE, P
variance_mat[,1] <- c(tissues,"AllTissues")
variance_mat[1:length(tissues),2] <- get_cov_alpha1alpha1_multitissue(alpha_z,X,N,weights1,weights2,weights3) * expected_true_cish2_genes
variance_mat[(length(tissues)+1),2] <- sum(as.numeric(variance_mat[1:length(tissues),2]))

jk <- matrix(0,nrow = chunks,ncol = length(tissues))
jk_weights <- matrix(0,nrow = chunks,ncol = 1)
jk_sum <- matrix(0,nrow = chunks,ncol = 1)
for (chunk in 1:chunks){
print(chunk)
remove_genes <- which(a[,5] == chunk)
print(length(remove_genes))
tab <- table(a[remove_genes,4]) #freq of tissues
subtract_genes <- rep(0,length(tissues))
m <- match(1:length(tissues), names(tab))
w <- which(!is.na(m))
subtract_genes[w] <- as.numeric(tab)[m[w]]
N_tissuespecific_jk <- sapply(1:length(N_tissuespecific), function(x) N_tissuespecific[x] - subtract_genes[x])
alpha_z_jk <- alpha_z[-remove_genes]
mean_chisq <- mean(alpha_z_jk^2, na.rm = T)
X_jk <- X[-remove_genes,]
totcoreg <- rowSums(X_jk)
crude_h2_est <- (mean_chisq - 1)/(N*mean(totcoreg))
weights1_jk <- (1 + N*crude_h2_est*totcoreg)^2
weights2_jk <- weights2[-remove_genes] 
weights3_jk <- weights3[-remove_genes] 
jk[chunk,] <- get_cov_alpha1alpha1_multitissue(alpha_z_jk,X_jk,N,weights1_jk,weights2_jk,weights3_jk)* N_tissuespecific_jk 
jk_sum[chunk,1] <- sum(as.numeric(jk[chunk,]))
print(jk[chunk,1])
jk_weights[chunk,1] <- sum(1/(weights1_jk*weights2_jk*weights3_jk))
} 

variance_mat[1:ncol(jk),3] <- sapply(1:ncol(jk), function(x) sqrt(wtd.var(jk[,x], jk_weights[,1]))*sqrt(length(which(jk_weights[,1] != 0)))) 
variance_mat[nrow(variance_mat),3] <- sqrt(wtd.var(jk_sum[,1],jk_weights[,1]))*sqrt(length(which(jk_weights[,1] != 0)))

variance_mat[,4] <- pnorm(q = 0, mean = as.numeric(variance_mat[,2]), sd = as.numeric(variance_mat[,3]))
variance_mat[1:ncol(jk),5] <- p.adjust(as.numeric(variance_mat[1:ncol(jk),4]), method = "fdr")
variance_mat[nrow(variance_mat),5] <- NA
variance_mat[,6] <- as.numeric(variance_mat[,2])/sum(as.numeric(variance_mat[1:ncol(jk),2]))

colnames(variance_mat) <- c("Tissue","Variance","JK_SE","P","FDRP","Variance_divideby_Total")
write.table(variance_mat, file = paste0("TCSC_",trait,".txt"), row.names = F, col.names = T, sep = "\t", quote = F)
