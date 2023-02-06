trait1_trait2 = commandArgs(trailingOnly=TRUE) #format: example: UKB_460K.body_BMIz,UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED 
s <- strsplit(trait1_trait2,split = ",")[[1]]
trait1 <- s[1]
trait2 <- s[2]

library(data.table)
library(Hmisc)

traits <- fread("TCSC/sumstats/alltraits.txt", header=T)
N1 <- traits$N[match(trait1,traits$Trait_Identifier)]
N2 <- traits$N[match(trait2,traits$Trait_Identifier)]

#read in covariance estimates for each trait pair 
covariance <- fread("TCSC/sumstats/AllPairs_Cov.txt", header=F, sep = "\t")
m1 <- match(paste0(trait1,",",trait2),covariance$V1)
m2 <- match(paste0(trait2,",",trait1),covariance$V1)
cov_ge <- ifelse(is.na(m1), covariance$V2[m2], covariance$V2[m1])

y <- fread("TCSC/analysis/TissueGroups.txt", header = T) #one row per GTEx tissue analyzed in Amariuta et al 2022 bioRxiv
tissues <- unique(y$MetaTissue) #tissues with eQTL sample size < 320 are grouped together via meta-analysis if the correlation of marginal eQTL effects is > 0.93. 

n_eqtl <- sapply(1:length(tissues), function(x) sum(y$N_EUR[which(y$MetaTissue == tissues[x])])) #find total available sample size from GTEx data after accounting for meta-analysis of tissues
small_tissues <- which(n_eqtl < 320)
normal_tissues <- c(1:length(tissues))[-small_tissues] #in primary analysis, these tissues are subsampled such that eQTL sample size = 320

get_cov_alpha1alpha1_multitissue <- function(alpha_z,X,nGWAS,weights1,weights2,weights3){ #single trait
y <- (alpha_z^2)/nGWAS
w1 <- 1/weights1
w2 <- 1/weights2
w3 <- 1/weights3
mod <- summary(glm(y ~ X, weights = w1*w2*w3))
cov_b1b1 <- coef(mod)[-1,1]  
cov_b1b1 <- cov_b1b1[1:length(tissues)]
return(cov_b1b1)
}

get_cov_alpha1alpha2_multitissue <- function(alpha_z1,alpha_z2,X,N1,N2,weights3,weights2){ #cross-trait
y <- (alpha_z1*alpha_z2)/(sqrt(N1)*sqrt(N2))
mean_chisq1 <- mean(alpha_z1^2, na.rm = T)
mean_chisq2 <- mean(alpha_z2^2, na.rm = T)
mean_coreg <- mean(rowSums(X), na.rm =T)
crude_h2_est_1 <- (mean_chisq1 - 1)/(N1*mean_coreg)
crude_h2_est_2 <- (mean_chisq2 - 1)/(N2*mean_coreg)

myweights <- weights_fun(alpha_z1,alpha_z2,X,N1,N2,weights3,weights2,crude_h2_est_1,crude_h2_est_2)
mod <- summary(glm(y ~ X, weights = myweights)) 
cov_a1a2 <- coef(mod)[-1,1]
cov_a1a2 <- cov_a1a2[1:length(tissues)] #then multiply by N_tissuespecific
return(list("covariance" = cov_a1a2, "weights" = myweights))
}

weights_fun <- function(alpha_z1,alpha_z2,X,N1,N2,weights3,weights2,crude_h2_est_1,crude_h2_est_2){
y <- (alpha_z1*alpha_z2)/(sqrt(N1)*sqrt(N2))
rho_g_est <- sum(alpha_z1*alpha_z2, na.rm = T)/(sqrt(N1)*sqrt(N2))
rs <- rowSums(X)
het_prod1 <- 1 + N1*crude_h2_est_1*rs 
het_prod2 <- 1 + N2*crude_h2_est_2*rs
het_prod3 <- sqrt(N1)*sqrt(N2)*rho_g_est*rs/mean(expected_true_cish2_genes)
het_prod4 <- 0
heteroscedasticity <- het_prod1*het_prod2 + (het_prod3 + het_prod4)^2
w1 <- 1/(heteroscedasticity)
w2 <- 1/weights2
w3 <- 1/weights3
mod <- summary(glm(y ~ rs, weights = w1*w2*w3))
intercept <- coef(mod)[1,1] 
rhop_NS <- intercept*sqrt(N1)*sqrt(N2)
rho_g_est2 <- coef(mod)[2,1]
het_prod3 <- sqrt(N1)*sqrt(N2)*rho_g_est2*rs/mean(expected_true_cish2_genes)
het_prod4 <- rhop_NS/(sqrt(N1)*sqrt(N2))
heteroscedasticity <- het_prod1*het_prod2 + (het_prod3 + het_prod4)^2
w1 <- 1/(heteroscedasticity) 
weights <- w1*w2*w3
return(weights)
}

load("TCSC/analysis/InputCoreg_TCSC.RData")
gtex <- fread("TCSC/analysis/gene_annotation.txt.gz", header = F, sep = "\t")

#### trait specific analysis #### 

outliers <- c() #union of both trait's outliers
for (k in 1:2){
alpha_z <- c()
trait <- get(paste0("trait",k))
for (i in 1:length(tissues)){
if(i %in% small_tissues){
transcript_key <- fread(paste0("TCSC/weights/heritablegenes/Nall/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1
keep <- fread(paste0("TCSC/weights/heritablegenes/Nall/TranscriptsIn",tissues[i],"Model_keep.txt"), header = F)$V1
z <- fread(paste0("TCSC/twas_statistics/allEUR_GTEx/",trait,"/Marginal_alphas_",trait,"_",tissues[i],".txt.gz"), header = F)$V1
}else{
transcript_key <- fread(paste0("TCSC/weights/heritablegenes/N320/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1
keep <- fread(paste0("TCSC/weights/heritablegenes/N320/TranscriptsIn",tissues[i],"Model_keep.txt"), header = F)$V1
z <- fread(paste0("TCSC/twas_statistics/320EUR_GTEx/",trait,"/Marginal_alphas_",trait,"_",tissues[i],".txt.gz"), header = F)$V1
}
w <- which(transcript_key %in% keep)
transcript_key <- transcript_key[w]
z <- z[w] 
m <- match(transcript_key,gtex$V7)
genetype <- gtex$V8[m]
alpha_z <- c(alpha_z,z[which(genetype == "protein_coding")])
}
N <- get(paste0("N",k))
alpha_z <- alpha_z[-remove_genes] #need to remove outliers from both traits. take union. 
w <- which(alpha_z^2 > max(80,0.001*N)) #trait specific qc 
if(length(w)>0){outliers <- c(outliers,w)}
assign(paste0("alpha_z",k),alpha_z)
}

w <- unique(outliers)
if(length(w) > 0){
alpha_z1 <- alpha_z1[-w]
alpha_z2 <- alpha_z2[-w]
transcripts <- transcripts[-w]
starts <- starts[-w]
chrs <- chrs[-w]
tissueassign <- tissueassign[-w] #add this to save
X <- X[-w,]
}

N_tissuespecific <- as.numeric(table(tissueassign))
a_transcripts <- table(transcripts)
weights3 <- sapply(1:length(transcripts), function(x) as.numeric(a_transcripts)[match(transcripts[x], names(a_transcripts))])
weights2 <- sapply(1:nrow(X), function(x) sum(X[x,-tissueassign[x]])) + 1
expected_true_cish2_genes <- N_tissuespecific

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
group_assignment <- c(group_assignment, rep(j,size_groups))}}
a <- cbind(a, group_assignment)
a <- a[order(a[,1], decreasing = F),]
a <- cbind(a,transcripts)

variance_mat <- matrix(0,1+length(tissues),7) #tissue, omega_ge_t, jackknife SE, P, FDR across tissues 
variance_mat[,1] <- c(tissues,"AllTissues")
variance_mat[1:length(tissues),2] <- get_cov_alpha1alpha2_multitissue(alpha_z1,alpha_z2,X,N1,N2,weights3,weights2)$covariance * expected_true_cish2_genes 
variance_mat[(length(tissues)+1),2] <- sum(as.numeric(variance_mat[1:length(tissues),2]))

chunks <- 200 
jk <- matrix(0,nrow = chunks,ncol = length(tissues))
jk_sum <- matrix(0,nrow = chunks,ncol = 1)
jk_weights <- matrix(0,nrow = chunks,ncol =1)
for (chunk in 1:chunks){
print(paste0("processing jackknife chunk ",chunk," out of 200"))
remove_genes <- which(a[,5] == chunk)
tab <- table(a[remove_genes,4]) #freq of tissues
subtract_genes <- rep(0,length(tissues))
m <- match(1:length(tissues), names(tab))
w <- which(!is.na(m))
subtract_genes[w] <- as.numeric(tab)[m[w]]
N_tissuespecific_jk <- sapply(1:length(N_tissuespecific), function(x) N_tissuespecific[x] - subtract_genes[x])
bb <- get_cov_alpha1alpha2_multitissue(alpha_z1[-remove_genes],alpha_z2[-remove_genes],X[-remove_genes,],N1,N2,weights3[-remove_genes],weights2[-remove_genes]) 
jk[chunk,] <- bb$covariance * N_tissuespecific_jk 
jk_sum[chunk,1] <- sum(as.numeric(jk[chunk,]))
jk_weights[,1] <- sum(bb$weights) 
} 

variance_mat[1:length(tissues),3] <- sapply(1:length(tissues), function(x) sqrt(wtd.var(jk[,x], jk_weights[,1]))*sqrt(length(which(jk_weights[,1] != 0))))
variance_mat[(length(tissues)+1),3] <- sqrt(wtd.var(jk_sum[,1],jk_weights[,1]))*sqrt(length(which(jk_weights[,1] != 0)))
variance_mat[,4] <- pnorm(q = 0, mean = as.numeric(variance_mat[,2]), sd = as.numeric(variance_mat[,3]))
twosidedP <- ifelse(as.numeric(variance_mat[1:length(tissues),2]) > 0, as.numeric(variance_mat[1:length(tissues),4]), 1-as.numeric(variance_mat[1:length(tissues),4]))
twosidedP <- twosidedP*2  
fdradjusted_twosidedp <- p.adjust(twosidedP, method = "fdr")
variance_mat[1:length(tissues),5] <- fdradjusted_twosidedp 
variance_mat[(length(tissues)+1),5] <- NA
variance_mat[1:nrow(variance_mat),6] <- as.numeric(variance_mat[,2]) / cov_ge
variance_mat[1:nrow(variance_mat),7] <- as.numeric(variance_mat[,3]) / cov_ge

colnames(variance_mat) <- c("Tissue","omega_ge_t","omega_ge_t_se","NomP","FDRP","PropCov","PropCov_se")
write.table(variance_mat, file = paste0("TCSC/results/CrossTraitTCSC_",trait1,"_",trait2,".txt"), row.names = F, col.names = T, sep = "\t", quote = F)
