
library(data.table)
#sims <- 50
tissues <- 0:9
samps <- c(100,200,300,500,1000,1500)
ve <- 0.05
res <- matrix(0,0,5)
setwd("TCSC/simulation_analysis/RTC_Coloc/output/Ongen_files/res/") #need to make this directory
files <- list.files(pattern = "ve0.1.txt") #"Ongen_files/res/")
for (sim in files){
    myres <- fread(sim, header = F)
    res <- rbind(res,myres)
}

power <- c()
fpr <- c()
for (i in 1:length(samps)){
    w <- which(res[,2] == i & res[,3] == 0)
    power[i] <- length(which(res[w,4] > 1 & res[w,5] < 0.05))/length(w)*100
    w <- which(res[,2] == i & res[,3] != 0)
    fpr[i] <- length(which(res[w,4] > 1 & res[w,5] < 0.05))/length(w)*100
}
print(power)
print(fpr)

x <- 10^(-rev(seq(1,100,1)))
pvals <- c(0,x,5*x,seq(0.01,1,0.00124)) #1000 points

simple_auc <- function(TPR, FPR){
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

for (i in 1:6){
sensitivity <- c()
oneminus_specificity <- c()
#i <- 3
w_causal <- which(res[,2] == i & res[,3] == 0)
w_tagging <- which(res[,2] == i & res[,3] != 0)
for (p in 1:length(pvals)){   
    sensitivity[p] <- length(which(res[w_causal,4] > 1 & res[w_causal,5] < pvals[p]))/length(w_causal)
    oneminus_specificity[p] <- length(which(res[w_tagging,4] > 1 & res[w_tagging,5] < pvals[p]))/length(w_tagging)
}
a <- cbind(oneminus_specificity,sensitivity)
a <- rbind(a, c(1,1))
a <- a[order(a[,1],decreasing = F),]
print(simple_auc(a[,2],a[,1]))
colnames(a) <- c("oneminus_specificity","sensitivity")
write.table(a, file = paste0("RTCColoc_ROCcurve_samp",samp,".txt"), row.names = F, col.names = T, sep = "\t", quote = F)
}

