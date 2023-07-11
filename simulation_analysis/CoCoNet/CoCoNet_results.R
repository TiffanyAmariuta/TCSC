library(data.table)
sims <- 1000
res_mat <- data.frame()
for (sim in 1:sims){
  tryCatch({
  res <- as.matrix(fread(paste0("TCSC/simulation_analysis/CoCoNet/results/res_",sim,".txt.gz")))
  res_mat <- rbind(res_mat, cbind(1:6,res))
  },error=function(cond){message(paste("Didn't finish ", sim))})
}
power <- c()
fpr <- c()
for (i in 1:6){
w <- which(res_mat[,1] == i)
res_mat2 <- res_mat[w,-1]
performance <- sapply(1:nrow(res_mat2), function(x) which.max(as.numeric(res_mat2[x,]))) #1 1 1 2
power[i] <- length(which(performance == 1))/length(performance)*100
fpr[i] <- 100-power[i]
}

