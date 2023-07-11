library(data.table)
set.seed(12345)
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

ngenes_tot <- 1000
ngenes_causal <- 100
ngenes_directsnpeffect <- 100 
rg <- 0.5 
nsims <- 2000 

y1 <- matrix(0,ngenes_tot,nsims)
y2 <- y1
y3 <- y1


for (sim in 1:nsims){

causal <- sample(ngenes_tot,size = ngenes_causal, replace = F)
direct <- sample(c(1:ngenes_tot)[-causal], size = ngenes_directsnpeffect, replace = F)
initialize_alphas <- 0 #is there the assumption that alpha is normally distributed? beta is usually normally distributed. if effect sizes are too small we might not have power to get corr. 
while(length(which(initialize_alphas == 0) > 0)){
initialize_alphas <- rnorm(n = ngenes_causal,0,1) #March 23
}
y1[causal,sim] <- initialize_alphas
y2[causal,sim] <- complement(y = initialize_alphas, rho = rg)
y3[direct,sim] <- 1

print(cor(y1[,sim],y2[,sim]))

}
write.table(y1, file = paste0("TCSC/simulation_analysis/Nov_Y1_alpha_nCausalGenes",ngenes_causal,".txt"), row.names = F, col.names = F, sep = "\t", quote = F)
write.table(y2, file = paste0("TCSC/simulation_analysis/Nov_Y2_alpha_nCausalGenes",ngenes_causal,".txt"), row.names = F, col.names = F, sep = "\t", quote = F)
system(paste0("gzip TCSC/simulation_analysis/Nov_Y1_alpha_nCausalGenes",ngenes_causal,".txt"))
system(paste0("gzip TCSC/simulation_analysis/Nov_Y2_alpha_nCausalGenes",ngenes_causal,".txt"))

