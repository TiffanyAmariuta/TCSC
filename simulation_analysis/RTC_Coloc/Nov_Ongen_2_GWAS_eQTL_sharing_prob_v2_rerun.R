sim_b_samp <- commandArgs(trailingOnly=TRUE)
s <- strsplit(sim_b_samp, split = "_")[[1]]
sim <- as.numeric(s[1])
b <- as.numeric(s[2])
samp <- as.numeric(s[3]) #ran within 5 hours for GWAS coloc over samps, but not for eqtl-eqtl coloc. can be NA for b = 1, b = 2 
ve <- as.numeric(s[4])

library(data.table)
source("TCSC/simulation_analysis/RTC_Coloc/Ongen_fxn.R")
ngenes <- 1000
genes <- fread("TCSC/simulation_analysis/1000genes_Nov.txt", header = T)
#tsses <- as.numeric(genes$START)
tissues <- 0:9
samplesizes <- c(100,200,300,500,1000,1500)
snpmat <- fread("TCSC/simulation_analysis/1KG_HM3_chr1.bim", header = F) #GRCh37
QTL_geno_all <- as.matrix(fread(paste0("TCSC/simulation_analysis/Simulated_eQTL_Cohort_for_ImputedExpression.txt.gz"), header = F))
colnames(QTL_geno_all) <- snpmat$V2
#coldspots <- fread("Ongen_files/coldspots_hg19.txt.gz", header = F) #too few snps in these coldspots
#want to take recombination hotspots (bed file) and define cold spots (regions between these hotpots)
#y <- fread("Ongen_files/hotspots_hg19.txt.gz", header = F)
#y <- y[which(y$V1 == "chr1"),]
#coldspots <- matrix(0,nrow(y)+1,2)
#coldspots[1,] <- c(1,y$V2[1])
#for (i in 1:nrow(y)){
#coldspots[i+1,] <- c(y$V3[i],y$V2[i+1])
#}
#coldspots[i+1,2] <- snpmat$V4[nrow(snpmat)]
#write.table(coldspots, file = "Ongen_files/coldspots_hg19.txt", row.names = F, col.names = F, sep = "\t", quote = F)

snplists <- c("gwas_snps","gwas_null_snps",paste0("eqtl_tissue",tissues))
outputs <- c("Ongen_GWASeQTL_res","Ongen_GWASeQTL_res_null",paste0("Ongen_eQTLeQTL_res_tissue",tissues))
boundary <- 500000

if(b > 2){samplesizes <- samplesizes[samp]}

for (samp in samplesizes){ #if eqtl-eqtl coloc, just iterate over chosen sample size. 
    print(samp)
    if(b < 3){snp_df <- fread(paste0("TCSC/simulation_analysis/RTC_Coloc/output/Ongen_files/",snplists[b],"/Nov_sim",sim,"_ve",ve,".txt.gz"), header = F)
    othertissues <- tissues
    } #snp index. 
    if(b > 2){#snp_df <- fread(paste0("Ongen_files/",snplists[b],"/SigeQTLs_samp",samp,"_sim",sim,".txt.gz"), header = F) #gene,snp index, move columns around so just use snp.
    #snp_df <- snp_df[,c(2,1)] #now snp index, gene. 
    #colnames(snp_df) <- c("V1","V2")
    #only do eqtl-eqtl coloc with eqtls of tissueA that are coloc with GWAS variant (and null variants). 
    #makes the most difference at large sample size! #at 300, basically the same. so still need to increase job time. 
    gwas_snp_df <- fread(paste0("TCSC/simulation_analysis/RTC_Coloc/output/Ongen_files/Ongen_GWASeQTL_res/GWAS_SNP_RTC_Tissue",b-3,"_sim",sim,"_ve",ve,"_samp",samp,"_Nov.txt"), header=F)
    #gwas_null_df <- fread(paste0("Ongen_files/Ongen_GWASeQTL_res_null/GWAS_SNP_RTC_Tissue",b-3,"_sim",sim,"_ve",ve,"_samp",samp,"_Nov.txt"), header=F)
    gwas_colo_eqtls <- unique(gwas_snp_df$V3) #,gwas_null_df$V3)
    snp_df <- data.frame(match(gwas_colo_eqtls,snpmat$V2))
    colnames(snp_df) <- "V1"
    othertissues <- tissues[-(b-2)] #second mod; only other tissues. 
    } 

    for (ti in othertissues){ 
        checkfile <- paste0("TCSC/simulation_analysis/RTC_Coloc/output/Ongen_files/",outputs[b],"/GWAS_SNP_RTC_Tissue",ti,"_sim",sim,"_ve",ve,"_samp",samp,"_Nov.txt")
        if(!file.exists(checkfile)){ #run code. 
        print(ti)
        tot_GE <- as.matrix(fread(paste0("TCSC/simulation_analysis/totalexpression/Nov_0.75_ggcoreg/",samp,"/TotExp_Nov_Group1_Tissue",ti,"_sim",sim,".txt.gz"),header = F))
        tissue_qtls <- fread(paste0("TCSC/simulation_analysis/RTC_Coloc/output/Ongen_files/eqtl_tissue",ti,"/SigeQTLs_samp",samp,"_sim",sim,".txt.gz"), header = F)
        GWAS_SNP_RTC <- matrix(0,0,8) #data.frame()
        eQTLP_forGWASvar <- c()
        for (j in 1:nrow(snp_df)){
            #print(paste0("working on GWAS SNP ",j," out of ",nrow(snp_df)))
            #find which eQTLs are within 500000 (recombination hotspot) of this GWAS variant. 
            bp <- snpmat$V4[snp_df$V1[j]] 
            rsid <- snpmat$V2[snp_df$V1[j]]
            snps_in_range <- which(snpmat$V4 >= bp-boundary & snpmat$V4 <= bp+boundary)
            eqtls <- intersect(snps_in_range,tissue_qtls$V2) #will be for different genes. eqtls that are nearby the gwas variant. WHAT HAPPENS WHEN NO EQTLS?
            if(length(eqtls) > 0){
                egenes <- unique(tissue_qtls$V1[which(tissue_qtls$V2 %in% eqtls)]) #extend boundary for 500kb +/- egene body. update snps_in_range and eqtls. 
                left_boundary <- min(bp-boundary,genes$START[egenes]-boundary) 
                right_boundary <- max(bp+boundary,genes$STOP[egenes]+boundary)
                snps_in_range <- which(snpmat$V4 >= left_boundary & snpmat$V4 <= right_boundary) #new cis window
                eqtls <- intersect(snps_in_range,tissue_qtls$V2)
                QTL_geno <- QTL_geno_all[1:samp,snps_in_range] #otherwise, gwas variant wont be there. 
                pos_gwassnp_inQTL_genomatrix <- match(rsid,colnames(QTL_geno))
                for (g in egenes){
                    tot_GE_gene <- tot_GE[g,]
                    independent_signals <- ConditionalAnalysis(QTL_geno, tot_GE_gene) #find eQTLs.
                    corP <- cor.test(tot_GE_gene, QTL_geno[,pos_gwassnp_inQTL_genomatrix])$p.val
                    eQTLP_forGWASvar <- c(eQTLP_forGWASvar, corP)
                    if(length(independent_signals) > 0){
                        for (is in 1:length(independent_signals)){
                            rtc <- RTC(tot_GE_gene, independent_signals, QTL_geno, pos_gwassnp_inQTL_genomatrix, is) #rank of snp
                            dump <- c(rtc[1], rsid, colnames(QTL_geno)[independent_signals[is]], g, rtc[3], rtc[2], left_boundary, right_boundary) #beta, int
                            GWAS_SNP_RTC <- rbind(GWAS_SNP_RTC,dump) #rtc, gwas snp, eqtl snp, nearby gene, beta, int, coldspot used
                         } #independent signals
                     } #if 
                 } #genes 
             } #if eqtls and gwas variant in coldspot. 
         } #over sig snps

        ## Pi statistics: pi1 stat = proportion of true positives based on a distribution of p values. 
        y <- hist(eQTLP_forGWASvar, probability = T, breaks = 20)
        pi0 <- min(1,mean(y$density[c(floor(length(y$density)/2):length(y$density))])) #proportion of null p values  
        pi_tp <- 1-pi0

        if(pi_tp > 0){
        ## Null and Alternative Hypotheses (need this for bayes step) (every gene is a coldspot, flanked by recomb hotspots)
        allcoldspots <- sapply(1:nrow(GWAS_SNP_RTC), function(x) paste0(GWAS_SNP_RTC[x,7:8],collapse = "_"))
        coldspots = unique(allcoldspots)
        trials <- 100
        H0_200_overchunks <- matrix(0,trials,length(coldspots)) #null: two variants in cold spot tag different functional effects 
        H1_200_overchunks <- matrix(0,trials,length(coldspots))
        for (j in 1:length(coldspots)){ 
            #print(paste0("working on coldspot ",j," out of ",length(coldspots)))
            s <- strsplit(coldspots[j],split = "_")[[1]]
            snps_in_range <- which(snpmat$V4 >= as.numeric(s[1]) & snpmat$V4 <= as.numeric(s[2]))
            rows_in_GWAS_SNP_RTC <- which(allcoldspots == coldspots[j])
            QTL_geno <- QTL_geno_all[1:samp,snps_in_range]
            immune <- GWAS_SNP_RTC[rows_in_GWAS_SNP_RTC,3] #supposed to handle multiple causal variants but poorly described.  
            betas <- as.numeric(GWAS_SNP_RTC[rows_in_GWAS_SNP_RTC,5])
            intercepts <- as.numeric(GWAS_SNP_RTC[rows_in_GWAS_SNP_RTC,6])
            possiblesnps <- c(1:ncol(QTL_geno))[-match(immune,colnames(QTL_geno))] #numeric value
            psuedo_pheno <- QTL_geno[,immune[1]] * betas[1] + intercepts[1]
            for (tries in 1:trials){
                #for null: pick two hidden causal variants
                #pick two random causal var (gwas, eqtl), then pick two variants in linkage (for our purposes, this is equiv to pick 2 random variants in gene locus)
                gwas_eqtl_hidden_linked <- sample(possiblesnps,size = 2, replace = F) #probably only ned 2 since we can only "see" the linked ones. 
                rtc <- RTC(psuedo_pheno,gwas_eqtl_hidden_linked[1],QTL_geno,gwas_eqtl_hidden_linked[2],1) #pretend there are 2 independent eqtl signals: eqtl real and causal
                H0_200_overchunks[tries,j] <- rtc[1] #RTC_null_nolinkage(QTL_geno)
                #pick one random causal var (gwas also is eqtl), then pick two vars in linkage (for our purposes, this is equiv to pick 2 neighbor variants in gene locus)
                s <- sample(c(1:length(possiblesnps))[-length(possiblesnps)],size = 1)
                gwas_eqtl_hidden_linked <- c(possiblesnps[s],possiblesnps[s+1]) #need random two neighboring snps. They are neighbors b/c in order. 
                rtc <- RTC(psuedo_pheno,gwas_eqtl_hidden_linked[1],QTL_geno,gwas_eqtl_hidden_linked[2],1)
                H1_200_overchunks[tries,j] <- rtc[1] #RTC_alt_nolinkage(QTL_geno) this value should be higher on average than H0
            }
        }

        ## P(share) with Bayes 
        pshare <- c()
        for (j in 1:nrow(GWAS_SNP_RTC)){
            #g <- as.numeric(GWAS_SNP_RTC[j,4]) #coldspot
            #gene <- match(g, coldspots)
            spot <- match(allcoldspots[j],coldspots) #match to unique list
            h0h1 <- c(H0_200_overchunks[,spot],H1_200_overchunks[,spot])  #200 values
            h0h1_ID <- c(rep("h0",length(H0_200_overchunks[,spot])), rep("h1",length(H1_200_overchunks[,spot])))
            h0h1 <- cbind(h0h1_ID,h0h1)
            h0h1 <- h0h1[order(h0h1[,2]),]
            matchh0h1 <- which.min(abs(as.numeric(GWAS_SNP_RTC[j,1]) - as.numeric(h0h1[,2]))) #closest match
            rtc_sim_range <- seq(matchh0h1-0.1*2*trials,matchh0h1+0.1*2*trials,1) #define range of flanking values
            rtc_sim_range <- rtc_sim_range[rtc_sim_range > 0]
            rtc_sim_range <- rtc_sim_range[rtc_sim_range <= nrow(h0h1)]
            h0prop <- length(which(h0h1[rtc_sim_range,1] == "h0"))/length(rtc_sim_range) #how many in this range were from h0? = p(rtc = rtc | not shared)
            h1prop <- length(which(h0h1[rtc_sim_range,1] == "h1"))/length(rtc_sim_range) #how many are from h1? p(RTC = rtc | shared)
            pshare[j] <- h1prop*pi_tp / (h0prop*pi0 + h1prop*pi_tp) #0.83 GWAS variant - eQTL sharing probability. 
        }
        }else{pshare <- rep(0,nrow(GWAS_SNP_RTC))}

        GWAS_SNP_RTC <- cbind(GWAS_SNP_RTC,pshare) #rtc, gwas snp (or eqtl snp), eqtl snp, nearby gene, beta, int, boundary used, pshare
        write.table(GWAS_SNP_RTC, file = paste0("TCSC/simulation_analysis/RTC_Coloc/output/Ongen_files/",outputs[b],"/GWAS_SNP_RTC_Tissue",ti,"_sim",sim,"_ve",ve,"_samp",samp,"_Nov.txt"), row.names = F, col.names = F, quote= F, sep = "\t")
        } #if file does not exist 
    } #tissues
} #samps



