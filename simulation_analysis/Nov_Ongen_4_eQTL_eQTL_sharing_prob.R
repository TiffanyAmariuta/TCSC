sim <- as.numeric(commandArgs(trailingOnly=TRUE))
library(data.table)

tissues <- 0:9
samps <- c(100,200,300,500,1000,1500)
inputs <- c("Ongen_GWASeQTL_res","Ongen_GWASeQTL_res_null",paste0("Ongen_eQTLeQTL_res_tissue",tissues))
ve <- 0.1

res <- matrix(0,0,5)
    print(sim)
    num_gwasvariants <- nrow(fread(paste0("Ongen_files/gwas_snps/Nov_sim",sim,"_ve",ve,".txt.gz"), header = F))
           
    for (samp in 1:length(samps)){
        print(samp)
        sigeQTLs <- c()
        for (ti in tissues){
            yy <- fread(paste0("Ongen_files/eqtl_tissue",ti,"/SigeQTLs_samp",samps[samp],"_sim",sim,".txt.gz"), header = F)
            sigeQTLs <- c(sigeQTLs,nrow(yy))
        }
        for (ti in tissues){
            othertissues <- tissues[-match(ti,tissues)]
            #read in GWAS-eQTL p_share for tissue ti
            #1.rtc, 2.gwas snp (or eqtl snp), 3.eqtl snp, 4.nearby gene, 5.beta, 6.int, 7.boundary used, 8.pshare
            y_real <- as.matrix(fread(paste0("Ongen_files/",inputs[1],"/GWAS_SNP_RTC_Tissue",ti,"_sim",sim,"_ve",ve,"_samp",samps[samp],"_Nov.txt"), header = F))
            y_null <- as.matrix(fread(paste0("Ongen_files/",inputs[2],"/GWAS_SNP_RTC_Tissue",ti,"_sim",sim,"_ve",ve,"_samp",samps[samp],"_Nov.txt"), header = F))
            #problem that mean(as.numeric(y_null[,9])) is greater than mean(as.numeric(y_real[,9]))
            #iterate through eqtl snp in col 3, divide col 8 by the pshare of that eqtl in the other tissues. 
            #note: these eqtls that were found here are at a very nonstrict threshold. i
            #should we factor in the zeros from the null and the real? 
            #maybe want each variant in y_real and y_null to only count once?  
            #sum(as.numeric(y_real[,9]))/sum(as.numeric(y_null[,9])), good when this is > 1, normalizing by tissue sharing shouldn't change much. 
            types <- c("real","null")
            for (j in 1:length(types)){
                y <- get(paste0("y_",types[j]))
                #alt soln: take the top pshare eqtl coloc for each variant. 
                #read in data first. 
                for (tis in othertissues){
                if(j == 1){
                yy <- as.matrix(fread(paste0("Ongen_files/",inputs[3+ti],"/GWAS_SNP_RTC_Tissue",tis,"_sim",sim,"_ve",ve,"_samp",samps[samp],"_Nov.txt"), header = F))
                }else{
                yy <- as.matrix(fread(paste0("Ongen_files/",inputs[3+ti],"/NULLGWAS_SNP_RTC_Tissue",tis,"_sim",sim,"_ve",ve,"_samp",samps[samp],"_Nov.txt"), header = F))
                }
                assign(paste0("yy_",tis),yy)
                } #true that pshare of gwas-eqtl  
                tissueshare <- c()
                for (i in 1:nrow(y)){
                    pshare <- 0
                    for (tis in othertissues){
                        yy <- get(paste0("yy_",tis))
                        w <- which(yy[,2] == y[i,3]) #which focal variants in eqtl-eqtl file match eqtl that colocalizes with the gwas variant?
                        if(length(w) > 0){pshare <- pshare + mean(as.numeric(yy[w,9]))}
                    }
                    tissueshare <- c(tissueshare,pshare)
                    #y[i,9] <- as.numeric(y[i,9])/(pshare) #update the value. sometimes pshare is 0. 
                }
                y <- cbind(y,tissueshare)
                assign(paste0("y_",types[j]),y)
                #ntcs <- as.numeric(y[,9])/sigeQTLs[ti+1]
                #ntcs <- c(ntcs,rep(0,num_gwasvariants-length(unique(y[,2]))))
                #assign(paste0("ntcs_",types[j]),ntcs)
            }
            #compare ntcs
            ntcs_real <- as.numeric(y_real[,9])/as.numeric(y_real[,10])/sigeQTLs[ti+1]
            ntcs_null <- as.numeric(y_null[,9])/as.numeric(y_null[,10])/sigeQTLs[ti+1]
            #if(length(ntcs_real) > length(ntcs_null)){ntcs_null <- c(ntcs_null,rep(0,length(ntcs_real)-length(ntcs_null)))}
            enrichment_p <- wilcox.test(ntcs_real, ntcs_null)$p.val
            #enrichment <- sum(ntcs_real)/(sum(ntcs_null)*length(ntcs_real)/length(ntcs_null))
            enrichment <- sum(ntcs_real)/sum(ntcs_null)
            dump <- c(sim, samp, ti, enrichment, enrichment_p)
            res <- rbind(res, dump)
        } #tissues
    } #samps
    write.table(res,file = paste0("Ongen_files/res/res_",sim,"_ve",ve,".txt"),row.names = F, col.names = F, sep = "\t", quote = F)


