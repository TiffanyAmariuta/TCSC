trait_tissue_type <- commandArgs(trailingOnly=TRUE)
library(data.table)
trait_tissue_type <- strsplit(trait_tissue,split = ",")[[1]]
trait <- trait_tissue_type[1]
tissue <- trait_tissue_type[2]
type <- trait_tissue_type[3]

if(type == "large"){

  transcript_key <- fread(paste0("TCSC/weights/heritablegenes/N320/TranscriptsIn",tissue,"Model.txt"), header = F)$V1 
  twas_sumstats_z <- c()
  twas_sumstats_transcripts <- c()
  for (chr in 1:22){
  twas_sumstats <- fread(paste0("results_320/v8_320EUR.",trait,"/v8_320EUR.",trait,".",tissue,".",chr,".dat"), header = T)
  twas_sumstats_z <- c(twas_sumstats_z, as.numeric(twas_sumstats$TWAS.Z))
  twas_sumstats_transcripts <- c(twas_sumstats_transcripts, twas_sumstats$ID)
  }

  m <- match(transcript_key,twas_sumstats_transcripts) 
  newalphaz <- twas_sumstats_z[m]
  write.table(newalphaz, file = paste0("TCSC/twas_statistics/320EUR_GTEx/",trait,"/Marginal_alphas_",trait,"_",tissue,".txt"), row.names = F, col.names = F, sep = "\t", quote = F)
  system(paste0("gzip TCSC/twas_statistics/320EUR_GTEx/",trait,"/Marginal_alphas_",trait,"_",tissue,".txt"))

}else{

  transcript_key <- fread(paste0("TCSC/weights/heritablegenes/Nall/TranscriptsIn",tissue,"Model.txt"), header = F)$V1 
  twas_sumstats_z <- c()
  twas_sumstats_transcripts <- c()
  for (chr in 1:22){
  twas_sumstats <- fread(paste0("results_all/v8_allEUR.",trait,"/v8_allEUR.",trait,".",tissue,".",chr,".dat"), header = T)
  twas_sumstats_z <- c(twas_sumstats_z, as.numeric(twas_sumstats$TWAS.Z))
  twas_sumstats_transcripts <- c(twas_sumstats_transcripts, twas_sumstats$ID)
  }

  m <- match(transcript_key,twas_sumstats_transcripts) 
  newalphaz <- twas_sumstats_z[m]
  write.table(newalphaz, file = paste0("TCSC/twas_statistics/allEUR_GTEx/",trait,"/Marginal_alphas_",trait,"_",tissue,".txt"), row.names = F, col.names = F, sep = "\t", quote = F)
  system(paste0("gzip TCSC/twas_statistics/allEUR_GTEx/",trait,"/Marginal_alphas_",trait,"_",tissue,".txt")) 
  
}


