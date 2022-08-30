#please set working directory with setwd() to your local path to TCSC/

library(data.table)
library(Hmisc)

gtex <- fread("analysis/gene_annotation.txt.gz", header = F, sep = "\t")

y <- fread("analysis/TissueGroups.txt", header = T)
tissues <- unique(y$MetaTissue)

n_eqtl <- sapply(1:length(tissues), function(x) sum(y$N_EUR[which(y$MetaTissue == tissues[x])])) #find total available sample size from GTEx data after accounting for meta-analysis of tissues
small_tissues <- which(n_eqtl < 320)
normal_tissues <- c(1:length(tissues))[-small_tissues] #in primary analysis, these tissues are subsampled such that eQTL sample size = 320

X <- as.matrix(fread("coregulation_scores/CoregulationMatrix_320orlessGTEx_062122.txt.gz", header = F))

remove_genes1 <- unique(which(is.na(rowSums(X)))) 
remove_genes2 <- unique(which(rowSums(X) == 0))
remove_genes <- unique(c(remove_genes1,remove_genes2))
transcripts <- c()
count <- 0
starts <- c()
chrs <- c()
tissueassign <- c()
for (i in 1:length(tissues)){
tissue <- tissues[i]
if(i %in% small_tissues){
transcript_key <- fread(paste0("weights/heritablegenes/Nall/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1
keep <- fread(paste0("weights/heritablegenes/Nall/TranscriptsIn",tissues[i],"Model_keep.txt"), header = F)$V1
}else{
transcript_key <- fread(paste0("weights/heritablegenes/N320/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1
keep <- fread(paste0("weights/heritablegenes/N320/TranscriptsIn",tissues[i],"Model_keep.txt"), header = F)$V1
}
w <- which(transcript_key %in% keep)
transcript_key <- transcript_key[w]
m <- match(transcript_key,gtex$V7)
genetype <- gtex$V8[m]
transcript_key <- transcript_key[which(genetype == "protein_coding")]
transcripts <- c(transcripts,transcript_key)
tissueassign <- c(tissueassign, rep(i,length(transcript_key)))

if(i %in% small_tissues){ #clean up by using gtex file? 
gene_annot <- fread(paste0("weights/allEUR_tissues/v8_allEUR_",tissue,"_blup.pos"), header = T)}else{
gene_annot <- fread(paste0("weights/320EUR_metatissues/",tissue,".pos"), header = T)}
gene_annot_transcript <- gene_annot$ID
m <- match(transcript_key,gene_annot_transcript) 
starts <- c(starts,gene_annot$P0[m])
chrs <- c(chrs,gene_annot$CHR[m])
count <- count + length(transcript_key)
}

mhc <- which(chrs == 6 & as.numeric(starts) > 29000000 & as.numeric(starts) < 33000000)
remove_genes <- c(remove_genes, mhc)

if(length(remove_genes)>0){
transcripts <- transcripts[-remove_genes]
starts <- starts[-remove_genes]
chrs <- chrs[-remove_genes]
tissueassign <- tissueassign[-remove_genes]
X <- X[-remove_genes,]
}

save(list=c("X","transcripts","starts","chrs","remove_genes","tissueassign"), file = "analysis/InputCoreg_TCSC.RData")

