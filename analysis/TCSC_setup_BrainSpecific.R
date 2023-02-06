#please set working directory with setwd() to your local path to TCSC/

library(data.table)
library(Hmisc)

gtex <- fread("TCSC/analysis/gene_annotation.txt.gz", header = F, sep = "\t")

y <- fread("TCSC/analysis/TissueGroups.txt", header = T)
tissues <- y$Tissues

X <- as.matrix(fread("TCSC/coregulation_scores/CoregulationMatrix_allGTEx_030422.txt.gz", header = F))

transcripts <- c()
count <- 0
starts <- c()
chrs <- c()
tissueassign <- c()
for (i in 1:length(tissues)){
tissue <- tissues[i]
transcript_key <- fread(paste0("TCSC/weights/heritablegenes/Nall/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1
m <- match(transcript_key,gtex$V7)
genetype <- gtex$V8[m]
transcript_key <- transcript_key[which(genetype == "protein_coding")]
transcripts <- c(transcripts,transcript_key)
tissueassign <- c(tissueassign, rep(i,length(transcript_key)))

gene_annot <- fread(paste0("TCSC/weights/allEUR_tissues/v8_allEUR_",tissue,"_blup.pos"), header = T)
gene_annot_transcript <- gene_annot$ID
m <- match(transcript_key,gene_annot_transcript)
starts <- c(starts,gene_annot$P0[m])
chrs <- c(chrs,gene_annot$CHR[m])
count <- count + length(transcript_key) 
}

#brain tissue filter
g <- grep("Brain",tissues)
tissues <- tissues[g] 
X <- X[,g]
wkeepX <- which(tissueassign %in% g)
X <- X[wkeepX,]
tissueassign <- tissueassign[wkeepX]
transcripts <- transcripts[wkeepX]
chrs <- chrs[wkeepX]
starts <- starts[wkeepX]

save(list=c("wkeepX","X","transcripts","starts","chrs","tissueassign"), file = "TCSC/analysis/InputCoreg_BrainSpecific.RData")

