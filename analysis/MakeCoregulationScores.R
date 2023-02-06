i = as.numeric(commandArgs(trailingOnly=TRUE)) #focal tissue
ciswindow <- 2e6 

library(data.table)

gene_ref <- fread("TCSC/analysis/AllGenesRef_in_v8320smalltoo_analysis.txt.gz") #might not have all genes including genes from small tissues (missing 2) from tissue 10
gtex <- fread("TCSC/analysis/gene_annotation.txt.gz", header = F, sep = "\t")
m <- match(gene_ref$gene_id, gtex$V7)
length(which(is.na(m))) #0
genetype <- gtex$V8[m]
gene_ref <- cbind(gene_ref,genetype)

y <- fread("TCSC/analysis/TissueGroups.txt", header = T)
tissues <- unique(y$MetaTissue)

small_tissues <- c(3,5,11,12,14,22,24,27:30,33:35,37:38) #add 1 to values 10 and above; update June 21st for artery aorta
normal_tissues <- c(1:length(tissues))[-small_tissues]

if(i %in% small_tissues){load(paste0("TCSC/predicted_expression/Nall/designmat_",tissues[i],"_v8_allEUR_double.RData"))}else{
load(paste0("TCSC/predicted_expression/N320/designmat_",tissues[i],"_v8_320EUR_double.RData"))}
print(dim(df1))
df1_focal <- df1
rm("df1")
if(i %in% small_tissues){transcript_key <- fread(paste0("TCSC/weights/heritablegenes/Nall/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1
keeptranscripts <- fread(paste0("TCSC/weights/heritablegenes/Nall/TranscriptsIn",tissues[i],"Model_keep.txt"), header = F)$V1
}else{
transcript_key <- fread(paste0("TCSC/weights/heritablegenes/N320/TranscriptsIn",tissues[i],"Model.txt"), header = F)$V1
keeptranscripts <- fread(paste0("TCSC/weights/heritablegenes/N320/TranscriptsIn",tissues[i],"Model_keep.txt"), header = F)$V1
}

#keep genes with gcta p < 0.01 aaaaand h2 estimate > 0
print(length(transcript_key) == ncol(df1_focal))
w <- which(transcript_key %in% keeptranscripts)
df1_focal <- df1_focal[,w]
transcript_key <- transcript_key[w]

#keep protein coding genes only 
m <- match(transcript_key, gene_ref$gene_id)
focal_tissue_gene_ref <- gene_ref[m,]
keepgenes_focal <- which(focal_tissue_gene_ref[,5] == "protein_coding")
transcript_key <- transcript_key[keepgenes_focal]
df1_focal <- df1_focal[,keepgenes_focal]
focal_tissue_gene_ref <- focal_tissue_gene_ref[keepgenes_focal,]

for (j in 1:length(tissues)){
if(j == i){ #apply bias correction
acc <- fread(paste0("TCSC/weights/accuracies/N320_smalltoo_accuracies_cish2_fromall_",tissues[i],".txt.gz"), header = F)
acc <- acc$V2[match(transcript_key,acc$V1)]

genegenecoreg_pc <- c() #vector of length g in focal tissue 
for (k in 1:ncol(df1_focal)){
#gene k is first gene in pair, therefore gene x must be only cis to gene k. 
gene_k_chr <- focal_tissue_gene_ref$chr[k]
gene_k_pos <- focal_tissue_gene_ref$start[k]
ww <- which(focal_tissue_gene_ref$chr == gene_k_chr)
focal_tissue_gene_ref_k_chrom <- focal_tissue_gene_ref[ww,]
ww <- which(focal_tissue_gene_ref_k_chrom$start > gene_k_pos - ciswindow & focal_tissue_gene_ref_k_chrom$start < gene_k_pos + ciswindow)
focal_tissue_gene_ref_k_cislocus <- focal_tissue_gene_ref_k_chrom[ww,]

cis_genes_pc <- match(focal_tissue_gene_ref_k_cislocus$gene_id, focal_tissue_gene_ref$gene_id)
if(length(cis_genes_pc) == 0){cor_est_pc <- 0}else{cor_est_pc <- sapply(cis_genes_pc, function(x) cor(as.numeric(df1_focal[,k]),as.numeric(df1_focal[,x])))}
#cor_est_pc[which(cor_est_pc == 1)] <- acc[k] 
coreg_before_sum_pc <- cor_est_pc^2
mysum <- sum(coreg_before_sum_pc) - 1 + acc[k] 
genegenecoreg_pc <- c(genegenecoreg_pc, mysum)
}
file=paste0("TCSC/coregulation_scores/CoRegScores_v8_320EUR_smalltoo_",tissues[i],"_with_",tissues[j],".txt")
write.table(genegenecoreg_pc, file=file, row.names = F, col.names = F, sep = "\t", quote = F)
system(paste0("gzip ", file,))
}else{ #target tissue is different from focal; normal co-regulation w/o bias correction


if(j %in% small_tissues){load(paste0("TCSC/predicted_expression/Nall/designmat_",tissues[j],"_v8_allEUR_double.RData"))}else{
load(paste0("TCSC/predicted_expression/N320/designmat_",tissues[j],"_v8_320EUR_double.RData"))}

if(j %in% small_tissues){transcript_key_other <- fread(paste0("TCSC/weights/heritablegenes/Nall/TranscriptsIn",tissues[j],"Model.txt"), header = F)$V1
keeptranscripts <- fread(paste0("TCSC/weights/heritablegenes/Nall/TranscriptsIn",tissues[j],"Model_keep.txt"), header = F)$V1
}else{
transcript_key_other <- fread(paste0("TCSC/weights/heritablegenes/N320/TranscriptsIn",tissues[j],"Model.txt"), header = F)$V1
keeptranscripts <- fread(paste0("TCSC/weights/heritablegenes/N320/TranscriptsIn",tissues[j],"Model_keep.txt"), header = F)$V1
}

w <- which(transcript_key_other %in% keeptranscripts)
df1 <- df1[,w]
transcript_key_other <- transcript_key_other[w]

m <- match(transcript_key_other, gene_ref$gene_id)
other_tissue_gene_ref <- gene_ref[m,]
keepgenes_other <- which(other_tissue_gene_ref[,5] == "protein_coding")
transcript_key_other <- transcript_key_other[keepgenes_other]
df1 <- df1[,keepgenes_other]
other_tissue_gene_ref <- other_tissue_gene_ref[keepgenes_other,]

genegenecoreg_pc <- c() #vector of length g in focal tissue 
for (k in 1:ncol(df1_focal)){
#gene k is first gene in pair, therefore gene x must be only cis to gene k. 
gene_k_chr <- focal_tissue_gene_ref$chr[k]
gene_k_pos <- focal_tissue_gene_ref$start[k]
ww <- which(other_tissue_gene_ref$chr == gene_k_chr)
other_tissue_gene_ref_k_chrom <- other_tissue_gene_ref[ww,]
ww <- which(other_tissue_gene_ref_k_chrom$start > gene_k_pos - ciswindow & other_tissue_gene_ref_k_chrom$start < gene_k_pos + ciswindow)
other_tissue_gene_ref_k_cislocus <- other_tissue_gene_ref_k_chrom[ww,]

cis_genes_pc <- match(other_tissue_gene_ref_k_cislocus$gene_id, other_tissue_gene_ref$gene_id)
if(length(cis_genes_pc) == 0){cor_est_pc <- 0}else{cor_est_pc <- sapply(cis_genes_pc, function(x) cor(as.numeric(df1_focal[,k]),as.numeric(df1[,x])))}
coreg_before_sum_pc <- cor_est_pc^2
genegenecoreg_pc <- c(genegenecoreg_pc, sum(coreg_before_sum_pc))
}

file=paste0("TCSC/coregulation_scores/CoRegScores_v8_320EUR_smalltoo_",tissues[i],"_with_",tissues[j],".txt")
write.table(genegenecoreg_pc, file=file, row.names = F, col.names = F, sep = "\t", quote = F)
system(paste0("gzip ", file,))

}
