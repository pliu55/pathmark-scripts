#! /usr/bin/Rscript --vanilla 
usage <- "\nusage: limma_ng input contrast top_count output.tab top_output.tab out.pdf"
library(edgeR)
library(limma)
args <- commandArgs(TRUE)

x = read.table(args[1], header = TRUE, row.names = 1, check.names=FALSE, stringsAsFactors=FALSE)
contrast = read.delim(args[2], header = TRUE, row.names = 1, check.names=FALSE, stringsAsFactors=FALSE)
u = intersect(colnames(x), rownames(contrast))
length(u)

x = x[, u]
dim(x)
rownames(contrast)
contrast = contrast[u, ]
length(contrast)
colnames(x)

dge = DGEList(counts=x)
isexpr = rowSums(cpm(dge) > 10) >= 2
flt = dge[isexpr,]
tmm = calcNormFactors(flt)
design = model.matrix(~ contrast)
y = voom(tmm, design, plot=TRUE)
pdf(args[6])
plotMDS(y,top=50,labels=contrast, col=ifelse(contrast=="SmallCell","blue","red"),gene.selection="common")
dev.off()

fit = eBayes(lmFit(y, design))
cn = sprintf("contrast%s", as.character(levels(as.factor(contrast))[2]))
cn
tt = topTable(fit, coef=cn, number=as.numeric(args[3]))
limma_out = list(design=design, y=y, fit=fit, tt=tt)
write.table( limma_out$fit, file=args[4],  sep="\t" , quote=FALSE, col.names=NA)
write.table( tt, file=args[5],  sep="\t" , quote=FALSE, col.names=NA)