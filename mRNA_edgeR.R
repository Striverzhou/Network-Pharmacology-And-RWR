#!/usr/bin/env Rscript
############Differential Expression Analysis of mRNA#############

if (!requireNamespace("BiocManager"))
	install.packages("BiocManager")
if (require("edgeR")){
	print("edgeR installed")}else{
		BiocManager::install("edgeR")
}

if (require("gplots")){
	print("gplots installed")}else{
		install.packages("gplots")
}


foldChange=1.5
padj=0.05

setwd(".")                    #set your workspace dir
library("edgeR")
rt=read.table("mRNA_symbol.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

group=c(rep("normal",normal_count),rep("tumor",tumor_count))
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y) #estimate the mean value of Empirical Bayes robust discrete
y <- estimateTagwiseDisp(y) #estimate the value of Empirical Bayes robust discrete for all genes
et <- exactTest(y,pair = c("normal","tumor"))
#topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F) 
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)

heatmapData <- newData[rownames(diffSig),]

#volcano
pdf(file="vol.pdf",width=5,height=5)
#tiff(file="vol.tiff",width =12,height =12,units ="cm",compression="lzw",bg="white",res=400)
if (max(-log10(allDiff$FDR)) < 350){
  xMax=max(-log10(allDiff$FDR))+1 }else{
    xMax=350
  }
yMax=12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",
     main="mRNA_Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC>foldChange,]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC<(-foldChange),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="blue",cex=0.4)
abline(h=0,lty=2,lwd=3)
while (!is.null(dev.list()))  dev.off()

