#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\独立验证集GEO\\GSE72922")

list.files()
expr<-read.delim("./geneMatrix.healthy.txt",header=T,row.names = 1)
#expr<-read.delim("./exp.tpm.for_tissue_origin.txt",header=T)
head(expr)

# 检查第一列
first_col <- expr[, 1]

# 找出重复的基因名
dup_genes <- unique(first_col[duplicated(first_col)])

if(length(dup_genes) > 0) {
  cat("重复的基因名:\n")
  print(dup_genes)
  
  # 统计每个重复基因出现的次数
  cat("\n每个重复基因的出现次数:\n")
  dup_counts <- table(first_col[first_col %in% dup_genes])
  print(dup_counts)
} else {
  cat("没有重复的行名\n")
}

#exp<-test1
####ref
mat<-read.delim("gtex.ref.16tissue.exp.matrix.txt",header=T,row.names = 1)
head(mat)
#tmp<-expr[rownames(mat),1:4]
#rownames(mat)<-tmp$symbol

#inputFile="exp.tpm.for_tissue_origin.txt"    

source("geoML24.CIBERSORT.R")      
exp<-expr
outTab=CIBERSORT("gtex.ref.16tissue.exp.matrix.txt", "geneMatrix.healthy.txt", perm=1000)
head(outTab)
write.table(outTab, file="GSE72922_CIBERSORT-Results.txt", sep="\t", quote=F, col.names=T)

outTab=CIBERSORT("gtex.ref.16tissue.exp.matrix.txt", "geneMatrix.CRC.txt", perm=1000)
head(outTab)
write.table(outTab, file="GSE72922_CIBERSORT-Results.CRC.txt", sep="\t", quote=F, col.names=T)

#
outTab=outTab[outTab[,"P-value"]<1,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="uEV_RNA_30sample_CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


