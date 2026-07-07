#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\自己构建组织特异性表达矩阵\\GTEx\\CIBERSORT.modify.nu")

list.files()
expr<-read.delim("./exp.tpm.protein_coding.txt",header=T,row.names = 1)
head(expr)

#exp<-test1
####ref
mat<-read.delim("gtex.ref.16tissue.exp.matrix.txt",header=T,row.names = 1)
head(mat)
#tmp<-expr[rownames(mat),1:4]
#rownames(mat)<-tmp$symbol

#inputFile="exp.tpm.protein_coding.txt"    

source("geoML24.CIBERSORT.modify.nu.R")      
exp<-expr
outTab=CIBERSORT("gtex.ref.16tissue.exp.matrix.txt", "exp.tpm.protein_coding.txt", perm=1000)
head(outTab$results)
write.table(outTab$results, file="uEV_RNA_30sample_CIBERSORT-Results.txt", sep="\t", quote=F, col.names=T)
#
outTab=outTab[outTab[,"P-value"]<1,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="uEV_RNA_30sample_CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


