#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\自己用CIBERSORT预测-参考矩阵还是EV-origin的")

list.files()
expr<-read.delim("exp.tpm.addbiotype.txt",header=T,row.names = 1)
head(expr)

expr.pro<-subset(expr,gene_biotype %in% c("protein_coding"))
expr.pro<-expr.pro[,c(1:31)]


test1<-aggregate(x=expr.pro[,2:31],by=list(expr.pro$symbol),FUN=mean,na.rm=T) 
dim(test1)
rownames(test1)<-test1$Group.1
test1<-test1[,-1]
head(test1)

write.table(test1,"exp.tpm.protein_coding.txt", 
            sep="\t", col.names=T,quote = F)

exp<-expr[,c(2:31)]

write.table(exp,"exp.tpm.protein_coding.ensembl.txt", 
            sep="\t", col.names=T,quote = F)

#exp<-test1
####ref
mat<-read.csv("Matrix_tissue.csv",header=T,row.names = 1)
#tmp<-expr[rownames(mat),1:4]
#rownames(mat)<-tmp$symbol
write.table(mat,"Matrix_tissue.ref.txt", 
            sep="\t", col.names=T,quote = F)

#inputFile="exp.tpm.protein_coding.txt"    

source("geoML24.CIBERSORT.R")      
exp<-expr[,c(2:31)]
outTab=CIBERSORT("Matrix_tissue.ref.txt", "exp.tpm.protein_coding.ensembl.txt", perm=1000)
write.table(outTab, file="uEV_RNA_30sample_CIBERSORT-Results.txt", sep="\t", quote=F, col.names=T)
#
outTab=outTab[outTab[,"P-value"]<1,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="uEV_RNA_30sample_CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


