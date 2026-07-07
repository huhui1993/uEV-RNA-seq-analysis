#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\自己构建组织特异性表达矩阵\\GTEx\\exoRBase3.0.EVP_urine")

list.files()
expr<-read.delim("./EVPs_in_Urine_longRNAs.txt",header=T,row.names = NULL)
head(expr);dim(expr)

test1<-aggregate(x=expr[,2:126],by=list(expr$Gene_name),FUN=mean,na.rm=T) 
dim(test1);head(test1)

rownames(test1)<-test1$Group.1
test1<-test1[,-1]
head(test1)

write.table(test1,"EVPs_in_Urine_longRNAs.final.txt", 
            sep="\t", col.names=T,quote = F)

#exp<-test1
####ref
mat<-read.delim("gtex.ref.16tissue.exp.matrix.txt",header=T,row.names = 1)
head(mat)
#tmp<-expr[rownames(mat),1:4]
#rownames(mat)<-tmp$symbol

#inputFile="exp.tpm.protein_coding.txt"    

source("geoML24.CIBERSORT.R")      
exp<-test1
outTab=CIBERSORT("gtex.ref.16tissue.exp.matrix.txt", "EVPs_in_Urine_longRNAs.final.txt", perm=1000)
head(outTab)
write.table(outTab, file="exoRBase3.0.EVP_urine.CIBERSORT-Results.txt", sep="\t", quote=F, col.names=T)

#
outTab=outTab[outTab[,"P-value"]<1,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="uEV_RNA_30sample_CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


