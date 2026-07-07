setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后")
list.files()

setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后")
setwd("/public/home/huhui/project/rna_urine_EV/result/STAR/DEG")
list.files()
expr<-read.delim("all_sample_count.addAnotation.txt",header=T,row.names = 1)
head(expr);dim(expr)
expr<-subset(expr,gene_biotype %in% c("protein_coding","lncRNA"))
expr[1:4,1:4]
dim(expr)
expr.pro<-subset(expr,gene_biotype %in% c("protein_coding"))
expr.lnc<-subset(expr,gene_biotype %in% c("lncRNA"))

write.table(expr.pro,
            "all_sample_count.addAnotation.protein_coding.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)
write.table(expr.lnc,
            "all_sample_count.addAnotation.lncRNA.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)

phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)

head(expr.pro);dim(expr.pro)

head(expr.lnc);dim(expr.lnc)

data.exp<-expr.lnc[,c(phenotype$sample)]
####expr.pro  condition只提取需要比较的两组,批量做差异
condition <- factor(phenotype$Time, 
                    levels = c("D1","D15","D29","D43"))
list<-c("D1","D15","D29","D43")

combn(list, 2)
combn(list, 2)[,1]
combn(list, 2)[,1][1]
combn(list, 2)[,1][2]
dim(combn(list, 2))
#[1] 2 6

colnames(data.exp)<-phenotype$Time
library(dplyr)
for(i in 1:dim(combn(list, 2))[2]){
gene.count.1<-data.exp[colnames(data.exp) %in% c(combn(list, 2)[,i][1],combn(list, 2)[,i][2])]

head(gene.count.1)
condition.1 <- condition[condition %in% c(combn(list, 2)[,i][1],combn(list, 2)[,i][2])]

library(DESeq2)

colData <- data.frame(row.names=colnames(gene.count.1), condition.1)
colData
dds <- DESeqDataSetFromMatrix(gene.count.1, colData, design= ~ condition.1)
dds <- DESeq(dds)
# 查看一下dds的内容
dds
#接下来，我们要查看case versus control的总体结果，并根据padj进行重新排序。利用summary命令统计显示一共多少个genes上调和下调（FDR0.1）
#res = results(dds, contrast=c("condition", "PM", "PN"))
#或下面命令
res= results(dds)
res = res[order(res$padj),]
head(res)
#summary(res)
#所有结果先进行输出,未筛选
#write.csv(res,file="All_results.csv")

####将表达量矩阵按照差异结果进行排序
exp_data<-gene.count.1[rownames(res),]
head(exp_data)
colnames(exp_data)<-phenotype$Time[phenotype$Time %in% c(combn(list, 2)[,i][1],combn(list, 2)[,i][2])]

####计算正常和疾病样本表达均值(这里是成对的正常样本和癌旁)
colnum<-ncol(exp_data)
test_1<-exp_data[colnames(exp_data) %in% c(combn(list, 2)[,i][1])]
test_2<-exp_data[colnames(exp_data) %in% c(combn(list, 2)[,i][2])]

exp_data$mean_2<-apply(test_2,1,mean)
exp_data$mean_1<-apply(test_1,1,mean)

###合并计算出来的均值和差异结果
res_addMean<-data.frame(exp_data,res)
head(res_addMean)
#####输出未筛选的差异结果
write.table(res_addMean,file=paste("lnc.DESeq2.",combn(list, 2)[,i][2],"vs",combn(list, 2)[,i][1],".compare.nofilter.xls",sep = ""),
            quote = F,sep='\t',row.names=T,col.names =NA )
#####提取筛选后的差异表达基因
##1.去掉低表达
deseq_gene<-res_addMean[apply(res_addMean[,c("mean_2","mean_1")],1,function(x){XXXX<-FALSE;if(any(x>5)){XXXX<-TRUE};return(XXXX)}),]
dim(deseq_gene)
#write.table(deseq_gene,file="gene.DESeq2.PMvsPN.compare.higher5.xls",
 #           quote = F,sep='\t',row.names=T,col.names =NA )
##2.根据padj和FC进一步筛选
#diff_gene_deseq2 <-subset(deseq_gene, padj < 0.05 & abs(log2FoldChange) > 0.585)
diff_gene_deseq2 <-subset(deseq_gene, pvalue < 0.05 & abs(log2FoldChange) > 0.585)
dim(diff_gene_deseq2)
##输出筛选后的差异结果
write.table(diff_gene_deseq2,file=paste("lnc.DESeq2.",combn(list, 2)[,i][2],"vs",combn(list, 2)[,i][1],".compare.p0.05.FC1.5.Higher5.xls",sep = ""),
            quote = F,sep='\t',row.names=T,col.names =NA )

diff_gene_deseq2 <-subset(deseq_gene, pvalue < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq2)
##输出筛选后的差异结果
write.table(diff_gene_deseq2,file=paste("lnc.DESeq2.",combn(list, 2)[,i][2],"vs",combn(list, 2)[,i][1],".compare.p0.05.FC2.Higher5.xls",sep = ""),
            quote = F,sep='\t',row.names=T,col.names =NA )
}


##################################################################################################
####condition只提取需要比较的两组PM vs H
#condition <- factor(c(rep("PM",3),rep("PN",3),rep("H",3)), 
 #                   levels = c("PM","PN","H"))
head(gene.count)
gene.count.1<-gene.count[,c(1:3,7:9)]
head(gene.count.1)
condition <- factor(c(rep("PM",3),rep("H",3)), 
                    levels = c("PM","H"))

library(DESeq2)

colData <- data.frame(row.names=colnames(gene.count.1), condition)
dds <- DESeqDataSetFromMatrix(gene.count.1, colData, design= ~ condition)
dds <- DESeq(dds)
# 查看一下dds的内容
dds
#接下来，我们要查看case versus control的总体结果，并根据padj进行重新排序。利用summary命令统计显示一共多少个genes上调和下调（FDR0.1）
res = results(dds, contrast=c("condition", "PM", "H"))
#或下面命令
res= results(dds)
res = res[order(res$padj),]
head(res)
summary(res)
#所有结果先进行输出,未筛选
#write.csv(res,file="All_results.csv")

####将表达量矩阵按照差异结果进行排序
exp_data<-gene.count.1[rownames(res),]
head(exp_data)
####计算正常和疾病样本表达均值(这里是成对的正常样本和癌旁)
colnum<-ncol(exp_data)
test_PM<-exp_data[,1:3]
test_H<-exp_data[,4:6]

exp_data$mean_PM<-apply(exp_data[,1:3],1,mean)
exp_data$mean_H<-apply(exp_data[,4:6],1,mean)

###合并计算出来的均值和差异结果
res_addMean<-data.frame(exp_data,res)
#####输出未筛选的差异结果
write.table(res_addMean,file="gene.DESeq2.PMvsH.compare.nofilter.xls",quote = F,sep='\t',row.names=T,col.names =NA )
#####提取筛选后的差异表达基因
##1.去掉低表达
deseq_gene<-res_addMean[apply(res_addMean[,1:6],1,function(x){XXXX<-FALSE;if(any(x>5)){XXXX<-TRUE};return(XXXX)}),]
dim(deseq_gene)
#write.table(deseq_gene,file="DESeq2.PMvsH.compare.higher5.xls",quote = F,sep='\t',row.names=T,col.names =NA )
##2.根据padj和FC进一步筛选
diff_gene_deseq2 <-subset(deseq_gene, padj < 0.05 & abs(log2FoldChange) > 0.585)
diff_gene_deseq2 <-subset(deseq_gene, pvalue < 0.05 & abs(log2FoldChange) > 0.585)
dim(diff_gene_deseq2)
##输出筛选后的差异结果
write.table(diff_gene_deseq2,file="gene.DESeq2.PMvsH.compare.p0.05.FC1.5.Higher5.xls",
            quote = F,sep='\t',row.names=T,col.names =NA )

diff_gene_deseq2 <-subset(deseq_gene, pvalue < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq2)
##输出筛选后的差异结果
write.table(diff_gene_deseq2,file="gene.DESeq2.PMvsH.compare.p0.05.FC2.Higher5.xls",
            quote = F,sep='\t',row.names=T,col.names =NA )

####condition只提取需要比较的两组PN vs H
#condition <- factor(c(rep("PN",3),rep("PN",3),rep("H",3)), 
#                   levels = c("PN","PN","H"))
head(gene.count)
gene.count.1<-gene.count[,c(4:9)]
head(gene.count.1);dim(gene.count.1)
gene.count.1<-gene.count.1[apply(gene.count.1,1,function(x){XXXX<-FALSE;if(any(x>=1)){XXXX<-TRUE};return(XXXX)}),]
dim(gene.count.1)

condition <- factor(c(rep("PN",3),rep("H",3)), 
                    levels = c("PN","H"))

library(DESeq2)

colData <- data.frame(row.names=colnames(gene.count.1), condition)
dds <- DESeqDataSetFromMatrix(gene.count.1, colData, design= ~ condition)
dds <- DESeq(dds)
# 查看一下dds的内容
dds
#接下来，我们要查看case versus control的总体结果，并根据padj进行重新排序。利用summary命令统计显示一共多少个genes上调和下调（FDR0.1）
res = results(dds, contrast=c("condition", "PN", "H"))
#或下面命令
res= results(dds)
res = res[order(res$padj),]
head(res)
summary(res)
#所有结果先进行输出,未筛选
#write.csv(res,file="All_results.csv")

####将表达量矩阵按照差异结果进行排序
exp_data<-gene.count.1[rownames(res),]
head(exp_data)
####计算正常和疾病样本表达均值(这里是成对的正常样本和癌旁)
colnum<-ncol(exp_data)
test_PN<-exp_data[,1:3]
test_H<-exp_data[,4:6]

exp_data$mean_PN<-apply(exp_data[,1:3],1,mean)
exp_data$mean_H<-apply(exp_data[,4:6],1,mean)

###合并计算出来的均值和差异结果
res_addMean<-data.frame(exp_data,res)
#####输出未筛选的差异结果
write.table(res_addMean,file="gene.DESeq2.PNvsH.compare.nofilter.xls",quote = F,sep='\t',row.names=T,col.names =NA )
#####提取筛选后的差异表达基因
##1.去掉低表达
head(res_addMean)
deseq_gene<-res_addMean[apply(res_addMean[,1:6],1,function(x){XXXX<-FALSE;if(any(x>5)){XXXX<-TRUE};return(XXXX)}),]
dim(deseq_gene)
#write.table(deseq_gene,file="DESeq2.PNvsH.compare.higher5.xls",quote = F,sep='\t',row.names=T,col.names =NA )
##2.根据padj和FC进一步筛选
diff_gene_deseq2 <-subset(deseq_gene, padj < 0.05 & abs(log2FoldChange) > 0.585)
diff_gene_deseq2 <-subset(deseq_gene, pvalue < 0.05 & abs(log2FoldChange) > 0.585)
dim(diff_gene_deseq2)
##输出筛选后的差异结果
write.table(diff_gene_deseq2,file="gene.DESeq2.PNvsH.compare.p0.05.FC1.5.Higher5.xls",
            quote = F,sep='\t',row.names=T,col.names =NA )

diff_gene_deseq2 <-subset(deseq_gene, pvalue < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq2)
##输出筛选后的差异结果
write.table(diff_gene_deseq2,file="gene.DESeq2.PNvsH.compare.p0.05.FC2.Higher5.xls",
            quote = F,sep='\t',row.names=T,col.names =NA )
#######################################################################
