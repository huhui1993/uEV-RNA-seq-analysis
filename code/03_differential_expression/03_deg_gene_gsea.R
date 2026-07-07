setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\DEG")
#教程：https://mp.weixin.qq.com/s/7DxD7JOa9ykAeFT1F9t1HA
BiocManager::install("ReactomePA")
library(ReactomePA) #加载所需R包
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(DOSE)
devtools::install_github("ToledoEM/msigdf")
library(msigdf)
#

list.files()

data<-read.delim("all.gene.FC.final.continue.up_and_down.final.xls",
                 header=T,sep='\t',row.names = NULL)
head(data)
#rownames(data)<-data$gene_name
#data$log2FoldChange.new<--data$log2FoldChange
#head(data)

genelist_input<-data[,c("geneid","log2FoldChange")]
head(genelist_input)
genename <- as.character(genelist_input[,1]) #提取第一列基因名
gene_map <- select(org.Hs.eg.db, keys=genename, keytype="ENSEMBL", columns=c("SYMBOL")) #将SYMBOL格式的ID换成ENTREZ格式的ID。

non_duplicates_idx <- which(duplicated(gene_map$SYMBOL) == FALSE)
gene_map <- gene_map[non_duplicates_idx, ] #去除重复值
colnames(gene_map)[1]<-"Gene" #输出mapping结果
head(gene_map)

#将ENTREZID与logFC对应起来，并根据最后一个时间点相较于第一个时间点的logFC的值降序排列，最终生成结果如图所示。
colnames(genelist_input)[1]<-"Gene"
temp<-inner_join(gene_map,genelist_input,by = "Gene")
head(temp)
temp<-temp[,-1]
#temp<-temp[,-2]
temp<-na.omit(temp)
temp$log2FoldChange<-sort(temp$log2FoldChange,decreasing = T)

#最后我们将文件内容整理成GSEA分析所需的格式（如图所示），我们就可以开始GSEA分析了！
head(temp)
geneList = temp[,2]
names(geneList) = as.character(temp[,1])
geneList
#352999        86     10900    644353     10235    199786       177    388849      5026     11189     51412 
#4.2745453 4.0051208 3.6282025 3.1360233 3.1177622 3.0793718 3.0790690 2.9663838 2.9397800 2.7978247 2.7736075 

#GSEA富集人的c2通路
c2 <- msigdf.human %>% 
  filter(category_code == "c2") %>% dplyr::select(geneset, symbol) %>% as.data.frame
head(c2)
colnames(c2)<-c("ont","gene")
head(geneList)
gsea_c2_human <- GSEA(geneList, TERM2GENE = c2, verbose=FALSE, pvalueCutoff =0.05)
length(gsea_c2_human@result$ID)
write.table(gsea_c2_human, 'gsea_human_KEGG.c2.continue.up.and.down.txt', 
            sep = '\t', row.names = FALSE, quote = FALSE)
library(DOSE)
gseaplot(gsea_c2_human, 1)

pdf(file="GSEA_KEGG.continue.up.and.down.all.pdf",4.5,5.2)
gseaplot2(gsea_c2_human, 1:5, title = "KEGG pathway in continue.gene", 
          pvalue_table = FALSE) 
dev.off()

#GSEA富集人的c5通路
c5 <- msigdf.human %>% 
  filter(category_code == "c5") %>% dplyr::select(geneset, symbol) %>% as.data.frame
head(c5)
colnames(c5)<-c("ont","gene")
head(geneList)
gsea_c5_human <- GSEA(geneList, TERM2GENE = c5, verbose=FALSE, pvalueCutoff = 0.05)
#gsea_c5_human@result$ID
length(gsea_c5_human@result$ID)
write.table(gsea_c5_human, 'gsea_human_GO.c5.continue.up.and.down.txt', 
            sep = '\t', row.names = FALSE, quote = FALSE)

library(DOSE)
gseaplot(gsea_c5_human, 1)

pdf(file="GSEA_GO.continue.up.and.down.all.pdf",4.5,5.2)
gseaplot2(gsea_c5_human, 1:7, title = "GO term in continue.gene", 
          pvalue_table = FALSE) 
dev.off()
#第四步：进行GSEA分析，并保存结果。
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "SYMBOL", ont="BP", pvalueCutoff=0.1)   #使用GSEA进行GO富集分析（'org.Hs.eg.db'：对应物种的数据库；ont：选择输出条目，可选“BP,MF,CC或者ALL”，pvalueCutoff：设置P的阈值）
KEGG_gseresult <- gseKEGG(geneList, pvalueCutoff=0.1) #使用GSEA进行KEGG富集分析
#保存富集分析结果
go_results<-as.data.frame(Go_gseresult)
kegg_results<-as.data.frame(KEGG_gseresult)
write.csv (go_results, file ="Go_gseresult.csv")
write.csv (kegg_results, file ="KEGG_gseresult.csv")

#第五步：结果可视化（以GO富集结果为例）
pdf(file="GSEA_KEGG.lowVShigh.group.pdf",4.5,5.2)
gseaplot2(KEGG_gseresult, 1:7, title = "KEGG pathway in low-risk VS high-risk group", 
          pvalue_table = FALSE)  #1:3：这表示在图上显示前3条富集结果，也可以根据自己分析需要指定输出某一条结果；Go_gseresult：GO富集分析结果；title：加上标题；pvalue_table：是否在图上显示P值列表。
dev.off()

gseaplot2(Go_gseresult, 1:6, title = "GO Biological Process in low-risk VS high-risk group", 
          pvalue_table = FALSE)

pdf(file="GSEA_KEGG.PMVSPN.group.pdf",4.5,5.2)
gseaplot2(KEGG_gseresult, 1, title = "Oocyte meiosis", 
          pvalue_table = FALSE) #1:3：这表示在图上显示前3条富集结果，也可以根据自己分析需要指定输出某一条结果；Go_gseresult：GO富集分析结果；title：加上标题；pvalue_table：是否在图上显示P值列表。
dev.off()

############################PMvsH
data<-read.delim("gene.DESeq2.PMvsH.compare.nofilter.symbol.xls",
                 header=T,sep='\t',row.names = 1)
head(data)
#rownames(data)<-data$gene_name
data$log2FoldChange.new<--data$log2FoldChange
head(data)

genelist_input<-data[,c("gene_name","log2FoldChange.new")]
head(genelist_input)
genename <- as.character(genelist_input[,1]) #提取第一列基因名
gene_map <- select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID")) #将SYMBOL格式的ID换成ENTREZ格式的ID。

non_duplicates_idx <- which(duplicated(gene_map$SYMBOL) == FALSE)
gene_map <- gene_map[non_duplicates_idx, ] #去除重复值
colnames(gene_map)[1]<-"Gene" #输出mapping结果
head(gene_map)

#将ENTREZID与logFC对应起来，并根据logFC的值降序排列，最终生成结果如图所示。
colnames(genelist_input)[1]<-"Gene"
temp<-inner_join(gene_map,genelist_input,by = "Gene")
head(temp)
#temp<-temp[,-1]
temp<-temp[,-2]
temp<-na.omit(temp)
temp$log2FoldChange.new<-sort(temp$log2FoldChange.new,decreasing = T)

#最后我们将文件内容整理成GSEA分析所需的格式（如图所示），我们就可以开始GSEA分析了！
head(temp)
geneList = temp[,2]
names(geneList) = as.character(temp[,1])
geneList
#352999        86     10900    644353     10235    199786       177    388849      5026     11189     51412 
#4.2745453 4.0051208 3.6282025 3.1360233 3.1177622 3.0793718 3.0790690 2.9663838 2.9397800 2.7978247 2.7736075 

#GSEA富集人的c2通路
c2 <- msigdf.human %>% 
  filter(category_code == "c2") %>% dplyr::select(geneset, symbol) %>% as.data.frame
head(c2)
colnames(c2)<-c("ont","gene")
head(geneList)
gsea_c2_human <- GSEA(geneList, TERM2GENE = c2, verbose=FALSE, pvalueCutoff = 0.1)
library(DOSE)
gseaplot(gsea_c2_human, 1)

pdf(file="GSEA_KEGG.PMVSH.group.pdf",4.5,5.2)
gseaplot2(gsea_c2_human, 1:4, title = "KEGG pathway in PM VS H", 
          pvalue_table = FALSE) 
dev.off()

#GSEA富集人的c5通路
c5 <- msigdf.human %>% 
  filter(category_code == "c5") %>% dplyr::select(geneset, symbol) %>% as.data.frame
head(c5)
colnames(c5)<-c("ont","gene")
head(geneList)
gsea_c5_human <- GSEA(geneList, TERM2GENE = c5, verbose=FALSE, pvalueCutoff = 0.05)
gsea_c5_human@result$ID
library(DOSE)
gseaplot(gsea_c5_human, 1)

pdf(file="GSEA_GO.PMVSPN.group.pdf",4.5,5.2)
gseaplot2(gsea_c5_human, 1:2, title = "GO term in PM VS PN", 
          pvalue_table = FALSE) 
dev.off()