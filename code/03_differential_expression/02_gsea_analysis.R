# =============================================================================
# 06. Gene Set Enrichment Analysis (GSEA)
# Description: GSEA using MSigDB Hallmark, GO, and KEGG gene sets
# Input: all_sample_tpm.protein_coding.txt
# Output: GSEA enrichment score matrices and heatmaps
# =============================================================================

list.files()
library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  


library(Seurat)
library(cowplot)
library(scater)
library(Matrix)
library(sva)
library(scran)
library(ggplot2)
library(dplyr)
library(GSVA) 
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(dplyr)
library(pheatmap)
list.files()

expr<-read.delim("all_sample_tpm.protein_coding.txt",header=T,row.names = 1)
head(expr);dim(expr)

#table(Idents(immune.combined.new2.filter))
#DimPlot(immune.combined.new2.filter,label = T)

#celltype.list<-unique(immune.combined.new2.filter@meta.data$celltye)
###Plasmablasts
#for(i in 1:length(celltype.list)){
#immune.combined.new2.1<- immune.combined.new2.filter %>% subset(celltye == celltype.list[i])
#Idents(immune.combined.new2.1)<-'orig.ident'
#av <-AverageExpression(immune.combined.new2.1 , 
 #                      assays = "RNA")
av[[1]] 
av=av[[1]]
expr.1<-expr[,c(1:30)]
cg=names(tail(sort(apply(expr.1, 1, sd)),1000)) 
#可以得到如下所示的各个单细胞亚群的表达量相关性
#pdf("celltype.subsets.correlation.pdf",width = 4,height = 4)
#pheatmap::pheatmap(cor(as.data.frame(av[cg,]))) 
#dev.off()

######################
head(av[cg,])
#Tex-prog   Tex-eff  Tex-term
#HMGA1  1.2481908 1.0714660 0.6828476
#GINS2  0.6956570 0.2665502 0.1449244

head(expr)
expr<-expr[,c(1:30)]
#####MF
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category='H') 
msigdbr_collections()
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category = "C5", subcategory =   "GO:MF"  ) 
# geneList = av[,1] 
gl = apply(expr, 2, function(geneList){
  geneList=sort(geneList,decreasing = T)
  print(head(geneList))
  print(tail(geneList))
  geneList=geneList[geneList>0.1]
  egmt <- GSEA(geneList, TERM2GENE= all_gene_sets[,c('gs_name','ensembl_gene')] , 
               minGSSize = 20, 
               pvalueCutoff = 1,
               verbose=FALSE)
  head(egmt)
  egmt@result 
  gsea_results_df <- egmt@result 
  return(gsea_results_df) 
  
})
path = unique(all_gene_sets$gs_name)
es.max <- do.call(cbind,
                  lapply(gl, function(x){
                    x[path,'enrichmentScore']
                  }))
rownames(es.max) = path
head(es.max)  
es.max=na.omit(es.max)
#pheatmap::pheatmap(es.max,show_colnames =T,show_rownames = F) 
###简单挑选了一下各个单细胞亚群特异性的结果，代码如下所示：
#每个单细胞亚群的特异性top5基因集的 富集分析结果

df = do.call(rbind,
             lapply(1:ncol(es.max), function(i){
               dat= data.frame(
                 path  = rownames(es.max),
                 cluster =   colnames(es.max)[i],
                 sd.1 = es.max[,i], #每一列原值
                 sd.2 = apply(es.max[,-i], 1, median)  #除当列以外每行的中位值
               )
             })) 
#df$fc = (df$sd.1)/(df$sd.2)#两值相减，变化越大说明越有意义（从中挑出top5）
#df$logfc<-log2(df$fc)

df$fc = df$sd.1 - df$sd.2#两值相减，变化越大说明越有意义（从中挑出top5）
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)#找出每个细胞类型的前五个通路

head(df)
#top5 <- df %>% group_by(cluster) %>% top_n(20, abs(logfc))#找出每个细胞类型的前五个通路
n=es.max[unique(top5$path),]

write.table(n,
            paste(celltype.list[i],".MF.txt",sep=""), 
            sep="\t", col.names=T,quote = F)

rownames(n)
rownames(n)[grepl('FGFR',rownames(n))]
rownames(n)=gsub('GOMF_','',rownames(n))
rownames(n)=substring(rownames(n),1,40)
pdf(paste(celltype.list[i],"GSEA.MF.pdf",sep=""),width = 6,height = 5)
pheatmap(n,show_rownames = T,cluster_cols = FALSE,scale = "row") 
dev.off()

#####BP
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category='H') 
msigdbr_collections()
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category = "C5", subcategory =   "GO:BP"  ) 
# geneList = av[,1] 
gl = apply(av, 2, function(geneList){
  geneList=sort(geneList,decreasing = T)
  print(head(geneList))
  print(tail(geneList))
  geneList=geneList[geneList>0.1]
  egmt <- GSEA(geneList, TERM2GENE= all_gene_sets[,c('gs_name','gene_symbol')] , 
               minGSSize = 20, 
               pvalueCutoff = 1,
               verbose=FALSE)
  head(egmt)
  egmt@result 
  gsea_results_df <- egmt@result 
  return(gsea_results_df) 
  
})
path = unique(all_gene_sets$gs_name)
es.max <- do.call(cbind,
                  lapply(gl, function(x){
                    x[path,'enrichmentScore']
                  }))
rownames(es.max) = path
head(es.max)  
es.max=na.omit(es.max)
#pheatmap::pheatmap(es.max,show_colnames =T,show_rownames = F) 
###简单挑选了一下各个单细胞亚群特异性的结果，代码如下所示：
#每个单细胞亚群的特异性top5基因集的 富集分析结果
df = do.call(rbind,
             lapply(1:ncol(es.max), function(i){
               dat= data.frame(
                 path  = rownames(es.max),
                 cluster =   colnames(es.max)[i],
                 sd.1 = apply(es.max[,c(1:3)], 1, median), #每一列原值,改为AD均值
                 sd.2 = apply(es.max[,c(4,5)], 1, median)  #除当列以外每行的中位值，改为AS均值
               )
             })) 
df$fc = (df$sd.1)/(df$sd.2)#两值相减，变化越大说明越有意义（从中挑出top5）
df$logfc<-log2(df$fc)
head(df)
top5 <- df %>% group_by(cluster) %>% top_n(20, abs(logfc))#找出每个细胞类型的前五个通路
n=es.max[unique(top5$path),]

write.table(n,
            paste(celltype.list[i],".BP.txt",sep=""), 
            sep="\t", col.names=T,quote = F)

rownames(n)
rownames(n)[grepl('FGFR',rownames(n))]
rownames(n)=gsub('GOBP_','',rownames(n))
rownames(n)=substring(rownames(n),1,40)
pdf(paste(celltype.list[i],"GSEA.BP.pdf",sep=""),width = 6,height = 6)
pheatmap(n,show_rownames = T,cluster_cols = FALSE,scale = "row")  
dev.off()


#####KEGG
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category='H') 
msigdbr_collections()
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category = "C2", subcategory =   "CP:KEGG"  ) 
# geneList = av[,1] 
gl = apply(av, 2, function(geneList){
  geneList=sort(geneList,decreasing = T)
  print(head(geneList))
  print(tail(geneList))
  geneList=geneList[geneList>0.1]
  egmt <- GSEA(geneList, TERM2GENE= all_gene_sets[,c('gs_name','gene_symbol')] , 
               minGSSize = 20, 
               pvalueCutoff = 1,
               verbose=FALSE)
  head(egmt)
  egmt@result 
  gsea_results_df <- egmt@result 
  return(gsea_results_df) 
  
})
path = unique(all_gene_sets$gs_name)
es.max <- do.call(cbind,
                  lapply(gl, function(x){
                    x[path,'enrichmentScore']
                  }))
rownames(es.max) = path
head(es.max)  
es.max=na.omit(es.max)
#pheatmap::pheatmap(es.max,show_colnames =T,show_rownames = F) 
###简单挑选了一下各个单细胞亚群特异性的结果，代码如下所示：
#每个单细胞亚群的特异性top5基因集的 富集分析结果
df = do.call(rbind,
             lapply(1:ncol(es.max), function(i){
               dat= data.frame(
                 path  = rownames(es.max),
                 cluster =   colnames(es.max)[i],
                 sd.1 = apply(es.max[,c(1:3)], 1, median), #每一列原值,改为AD均值
                 sd.2 = apply(es.max[,c(4,5)], 1, median)  #除当列以外每行的中位值，改为AS均值
               )
             })) 
df$fc = (df$sd.1)/(df$sd.2)#两值相减，变化越大说明越有意义（从中挑出top5）
df$logfc<-log2(df$fc)
head(df)
top5 <- df %>% group_by(cluster) %>% top_n(20, abs(logfc))#找出每个细胞类型的前五个通路
n=es.max[unique(top5$path),]

write.table(n,
            paste(celltype.list[i],".KEGG.txt",sep=""), 
            sep="\t", col.names=T,quote = F)

rownames(n)
rownames(n)[grepl('FGFR',rownames(n))]
rownames(n)=gsub('KEGG_','',rownames(n))
rownames(n)=substring(rownames(n),1,40)
pdf(paste(celltype.list[i],"GSEA.KEGG.pdf",sep=""),width = 6,height = 6)
pheatmap(n,show_rownames = T,cluster_cols = FALSE,scale = "row") 
dev.off()
}
