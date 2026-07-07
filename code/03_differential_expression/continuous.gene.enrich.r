
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后")
list.files()
expr<-read.delim("all_sample_count.addAnotation.txt",header=T,row.names = 1)
head(expr);dim(expr)
expr<-subset(expr,gene_biotype %in% c("protein_coding","lncRNA"))
expr[1:4,1:4]
dim(expr)
expr.pro<-subset(expr,gene_biotype %in% c("protein_coding"))
expr.lnc<-subset(expr,gene_biotype %in% c("lncRNA"))
##########continue up vs down进行比较
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library("ggplot2")
data(geneList, package="DOSE")

setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\DEG")
list.files(pattern = "*.xls")

data.up<-read.delim("all.gene.FC.final.continue.up.final.xls" ,header=T,row.names = 1)
data.down<-read.delim("all.gene.FC.final.continue.down.final.xls" ,header=T,row.names = 1)

data.pro.up<-expr.pro[rownames(data.up),]
data.pro.down<-expr.pro[rownames(data.down),]

head(data.pro.up)
head(data.pro.down)


####Up and down 2 files
keytypes(org.Hs.eg.db)

Gene_list1<- bitr(rownames(data.pro.up),fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #up
Gene_list2<- bitr(rownames(data.pro.up),fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #down

cp = list(Up=Gene_list1$ENTREZID,
          Down=Gene_list2$ENTREZID)  # 合并两个数据集，并转换为列表
## 最新版本需要加上这个命令
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

xx <- compareCluster(cp, fun="enrichKEGG", organism="hsa", 
                     pvalueCutoff=0.1,pAdjustMethod = "BH",qvalueCutoff = 0.5) 

dotplot(xx,showCategory=25,label_format=40,color = "pvalue")
kegg <- setReadable(xx, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
write.table(kegg,file="upAndDown.enrich.KEGG.clusterProfiler.txt",
            sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)


p1<-dotplot(xx,showCategory=16,label_format=40,includeAll=TRUE,
            title="KEGG Enrichment Comparison")
ggsave(filename = "KEGG.upAndDown.compare.clusterProfiler.pdf", 
       plot =p1,width = 15,
       height = 15, units = 'cm')

#
xx <- compareCluster(cp, fun="enrichGO", OrgDb='org.Hs.eg.db',ont= "BP",
                     pvalueCutoff=0.05, pAdjustMethod = "BH",qvalueCutoff = 0.1) 
dotplot(xx,showCategory=25,label_format=40,includeAll=TRUE,
        title="BP Enrichment Comparison")
p1<-dotplot(xx,showCategory=15,label_format=40,
            includeAll=TRUE,title="BP Enrichment Comparison")
ggsave(filename = "GOBP.upAndDown.compare.clusterProfiler.pdf", plot =p1,
       width = 15, height = 15, units = 'cm')

bp <- setReadable(xx, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
write.table(bp,file="enrich.upAndDown.clusterProfiler.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)

#
xx1 <- compareCluster(cp, fun="enrichGO", OrgDb='org.Hs.eg.db',ont= "MF",
                      pvalueCutoff=0.05, pAdjustMethod = "BH",qvalueCutoff = 0.05) 
xx2 <- compareCluster(cp, fun="enrichGO", OrgDb='org.Hs.eg.db',ont= "CC",
                      pvalueCutoff=0.05, pAdjustMethod = "BH",qvalueCutoff = 0.05) 

dotplot(xx1,showCategory=20,label_format=40,includeAll=TRUE)
p2<-dotplot(xx1,showCategory=20,label_format=40,includeAll=TRUE,
            title="MF Enrichment Comparison")
ggsave(filename = "GOMF.5cluster.compare.pdf", plot =p2,
       width = 15, height = 15, units = 'cm')

dotplot(xx2,showCategory=20,label_format=40,includeAll=TRUE)
p2<-dotplot(xx2,showCategory=20,label_format=40,includeAll=TRUE,title="CC Enrichment Comparison")
ggsave(filename = "GOCC.5cluster.compare.pdf", plot =p2,
       width = 15, height = 15, units = 'cm')

mf <- setReadable(xx1, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
write.table(mf,file="enrich.MF.5cluster.txt",sep='\t',
            quote=FALSE,row.names=FALSE,col.names=TRUE)

cc <- setReadable(xx2, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
write.table(cc,file="enrich.CC.5cluster.txt",sep='\t',
            quote=FALSE,row.names=FALSE,col.names=TRUE)

################################用自己代码绘制
#data<-read.delim("Exhausted_T.CD19VSCD26.CD19VSCD44.consistent_up_and_down.union.enrichKEGG.filter2.txt",header=T,sep='\t')
data<-read.delim('enrich.DEG_downAandUp.KEGG.new.txt', stringsAsFactors = FALSE)
head(data)
data
data<-data[which(data$Count>=4),]
library(ggplot2)
head(data)
data$LogFDR<--log10(data$pvalue)
data <- data[order(data$Cluster,data$LogFDR,decreasing = c(TRUE, FALSE),method = "radix"), ]
data$Description <- factor(data$Description, levels = data$Description)

pp = ggplot(data,aes(Cluster,Description))
p.final<-pp + geom_point(aes(color=(LogFDR),size=Count))+theme_bw()+
  labs(x="Cluster",y="",title="KEGG enrichment")+
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=11,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=12,colour="black"),
        axis.title.y = element_text(size=12,colour="black"),
        plot.title = element_text(size=13,colour="black"),
        legend.key.size = unit(5,'mm'))+
  scale_colour_gradient(low="DodgerBlue",high="NavyBlue",name="-log P-Value")

p.final
ggsave(filename="KEGG_upAndDown.pdf",plot=p.final,
       device='pdf',path=".",width=6.5,height=5)
