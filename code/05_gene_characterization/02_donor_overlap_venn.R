# =============================================================================
# 02_donor_overlap_venn
# Description: Donor Overlap Venn Diagram
# =============================================================================

##
df1<-df.pro[,phenotype.d4$sample]
df1$geneid<-rownames(df1)
head(df1)
##多数据集求交集
gene4 <- Reduce(intersect,list(
  df1[df1[,1]>=100,5], df1[df1[,2]>=100,5], 
  df1[df1[,3]>=100,5], df1[df1[,4]>=100,5]))
length(gene4)
gene4
##
df1<-df.pro[,phenotype.d5$sample]
df1$geneid<-rownames(df1)
head(df1)
##多数据集求交集
gene5 <- Reduce(intersect,list(
  df1[df1[,1]>=100,5], df1[df1[,2]>=100,5], 
  df1[df1[,3]>=100,5], df1[df1[,4]>=100,5]))
length(gene5)
gene5
##
df1<-df.pro[,phenotype.d6$sample]
df1$geneid<-rownames(df1)
head(df1)
##多数据集求交集
gene6 <- Reduce(intersect,list(
  df1[df1[,1]>=100,5], df1[df1[,2]>=100,5], 
  df1[df1[,3]>=100,5], df1[df1[,4]>=100,5]))
length(gene6)
gene6
##
df1<-df.pro[,phenotype.d11$sample]
df1$geneid<-rownames(df1)
head(df1)
##多数据集求交集
gene11 <- Reduce(intersect,list(
  df1[df1[,1]>=100,5], df1[df1[,2]>=100,5], 
  df1[df1[,3]>=100,5], df1[df1[,4]>=100,5]))
length(gene11)
gene11
##
df1<-df.pro[,phenotype.d12$sample]
df1$geneid<-rownames(df1)
head(df1)
##多数据集求交集
gene12 <- Reduce(intersect,list(
  df1[df1[,1]>=100,5], df1[df1[,2]>=100,5], 
  df1[df1[,3]>=100,5], df1[df1[,4]>=100,5]))
length(gene12)
gene12
##
df1<-df.pro[,phenotype.d13$sample]
df1$geneid<-rownames(df1)
head(df1)
##多数据集求交集
gene13 <- Reduce(intersect,list(
  df1[df1[,1]>=100,5], df1[df1[,2]>=100,5], 
  df1[df1[,3]>=100,5], df1[df1[,4]>=100,5]))
length(gene13)
gene13

####venn
library(venn)
library(VennDiagram)

venn_list <- list(gene4, gene5, gene6,
                  gene11, gene12,gene13)   # 制作韦恩图搜所需要的列表文件
names(venn_list) <- c("donor4","donor5","donor6",
                      "donor11","donor12","donor13")    # 把列名赋值给列表的key值
venn_list = purrr::map(venn_list,na.omit)      # 删除列表中每个向量中的NA

#作图
venn(venn_list,
     zcolor='style', # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
     opacity = 0.3,  # 调整颜色透明度
     box = F,        # 是否添加边框
     ilcs = 1.4,     # 数字大小
     sncs = 1.5,        # 组名字体大小
     ilabels = "counts"
)

#save as donor12.high.exp.gene.overlap
# 更多参数 ?venn查看

# 查看交集详情,并导出结果
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )
write.table(inter, "./6donors.high.exp.gene.venn.result.csv", 
            row.names = FALSE, sep = ',', quote = FALSE)

gene.inter <- Reduce(intersect,list(
  gene4, gene5, gene6,
  gene11, gene12,gene13))
length(gene.inter)
gene.inter

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library("ggplot2")
data(geneList, package="DOSE")

Gene_list.1<- bitr(gene.inter,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
Gene_list.2<- bitr(gene.inter,fromType="ENSEMBL",toType="SYMBOL",OrgDb="org.Hs.eg.db")
gene376<-cbind(Gene_list.1,Gene_list.2)
write.table(gene376,
            "gene376.info.6donors.overlap.txt", 
            sep="\t", col.names=T,quote = F)

###D1 time samples high exp gene
#####高表达并集
df1<-df.pro[,phenotype.d1$sample]
df1$geneid<-rownames(df1)
##多数据集求并集
genes <- Reduce(intersect,list(
  df1[df1[,1]>=100,13], df1[df1[,2]>=100,13], 
  df1[df1[,3]>=100,13], df1[df1[,4]>=100,13],
  df1[df1[,5]>=100,13],df1[df1[,6]>=100,13],
  df1[df1[,7]>=100,13],df1[df1[,8]>=100,13],
  df1[df1[,9]>=100,13],df1[df1[,10]>=100,13],
  df1[df1[,11]>=100,13],df1[df1[,12]>=100,13]))
length(genes)
gene420<-genes
length(intersect(gene376$ENSEMBL,gene420))
gene<-intersect(gene376$ENSEMBL,gene420)

Gene_list.1<- bitr(gene,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
Gene_list.2<- bitr(gene,fromType="ENSEMBL",toType="SYMBOL",OrgDb="org.Hs.eg.db")
gene335<-cbind(Gene_list.1,Gene_list.2)
write.table(gene335,
            "gene335.info.6donors.and.D1.samples.overlap.txt", 
            sep="\t", col.names=T,quote = F)
####venn
library(venn)
library(VennDiagram)

venn_list <- list(gene376$ENSEMBL,
                  gene420)   # 制作韦恩图搜所需要的列表文件
names(venn_list) <- c("6donors.high.exp.genes",
                      "D1.samples.high.exp.genes")    # 把列名赋值给列表的key值
venn_list = purrr::map(venn_list,na.omit)      # 删除列表中每个向量中的NA

#作图
venn(venn_list,
     zcolor='style', # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
     opacity = 0.3,  # 调整颜色透明度
     box = F,        # 是否添加边框
     ilcs = 1.4,     # 数字大小
     sncs = 1.5,        # 组名字体大小
     ilabels = "counts"
)

#save as donor12.high.exp.gene.overlap
# 更多参数 ?venn查看

##################################################
####富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library("ggplot2")
data(geneList, package="DOSE")

#keygene=df.c[which(df.c$annotation1=="Promoter (<=1kb)"),]
head(gene335)
Gene_list<- bitr(gene335$ENSEMBL,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

ego_KEGG <- enrichKEGG(gene         = Gene_list$ENTREZID,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.2)
ego_KEGG <- setReadable(ego_KEGG, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
write.table(ego_KEGG,file="./个体内和个体间都高表达的基因/all.highexp.gene.enrichKEGG.txt",
            sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
dotplot(ego_KEGG,showCategory=40,title="KEGG Enrichment")
p1<-dotplot(ego_KEGG,showCategory=60,title="KEGG Enrichment")
ggsave(filename = "./个体内和个体间都高表达的基因/all.highexp.gene.enrichKEGG.pdf", plot =p1,
       width = 13, height = 17, units = 'cm')

p2<-cnetplot(ego_KEGG, categorySize="pvalue", circular = TRUE, colorEdge = TRUE)
p2

categorys <- c("Ribosome", "Oxidative phosphorylation",
               "Thermogenesis", "Diabetic cardiomyopathy")
p2<-cnetplot(ego_KEGG, categorySize="pvalue", circular = TRUE, 
         colorEdge = TRUE, showCategory = categorys)

ggsave(filename = "./个体内和个体间都高表达的基因/335gene.enrichKEGG.cnetplot.pdf",
       plot =p2,width = 20, 
       height = 17, units = 'cm')

ego_bp <- enrichGO(gene         = Gene_list$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
ego_bp <- setReadable(ego_bp, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

write.table(ego_bp,file="./个体内和个体间都高表达的基因/all.highexp.gene.enrichGO.txt",
            sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
dotplot(ego_bp,showCategory=60,title="GO Enrichment")
#p1<-dotplot(ego_bp,showCategory=40,title="BP Enrichment")
#ggsave(filename = "enrichBP.DEG.final.pdf", 
#      plot =p1,width = 19, height = 19, units = 'cm')
#GO的三种条目分开画
go<-read.delim('./个体内和个体间都高表达的基因/all.highexp.gene.enrichGO.filter.metascape.txt', 
               stringsAsFactors = FALSE)
dim(go)
#go$term <- paste(go$ID, go$Description, sep = ': ')
go$term<-go$Description

go <- go[order(go$ONTOLOGY, go$p.adjust, decreasing = c(TRUE, TRUE),method = "radix"), ]
go$term <- factor(go$term, levels = go$term)

p<-ggplot(go, aes(term, -log10(p.adjust))) +
  geom_col(aes(fill = ONTOLOGY), width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c('#D06660', '#5AAD36', '#6C85F5')) +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x = '', y = '-Log10 P-Value\n',title="GO enrichment")
p

p<-ggplot(go, aes(term, -log10(p.adjust))) +
  geom_col(aes(fill = ONTOLOGY), width = 0.5) +
  scale_fill_manual(values = c('#D06660', '#5AAD36', '#6C85F5')) +
  theme_bw()+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 60, hjust = 1,vjust =1),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  labs(x = '', y = '-Log10 P-Value\n',title="GO enrichment")
p
ggsave(filename="./个体内和个体间都高表达的基因/all.highexp.gene.enrichGO.barplot.pdf",plot=p,
       device='pdf',path=".",width=11,height=6)

#########box plot
exp.tpm.1<-exp.tpm[c(gene335$ENSEMBL),]
exp.tpm.1.pro<-exp.tpm.1[,c(1:30)]
dim(exp.tpm.1.pro)

data1<-exp.tpm.1.pro
normalized_data<-as.data.frame(data1)
normalized_data$geneid<-rownames(normalized_data)
normalized_data

library("reshape")
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","tpm")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype,by = "sample",all=T)
head(exp_melt.final)
#exp_melt.final.1<-exp_melt.final %>% subset(scaled_count>=-0.3 & scaled_count<=0.2) 
dim(exp_melt.final);
#dim(exp_melt.final.1)

head(exp_melt.final)

library("ggpubr")
library("reshape")
library("ggplot2")

####by sample
exp_melt.final$logtpm<-log2(exp_melt.final$tpm)
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt.final, x = "sample", y = "logtpm",
               fill = "DarkTurquoise",
               outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  xlab("") + ylab("Expression (log2TPM)")+
  theme(axis.text.x = 
          element_text(size = 14, face = "bold",angle = 45,
                       hjust = 1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=14,colour="black",face = "bold"),
        legend.title=element_text(size=12,colour="black",face = "bold"),
        legend.text=element_text(size=12),
        axis.title.x = element_text(size=17,face = "bold"),
        axis.title.y = element_text(size=15,face = "bold"),
        plot.title = element_text(size=16),legend.key.size = unit(4,'mm'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
p
#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="./个体内和个体间都高表达的基因/335gene.boxplot.logTPM.pdf",plot=p,
       device='pdf',path=".",width=8,height=3.2)

