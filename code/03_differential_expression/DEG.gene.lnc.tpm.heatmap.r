
###根据行注释提取出基因group，可以去做富集分析，后期加到热图的旁边
gene_group=pht2$annotation_row
gene_group$gene=rownames(gene_group)
gene_group=pht2$annotation_row
gene_group$gene=rownames(gene_group)
head(gene_group)
write.table(gene_group,
            "gene.cluster.txt", 
            sep="\t", col.names=T,quote = F)

library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])
save(allcluster_go,file="allcluster_go.RData")
write.table(allcluster_go,
            "BP.genes_branched.cluster.txt", 
            sep="\t", col.names=T,quote = F)
###kegg
library(clusterProfiler)
library(org.Hs.eg.db)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", 
               toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichKEGG(gene         = unique(df_name$ENTREZID),
                   organism     = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.2)
  #go_res=go@result
  kegg <- setReadable(go, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
  go_res=kegg@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])

save(allcluster_go,file="allcluster_kegg.RData")
write.table(allcluster_go,
            "KEGG.genes_branched.cluster.txt", 
            sep="\t", col.names=T,quote = F)

#####################################
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后")
list.files()
expr<-read.delim("all_sample_count.addAnotation.txt",header=T,row.names = 1)
head(expr);dim(expr)
expr<-subset(expr,gene_biotype %in% c("protein_coding","lncRNA"))
expr[1:4,1:4]
dim(expr)
expr.pro<-subset(expr,gene_biotype %in% c("protein_coding"))
expr.lnc<-subset(expr,gene_biotype %in% c("lncRNA"))

#####DEG gene
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\DEG")
list.files(pattern = "*.txt")

data<-read.delim("all.gene.DEG.union.7840.list.addsymbol.txt" ,header=T,row.names = 1)
data.pro<-df[rownames(data),]

library(dplyr)
#data.1<-data %>% subset(core >=0.5)
head(data.pro);dim(data.pro)

#exp.tpm.1.pro<- exp.tpm.1 %>% subset(gene_biotype=="protein_coding")
#exp.tpm.1.lnc<- exp.tpm.1 %>% subset(gene_biotype=="lncRNA")

#exp.tpm.1.pro<-data.pro[,c(3:32)]

#phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",header=T,row.names = NULL)
#head(phenotype);dim(phenotype)
##########heatmap protein_coding
final_data<-data.pro
final_data<-apply(final_data,2,as.numeric)
rownames(final_data)<-rownames(data.pro)
final_data[1:4,1:4]
str(final_data)
final_data<-as.data.frame(final_data)
str(final_data)

data.m<-apply(final_data,1,scale)
data.m<-t(data.m)
colnames(data.m)<-colnames(final_data)
library("ComplexHeatmap")
library(circlize)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

#################
phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
head(phenotype)
annot_df<-data.frame(type=phenotype$Time)
#info$stim<-factor(info$stim,levels=c("NC","D1.5","D7"))
rownames(phenotype)<-phenotype$sample

ha2 = HeatmapAnnotation(bar = phenotype$Sex,
                        col = list(bar = c("female" = "Salmon", "male" = "DodgerBlue1")))
data.m<-na.omit(data.m)
#data<-data[rownames(data.m),]
#dim(data);dim(data.m)
#ha = rowAnnotation(bar = data$Type,
#                  col = list(bar = c("Cytotoxic" = "Orange1", "Exhausted" = "LightYellow4")))

#ha = HeatmapAnnotation(bar = phenotype$Time,
#                  col = list(type = c(rainbow(length(unique(phenotype$Time)))))
#)

#annot_df <- data.frame(Sex = phenotype$Sex, Time = phenotype$Time)
annot_df <- data.frame(Sex = phenotype$Sex)
rownames(annot_df)<-rownames(phenotype)
col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
           Time = c("D1" = "Firebrick1","D15" = "RoyalBlue1",
                    "D29" = "orange","D43" = "LimeGreen") )
col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"))
ha <- HeatmapAnnotation(Sex = phenotype$Sex,  
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1")))

pdf(file = "heatmap.7839DEG.gene.up.expTPM.scaled.addlegend.pdf", 6.2, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()
#######DEG lnc
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\DEG")
list.files(pattern = "*.txt")

data<-read.delim("all.lnc.DEG.union.4802.list.addsymbol.txt" ,header=T,row.names = 1)
data.pro<-df[rownames(data),]

library(dplyr)
#data.1<-data %>% subset(core >=0.5)
head(data.pro);dim(data.pro)

#exp.tpm.1.pro<- exp.tpm.1 %>% subset(gene_biotype=="protein_coding")
#exp.tpm.1.lnc<- exp.tpm.1 %>% subset(gene_biotype=="lncRNA")
####heatmap lnc
final_data<-data.pro
final_data<-apply(final_data,2,as.numeric)
rownames(final_data)<-rownames(data.pro)
final_data[1:4,1:4]
str(final_data)
final_data<-as.data.frame(final_data)
str(final_data)

data.m<-apply(final_data,1,scale)
data.m<-t(data.m)
colnames(data.m)<-colnames(final_data)
library("ComplexHeatmap")
library(circlize)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

###
phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
head(phenotype)
annot_df<-data.frame(type=phenotype$Time)
#info$stim<-factor(info$stim,levels=c("NC","D1.5","D7"))
rownames(phenotype)<-phenotype$sample

ha2 = HeatmapAnnotation(bar = phenotype$Sex,
                        col = list(bar = c("female" = "Salmon", "male" = "DodgerBlue1")))
data.m<-na.omit(data.m)
#data<-data[rownames(data.m),]
#dim(data);dim(data.m)
annot_df <- data.frame(Sex = phenotype$Sex)
rownames(annot_df)<-rownames(phenotype)
col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
           Time = c("D1" = "Firebrick1","D15" = "RoyalBlue1",
                    "D29" = "orange","D43" = "LimeGreen") )
col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"))
ha <- HeatmapAnnotation(Sex = phenotype$Sex,  
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1")))

pdf(file = "heatmap.4801DEG.lnc.expTPM.scaled.addlegend.pdf", 6.2, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()

############################富集分析
###########富集分析
data<-read.delim("all.gene.DEG.union.7840.list.addsymbol.txt" ,header=T,row.names = 1)
data$gene_name
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library("ggplot2")
data(geneList, package="DOSE")

#keygene=df.c[which(df.c$annotation1=="Promoter (<=1kb)"),]

Gene_list<- bitr(rownames(data),fromType="ENSEMBL",toType="ENTREZID",
                 OrgDb="org.Hs.eg.db")

library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

ego_KEGG <- enrichKEGG(gene         = Gene_list$ENTREZID,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.2)
ego_KEGG <- setReadable(ego_KEGG, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
write.table(ego_KEGG,file="DEGgene.enrichKEGG.txt",
            sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
dotplot(ego_KEGG,showCategory=40,title="KEGG Enrichment")
p1<-dotplot(ego_KEGG,showCategory=40,title="KEGG Enrichment")
ggsave(filename = "DEGgene.enrichKEGG.pdf", plot =p1,
       width = 14, height = 12, units = 'cm')


ego_bp <- enrichGO(gene         = Gene_list$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
ego_bp <- setReadable(ego_bp, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

write.table(ego_bp,file="DEGgene.enrichGO.txt",
            sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
dotplot(ego_bp,showCategory=40,title="GO Enrichment")
#p1<-dotplot(ego_bp,showCategory=40,title="BP Enrichment")
#ggsave(filename = "enrichBP.DEG.final.pdf", 
#      plot =p1,width = 19, height = 19, units = 'cm')
#GO的三种条目分开画
go<-read.delim('DEGgene.enrichGO.txt', stringsAsFactors = FALSE)
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
ggsave(filename="DEGgene.enrichGO.barplot.pdf",plot=p,
       device='pdf',path=".",width=5.5,height=6.5)

