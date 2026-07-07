# =============================================================================
# 03_temporal_heatmap
# Description: Temporal Expression Heatmap
# =============================================================================

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
list.files()
expr<-read.delim("all_sample_count.addAnotation.txt",header=T,row.names = 1)
head(expr);dim(expr)
expr<-subset(expr,gene_biotype %in% c("protein_coding","lncRNA"))
expr[1:4,1:4]
dim(expr)
exp.tpm<-read.delim("all_sample_tpm.txt",header=T,row.names = 1)
dim(exp.tpm)
exp.tpm[1:4,1:4]
exp.tpm<-exp.tpm[rownames(expr),]
dim(exp.tpm)
identical(rownames(expr),rownames(exp.tpm))
exp.tpm$gene_biotype<-expr$gene_biotype
head(exp.tpm)
write.table(exp.tpm,
            "exp.tpm.addbiotype.txt", 
            sep="\t", col.names=T,quote = F)

list.files()
data<-read.csv("mfuzz_male_3.csv",header = TRUE,row.names = 1)
colnames(data)<-c("core")
library(dplyr)
data.1<-data %>% subset(core >=0.5)
head(data.1);dim(data.1)
exp.tpm.1<-exp.tpm[rownames(data.1),]
exp.tpm.1.pro<- exp.tpm.1 %>% subset(gene_biotype=="protein_coding")
exp.tpm.1.lnc<- exp.tpm.1 %>% subset(gene_biotype=="lncRNA")
head(exp.tpm.1.pro)
dim(exp.tpm.1.pro)
exp.tpm.1.pro<-exp.tpm.1.pro[,c(2:31)]

phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)

##########heatmap protein_coding
final_data<-exp.tpm.1.pro
final_data<-apply(final_data,2,as.numeric)
rownames(final_data)<-rownames(exp.tpm.1.pro)
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

pdf(file = "heatmap.male_3.expTPM.scaled.pdf", 5.5, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()

######################################################################################3
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

annot_df <- data.frame(Sex = phenotype$Sex, Time = phenotype$Time)
rownames(annot_df)<-rownames(phenotype)
col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
           Time = c("D1" = "Firebrick1","D15" = "RoyalBlue1",
                    "D29" = "orange","D43" = "LimeGreen") )

ha <- HeatmapAnnotation(Sex = phenotype$Sex, Time = phenotype$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "heatmap.male_3.expTPM.scaled.addlegend.pdf", 6.2, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()

##########heatmap lnc
exp.tpm.1.lnc<-exp.tpm.1.lnc[,c(2:31)]
final_data<-exp.tpm.1.lnc
final_data<-apply(final_data,2,as.numeric)
rownames(final_data)<-rownames(exp.tpm.1.lnc)
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
head(phenotype)
annot_df<-data.frame(type=phenotype$Time)
rownames(phenotype)<-phenotype$sample

ha2 = HeatmapAnnotation(bar = phenotype$Sex,
                        col = list(bar = c("female" = "Salmon", "male" = "DodgerBlue1")))
data.m<-na.omit(data.m)

ha <- HeatmapAnnotation(Sex = phenotype$Sex, Time = phenotype$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "heatmap.lnc.male_3.expTPM.scaled.addlegend.pdf", 6.2, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()
