# =============================================================================
# 03_common_genes_heatmap
# Description: Common Genes Heatmap
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
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"),
               OrgDb="org.Hs.eg.db")
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
###################################################################################################
###################################################################################################
############################################################
exp.tpm<-read.delim("exp.tpm.addbiotype.txt",header=T,row.names = 1)
phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",
                      header=T,row.names = NULL)
head(phenotype);dim(phenotype)
library(dplyr)
###提取特定供体
phenotype.1<-phenotype %>% subset(Donor %in% c("4","5","6","11","12","13"))
head(phenotype.1);dim(phenotype.1)

exp.tpm.1<-exp.tpm[c("ENSG00000225528","ENSG00000293064"),]
exp.tpm.1.pro<-exp.tpm.1[,phenotype.1$sample]
dim(exp.tpm.1.pro)
#[1]  2 24

gene.tmp<-c("ENSG00000225528","ENSG00000293064")
list.1<- bitr(gene.tmp,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

list.2<- bitr(gene.tmp,fromType="ENSEMBL",toType="SYMBOL",OrgDb="org.Hs.eg.db")

list.tmp<-cbind(Gene_list,Gene_list.tmp)
write.table(list.tmp,
            "2888gene.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)

##########heatmap 2 gene
final_data<-exp.tpm.1.pro
final_data<-apply(final_data,2,as.numeric)
rownames(final_data)<-rownames(exp.tpm.1.pro)
final_data[1:4,1:4]
final_data
str(final_data)
final_data<-as.data.frame(final_data)
str(final_data)

data.m<-apply(final_data,1,scale)
data.m<-t(data.m)
colnames(data.m)<-colnames(final_data)
library("ComplexHeatmap")
library(circlize)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf(file = "./2个共有基因表达/heatmap.2gene.expTPM.scaled.pdf", 6.5, 1.6)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = TRUE)
dev.off()

pdf(file = "./2个共有基因表达/heatmap.2gene.expTPM.noscaled.pdf", 6.5, 1.6)
mycol <- colorRamp2(c(0, 2, 4), c("blue", "white", "red"))
Heatmap(as.matrix(exp.tpm.1.pro),cluster_columns = FALSE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = TRUE)
dev.off()
######################################################################################3
head(phenotype.1)
annot_df<-data.frame(type=phenotype.1$Time)
#info$stim<-factor(info$stim,levels=c("NC","D1.5","D7"))
rownames(phenotype.1)<-phenotype.1$sample

ha2 = HeatmapAnnotation(bar = phenotype.1$Sex,
                        col = list(bar = c("female" = "Salmon", "male" = "DodgerBlue1")))
data.m<-na.omit(data.m)
#data<-data[rownames(data.m),]
#dim(data);dim(data.m)
#ha = rowAnnotation(bar = data$Type,
 #                  col = list(bar = c("Cytotoxic" = "Orange1", "Exhausted" = "LightYellow4")))

#ha = HeatmapAnnotation(bar = phenotype$Time,
 #                  col = list(type = c(rainbow(length(unique(phenotype$Time)))))
#)

annot_df <- data.frame(Sex = phenotype.1$Sex, Time = phenotype.1$Time)
rownames(annot_df)<-rownames(phenotype.1)
col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
           Time = c("D1" = "Firebrick1","D15" = "RoyalBlue1",
                    "D29" = "orange","D43" = "LimeGreen") )

ha <- HeatmapAnnotation(Sex = phenotype.1$Sex, Time = phenotype.1$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "./2个共有基因表达/heatmap.2gene.expTPM.scaled.addlegend.pdf", 7, 2.2)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = TRUE,top_annotation = ha)

dev.off()
#
pdf(file = "./2个共有基因表达/heatmap.2gene.expTPM.noscaled.addlegend.pdf", 7, 2.2)
mycol <- colorRamp2(c(0, 2, 4), c("blue", "white", "red"))
Heatmap(as.matrix(exp.tpm.1.pro),cluster_columns = FALSE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = TRUE,top_annotation = ha)
dev.off()
############################################################################
######box plot
#data1<-df[apply(df,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(data1)
#log(TPM+1)
#data1<-log2(data1+1)
data1<-exp.tpm.1.pro
normalized_data<-as.data.frame(data1)
normalized_data$geneid<-rownames(normalized_data)
normalized_data

library("reshape")
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","tpm")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype.1,by = "sample",all=T)
head(exp_melt.final)
#exp_melt.final.1<-exp_melt.final %>% subset(scaled_count>=-0.3 & scaled_count<=0.2) 
dim(exp_melt.final);
#dim(exp_melt.final.1)

head(exp_melt.final)

library("ggpubr")
library("reshape")
library("ggplot2")

####by time
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt.final, x = "geneid", y = "tpm",
               fill="Time",
               palette =c("jco"),
               outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  xlab("") + ylab("Expression (log2(TPM+1))")+
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
ggsave(filename="./2个共有基因表达/boxplot.bytime.2gene.log2TPM.pdf",plot=p,
       device='pdf',path=".",width=6,height=4.5)

################add p
p <- ggboxplot(exp_melt.final, x = "geneid", y = "tpm",
               color ="Time",
               palette = "nejm",show.legend = FALSE)+
  theme_classic()+
  xlab("") + ylab("Expression (log2(TPM+1))")+
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
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p +
  stat_compare_means(comparisons = my_comparisons,label = "p.format",
                     size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2
ggsave(filename="./2个共有基因表达/boxplot.bytime.2gene.log2TPM.pdf",plot=p2,
       device='pdf',path=".",width=6,height=4.6)
#
p <- ggboxplot(exp_melt.final, x = "geneid", y = "tpm",
               fill="Time",
               palette =c("jco"),
               outlier.colour = NA,show.legend = FALSE)+
  scale_color_manual(values = c("#E31A1C","#1F78B4","#FF7F00","#33A02C"))+
  theme_classic()+
  xlab("") + ylab("Expression (log2(TPM+1))")+
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

head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")

library(dplyr)
#exp_melt.final.1<-exp_melt.final %>% subset(scaled_count>=-0.3 & scaled_count<=0.2) 
dim(exp_melt.final);
#dim(exp_melt.final.1)

p<-ggplot(exp_melt.final, aes(x=Time, y=tpm,color=Time) )+
  geom_violin(trim=FALSE,position = "dodge",scale="width") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.shape = NA)+ #绘制箱线图
  scale_color_manual(values = c("#E31A1C","#1F78B4","#FF7F00","#33A02C"))+
  theme_classic()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1,colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(fill="",x="",y="Expression (log2(TPM+1))",title="")
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p +
  stat_compare_means(comparisons = my_comparisons,label = "p.format",
                     size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2
#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
#ggsave(filename="ggplot2.violinplot.bytime.log2TPM_exp.pdf",plot=p2,
 #      device='pdf',path=".",width=4,height=3.5)
#
head(exp_melt.final)
exp_melt.final.1.a <- exp_melt.final %>% subset(geneid=="ENSG00000225528")
exp_melt.final.1.b <- exp_melt.final %>% subset(geneid=="ENSG00000293064")

p <- ggboxplot(exp_melt.final.1.a, x = "Time", y = "tpm",
               color = "Time", palette = "nejm",
               outlier.colour = NA)+
  theme_classic()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 45, hjust = 1,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(fill="type",x="",y="Expression (log2(TPM+1))",title="") 
#p
# Use only p.format as label. Remove method name.
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p + stat_compare_means(comparisons = my_comparisons,label = "p.format",
                           size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="./2个共有基因表达/ENSG00000225528.boxplot.bytime.log2TPM.addP.pdf",plot=p2,
       device='pdf',path=".",width=3.2,height=3.2)

p <- ggboxplot(exp_melt.final.1.a, x = "Donor", y = "tpm",
               color = "Donor", 
               outlier.colour = NA)+
  theme_classic()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 45, hjust = 1,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(fill="type",x="",y="Expression (log2(TPM+1))",title="") 
#p
# Use only p.format as label. Remove method name.
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p + stat_compare_means(comparisons = my_comparisons,label = "p.format",
                           size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

#
exp_melt.final.1.b <- exp_melt.final %>% subset(geneid=="ENSG00000293064")

p <- ggboxplot(exp_melt.final.1.b, x = "Time", y = "tpm",
               color = "Time", palette = "nejm",
               outlier.colour = NA)+
  theme_classic()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 45, hjust = 1,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(fill="type",x="",y="Expression (log2(TPM+1))",title="") 
#p
# Use only p.format as label. Remove method name.
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p + stat_compare_means(comparisons = my_comparisons,label = "p.format",
                           size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="./2个共有基因表达/ENSG00000293064.boxplot.bytime.log2TPM.addP.pdf",plot=p2,
       device='pdf',path=".",width=3.2,height=3.2)


p <- ggboxplot(exp_melt.final.1.b, x = "Donor", y = "tpm",
               color = "Donor", 
               outlier.colour = NA)+
  theme_classic()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 45, hjust = 1,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(fill="type",x="",y="Expression (log2(TPM+1))",title="") 
#p
# Use only p.format as label. Remove method name.
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p + stat_compare_means(comparisons = my_comparisons,label = "p.format",
                           size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2
