
###################################################################################################
############################################################
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后")
exp.tpm<-read.delim("exp.tpm.addbiotype.txt",header=T,row.names = 1)
phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",
                      header=T,row.names = NULL)
head(phenotype);dim(phenotype)
library(dplyr)
###提取特定供体
phenotype.1<-phenotype %>% subset(Donor %in% c("4","5","6","11","12","13"))
head(phenotype.1);dim(phenotype.1)

phenotype.d1<-phenotype %>% subset(Time %in% c("D1"))
head(phenotype.d1);dim(phenotype.d1)

exp.tpm.1<-exp.tpm[c("ENSG00000225528","ENSG00000293064"),]
exp.tpm.1.pro<-exp.tpm.1[,phenotype.1$sample]
dim(exp.tpm.1.pro)
#[1]  2 24

##########heatmap D1 sample
df1<-df[,phenotype.d1$sample]
final_data<-df1
final_data<-apply(final_data,2,as.numeric)
rownames(final_data)<-rownames(df1)
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

pdf(file = "./D1/heatmap.D1sample.expTPM.scaled.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()

ha <- HeatmapAnnotation(Sex = phenotype.d1$Sex, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1") ))

data.m<-na.omit(data.m)
pdf(file = "./D1/heatmap.D1sample.expTPM.scaled.addlegend.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()

#pdf(file = "./D1/heatmap.D1sample.expTPM.noscaled.pdf", 6.5, 1.6)
#mycol <- colorRamp2(c(0, 2, 4), c("blue", "white", "red"))
#Heatmap(as.matrix(df1),cluster_columns = TRUE,cluster_rows = TRUE,
#        col=mycol,name = "expression", 
#        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
#        row_names_gp=gpar(fontsize = 10),
#        row_title_gp = gpar(col = "black",fontsize = 12),
#        show_row_names = FALSE)
#dev.off()
#########################D1样本低中高表达基因的表达量热图
#####高表达并集
df1<-df.pro[,phenotype.d1$sample]
df1$geneid<-rownames(df1)
#S0102=as.factor(df1[df1[,1]>=100,13])
#S0104=as.factor(df1[df1[,2]>=100,13])
#S0105=as.factor(df1[df1[,3]>=100,13])
#S0106=as.factor(df1[df1[,4]>=100,13])
#S0107=as.factor(df1[df1[,5]>=100,13])
#S0108=as.factor(df1[df1[,6]>=100,13])
#S0111=as.factor(df1[df1[,7]>=100,13])
#S0112=as.factor(df1[df1[,8]>=100,13])
#S0113=as.factor(df1[df1[,9]>=100,13])
#S0115=as.factor(df1[df1[,10]>=100,13])
#S0116=as.factor(df1[df1[,11]>=100,13])
#S0117=as.factor(df1[df1[,12]>=100,13])
##多数据集求并集
genes <- Reduce(union,list(
df1[df1[,1]>=100,13], df1[df1[,2]>=100,13], 
df1[df1[,3]>=100,13], df1[df1[,4]>=100,13],
df1[df1[,5]>=100,13],df1[df1[,6]>=100,13],
df1[df1[,7]>=100,13],df1[df1[,8]>=100,13],
df1[df1[,9]>=100,13],df1[df1[,10]>=100,13],
df1[df1[,11]>=100,13],df1[df1[,12]>=100,13]))
length(genes)
dim(df)
df2<-df.pro[genes,phenotype.d1$sample]
dim(df2)
final_data<-df2
final_data<-apply(final_data,2,as.numeric)
rownames(final_data)<-rownames(df2)
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
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)

pdf(file = "./D1/heatmap.pro.D1sample.high.expTPM.scaled.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()

####
#####高表达交集
df1<-df.pro[,phenotype.d1$sample]
df1$geneid<-rownames(df1)

##多数据集求交集
genes <- Reduce(intersect,list(
  df1[df1[,1]>=100,13], df1[df1[,2]>=100,13], 
  df1[df1[,3]>=100,13], df1[df1[,4]>=100,13],
  df1[df1[,5]>=100,13],df1[df1[,6]>=100,13],
  df1[df1[,7]>=100,13],df1[df1[,8]>=100,13],
  df1[df1[,9]>=100,13],df1[df1[,10]>=100,13],
  df1[df1[,11]>=100,13],df1[df1[,12]>=100,13]))
length(genes)
dim(df)
df2<-df.pro[genes,phenotype.d1$sample]
dim(df2)
final_data<-df2
final_data<-apply(final_data,2,as.numeric)
rownames(final_data)<-rownames(df2)
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
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)

pdf(file = "./D1/heatmap.pro.D1sample.high.expTPM.intersect.scaled.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()
##
##############################################
######box plot
#data1<-df[apply(df,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(data1)

data1<-df.pro[genes,phenotype.d1$sample]
#log(TPM+1)
data1<-log2(data1+1)
normalized_data<-as.data.frame(data1)
normalized_data$geneid<-rownames(normalized_data)
normalized_data

library("reshape")
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","tpm")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype.d1,by = "sample",all=T)
head(exp_melt.final)
#exp_melt.final.1<-exp_melt.final %>% subset(Donor %in% c("4","5","6","11","12","13")) 
dim(exp_melt.final);
dim(exp_melt.final.1)

head(exp_melt.final)

library("ggpubr")
library("reshape")
library("ggplot2")

####by time
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt.final, x = "Donor", y = "tpm",
               fill="Donor",
               
               outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  xlab("") + ylab("Expression (log2(TPM+1))")+
  theme(axis.text.x = 
          element_text(size = 14, face = "bold",angle = 0,
                       hjust = 0,colour = "black",lineheight=100),
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
ggsave(filename="./D1/boxplot.bydonor.420gene.log2TPM.pdf",plot=p,
       device='pdf',path=".",width=6,height=4.5)

################add p
p <- ggboxplot(exp_melt.final, x = "Donor", y = "tpm",
               color ="Donor",
              show.legend = FALSE)+
  theme_classic()+
  xlab("") + ylab("Expression (log2(TPM+1))")+
  theme(axis.text.x = 
          element_text(size = 14, face = "bold",angle = 0,
                       hjust = 0,colour = "black",lineheight=100),
        axis.text.y = element_text(size=14,colour="black",face = "bold"),
        legend.title=element_text(size=12,colour="black",face = "bold"),
        legend.text=element_text(size=12),
        axis.title.x = element_text(size=17,face = "bold"),
        axis.title.y = element_text(size=15,face = "bold"),
        plot.title = element_text(size=16),legend.key.size = unit(4,'mm'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
p
p+stat_compare_means()
p1<-p+stat_compare_means(method = "anova", label.y = 18)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "2")                    # Pairwise comparison against reference
p1
ggsave(filename="./D1/boxplot.bydonor.420gene.log2TPM.addP.d2.pdf",plot=p1,
       device='pdf',path=".",width=6,height=3.6)

#p1<-p+stat_compare_means(method = "anova", label.y = 18)+      # Add global p-value
 # stat_compare_means(label = "p.signif", method = "t.test",
#                     ref.group = "13")                    # Pairwise comparison against reference
#p1
#ggsave(filename="./boxplot.bydonor.2888gene.log2TPM.addP.13.pdf",plot=p1,
#       device='pdf',path=".",width=6,height=3.6)
######共同高表达基因功能分析
###########富集分析
genes
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library("ggplot2")
data(geneList, package="DOSE")

#keygene=df.c[which(df.c$annotation1=="Promoter (<=1kb)"),]

Gene_list<- bitr(genes,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

ego_KEGG <- enrichKEGG(gene         = Gene_list$ENTREZID,
                       organism     = 'hsa',
                       pvalueCutoff = 0.1,
                       qvalueCutoff  = 1)
ego_KEGG <- setReadable(ego_KEGG, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
write.table(ego_KEGG,file="./D1/420gene.enrichKEGG.txt",
            sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
dotplot(ego_KEGG,showCategory=40,title="KEGG Enrichment")
p1<-dotplot(ego_KEGG,showCategory=40,title="KEGG Enrichment")
ggsave(filename = "./D1/420gene.enrichKEGG.pdf", plot =p1,
       width = 16, height = 16, units = 'cm')
#########
#自己的代码画KEGG
kegg<-read.delim('./D1/420gene.enrichKEGG.filter.txt', stringsAsFactors = FALSE)
dim(kegg)
head(kegg)

kegg<-kegg[which(kegg$Count>=5),]
library(ggplot2)
data<-kegg
head(data)
data$LogFDR<--log10(data$pvalue)
data <- data[order(data$Count,data$LogFDR,decreasing = c(TRUE, FALSE),method = "radix"), ]
data$Description <- factor(data$Description, levels = data$Description)

pp = ggplot(data,aes(Count,Description))
p.final<-pp + geom_point(aes(color=(LogFDR),size=Count))+theme_bw()+
  labs(x="Count",y="",title="KEGG enrichment")+
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=11,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=12,colour="black"),
        axis.title.y = element_text(size=12,colour="black"),
        plot.title = element_text(size=13,colour="black"),
        legend.key.size = unit(5,'mm'))+
  scale_colour_gradient(low="DodgerBlue",high="NavyBlue",name="-log P")+
  scale_y_discrete(labels=c("Chemical carcinogenesis - reactive oxygen species"="Chemical carcinogenesis-\nreactive oxygen species"))

p.final
ggsave(filename="./D1/420gene.enrichKEGG.new.pdf",plot=p.final,
       device='pdf',path=".",width=5.7,height=5.4)

####
ego_bp <- enrichGO(gene         = Gene_list$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
ego_bp <- setReadable(ego_bp, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

write.table(ego_bp,file="./D1/420gene.enrichGO.txt",
            sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
dotplot(ego_bp,showCategory=40,title="GO Enrichment")
#p1<-dotplot(ego_bp,showCategory=40,title="BP Enrichment")
#ggsave(filename = "enrichBP.DEG.final.pdf", 
#      plot =p1,width = 19, height = 19, units = 'cm')
#GO的三种条目分开画
go<-read.delim('./D1/420gene.enrichGO.filter.txt', stringsAsFactors = FALSE)
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
ggsave(filename="./D1/420gene.enrichGO.barplot.pdf",plot=p,
       device='pdf',path=".",width=8,height=7)



#####中表达
df1<-df.pro[,phenotype.d1$sample]
df1$geneid<-rownames(df1)

##多数据集求并集
genes <- Reduce(union,list(
df1[df1[,1]>1 & df1[,1]<100,13],
df1[df1[,2]>1 & df1[,2]<100,13],
df1[df1[,3]>1 & df1[,3]<100,13],
df1[df1[,4]>1 & df1[,4]<100,13],
df1[df1[,5]>1 & df1[,5]<100,13],
df1[df1[,6]>1 & df1[,6]<100,13],
df1[df1[,7]>1 & df1[,7]<100,13],
df1[df1[,8]>1 & df1[,8]<100,13],
df1[df1[,9]>1 & df1[,9]<100,13],
df1[df1[,10]>1 & df1[,10]<100,13],
df1[df1[,11]>1 & df1[,11]<100,13],
df1[df1[,12]>1 & df1[,12]<100,13]))

length(genes)
dim(df)
df2<-df.pro[genes,phenotype.d1$sample]
dim(df2)
final_data<-df2
final_data<-apply(final_data,2,as.numeric)
rownames(final_data)<-rownames(df2)
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
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)

pdf(file = "./D1/heatmap.pro.D1sample.medium.expTPM.scaled.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()
#########lnc
df1<-df.lnc[,phenotype.d1$sample]
dim(df1)
df1$geneid<-rownames(df1)
dim(df1)
##多数据集求并集
genes <- Reduce(union,list(
  df1[df1[,1]>=100,13], df1[df1[,2]>=100,13], 
  df1[df1[,3]>=100,13], df1[df1[,4]>=100,13],
  df1[df1[,5]>=100,13],df1[df1[,6]>=100,13],
  df1[df1[,7]>=100,13],df1[df1[,8]>=100,13],
  df1[df1[,9]>=100,13],df1[df1[,10]>=100,13],
  df1[df1[,11]>=100,13],df1[df1[,12]>=100,13]))
length(genes)
dim(df)
df2<-df.lnc[genes,phenotype.d1$sample]
dim(df2)
final_data<-df2
final_data<-apply(final_data,2,as.numeric)
rownames(final_data)<-rownames(df2)
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
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)

pdf(file = "./D1/heatmap.lnc.D1sample.high.expTPM.scaled.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
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
