# =============================================================================
# 02_mfuzz_clustering_tpm
# Description: Time Series Clustering with Mfuzz (TPM-based)
# =============================================================================

library(Mfuzz)
library(Biobase)
library(limma)
library(clusterProfiler)

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

#df<-expr[,c(3:32)]
df<-exp.tpm[,c(2:31)]
head(df)

df.1<-df
df.1<-apply(df.1,2,as.numeric)
rownames(df.1)<-rownames(df)
df.1[1:4,1:4]
str(df.1)
df.1<-as.data.frame(df.1)
str(df.1)

head(df.1);dim(df.1)

phenotype<-read.delim("sample.type.txt",header=T,row.names = 1)
head(phenotype);dim(phenotype)
#phenotype<-phenotype[rownames(phenotype)!="N6MVA",]
colnames(df.1)<-phenotype
identical(colnames(df.1),rownames(phenotype))
head(df.1);dim(df.1)

colnames(df.1)<-phenotype$Time
head(df.1);dim(df.1)

dat<-df.1

library(limma)
avereps_df  <- t(limma::avereps( t(dat) , ID = colnames(dat)))##对相同时间序列的表达值取平均
avereps_df[1:4,1:4]
colnames(avereps_df)
save(avereps_df,file = 'avereExp_TPM.Rdata')


#gene <- as.matrix(df.1)

## 2.1 Filtering----
# 去除表达量太低或者在不同时间点间变化太小的基因等步骤
# Mfuzz聚类时要求是一个ExpressionSet类型的对象，所以需要先用表达量构建这样一个对象。
eset <- new("ExpressionSet",exprs = avereps_df)

eset <- filter.NA(eset, thres = 0.25)

eset <- fill.NA(eset, mode = 'mean')

# 根据标准差去除样本间差异太小的基因
eset <- filter.std(eset,min.std=0)
# 804 genes excluded ，不同的数据集去除的基因数量不一样
eset


## 2.2 Standardisation----
# 聚类时需要用一个数值来表征不同基因间的距离，Mfuzz中采用的是欧式距离，
# 由于普通欧式距离的定义没有考虑不同维度间量纲的不同，所以需要先进行标准化
eset <- standardise(eset)

## 2.3 Setting of parameters for FCM clustering----
# Mfuzz中的聚类算法需要提供两个参数，
# 第一个参数为希望最终得到的聚类的个数，这个参数由我们直接指定
# 第二个参数称之为fuzzifier值，用小写字母m表示，可以通过函数评估一个最佳取值
c <- 5
m <- mestimate(eset) #  评估出最佳的m值
m
set.seed(2023)
cl <- mfuzz(eset, c = c, m = m) # 聚类

## 2.4 glimpse results----
# 在cl这个对象中就保存了聚类的完整结果，对于这个对象的常见操作如下
cl$size # 查看每个cluster中的基因个数
cl$cluster[cl$cluster == 1] # 提取某个cluster下的基因
## cluster cores
# membership values can also indicate the similarity of vectors to each other.
eset
cl.thres <- acore(eset,cl,min.acore=0.5)  ## extracts genes forming the alpha cores of soft clusters
head(cl.thres[[1]])
table(cl$cluster)
unlist(lapply(cl.thres, nrow))
#查看集群之间的耦合
# coupling between clusters
overlap.cl <- overlap(cl)
pdf('mfuzz_overlap_plot.tpm.pdf',height = 4,width = 5)
p.overlaps <- overlap.plot(cl, over = overlap.cl, thres = 0.05)
p.overlaps
dev.off()

## 2.5 visualise----
library(RColorBrewer)
color.2 <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
pdf('mfuzz_clusters_plot.tpm.pdf',height = 7,width = 12)
mfuzz.plot(eset,cl,mfrow=c(3,3),
           new.window= FALSE,
           time.labels= colnames(eset) ,
           colo = color.2)
dev.off()
pdf('mfuzz_clusters_plot01.tpm.pdf',height = 7,width = 12)
mfuzz.plot2(eset, cl, mfrow = c(3, 3),
            time.labels= colnames(eset) ,
            centre = T, x11 = F, centre.lwd = 0.2)
dev.off()

# 3. 批量输出聚类所含的基因 ----------------------------------------------------------
getwd()
dir.create(path = "mfuzzGenes.tpm", recursive = T)
for (i in 1:5) {
  potname = names(cl$cluster[unname(cl$cluster) == i])
  write.csv(cl[[4]][potname, i], paste0("mfuzzGenes.tpm", "/mfuzz_", i, ".csv"))
}
########################################################
#####male 男性
library(Mfuzz)
library(Biobase)
library(limma)
library(clusterProfiler)

list.files()
expr<-read.delim("all_sample_count.addAnotation.txt",header=T,row.names = 1)
head(expr);dim(expr)
expr<-subset(expr,gene_biotype %in% c("protein_coding","lncRNA"))
expr[1:4,1:4]
dim(expr)

df<-expr[,c(3:32)]

head(df)
head(df);dim(df)

phenotype<-read.delim("sample.type.txt",header=T,row.names = 1)
head(phenotype);dim(phenotype)
male.s<-rownames(phenotype[which(phenotype$Sex=="male"),])
male.s
female.s<-rownames(phenotype[which(phenotype$Sex=="female"),])
female.s
df.male<-df.1[,male.s]
df.female<-df.1[,female.s]
#phenotype<-phenotype[rownames(phenotype)!="N6MVA",]
colnames(df.male)<-phenotype[male.s,]$Time
identical(colnames(df.male),phenotype[male.s,]$Time)
head(df.male);dim(df.male)
colnames(df.male)
head(df.male)

dat<-df.male

library(limma)
avereps_df  <- t(limma::avereps( t(dat) , ID = colnames(dat)))##对相同时间序列的表达值取平均
avereps_df[1:4,1:4]
colnames(avereps_df)
save(avereps_df,file = 'avereExp_TPM.male.Rdata')

## 2.1 Filtering----
# 去除表达量太低或者在不同时间点间变化太小的基因等步骤
# Mfuzz聚类时要求是一个ExpressionSet类型的对象，所以需要先用表达量构建这样一个对象。
eset <- new("ExpressionSet",exprs = avereps_df)

eset <- filter.NA(eset, thres = 0.25)

eset <- fill.NA(eset, mode = 'mean')

# 根据标准差去除样本间差异太小的基因
eset <- filter.std(eset,min.std=0)
# 2768 genes excluded ，不同的数据集去除的基因数量不一样
eset


## 2.2 Standardisation----
# 聚类时需要用一个数值来表征不同基因间的距离，Mfuzz中采用的是欧式距离，
# 由于普通欧式距离的定义没有考虑不同维度间量纲的不同，所以需要先进行标准化
eset <- standardise(eset)

## 2.3 Setting of parameters for FCM clustering----
# Mfuzz中的聚类算法需要提供两个参数，
# 第一个参数为希望最终得到的聚类的个数，这个参数由我们直接指定
# 第二个参数称之为fuzzifier值，用小写字母m表示，可以通过函数评估一个最佳取值
c <- 5
m <- mestimate(eset) #  评估出最佳的m值
m
set.seed(2023)
cl <- mfuzz(eset, c = c, m = m) # 聚类

## 2.4 glimpse results----
# 在cl这个对象中就保存了聚类的完整结果，对于这个对象的常见操作如下
cl$size # 查看每个cluster中的基因个数
cl$cluster[cl$cluster == 1] # 提取某个cluster下的基因
## cluster cores
# membership values can also indicate the similarity of vectors to each other.
eset
cl.thres <- acore(eset,cl,min.acore=0.5)  ## extracts genes forming the alpha cores of soft clusters
head(cl.thres[[1]])
table(cl$cluster)
unlist(lapply(cl.thres, nrow))
#查看集群之间的耦合
# coupling between clusters
overlap.cl <- overlap(cl)
pdf('mfuzz_overlap_plot.male.tpm.pdf',height = 4,width = 5)
p.overlaps <- overlap.plot(cl, over = overlap.cl, thres = 0.05)
p.overlaps
dev.off()

## 2.5 visualise----
library(RColorBrewer)
color.2 <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
pdf('mfuzz_clusters_plot.male.tpm.pdf',height = 7,width = 12)
mfuzz.plot(eset,cl,mfrow=c(3,3),
           new.window= FALSE,
           time.labels= colnames(eset) ,
           colo = color.2)
dev.off()
pdf('mfuzz_clusters_plot01.male.tpm.pdf',height = 7,width = 12)
mfuzz.plot2(eset, cl, mfrow = c(3, 3),
            time.labels= colnames(eset) ,
            centre = T, x11 = F, centre.lwd = 0.2)
dev.off()

# 3. 批量输出聚类所含的基因 
dir.create(path = "mfuzzGenes.male.tpm", recursive = T)
for (i in 1:5) {
  potname = names(cl$cluster[unname(cl$cluster) == i])
  write.csv(cl[[4]][potname, i], paste0("mfuzzGenes.male.tpm", "/mfuzz_male_", i, ".csv"))
}

########################################################
#####female 女性
library(Mfuzz)
library(Biobase)
library(limma)
library(clusterProfiler)

list.files()
expr<-read.delim("all_sample_count.addAnotation.txt",header=T,row.names = 1)
head(expr);dim(expr)
expr<-subset(expr,gene_biotype %in% c("protein_coding","lncRNA"))
expr[1:4,1:4]
dim(expr)

df<-expr[,c(3:32)]

head(df)
head(df);dim(df)

phenotype<-read.delim("sample.type.txt",header=T,row.names = 1)
head(phenotype);dim(phenotype)
male.s<-rownames(phenotype[which(phenotype$Sex=="male"),])
male.s
female.s<-rownames(phenotype[which(phenotype$Sex=="female"),])
female.s
df.male<-df.1[,male.s]
df.female<-df.1[,female.s]
#phenotype<-phenotype[rownames(phenotype)!="N6MVA",]
colnames(df.female)<-phenotype[female.s,]$Time
identical(colnames(df.female),phenotype[female.s,]$Time)
head(df.female);dim(df.female)
colnames(df.female)

dat<-df.female

library(limma)
avereps_df  <- t(limma::avereps( t(dat) , ID = colnames(dat)))##对相同时间序列的表达值取平均
avereps_df[1:4,1:4]
colnames(avereps_df)
save(avereps_df,file = 'avereExp_TPM.female.Rdata')
## 2.1 Filtering----
# 去除表达量太低或者在不同时间点间变化太小的基因等步骤
# Mfuzz聚类时要求是一个ExpressionSet类型的对象，所以需要先用表达量构建这样一个对象。
eset <- new("ExpressionSet",exprs = avereps_df)

eset <- filter.NA(eset, thres = 0.25)

eset <- fill.NA(eset, mode = 'mean')

# 根据标准差去除样本间差异太小的基因
eset <- filter.std(eset,min.std=0)
# 2222 genes excluded ，不同的数据集去除的基因数量不一样
eset


## 2.2 Standardisation----
# 聚类时需要用一个数值来表征不同基因间的距离，Mfuzz中采用的是欧式距离，
# 由于普通欧式距离的定义没有考虑不同维度间量纲的不同，所以需要先进行标准化
eset <- standardise(eset)

## 2.3 Setting of parameters for FCM clustering----
# Mfuzz中的聚类算法需要提供两个参数，
# 第一个参数为希望最终得到的聚类的个数，这个参数由我们直接指定
# 第二个参数称之为fuzzifier值，用小写字母m表示，可以通过函数评估一个最佳取值
c <- 5
m <- mestimate(eset) #  评估出最佳的m值
set.seed(2023)
cl <- mfuzz(eset, c = c, m = m) # 聚类

## 2.4 glimpse results----
# 在cl这个对象中就保存了聚类的完整结果，对于这个对象的常见操作如下
cl$size # 查看每个cluster中的基因个数
cl$cluster[cl$cluster == 1] # 提取某个cluster下的基因
## cluster cores
# membership values can also indicate the similarity of vectors to each other.
eset
cl.thres <- acore(eset,cl,min.acore=0.5)  ## extracts genes forming the alpha cores of soft clusters
head(cl.thres[[1]])
table(cl$cluster)
unlist(lapply(cl.thres, nrow))
#查看集群之间的耦合
# coupling between clusters
overlap.cl <- overlap(cl)
pdf('mfuzz_overlap_plot.female.tpm.pdf',height = 4,width = 5)
p.overlaps <- overlap.plot(cl, over = overlap.cl, thres = 0.05)
p.overlaps
dev.off()

## 2.5 visualise----
library(RColorBrewer)
color.2 <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
pdf('mfuzz_clusters_plot.female.tpm.pdf',height = 7,width = 12)
mfuzz.plot(eset,cl,mfrow=c(3,3),
           new.window= FALSE,
           time.labels= colnames(eset) ,
           colo = color.2)
dev.off()
pdf('mfuzz_clusters_plot01.female.tpm.pdf',height = 7,width = 12)
mfuzz.plot2(eset, cl, mfrow = c(3, 3),
            time.labels= colnames(eset) ,
            centre = T, x11 = F, centre.lwd = 0.2)
dev.off()

# 3. 批量输出聚类所含的基因 
dir.create(path = "mfuzzGenes.female.tpm", recursive = T)
for (i in 1:5) {
  potname = names(cl$cluster[unname(cl$cluster) == i])
  write.csv(cl[[4]][potname, i], paste0("mfuzzGenes.female.tpm", "/mfuzz_female_", i, ".csv"))
}
######
####################数据基本分析
library(Mfuzz)
library(Biobase)
library(limma)
library(clusterProfiler)

list.files()
expr<-read.delim("all_sample_count.addAnotation.txt",header=T,row.names = 1)
head(expr);dim(expr)
expr<-subset(expr,gene_biotype %in% c("protein_coding","lncRNA"))
expr[1:4,1:4]
dim(expr)

df<-expr[,c(3:32)]

head(df);dim(df)

library(genefilter)
library(Biobase)
data1<-df[apply(df,1,function(x){XXXX<-FALSE;if(any(x>5)){XXXX<-TRUE};return(XXXX)}),];dim(data1)

# 在基因表达数据中，归一化是为了使不同基因的表达值在数量级上保持一致，便于比较。
normalized_data <- scale(data1)
head(normalized_data)
# Boxplot
pdf(file = "boxplot.all_sample_scale_count_exp.pdf",15,5)
boxplot(normalized_data, col = "lightblue", main = "Normalized Gene Count",
          xlab = "Gene", ylab = "Normalized Expression",ylim=c(-0.5,0.5))
dev.off()

# Scatter plot
pdf("pairs.plot",20,20)
pairs(normalized_data, pch = 19, cex = 0.8)#太大画不了
dev.off()

normalized_data<-as.data.frame(normalized_data)
normalized_data$geneid<-rownames(normalized_data)
library("reshape")
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","scaled_count")
head(exp_melt)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt, x = "sample", y = "scaled_count",
              palette =c("lightblue"),
               add = "jitter")+
  xlab("") + ylab("Gene expression")+
  theme(axis.text.x = 
          element_text(size = 14, face = "bold",angle = 45, hjust = 1,
                       vjust =0,colour = "black",lineheight=100),
        axis.text.y = element_text(size=14,colour="black",face = "bold"),
        legend.title=element_text(size=12,colour="black",face = "bold"),
        legend.text=element_text(size=12),
        axis.title.x = element_text(size=17,face = "bold"),
        axis.title.y = element_text(size=15,face = "bold"),
        plot.title = element_text(size=16),legend.key.size = unit(4,'mm'))+
  theme(legend.position='none')+
  theme(axis.line = element_line(size=1, colour = "black"))+ylim(-0.5,0.5) 
p


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
