# =============================================================================
# 04_remove_within_group_variation
# Description: Remove Within-group Variation for Time Series
# =============================================================================

stat<-exp_melt.final.2 %>% group_by(Donor,Time) %>% 
  dplyr::summarise(median = median(log2_count))
###########################################################################
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
dim(exp.tpm)

exp.tpm.pro<-subset(exp.tpm,gene_biotype %in% c("protein_coding"))
exp.tpm.lnc<-subset(exp.tpm,gene_biotype %in% c("lncRNA"))

write.table(exp.tpm.pro,
            "all_sample_tpm.protein_coding.txt", 
            sep="\t", col.names=T,quote = F)

write.table(exp.tpm.lnc,
            "all_sample_tpm.lnc.txt", 
            sep="\t", col.names=T,quote = F)


library(genefilter)
library(Biobase)
#data1<-df[apply(df,1,function(x){XXXX<-FALSE;if(any(x>5)){XXXX<-TRUE};return(XXXX)}),];dim(data1)

exp.tpm.test<-exp.tpm[,c(1:30)]
head(exp.tpm.test)
dim(exp.tpm.test)


final_data<-exp.tpm.test
final_data<-apply(final_data,2,as.numeric)
rownames(final_data)<-rownames(exp.tpm.test)
final_data[1:4,1:4]
str(final_data)
final_data<-as.data.frame(final_data)
str(final_data)
write.table(final_data,
            "all_sample_tpm.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)

data1<-final_data
data1$geneid<-rownames(data1)

library("reshape")
data_melt<-melt(data1,id=c("geneid"))
head(data_melt)
colnames(data_melt)<-c("geneid","sample","count")
head(data_melt);dim(data_melt)

phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)

data_melt.final<-merge(data_melt,phenotype,by = "sample",all=T)
head(data_melt.final);dim(data_melt.final)

library(dplyr)
data_sd<-data_melt.final %>% group_by(geneid,Time) %>% 
  dplyr::summarise(sd = sd(count))

data_mean<-data_melt.final %>% group_by(geneid,Time) %>% 
  dplyr::summarise(mean = mean(count))
data_mean<-cbind(data_mean,expr[data_mean$geneid,c("gene_name","gene_biotype")])
head(data_mean);dim(data_mean)

write.table(data_sd,
            "data_sd.TPM.4group.txt", 
            sep="\t", col.names=T,quote = F,row.names = F)

write.table(data_mean,
            "data_mean.TPM.4group.txt", 
            sep="\t", col.names=T,quote = F,row.names = F)

m1 <- merge(data_mean, data_sd, by = c("geneid","Time"))

write.table(m1,
            "data_mean.add_sd.TPM.4group.txt", 
            sep="\t", col.names=T,quote = F,row.names = F)

head(m1)
D1<-subset(m1,Time=="D1")
D15<-subset(m1,Time=="D15")
D29<-subset(m1,Time=="D29")
D43<-subset(m1,Time=="D43")

summary(D1$sd)
summary(D15$sd)
summary(D29$sd)
summary(D43$sd)

#threshold_sd <- 2  # 设置标准差的阈值
selected_genes.D1 <- D1[D1$sd < summary(D1$sd)[4], ]
selected_genes.D1 <- selected_genes.D1[selected_genes.D1$mean > 1, ]

selected_genes.D15 <-  D15[D15$sd < summary(D15$sd)[4], ]
selected_genes.D15 <- selected_genes.D15[selected_genes.D15$mean > 1, ]

selected_genes.D29 <-  D29[D29$sd < summary(D29$sd)[4], ]
selected_genes.D29 <- selected_genes.D29[selected_genes.D29$mean > 1, ]

selected_genes.D43 <-  D43[D43$sd < summary(D43$sd)[4], ]
selected_genes.D43 <- selected_genes.D43[selected_genes.D43$mean > 1, ]

dim(selected_genes.D1);dim(selected_genes.D15);dim(selected_genes.D29);dim(selected_genes.D43)

dim(selected_genes.D1);dim(selected_genes.D15);dim(selected_genes.D29);dim(selected_genes.D43)

in1<-intersect(unique(selected_genes.D1$geneid),unique(selected_genes.D15$geneid))
in2<-intersect(unique(selected_genes.D29$geneid),unique(selected_genes.D43$geneid))
in3<-intersect(in1,in2)
length(in3)

head(final_data)
sdfiltergene5534<-final_data[in3,]
head(sdfiltergene5534);dim(sdfiltergene5534)
write.table(sdfiltergene5534,
            "sdfiltergene5534_tpm.for_time_series.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)

head(expr);dim(expr)

sdfiltergene5534.annotaion<-cbind(expr[rownames(sdfiltergene5534),c("gene_name","gene_biotype")],sdfiltergene5534)
head(sdfiltergene5534.annotaion)
dim(sdfiltergene5534.annotaion)
write.table(sdfiltergene5534.annotaion,
            "sdfiltergene5534.addAnnotation_tpm.for_time_series.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)
#####################################################################
#####时序分析

library(Mfuzz)
library(Biobase)
library(limma)
library(clusterProfiler)

head(sdfiltergene5534)
phenotype<-read.delim("sample.type.txt",header=T,row.names = 1)
head(phenotype);dim(phenotype)
#phenotype<-phenotype[rownames(phenotype)!="N6MVA",]
#colnames(df.1)<-phenotype
identical(colnames(sdfiltergene5534),rownames(phenotype))
head(sdfiltergene5534);dim(sdfiltergene5534)

colnames(sdfiltergene5534)<-phenotype$Time
head(sdfiltergene5534);dim(sdfiltergene5534)

dat<-sdfiltergene5534

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

#####male 男性
library(Mfuzz)
library(Biobase)
library(limma)
library(clusterProfiler)

head(sdfiltergene5534)

head(phenotype);dim(phenotype)
male.s<-rownames(phenotype[which(phenotype$Sex=="male"),])
male.s
female.s<-rownames(phenotype[which(phenotype$Sex=="female"),])
female.s
df.male<-sdfiltergene5534[,male.s]
df.female<-sdfiltergene5534[,female.s]
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

#####################
#####female 女性
library(Mfuzz)
library(Biobase)
library(limma)
library(clusterProfiler)

head(phenotype);dim(phenotype)
male.s<-rownames(phenotype[which(phenotype$Sex=="male"),])
male.s
female.s<-rownames(phenotype[which(phenotype$Sex=="female"),])
female.s
df.male<-sdfiltergene5534[,male.s]
df.female<-sdfiltergene5534[,female.s] 

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
######热图
list.files()
data<-read.csv("mfuzz_male_4.csv",header = TRUE,row.names = 1)
colnames(data)<-c("core")
library(dplyr)
data.1<-data %>% subset(core >=0.5)
head(data.1);dim(data.1)

male.s<-rownames(phenotype[which(phenotype$Sex=="male"),])
male.s
female.s<-rownames(phenotype[which(phenotype$Sex=="female"),])
female.s

exp.tpm.1<-sdfiltergene5534.annotaion[rownames(data.1),]
exp.tpm.1.pro<- exp.tpm.1 %>% subset(gene_biotype=="protein_coding")
exp.tpm.1.lnc<- exp.tpm.1 %>% subset(gene_biotype=="lncRNA")
head(exp.tpm.1.pro)
dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)

write.table(exp.tpm.1.pro,
            "male_4_core0.5_protein_coding.tpm.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)
write.table(exp.tpm.1.lnc,
            "male_4_core0.5_lnc.tpm.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)

exp.tpm.1.pro<-exp.tpm.1.pro[,c(male.s)]
exp.tpm.1.lnc<-exp.tpm.1.lnc[,c(male.s)]
dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)

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

phenotype.male<-phenotype %>% subset(Sex=="male")
head(phenotype)
#info$stim<-factor(info$stim,levels=c("NC","D1.5","D7"))
#rownames(phenotype)<-phenotype$sample

data.m<-na.omit(data.m)


ha <- HeatmapAnnotation(Sex = phenotype.male$Sex, Time = phenotype.male$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "heatmap.onlymaleSample.male_4.expTPM.scaled.addlegend.pdf", 4.5, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.male$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()

##########heatmap lnc
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

rownames(phenotype)<-phenotype$sample

data.m<-na.omit(data.m)
phenotype.male<-phenotype %>% subset(Sex=="male")
ha <- HeatmapAnnotation(Sex = phenotype.male$Sex, Time = phenotype.male$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "heatmap.lnc.onlymaleSample.male_4.expTPM.scaled.addlegend.pdf", 4.5, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.male$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()
#########male_4 in female
exp.tpm.1<-sdfiltergene5534.annotaion[rownames(data.1),]
exp.tpm.1.pro<- exp.tpm.1 %>% subset(gene_biotype=="protein_coding")
exp.tpm.1.lnc<- exp.tpm.1 %>% subset(gene_biotype=="lncRNA")
head(exp.tpm.1.pro)
dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)

exp.tpm.1.pro<-exp.tpm.1.pro[,c(female.s)]
exp.tpm.1.lnc<-exp.tpm.1.lnc[,c(female.s)]
dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)

#phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",header=T,row.names = NULL)
#head(phenotype);dim(phenotype)

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

phenotype.female<-phenotype %>% subset(Sex=="female")
head(phenotype.female)
#info$stim<-factor(info$stim,levels=c("NC","D1.5","D7"))
#rownames(phenotype)<-phenotype$sample

data.m<-na.omit(data.m)


ha <- HeatmapAnnotation(Sex = phenotype.female$Sex, Time = phenotype.female$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "heatmap.onlyfemaleSample.male_4.expTPM.scaled.addlegend.pdf", 4.5, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.female$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()

##########heatmap lnc
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

#rownames(phenotype)<-phenotype$sample

data.m<-na.omit(data.m)
phenotype.female<-phenotype %>% subset(Sex=="female")
ha <- HeatmapAnnotation(Sex = phenotype.female$Sex, Time = phenotype.female$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "heatmap.lnc.onlyfemaleSample.male_4.expTPM.scaled.addlegend.pdf", 4.5, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.male$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()


#####male_1\2\3\5
list.files(pattern = "*.csv")

data<-read.csv("mfuzz_male_5.csv",header = TRUE,row.names = 1)
colnames(data)<-c("core")
library(dplyr)
data.1<-data %>% subset(core >=0.5)
head(data.1);dim(data.1)

exp.tpm.1<-sdfiltergene5534.annotaion[rownames(data.1),]
exp.tpm.1.pro<- exp.tpm.1 %>% subset(gene_biotype=="protein_coding")
exp.tpm.1.lnc<- exp.tpm.1 %>% subset(gene_biotype=="lncRNA")
head(exp.tpm.1.pro)
dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)
exp.tpm.1.pro<-exp.tpm.1.pro[,c(3:32)]
exp.tpm.1.lnc<-exp.tpm.1.lnc[,c(3:32)]


#phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",header=T,row.names = NULL)
#head(phenotype);dim(phenotype)

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

#################
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

pdf(file = "heatmap.male_5.expTPM.scaled.addlegend.pdf", 6.2, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()

##########heatmap lnc
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

pdf(file = "heatmap.lnc.male_5.expTPM.scaled.addlegend.pdf", 6.2, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()
###########################################################################################
###########################################################################################
###female
######热图
list.files()
data<-read.csv("mfuzz_female_2.csv",header = TRUE,row.names = 1)
colnames(data)<-c("core")
library(dplyr)
data.1<-data %>% subset(core >=0.5)
head(data.1);dim(data.1)

male.s<-rownames(phenotype[which(phenotype$Sex=="male"),])
male.s
female.s<-rownames(phenotype[which(phenotype$Sex=="female"),])
female.s

exp.tpm.1<-sdfiltergene5534.annotaion[rownames(data.1),]
exp.tpm.1.pro<- exp.tpm.1 %>% subset(gene_biotype=="protein_coding")
exp.tpm.1.lnc<- exp.tpm.1 %>% subset(gene_biotype=="lncRNA")
head(exp.tpm.1.pro)
dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)

dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)

write.table(exp.tpm.1.pro,
            "female_2_core0.5_protein_coding.tpm.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)
write.table(exp.tpm.1.lnc,
            "female_2_core0.5_lnc.tpm.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)

###female_2 in male
exp.tpm.1.pro<-exp.tpm.1.pro[,c(male.s)]
exp.tpm.1.lnc<-exp.tpm.1.lnc[,c(male.s)]
dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)

#phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",header=T,row.names = NULL)
#head(phenotype);dim(phenotype)

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

phenotype.male<-phenotype %>% subset(Sex=="male")
head(phenotype)
#info$stim<-factor(info$stim,levels=c("NC","D1.5","D7"))
#rownames(phenotype)<-phenotype$sample

data.m<-na.omit(data.m)


ha <- HeatmapAnnotation(Sex = phenotype.male$Sex, Time = phenotype.male$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "heatmap.onlymaleSample.female_2.expTPM.scaled.addlegend.pdf", 4.5, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.male$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()

##########heatmap lnc
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

rownames(phenotype)<-phenotype$sample

data.m<-na.omit(data.m)
phenotype.male<-phenotype %>% subset(Sex=="male")
ha <- HeatmapAnnotation(Sex = phenotype.male$Sex, Time = phenotype.male$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "heatmap.lnc.onlymaleSample.female_2.expTPM.scaled.addlegend.pdf", 4.5, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.male$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()
#########female_2 in female
exp.tpm.1<-sdfiltergene5534.annotaion[rownames(data.1),]
exp.tpm.1.pro<- exp.tpm.1 %>% subset(gene_biotype=="protein_coding")
exp.tpm.1.lnc<- exp.tpm.1 %>% subset(gene_biotype=="lncRNA")
head(exp.tpm.1.pro)
dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)

exp.tpm.1.pro<-exp.tpm.1.pro[,c(female.s)]
exp.tpm.1.lnc<-exp.tpm.1.lnc[,c(female.s)]
dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)

#phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",header=T,row.names = NULL)
#head(phenotype);dim(phenotype)

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

phenotype.female<-phenotype %>% subset(Sex=="female")
head(phenotype.female)
#info$stim<-factor(info$stim,levels=c("NC","D1.5","D7"))
#rownames(phenotype)<-phenotype$sample

data.m<-na.omit(data.m)


ha <- HeatmapAnnotation(Sex = phenotype.female$Sex, Time = phenotype.female$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "heatmap.onlyfemaleSample.female_2.expTPM.scaled.addlegend.pdf", 4.5, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.female$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()

##########heatmap lnc
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

#rownames(phenotype)<-phenotype$sample

data.m<-na.omit(data.m)
phenotype.female<-phenotype %>% subset(Sex=="female")
ha <- HeatmapAnnotation(Sex = phenotype.female$Sex, Time = phenotype.female$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "heatmap.lnc.onlyfemaleSample.female_2.expTPM.scaled.addlegend.pdf", 4.5, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.male$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()


#####male_1\2\3\5
list.files(pattern = "*.csv")

data<-read.csv("mfuzz_male_5.csv",header = TRUE,row.names = 1)
colnames(data)<-c("core")
library(dplyr)
data.1<-data %>% subset(core >=0.5)
head(data.1);dim(data.1)

exp.tpm.1<-sdfiltergene5534.annotaion[rownames(data.1),]
exp.tpm.1.pro<- exp.tpm.1 %>% subset(gene_biotype=="protein_coding")
exp.tpm.1.lnc<- exp.tpm.1 %>% subset(gene_biotype=="lncRNA")
head(exp.tpm.1.pro)
dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)
exp.tpm.1.pro<-exp.tpm.1.pro[,c(3:32)]
exp.tpm.1.lnc<-exp.tpm.1.lnc[,c(3:32)]

#phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",header=T,row.names = NULL)
#head(phenotype);dim(phenotype)

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

#################
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

pdf(file = "heatmap.male_5.expTPM.scaled.addlegend.pdf", 6.2, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()

##########heatmap lnc
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

pdf(file = "heatmap.lnc.male_5.expTPM.scaled.addlegend.pdf", 6.2, 6)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype$Time,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()

#################################################################
######表达模式相同的基因在男女中是否有重合？
getwd()

male.1<-read.csv("./mfuzzGenes.male.tpm/mfuzz_male_5.csv",header = TRUE,row.names = 1)
female.1<-read.csv("./mfuzzGenes.female.tpm/mfuzz_female_4.csv",header = TRUE,row.names = 1)

library(dplyr)
m.1<-male.1 %>% subset(x >=0.5)
f.1<-female.1 %>% subset(x >=0.5)
dim(m.1);dim(f.1)
length(intersect(rownames(m.1),rownames(f.1)))



head(data.1);dim(data.1)

exp.tpm.1<-sdfiltergene5534.annotaion[rownames(data.1),]
exp.tpm.1.pro<- exp.tpm.1 %>% subset(gene_biotype=="protein_coding")
exp.tpm.1.lnc<- exp.tpm.1 %>% subset(gene_biotype=="lncRNA")
head(exp.tpm.1.pro)
dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)

dim(exp.tpm.1.pro);dim(exp.tpm.1.lnc)

write.table(exp.tpm.1.pro,
            "female_2_core0.5_protein_coding.tpm.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)
write.table(exp.tpm.1.lnc,
            "female_2_core0.5_lnc.tpm.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)