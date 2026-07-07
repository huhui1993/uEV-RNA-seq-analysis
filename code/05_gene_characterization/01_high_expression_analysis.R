# =============================================================================
# 09. High Expression Gene Analysis
# Description: Identify highly expressed genes by group and donor
# Input: all_sample_count.addAnotation.txt, sample_info.txt
# Output: Mean/sd expression tables (.txt)
# =============================================================================

stat<-exp_melt.final.2 %>% group_by(Donor,Time) %>% 
  dplyr::summarise(median = median(log2_count))
###########################################################################
####以中位数表达作为代表进行表达量比较
list.files()
expr<-read.delim("all_sample_count.addAnotation.txt",header=T,row.names = 1)
head(expr);dim(expr)
expr<-subset(expr,gene_biotype %in% c("protein_coding","lncRNA"))
expr[1:4,1:4]
dim(expr)
expr.pro<-subset(expr,gene_biotype %in% c("protein_coding"))
expr.lnc<-subset(expr,gene_biotype %in% c("lncRNA"))

write.table(expr.pro,
            "all_sample_count.addAnotation.protein_coding.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)
write.table(expr.lnc,
            "all_sample_count.addAnotation.lncRNA.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)


df<-expr[,c(3:32)]
head(df);dim(df)

library(genefilter)
library(Biobase)
data1<-df[apply(df,1,function(x){XXXX<-FALSE;if(any(x>5)){XXXX<-TRUE};return(XXXX)}),];dim(data1)
head(data1)
dim(data1)
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
data_sd<-data_melt.final %>% group_by(geneid,type) %>% 
  dplyr::summarise(sd = sd(count))

data_mean<-data_melt.final %>% group_by(geneid,type) %>% 
  dplyr::summarise(mean = mean(count))
data_mean<-cbind(data_mean,expr[data_mean$geneid,c("gene_name","gene_biotype")])
head(data_mean);dim(data_mean)

write.table(data_sd,
            "data_sd.count.4group.txt", 
            sep="\t", col.names=T,quote = F,row.names = F)

write.table(data_mean,
            "data_mean.count.4group.txt", 
            sep="\t", col.names=T,quote = F,row.names = F)

m1 <- merge(data_mean, data_sd, by = c("geneid","type"))

write.table(m1,
            "data_mean.add_sd.count.4group.txt", 
            sep="\t", col.names=T,quote = F,row.names = F)

library(dplyr)
head(m1)
data_mean.top<-m1 %>% group_by(type,gene_biotype) %>% top_n(n = 200, wt = mean) 
write.table(data_mean.top,
            "data_mean.top200gene.txt", 
            sep="\t", col.names=T,quote = F,row.names = F)

data_mean.top<-m1 %>% group_by(type,gene_biotype) %>% top_n(n = 500, wt = mean) 
write.table(data_mean.top,
            "data_mean.top500gene.txt", 
            sep="\t", col.names=T,quote = F,row.names = F)

data_mean.top<-m1 %>% group_by(type,gene_biotype) %>% top_n(n = 700, wt = mean) 
write.table(data_mean.top,
            "data_mean.top700gene.txt", 
            sep="\t", col.names=T,quote = F,row.names = F)

phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
D1<-phenotype %>% subset(Time=="D1")
D15<-phenotype %>% subset(Time=="D15")
D29<-phenotype %>% subset(Time=="D29")
D43<-phenotype %>% subset(Time=="D43")

threshold_sd <- 2  # 设置标准差的阈值
selected_genes.D1 <- data1[rowSds(data1[,D1$sample]) < threshold_sd, ]
selected_genes.D15 <- data1[rowSds(data1[,D15$sample]) < threshold_sd, ]
selected_genes.D29 <- data1[rowSds(data1[,D29$sample]) < threshold_sd, ]
selected_genes.D43 <- data1[rowSds(data1[,D43$sample]) < threshold_sd, ]

dim(selected_genes.D1);dim(selected_genes.D15);dim(selected_genes.D29);dim(selected_genes.D43)
selected_genes.D1<- selected_genes.D1[,D1$sample]
selected_genes.D15<- selected_genes.D15[,D15$sample]
selected_genes.D29<- selected_genes.D29[,D29$sample]
selected_genes.D43<- selected_genes.D43[,D43$sample]

selected_genes.D1$geneid<-rownames(selected_genes.D1)
selected_genes.D15$geneid<-rownames(selected_genes.D15)
selected_genes.D29$geneid<-rownames(selected_genes.D29)
selected_genes.D43$geneid<-rownames(selected_genes.D43)

library("reshape")
data_D1<-melt(selected_genes.D1,id=c("geneid"))
colnames(data_D1)<-c("geneid","sample","count")

data_D15<-melt(selected_genes.D15,id=c("geneid"))
colnames(data_D15)<-c("geneid","sample","count")
data_D29<-melt(selected_genes.D29,id=c("geneid"))
colnames(data_D29)<-c("geneid","sample","count")
data_D43<-melt(selected_genes.D43,id=c("geneid"))
colnames(data_D43)<-c("geneid","sample","count")

data_D1<-data_D1 %>% group_by(geneid) %>% 
  dplyr::summarise(mean = mean(count))
data_D15<-data_D15 %>% group_by(geneid) %>% 
  dplyr::summarise(mean = mean(count))
data_D29<-data_D29 %>% group_by(geneid) %>% 
  dplyr::summarise(mean = mean(count))
data_D43<-data_D43 %>% group_by(geneid) %>% 
  dplyr::summarise(mean = mean(count))

data_D1.1<-cbind(data_D1,expr[data_D1$geneid,c("gene_name","gene_biotype")])
data_D15.1<-cbind(data_D15,expr[data_D15$geneid,c("gene_name","gene_biotype")])
data_D29.1<-cbind(data_D29,expr[data_D29$geneid,c("gene_name","gene_biotype")])
data_D43.1<-cbind(data_D43,expr[data_D43$geneid,c("gene_name","gene_biotype")])

library(dplyr)
D1.top<-data_D1.1 %>% group_by(gene_biotype) %>% top_n(n = 200, wt = mean) 
write.table(D1.top,
            "D1.top200gene.txt", 
            sep="\t", col.names=T,quote = F,row.names = F)

phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)

library("reshape")
data_melt<-melt(data1,id=c("geneid"))
head(data_melt)
colnames(data_melt)<-c("geneid","sample","count")
head(data_melt);dim(data_melt)
data_melt.final<-merge(data_melt,phenotype,by = "sample",all=T)
head(data_melt.final);dim(data_melt.final)

data_median<-data_melt.final %>% group_by(Donor,Time,Sex) %>% 
  dplyr::summarise(median = median(count))

data_median
####同一个人不同时间点比较,raw count median  line_plot
head(data_median);dim(data_median)
library("reshape")

unique(data_median$Donor)
data_median.2 <- data_median %>% subset(Donor %in% c("4","5","6","11","12","13")) 
#exp_melt.final.2.1<-exp_melt.final.2 %>% subset(scaled_count>=-0.3 & scaled_count<=0.2)
dim(data_median.2);dim(data_median)
head(data_median.2)
#exp_melt.final.1.a <- exp_melt.final.1 %>% subset(Sex=="female")
#exp_melt.final.1.b <- exp_melt.final.1 %>% subset(Sex=="male")
data_median.2$Donor<-as.factor(data_median.2$Donor)

