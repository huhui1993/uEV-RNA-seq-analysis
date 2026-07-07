
######
####################数据基本分析
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后")
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

dim(expr)
exp.tpm<-read.delim("all_sample_tpm.txt",header=T,row.names = 1)
dim(exp.tpm)
exp.tpm[1:4,1:4]
exp.tpm<-exp.tpm[rownames(expr),]
dim(exp.tpm)
identical(rownames(expr),rownames(exp.tpm))
exp.tpm$gene_biotype<-expr$gene_biotype
head(exp.tpm)

exp.tpm.pro<-subset(exp.tpm,gene_biotype %in% c("protein_coding"))
exp.tpm.lnc<-subset(exp.tpm,gene_biotype %in% c("lncRNA"))
dim(exp.tpm)
head(exp.tpm)
df<-exp.tpm[,c(2:31)]

head(df);dim(df)

library(genefilter)
library(Biobase)
data1<-df[apply(df,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(data1)
#log(TPM+1)
data1<-log2(data1+1)
# 在基因表达数据中，归一化是为了使不同基因的表达值在数量级上保持一致，便于比较。
#normalized_data <- scale(data1)
#head(normalized_data)
# Boxplot
#pdf(file = "boxplot.all_sample_scale_count_exp.pdf",15,5)
#boxplot(normalized_data, col = "lightblue", main = "Normalized Gene Count",
#          xlab = "Gene", ylab = "Normalized Expression",ylim=c(-0.5,0.5))
#dev.off()

# Scatter plot
#pdf("pairs.plot",20,20)
#pairs(normalized_data, pch = 19, cex = 0.8)#太大画不了
#dev.off()

normalized_data<-as.data.frame(data1)
normalized_data$geneid<-rownames(normalized_data)

phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)

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

###################################

exp_melt.final %>% group_by(Time,Sex) %>% 
  dplyr::summarise(median = median(log2_count))



####同一个人不同时间点
exp_melt.final.2 <- exp_melt.final %>% 
  subset(Donor %in% c("4","5","6","11","12","13")) 

stat<-exp_melt.final.2 %>% group_by(Donor,Time) %>% 
  dplyr::summarise(median = median(log2_count))
###########################################################################
####以中位数表达作为代表进行表达量比较

head(exp_melt.final);dim(exp_melt.final)

data_median<-exp_melt.final %>% group_by(Donor,Time,Sex) %>% 
  dplyr::summarise(median = median(tpm))

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


#####################################################################
############样本聚类热图D1
exp.tpm<-read.delim("all_sample_tpm.txt",header=T,row.names = 1)
dim(exp.tpm)
exp.tpm[1:4,1:4]
exp.tpm<-exp.tpm[rownames(expr),]
dim(exp.tpm)
identical(rownames(expr),rownames(exp.tpm))
exp.tpm$gene_biotype<-expr$gene_biotype
head(exp.tpm)

exp.tpm.pro<-subset(exp.tpm,gene_biotype %in% c("protein_coding"))
exp.tpm.lnc<-subset(exp.tpm,gene_biotype %in% c("lncRNA"))
dim(exp.tpm)
head(exp.tpm)
df<-exp.tpm[,c(2:31)]

head(df);dim(df)
phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
phenotype.d1<-phenotype %>% subset(Time=="D1")

####D1 exp matrix
df1<-df[,phenotype.d1$sample]

library(genefilter)
library(Biobase)
data1<-df1[apply(df1,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(data1)
#log(TPM+1)
#data1<-log2(data1+1)

normalized_data<-as.data.frame(data1)
normalized_data$geneid<-rownames(normalized_data)

phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
phenotype.d1<-phenotype %>% subset(Time=="D1")

head(normalized_data)
library("reshape")
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","tpm")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype.d1,by = "sample",all=T)
head(exp_melt.final)

# init 
decor <- matrix(1,length(colnames(data1)),length(colnames(data1)))
colnames(decor) <- colnames(data1)
rownames(decor) <- colnames(data1)
decor<-round(cor(data1),4)
# plot heatmap
library(latticeExtra)
library(lattice)
library(RColorBrewer)

x <- decor
dd.row <- as.dendrogram(hclust(dist(x)))
row.ord <- order.dendrogram(dd.row)

dd.col <- as.dendrogram(hclust(dist(t(x))))
col.ord <- order.dendrogram(dd.col)

p <- levelplot(x[row.ord, col.ord],
               aspect = "fill",
               scales = list(x = list(rot = 90)),
               colorkey = list(space = "left"),
               col.regions=colorRampPalette(brewer.pal(5,"YlOrRd"))(256),
               main="All D1 samples TPM correlation",
               legend =
                 list(right =
                        list(fun = dendrogramGrob,
                             args =
                               list(x = dd.col, ord = col.ord,
                                    side = "right",
                                    size = 10)),
                      top =
                        list(fun = dendrogramGrob,
                             args =
                               list(x = dd.row,
                                    side = "top",
                                    size = 10))), xlab="Sample",ylab="")

# output
file <- paste("./D1", "/all_sample_cluster.D1.png", sep="")
png(filename=file, height = 2200, width = 2200, res = 300, units = "px")
print(p)
dev.off()

file <- paste("./D1", "/all_sample_cluster.D1.pdf", sep="")
pdf(file=file, height = 8, width = 8)
print(p)
dev.off()

###############################################################################
#################################
####相关性画图
file <- paste("./D1", "/all_sample_pearson.cor.D1.xls", sep="")
write.table(decor, file=file, sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
###########################plot correlation ###########################################################
colors =colorRampPalette(brewer.pal(9,"YlOrRd")[3:7])((max(decor)-min(decor))*10000+1)
#colors=rev(colors)
#colors
pal=brewer.pal(9,"Set1")
panel.cor <- function(x, y, digits = 4, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  x=10^x-1
  
  y=10^y-1
  r <- round(cor(x, y),digits)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  #txt <- paste(prefix, r, sep = "")
  if (missing(cex.cor)) {
    cex.cor <- 0.8/strwidth(txt)
  }
  marg<-par("usr")
  mcol=floor((r-min(decor))*10000)+1
  #cat(paste(colors[mcol]),mcol,r,"\n")
  rect(marg[1],marg[3],marg[2],marg[4],col=colors[mcol])
  
  #rect(x[1],x[3],x[2],x[4],col="lightgray ")
  text(0.5, 0.5, txt, cex = cex.cor * (1 + r)/2,col="black")
}
# To show histograms of each variable along the diagonal
panel.hist <- function(x, ...) {
  x=log10(10^x-1)
  #x=x[!is.nan(x)]
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  #par(usr = c(usr[1:2], 0, 0.5))
  #  lines(density(x),col="red",lwd=2)
  h <- hist(x, breaks=30,plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="white",border="black", ...)
}

panel.lm <- function(x, y, col =pal[1], bg = NA, pch = par("pch"), cex = 1, 
                     col.smooth = pal[2], ...) {
  #x=log10(x+1)
  #y=log10(y+1)
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  rug(x,side=1,col=col,ticksize=0.02,lwd=0.2)
  rug(y,side=2,col=col,ticksize=0.02,lwd=0.2)
  #colormap(flipud(cool(256)))
  abline(stats::lm(y ~ x), col = col.smooth, ...)
}


file <- paste("./D1", "/all_D1sample_cor.png", sep="")
png(filename=file, height = 2200, width = 2200, res = 400, units = "px")

pairs(log10(data1+1), pch = ".", upper.panel = panel.cor, diag.panel = panel.hist, 
      lower.panel = panel.lm, main = "Pearson Correlation of all samples",gap=0)
dev.off()


file <- paste("./D1", "/all_D1sample_cor.pdf", sep="")
pdf(file=file, height = 4, width = 4)
pairs(log10(data1+1), pch = ".", upper.panel = panel.cor, diag.panel = panel.hist, 
      lower.panel = panel.lm, main = "Pearson Correlation of all D1 samples",gap=0)
dev.off()

#################################################################3
####density plot all gene
head(exp_melt);dim(exp_melt)
head(exp_melt.final)
exp_melt.final$log10tmp<-log10(exp_melt.final$tpm)
test<-ggplot(exp_melt.final, aes(x=log10tmp, color=sample)) +
  geom_density()+
  theme_bw()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
test
ggsave(filename="./D1/density.plot.all_D1sample.tpm.pdf",plot=test,
       device='pdf',path=".",width=4.5,height=2.8)
###by time
head(exp_melt);dim(exp_melt)
head(exp_melt.final)
#exp_melt.final$log10tmp<-log10(exp_melt.final$tpm)
head(exp_melt.final)
exp_melt.final$Donor<-as.factor(exp_melt.final$Donor)
test2<-ggplot(exp_melt.final, aes(x=log10tmp, color=Donor)) +
  geom_density()+
  theme_bw()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  facet_grid( ~ exp_melt.final$Time, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  
test2
ggsave(filename="./D1/density.byTime.plot.all_D1sample.tpm.pdf",plot=test2,
       device='pdf',path=".",width=4.5,height=2.8)
###by donor
head(exp_melt.final)
test3<-ggplot(exp_melt.final, aes(x=log10tmp, color=Time)) +
  geom_density()+
  theme_bw()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  facet_grid( ~ exp_melt.final$Donor, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  
test3
#ggsave(filename="density.byDonor.plot.all_sample.tpm.pdf",plot=test3,
 #      device='pdf',path=".",width=20,height=3)
######protein coding gene
####density plot protein coding gene
head(df.pro);dim(df.pro)
df1<-df.pro[,phenotype.d1$sample]

data1<-df1[apply(df1,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(data1)
#log(TPM+1)
#data1<-log2(data1+1)

normalized_data<-as.data.frame(data1)
normalized_data$geneid<-rownames(normalized_data)

head(normalized_data)
library("reshape")
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","tpm")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype.d1,by = "sample",all=T)
head(exp_melt.final)

head(exp_melt);dim(exp_melt)
head(exp_melt.final)
exp_melt.final$log10tmp<-log10(exp_melt.final$tpm)
test<-ggplot(exp_melt.final, aes(x=log10tmp, color=sample)) +
  geom_density()+
  theme_bw()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(title="Protein coding")
test
ggsave(filename="./D1/density.plot.protein_coding.all_D1sample.tpm.pdf",plot=test,
       device='pdf',path=".",width=4.5,height=2.8)
###by time
head(exp_melt);dim(exp_melt)
head(exp_melt.final)
#exp_melt.final$log10tmp<-log10(exp_melt.final$tpm)
head(exp_melt.final)
exp_melt.final$Donor<-as.factor(exp_melt.final$Donor)
test2<-ggplot(exp_melt.final, aes(x=log10tmp, color=Donor)) +
  geom_density()+
  theme_bw()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  facet_grid( ~ exp_melt.final$Time, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title="Protein coding")
test2
ggsave(filename="./D1/density.protein_coding.byTime.plot.all_D1sample.tpm.pdf",plot=test2,
       device='pdf',path=".",width=4.5,height=2.8)
######lncRNA gene
####density plot lncRNA gene
head(df.lnc);dim(df.lnc)
df1<-df.lnc[,phenotype.d1$sample]

data1<-df1[apply(df1,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(data1)
#log(TPM+1)
#data1<-log2(data1+1)

normalized_data<-as.data.frame(data1)
normalized_data$geneid<-rownames(normalized_data)

head(normalized_data)
library("reshape")
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","tpm")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype.d1,by = "sample",all=T)
head(exp_melt.final)

head(exp_melt);dim(exp_melt)
head(exp_melt.final)
exp_melt.final$log10tmp<-log10(exp_melt.final$tpm)
test<-ggplot(exp_melt.final, aes(x=log10tmp, color=sample)) +
  geom_density()+
  theme_bw()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(title="LncRNA")
test
ggsave(filename="./D1/density.plot.lncrna.all_D1sample.tpm.pdf",plot=test,
       device='pdf',path=".",width=4.5,height=2.8)
###by time
head(exp_melt);dim(exp_melt)
head(exp_melt.final)
#exp_melt.final$log10tmp<-log10(exp_melt.final$tpm)
head(exp_melt.final)
exp_melt.final$Donor<-as.factor(exp_melt.final$Donor)
test2<-ggplot(exp_melt.final, aes(x=log10tmp, color=Donor)) +
  geom_density()+
  theme_bw()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  facet_grid( ~ exp_melt.final$Time, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title="Lncrna")
test2
ggsave(filename="./D1/density.lncrna.byTime.plot.all_D1sample.tpm.pdf",plot=test2,
       device='pdf',path=".",width=4.5,height=2.8)
############################################################
############统计低中高表达基因
###所有基因（不区分编码与否）
head(df1)
dim(df1)
# TPM：低<1，中1~100，高>100
#length(df1$S0112[df1$S0112!=0])
#length(df1$S0112[df1$S0112<=1])
#length(df1$S0112[df1$S0112>1 & df1$S0112<100])
#length(df1$S0112[df1$S0112>=100])
list<-colnames(df1)
all_result=data.frame()
for(i in 1:length(list)){
  t<-df1[,i]
  a<-length(t[t>0 & t<=1])
  b<-length(t[t>1 & t<100])
  c<-length(t[t>=100])
  d<-cbind(a,b,c)
  all_result<-rbind(all_result,d)
}
head(all_result)
rownames(all_result)<-colnames(df1)
colnames(all_result)<-c("low","medium","high")
head(all_result)
save(all_result,file="./D1/all_result.RData")
write.table(all_result,
            "./D1/all_result.D1.sample.low.medium.high.gene.count.txt", 
            sep="\t", col.names=T,quote = F)

########protein coding
exp.tpm.pro<-subset(exp.tpm,gene_biotype %in% c("protein_coding"))
exp.tpm.lnc<-subset(exp.tpm,gene_biotype %in% c("lncRNA"))
dim(exp.tpm.pro)
head(exp.tpm.pro)
df.pro<-exp.tpm.pro[,c(2:31)]

head(df.pro);dim(df.pro)
phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
phenotype.d1<-phenotype %>% subset(Time=="D1")

####D1 exp matrix
df1<-df.pro[,phenotype.d1$sample]
head(df1)
dim(df1)
# TPM：低<1，中1~100，高>100
#length(df1$S0112[df1$S0112!=0])
#length(df1$S0112[df1$S0112<=1])
#length(df1$S0112[df1$S0112>1 & df1$S0112<100])
#length(df1$S0112[df1$S0112>=100])
list<-colnames(df1)
all_result=data.frame()
for(i in 1:length(list)){
  t<-df1[,i]
  a<-length(t[t>0 & t<=1])
  b<-length(t[t>1 & t<100])
  c<-length(t[t>=100])
  d<-cbind(a,b,c)
  all_result<-rbind(all_result,d)
}
head(all_result)
rownames(all_result)<-colnames(df1)
colnames(all_result)<-c("low","medium","high")
head(all_result)
save(all_result,file="./D1/all_result.pro.RData")
write.table(all_result,
            "./D1/all_result.protein_coding.D1.sample.low.medium.high.gene.count.txt", 
            sep="\t", col.names=T,quote = F)
################lnc
exp.tpm.pro<-subset(exp.tpm,gene_biotype %in% c("protein_coding"))
exp.tpm.lnc<-subset(exp.tpm,gene_biotype %in% c("lncRNA"))
dim(exp.tpm.lnc)
head(exp.tpm.lnc)
df.lnc<-exp.tpm.lnc[,c(2:31)]

head(df.lnc);dim(df.lnc)
phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
phenotype.d1<-phenotype %>% subset(Time=="D1")

####D1 exp matrix
df1<-df.lnc[,phenotype.d1$sample]
head(df1)
dim(df1)
# TPM：低<1，中1~100，高>100
#length(df1$S0112[df1$S0112!=0])
#length(df1$S0112[df1$S0112<=1])
#length(df1$S0112[df1$S0112>1 & df1$S0112<100])
#length(df1$S0112[df1$S0112>=100])
list<-colnames(df1)
all_result=data.frame()
for(i in 1:length(list)){
  t<-df1[,i]
  a<-length(t[t>0 & t<=1])
  b<-length(t[t>1 & t<100])
  c<-length(t[t>=100])
  d<-cbind(a,b,c)
  all_result<-rbind(all_result,d)
}
head(all_result)
rownames(all_result)<-colnames(df1)
colnames(all_result)<-c("low","medium","high")
head(all_result)
save(all_result,file="./D1/all_result.lnc.RData")
write.table(all_result,
            "./D1/all_result.lnc.D1.sample.low.medium.high.gene.count.txt", 
            sep="\t", col.names=T,quote = F)
##################################################################
##山脊线图

# library
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)

head(df.pro);dim(df.pro)
df1<-df.pro[,phenotype.d1$sample]

data1<-df1[apply(df1,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(data1)
#log(TPM+1)
#data1<-log2(data1+1)

normalized_data<-as.data.frame(data1)
normalized_data$geneid<-rownames(normalized_data)

head(normalized_data)
library("reshape")
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","tpm")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype.d1,by = "sample",all=T)
head(exp_melt.final)

head(exp_melt);dim(exp_melt)
head(exp_melt.final)
exp_melt.final$tpm<-as.numeric(exp_melt.final$tpm)
exp_melt.final$log10tmp<-log10(exp_melt.final$tpm)
head(exp_melt.final)
# Plot
p1<-ggplot(exp_melt.final, aes(x = tpm, y = sample, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_gradientn(
    colours = c("blue", "#CC4678FF", "#F0F921FF"),
    name = "TPM"
  )+
  labs(title = 'Gene expression in D1 samples') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )+xlim(100,1000)+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title="High expression--protein coding")
p1

p1<-ggplot(exp_melt.final, aes(x = tpm, y = sample, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_gradientn(
    colours = c("blue", "#CC4678FF", "#F0F921FF"),
    name = "TPM"
  )+
  labs(title = 'Gene expression in D1 samples') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )+xlim(1,100)+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title="Medium expression--protein coding")
p1

p1<-ggplot(exp_melt.final, aes(x = tpm, y = sample, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_gradientn(
    colours = c("blue", "#CC4678FF", "#F0F921FF"),
    name = "TPM"
  )+
  labs(title = 'Gene expression in D1 samples') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )+xlim(0,1)+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title="Low expression--protein coding")
p1
#pdf(file="./D1/ggridges.pro.gene.high.exp.D1sample.tpm.pdf",5,3)
#p1
#dev.off()

###山脊线图lnc
head(df.lnc);dim(df.lnc)
df1<-df.lnc[,phenotype.d1$sample]

data1<-df1[apply(df1,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(data1)
#log(TPM+1)
#data1<-log2(data1+1)

normalized_data<-as.data.frame(data1)
normalized_data$geneid<-rownames(normalized_data)

head(normalized_data)
library("reshape")
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","tpm")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype.d1,by = "sample",all=T)
head(exp_melt.final)

head(exp_melt);dim(exp_melt)
head(exp_melt.final)
exp_melt.final$log10tmp<-log10(exp_melt.final$tpm)
head(exp_melt.final)
# Plot
p1<-ggplot(exp_melt.final, aes(x = tpm, y = sample, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_gradientn(
    colours = c("blue", "#CC4678FF", "#F0F921FF"),
    name = "Temp"
  )+
  labs(title = 'Gene expression in D1 samples') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )+xlim(100,1000)+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title="High expression--lncrna")
p1

p1<-ggplot(exp_melt.final, aes(x = tpm, y = sample, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_gradientn(
    colours = c("blue", "#CC4678FF", "#F0F921FF"),
    name = "Temp"
  )+
  labs(title = 'Gene expression in D1 samples') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )+xlim(1,100)+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title="Medium expression--lncrna")
p1

p1<-ggplot(exp_melt.final, aes(x = tpm, y = sample, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_gradientn(
    colours = c("blue", "#CC4678FF", "#F0F921FF"),
    name = "Temp"
  )+
  labs(title = 'Gene expression in D1 samples') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )+xlim(0,1)+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold",
                                   angle = 0, hjust = 0,lineheight=100
  ),
  axis.text.y = element_text(size=10,colour="black",face="bold"),
  legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title="Low expression--lncrna")
p1
#pdf(file="./D1/ggridges.lnc.gene.high.exp.D1sample.tpm.pdf",5,3)
#p1
#dev.off()
###################################################################
###############################################################################
# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)  # 用于科学期刊配色
library(scales)  # 用于百分比格式化

# 1. 创建表达水平分类列
exp_melt.final <- exp_melt.final %>%
  mutate(
    expression_level = case_when(
      tpm < 1 ~ "Low (0-1 TPM)",
      tpm >= 1 & tpm < 100 ~ "Medium (1-100 TPM)",
      tpm >= 100 ~ "High (≥100 TPM)",
      TRUE ~ NA_character_
    ),
    expression_level = factor(expression_level, 
                              levels = c("High (≥100 TPM)", "Medium (1-100 TPM)", "Low (0-1 TPM)"))
  )

# 2. 计算每个样本的表达水平分布
expression_summary <- exp_melt.final %>%
  group_by(sample, expression_level, Time, Sex, Donor) %>%
  summarise(
    gene_count = n(),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  mutate(
    total_genes = sum(gene_count),
    percentage = gene_count / total_genes * 100
  ) %>%
  ungroup() %>%
  # 为排序添加信息
  mutate(
    sample_order = factor(sample, levels = unique(sample[order(Time, Sex, Donor)]))
  )

# 3. 创建分组热图式条形图
p1 <- ggplot(expression_summary, aes(x = sample_order, y = percentage, fill = expression_level)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", size = 0.2) +
  scale_fill_jco(  # JAMA期刊配色，对色盲友好
    name = "Expression Level",
    labels = c("High (≥100 TPM)", "Medium (1-100 TPM)", "Low (0-1 TPM)")
  ) +
  labs(
    title = "Distribution of Gene Expression Levels Across Samples",
    subtitle = "Classification based on TPM values: Low (0-1), Medium (1-100), High (≥100)",
    x = "Sample ID",
    y = "Percentage of Genes (%)",
    caption = paste("Total genes per sample:", unique(expression_summary$total_genes)[1])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 15)),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor.y = element_blank(),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0)
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05))
  )

print(p1)

# 4. 创建点线图展示具体数值
p2 <- expression_summary %>%
  ggplot(aes(x = sample_order, y = percentage, color = expression_level, group = expression_level)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_jco(
    name = "Expression Level",
    labels = c("High (≥100 TPM)", "Medium (1-100 TPM)", "Low (0-1 TPM)")
  ) +
  labs(
    title = "Expression Level Trends Across Samples",
    subtitle = "Line plot showing percentage distribution of low/medium/high expression genes",
    x = "Sample ID",
    y = "Percentage of Genes (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 15)),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, NA)
  )

print(p2)

# 5. 创建热图式矩阵图（更紧凑的展示）
# 计算每个样本的平均表达水平数值
sample_summary <- expression_summary %>%
  pivot_wider(
    id_cols = c(sample, Time, Sex, Donor),
    names_from = expression_level,
    values_from = percentage,
    values_fill = 0
  ) %>%
  mutate(
    sample_group = paste(Time, Sex, sep = "_")
  )

# 转换为长格式用于热图
sample_long <- sample_summary %>%
  pivot_longer(
    cols = c(`High (≥100 TPM)`, `Medium (1-100 TPM)`, `Low (0-1 TPM)`),
    names_to = "expression_level",
    values_to = "percentage"
  ) %>%
  mutate(
    expression_level = factor(expression_level, 
                              levels = c("High (≥100 TPM)", "Medium (1-100 TPM)", "Low (0-1 TPM)")),
    sample = factor(sample, levels = unique(sample[order(Time, Sex, Donor)]))
  )

p3 <- ggplot(sample_long, aes(x = expression_level, y = sample, fill = percentage)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradientn(
    colors = c("#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C"),
    name = "Percentage (%)",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0%", "25%", "50%", "75%", "100%")
  ) +
  labs(
    title = "Heatmap of Expression Level Distribution",
    subtitle = "Each cell shows the percentage of genes in each expression category",
    x = "Expression Level Category",
    y = "Sample ID"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black", face = "bold"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 15)),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  # 添加数值标签
  geom_text(
    aes(label = ifelse(percentage > 5, sprintf("%.1f%%", percentage), "")),
    color = "black",
    size = 3
  )

print(p3)

# 6. 按分组展示的小提琴图 + 箱线图
p4 <- expression_summary %>%
  ggplot(aes(x = expression_level, y = percentage, fill = expression_level)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  scale_fill_jco() +
  labs(
    title = "Distribution of Expression Level Categories Across All Samples",
    subtitle = "Violin plots show density, boxplots show quartiles",
    x = "Expression Level Category",
    y = "Percentage of Genes (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, color = "black", face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 15)),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, NA)
  ) +
  # 添加统计摘要
  stat_summary(
    fun.data = function(x) {
      data.frame(
        y = max(x) * 1.05,
        label = paste(
          "n =", length(x), "\n",
          "Mean =", round(mean(x), 1), "%\n",
          "SD =", round(sd(x), 1), "%"
        )
      )
    },
    geom = "text",
    size = 3,
    vjust = 0
  )

print(p4)

# 7. 保存图表
# ggsave("expression_level_distribution.pdf", p1, width = 12, height = 8, dpi = 300)
# ggsave("expression_level_trends.pdf", p2, width = 12, height = 8, dpi = 300)
# ggsave("expression_heatmap.pdf", p3, width = 10, height = 12, dpi = 300)
# ggsave("expression_violin_plots.pdf", p4, width = 10, height = 8, dpi = 300)
##########################################################3
###山脊线图
# 加载必要的包
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(viridis)

# 1. 创建表达水平分类列
exp_melt.final <- exp_melt.final %>%
  mutate(
    expression_level = case_when(
      tpm < 1 ~ "Low (0-1 TPM)",
      tpm >= 1 & tpm < 100 ~ "Medium (1-100 TPM)",
      tpm >= 100 ~ "High (≥100 TPM)",
      TRUE ~ NA_character_
    ),
    expression_level = factor(expression_level, 
                              levels = c("Low (0-1 TPM)", "Medium (1-100 TPM)", "High (≥100 TPM)"))
  )

# 2. 计算每个样本每个表达水平的log10(tpm+1)密度分布
# 添加1避免log10(0)的问题
exp_melt.final <- exp_melt.final %>%
  mutate(
    log10_tpm_plus1 = log10(tpm + 1)
  )

# 3. 创建山脊线图
p <- ggplot(exp_melt.final, aes(
  x = log10_tpm_plus1, 
  y = sample, 
  fill = expression_level
)) +
  geom_density_ridges(
    scale = 1.2,
    rel_min_height = 0.01,
    alpha = 0.7,
    size = 0.3,
    color = "white"
  ) +
  # 使用viridis配色，对色盲友好
  scale_fill_viridis_d(
    name = "Expression Level",
    option = "D",  # 使用viridis色系的D变体
    direction = -1,  # 反转颜色顺序
    labels = c("Low (0-1 TPM)", "Medium (1-100 TPM)", "High (≥100 TPM)")
  ) +
  # 添加垂直线标记TPM阈值
  geom_vline(xintercept = log10(1 + 1), linetype = "dashed", color = "gray30", size = 0.4) +
  geom_vline(xintercept = log10(100 + 1), linetype = "dashed", color = "gray30", size = 0.4) +
  # 添加阈值标签
  annotate("text", 
           x = log10(1 + 1), 
           y = length(unique(exp_melt.final$sample)) + 1.5, 
           label = "1 TPM", 
           hjust = -0.1, 
           vjust = 0.5, 
           size = 3.5,
           color = "gray30") +
  annotate("text", 
           x = log10(100 + 1), 
           y = length(unique(exp_melt.final$sample)) + 1.5, 
           label = "100 TPM", 
           hjust = -0.1, 
           vjust = 0.5, 
           size = 3.5,
           color = "gray30") +
  labs(
    title = "Expression Level Distribution Across Samples",
    subtitle = "Ridgeline plot showing low (0-1 TPM), medium (1-100 TPM), and high (≥100 TPM) expression genes",
    x = expression(log[10](TPM + 1)),
    y = "Sample ID",
    caption = paste("Total samples:", length(unique(exp_melt.final$sample)), 
                    "| Total genes:", nrow(exp_melt.final))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 15)),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.height = unit(1.2, "cm"),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_line(color = "gray95", size = 0.2),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0)
  ) +
  # 调整x轴范围，确保显示所有数据
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.05)),
    breaks = c(0, log10(1 + 1), 1, log10(100 + 1), 2, 3),
    labels = c("0", "1", "10", "100", "1000", "10000")
  ) +
  # 按样本分组排序
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.1)))

# 显示图表
print(p)

# 4. 可选：添加统计摘要信息到图表中
# 计算每个样本每个表达水平的基因数量
summary_stats <- exp_melt.final %>%
  group_by(sample, expression_level) %>%
  summarise(
    gene_count = n(),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  mutate(
    total_genes = sum(gene_count),
    percentage = round(gene_count / total_genes * 100, 1)
  ) %>%
  ungroup()

# 创建带有统计信息的版本
p_with_stats <- p +
  # 在右侧添加统计信息
  geom_text(
    data = summary_stats,
    aes(x = 3.5, 
        y = sample, 
        label = paste0(
          expression_level, ": ", 
          percentage, "% (n=", gene_count, ")"
        ),
        color = expression_level),
    hjust = 0,
    size = 2.8,
    position = position_dodge(width = 0.8),
    show.legend = FALSE
  ) +
  scale_color_viridis_d(
    option = "D",
    direction = -1
  ) +
  # 扩展x轴范围以容纳统计信息
  scale_x_continuous(
    limits = c(NA, 5),
    expand = expansion(mult = c(0.02, 0.15)),
    breaks = c(0, log10(1 + 1), 1, log10(100 + 1), 2, 3),
    labels = c("0", "1", "10", "100", "1000", "10000")
  ) +
  labs(
    subtitle = "Ridgeline plot with expression level distributions and percentage statistics"
  )

print(p_with_stats)

# 5. 可选：分面显示不同表达水平
p_faceted <- ggplot(exp_melt.final, aes(
  x = log10_tpm_plus1, 
  y = sample, 
  fill = expression_level
)) +
  geom_density_ridges(
    scale = 1.2,
    rel_min_height = 0.01,
    alpha = 0.8,
    size = 0.3,
    color = "white"
  ) +
  scale_fill_viridis_d(
    name = "Expression Level",
    option = "D",
    direction = -1,
    labels = c("Low (0-1 TPM)", "Medium (1-100 TPM)", "High (≥100 TPM)")
  ) +
  facet_wrap(~ expression_level, ncol = 1, scales = "free_y") +
  labs(
    title = "Expression Level Distribution by Category",
    subtitle = "Separate ridgeline plots for low, medium, and high expression genes",
    x = expression(log[10](TPM + 1)),
    y = "Sample ID"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 15)),
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold", margin = margin(b = 5)),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_line(color = "gray95", size = 0.2)
  ) +
  scale_x_continuous(
    breaks = c(0, log10(1 + 1), 1, log10(100 + 1), 2, 3),
    labels = c("0", "1", "10", "100", "1000", "10000")
  )

print(p_faceted)

# 6. 保存图表
# ggsave("expression_ridgeline_plot.pdf", p, width = 12, height = 10, dpi = 300)
# ggsave("expression_ridgeline_with_stats.pdf", p_with_stats, width = 14, height = 10, dpi = 300)
# ggsave("expression_ridgeline_faceted.pdf", p_faceted, width = 10, height = 14, dpi = 300)