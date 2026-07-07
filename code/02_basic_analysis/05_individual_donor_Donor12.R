
######
####################数据基本分析
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后")
library(Mfuzz)
library(Biobase)
library(limma)
library(clusterProfiler)


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
############样本聚类热图donor12
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后")
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
df<-exp.tpm[,c(1:30)]

head(df);dim(df)
phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
phenotype.d6<-phenotype %>% subset(Donor=="6")
phenotype.d13<-phenotype %>% subset(Donor=="13")
phenotype.d12<-phenotype %>% subset(Donor=="12")
####donor12 exp matrix
df1<-df[,phenotype.d12$sample]

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
exp_melt.final<-merge(exp_melt,phenotype.d12,by = "sample",all=T)
head(exp_melt.final)

# init 
decor <- matrix(1,length(colnames(data1)),length(colnames(data1)))
colnames(decor) <- colnames(data1)
rownames(decor) <- colnames(data1)
decor<-round(cor(data1),4)
head(decor)
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
               main="All donor12 samples TPM correlation",
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
file <- paste("./donor12", "/all_sample_cluster.donor12.png", sep="")
png(filename=file, height = 1800, width = 1800, res = 300, units = "px")
print(p)
dev.off()

file <- paste("./donor12", "/all_sample_cluster.donor12.pdf", sep="")
pdf(file=file, height = 6, width = 6)
print(p)
dev.off()

###############################################################################
#################################
####相关性画图
file <- paste("./donor12", "/all_sample_pearson.cor.donor12.xls", sep="")
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


file <- paste("./donor12", "/all_donor12sample_cor.png", sep="")
png(filename=file, height = 3000, width = 3000, res = 300, units = "px")

pairs(log10(data1+1), pch = ".", upper.panel = panel.cor, diag.panel = panel.hist, 
      lower.panel = panel.lm, main = "Pearson Correlation of donor12 samples",gap=0)
dev.off()


file <- paste("./donor12", "/all_donor12sample_cor.pdf", sep="")
pdf(file=file, height = 3, width = 3)
pairs(log10(data1+1), pch = ".", upper.panel = panel.cor, diag.panel = panel.hist, 
      lower.panel = panel.lm, main = "Pearson Correlation of all donor12 samples",gap=0)
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
ggsave(filename="./donor12/density.plot.all_donor12sample.tpm.pdf",plot=test,
       device='pdf',path=".",width=3.5,height=2)

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
df1<-df.pro[,phenotype.d12$sample]

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
exp_melt.final<-merge(exp_melt,phenotype.d12,by = "sample",all=T)
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
ggsave(filename="./donor12/density.plot.protein_coding.all_donor12sample.tpm.pdf",plot=test,
       device='pdf',path=".",width=3.5,height=2)
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
#ggsave(filename="./D1/density.protein_coding.byTime.plot.all_D1sample.tpm.pdf",plot=test2,
    #   device='pdf',path=".",width=4.5,height=2.8)
######lncRNA gene
####density plot lncRNA gene
head(df.lnc);dim(df.lnc)
df1<-df.lnc[,phenotype.d12$sample]

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
exp_melt.final<-merge(exp_melt,phenotype.d12,by = "sample",all=T)
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
ggsave(filename="./donor12/density.plot.lncrna.all_donor12sample.tpm.pdf",plot=test,
       device='pdf',path=".",width=3.5,height=2)
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
df1<-df[,phenotype.d12$sample]
head(df)
dim(df)
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
save(all_result,file="./donor12/all_result.RData")
write.table(all_result,
            "./donor12/all_result.donor12.sample.low.medium.high.gene.count.txt", 
            sep="\t", col.names=T,quote = F)

########protein coding
exp.tpm.pro<-subset(exp.tpm,gene_biotype %in% c("protein_coding"))
exp.tpm.lnc<-subset(exp.tpm,gene_biotype %in% c("lncRNA"))
dim(exp.tpm.pro)
head(exp.tpm.pro)
df.pro<-exp.tpm.pro[,c(1:30)]

head(df.pro);dim(df.pro)
phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
phenotype.d1<-phenotype %>% subset(Time=="D1")

####D1 exp matrix
df1<-df.pro[,phenotype.d12$sample]
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
save(all_result,file="./donor12/all_result.pro.RData")
write.table(all_result,
            "./donor12/all_result.protein_coding.donor12.sample.low.medium.high.gene.count.txt", 
            sep="\t", col.names=T,quote = F)
################lnc
exp.tpm.pro<-subset(exp.tpm,gene_biotype %in% c("protein_coding"))
exp.tpm.lnc<-subset(exp.tpm,gene_biotype %in% c("lncRNA"))
dim(exp.tpm.lnc)
head(exp.tpm.lnc)
df.lnc<-exp.tpm.lnc[,c(1:30)]

head(df.lnc);dim(df.lnc)
phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
phenotype.d1<-phenotype %>% subset(Time=="D1")

####D1 exp matrix
df1<-df.lnc[,phenotype.d12$sample]
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
save(all_result,file="./donor12/all_result.lnc.RData")
write.table(all_result,
            "./donor12/all_result.lnc.donor12.sample.low.medium.high.gene.count.txt", 
            sep="\t", col.names=T,quote = F)
##################################################################
##山脊线图

# library
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)

head(df.pro);dim(df.pro)
df1<-df.pro[,phenotype.d12$sample]

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
exp_melt.final<-merge(exp_melt,phenotype.d12,by = "sample",all=T)
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
  labs(title="High expression--protein coding")
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
  labs(title="Medium expression--protein coding")
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
  labs(title="Low expression--protein coding")
p1
#pdf(file="./D1/ggridges.pro.gene.high.exp.D1sample.tpm.pdf",5,3)
#p1
#dev.off()

###山脊线图lnc
head(df.lnc);dim(df.lnc)
df1<-df.lnc[,phenotype.d12$sample]

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
exp_melt.final<-merge(exp_melt,phenotype.d12,by = "sample",all=T)
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
