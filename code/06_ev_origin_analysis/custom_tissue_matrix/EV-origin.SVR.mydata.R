setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\自己构建组织特异性表达矩阵\\GTEx\\EV-origin")
library(e1071)
### nu-SVR
### read the reference matrix as "ref.m"
list.files()
expr<-read.delim("./exp.tpm.protein_coding.txt",header=T,row.names = 1)
head(expr)

exp<-expr

x <- data.matrix(exp)
x <- (x - mean(x)) / sd(as.vector(x))
avdata.m <- x#the expression profile normalization
nu.v = c(0.25, 0.5, 0.75)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do.exLR_origin <- function(avdata.m, ref.m, nu.v) {
  map.idx <- match(rownames(ref.m), rownames(avdata.m))
  rep.idx <- which(is.na(map.idx) == FALSE)
  data2.m <- avdata.m[map.idx[rep.idx], ]
  ref2.m <- ref.m[rep.idx, ]
  est.ab.lm <- list()
  est.lm <- list()
  nui <- 1
  for (nu in nu.v) {
    est.m <- matrix(NA, nrow = ncol(data2.m), ncol = ncol(ref2.m))
    colnames(est.m) <- colnames(ref2.m)
    rownames(est.m) <- colnames(data2.m)
    est.ab.m <- matrix(NA, nrow = ncol(data2.m), ncol = ncol(ref2.m))
    colnames(est.ab.m) <- colnames(ref2.m)
    rownames(est.ab.m) <- colnames(data2.m)
    
    for (s in seq_len(ncol(data2.m))) {
      svm.o <- svm(x = ref2.m, y = data2.m[, s], scale = TRUE, type = "nu-regression", kernel = "linear", nu = nu)
      coef.v <- t(svm.o$coefs) %*% svm.o$SV
      coef.v[which(coef.v < 0)] <- 1*10^-10
      est.ab.m[s,] <- coef.v
      total <- sum(coef.v)
      coef.v <- coef.v/total
      est.m[s, ] <- coef.v
    }
    est.lm[[nui]] <- est.m
    est.ab.lm[[nui]] <- est.ab.m
    nui <- nui + 1
  }
  
  #### select best nu using RMSE
  rmse.m <- matrix(NA, nrow = ncol(avdata.m), ncol = length(nu.v))
  for (nui in seq_along(nu.v)) {
    reconst.m <- ref2.m %*% t(est.lm[[nui]])
    s <- seq_len(ncol(avdata.m))
    rmse.m[s, nui] <- sqrt(colMeans((data2.m[, s] - reconst.m[, s])^2))
    message(nui)
  }
  colnames(rmse.m) <- nu.v
  nu.idx <- apply(rmse.m, 1, which.min)
  estF.m <- est.m
  for (s in seq_len(nrow(estF.m))) {
    estF.m[s, ] <- est.lm[[nu.idx[s]]][s, ]
  }
  estF.ab.m <- est.ab.m
  for (s in seq_len(nrow(estF.ab.m))) {
    estF.ab.m[s, ] <- est.ab.lm[[nu.idx[s]]][s, ]
  }
  #selecting min RMSE
  rmse.min.value <- as.data.frame(apply(rmse.m, 1, min))
  rownames(rmse.min.value) <- colnames(avdata.m)
  colnames(rmse.min.value)[1] <- 'RMSE'
  #caculating PCC
  pearson.corr.value <- c()
  for (i in 1:ncol(data2.m)) {
    cor.index <- cor.test(data2.m[, i], reconst.m[, i])
    cor.p <- as.numeric(cor.index$estimate)
    pearson.corr.value <- c(pearson.corr.value, cor.p)
  }
  estF.m <- cbind.data.frame(estF.m, pearson.corr.value, rmse.min.value)
  return(list(estF = estF.m, est.ab.sum = estF.ab.m, nu = nu.v[nu.idx]))
}
ref.m<-read.delim("gtex.ref.16tissue.exp.matrix.txt",header=T,row.names = 1)
ref.m<-as.matrix(ref.m)
exLR_origin.results <- do.exLR_origin(avdata.m, ref.m, nu.v)
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
save(exLR_origin.results,file="exLR_origin.results.tissue.new.Rdata")
###结果导出
exLR_origin.results$est.ab.sum
write.table(exLR_origin.results$est.ab.sum,
            "mydata.EV.origin.tissue.ab.sum.new.txt", 
            sep="\t", col.names=T,quote = F,row.names = T)
###Matrix_Hemopoietic.cell.csv
#ref.m<-read.csv("Matrix_Hemopoietic.cell.csv",row.names=1,
#                stringsAsFactors = TRUE,na.strings=c('NA',''))
#ref.m<-as.matrix(ref.m)
#exLR_origin.results.cell <- do.exLR_origin(avdata.m, ref.m, nu.v)
#save(exLR_origin.results.cell,file="exLR_origin.results.cell.Rdata")
###结果导出
#exLR_origin.results.cell$est.ab.sum
#write.table(exLR_origin.results.cell$est.ab.sum,
 #           "mydata.EV.origin.Hemopoietic.cell.ab.sum.txt", 
  #          sep="\t", col.names=T,quote = F,row.names = T)

####################################################################
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\自己构建组织特异性表达矩阵\\GTEx\\EV-origin")
load("exLR_origin.results.tissue.new.Rdata")

head(exLR_origin.results)

phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)

#########bar plot
library(RColorBrewer)
display.brewer.all()
col <- colorRampPalette(brewer.pal(8,'Accent'))(7)
mycol<-colorRampPalette(brewer.pal(12,'Paired'),alpha = TRUE)(16)

col1<-colorRampPalette(brewer.pal(12,'Set3'),alpha = TRUE)(12)
col2<-colorRampPalette(brewer.pal(4,'Set1'),alpha = TRUE)(4)
mycol<-c(col1,col2)

b<-exLR_origin.results$est.ab.sum
head(b);dim(b)

a<-as.data.frame(b)

a$sample<-rownames(a)

head(a)


library("reshape")
library("ggplot2")
library(dplyr)
data_melt<-melt(a,id=c("sample"))
head(data_melt)
colnames(data_melt)<-c("sample","Type","Abundance")
head(data_melt)

c<-data_melt %>% group_by(sample) %>% mutate(RA=Abundance/sum(Abundance)*100)
head(c)

#####百分比相对丰度
p<-ggplot(data=c,aes(x=sample,y=RA,fill=Type))+
  geom_bar(stat="identity")+theme_bw()
p
p.final<-p+
  labs(y="Relative Abundance(%)",x="",title="")+
  scale_fill_manual(values=c(mycol))+
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1,
                                   vjust =1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=13,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=13,colour="black"),
        axis.title.y = element_text(size=13,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  scale_x_discrete(limit = unique(a$sample))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.final

ggsave(filename="barplot_EV_origin_tissue_AB.pdf",plot=p.final,
       device='pdf',path=".",width=10,height=4.5)
####D1 sample
head(c)
c.1<-c %>% subset(sample %in% c(phenotype.d1$sample))
c.1
#百分比相对丰度
p<-ggplot(data=c.1,aes(x=sample,y=RA,fill=Type))+
  geom_bar(stat="identity")+theme_bw()
p
p.final<-p+
  labs(y="Relative Abundance(%)",x="",title="")+
  scale_fill_manual(values=c(mycol))+
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1,
                                   vjust =1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=13,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=13,colour="black"),
        axis.title.y = element_text(size=13,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  scale_x_discrete(limit = unique(c.1$sample))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.final

ggsave(filename="barplot_EV_origin_tissue_AB.D1sample.new.pdf",plot=p.final,
       device='pdf',path=".",width=5,height=5)
##########
####donor 4\5\6\11\12\13
head(c)
phenotype.t<-phenotype %>% subset(Donor=="11")
c.2<-c %>% subset(sample %in% c(phenotype.t$sample))
c.2
#百分比相对丰度
p<-ggplot(data=c.2,aes(x=sample,y=RA,fill=Type))+
  geom_bar(stat="identity")+theme_bw()
p
p.final<-p+
  labs(y="Relative Abundance(%)",x="",title="")+
  scale_fill_manual(values=c(mycol))+
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1,
                                   vjust =1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=13,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=13,colour="black"),
        axis.title.y = element_text(size=13,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  scale_x_discrete(limit = unique(c.2$sample))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.final

ggsave(filename="barplot_EV_origin_tissue_AB.Donor11sample.pdf",plot=p.final,
       device='pdf',path=".",width=3.5,height=4.5)
#####
###丰度绝对值
head(c)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(c, x = "Type", y = "Abundance",
               palette =c("lightblue"),
               outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  xlab("") + ylab("Abundance")+
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

ggsave(filename="boxplot_EV_origin_tissue_Aboundance.pdf",plot=p,
       device='pdf',path=".",width=8,height=3.5)

###by sex
phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
head(c)

exp_melt.final<-merge(c,phenotype,by = "sample",all=T)
dim(exp_melt.final);head(exp_melt.final)

head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")

p <- ggboxplot(exp_melt.final, x = "Type", y = "Abundance",
               color = "Sex", palette = "nejm",
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
  labs(fill="Sex",x="",y="Abundance",title="") 
#p
# Use only p.format as label. Remove method name.
p2<-p + stat_compare_means(aes(group = Sex),label = "p.format",size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="Tissue.Abundance.addp.boxplot.bysex.pdf",plot=p2,
       device='pdf',path=".",width=6,height=3.5)
####split sex
head(exp_melt.final)
co<-c("DeepSkyBlue","Tomato")
p<-ggplot(data=exp_melt.final, aes(x=Sex,y=Abundance))+geom_boxplot(aes(fill=Sex),outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  labs(fill="",y="",title="")+scale_fill_manual(values=co)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p
#调整P值的位置：
p2<-p+facet_wrap(~ Type, scales="free",nrow = 4)+
  xlab("") + ylab("Abundance")+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)),size=4,method="t.test",
                     vjust = 0.5,label.x = 1.2)+
  theme(strip.text.x = element_text(size=10, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")

p2

#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="Tissue.Abundance.addp.boxplot.bysex.split.pdf",plot=p2,
       device='pdf',path=".",width=7.5,height=6)

####split Time
head(exp_melt.final)
#co<-c("DeepSkyBlue","Tomato") scale_fill_manual(values=co)
p<-ggplot(data=exp_melt.final, aes(x=Time,y=Abundance))+
  geom_boxplot(aes(fill=Time),outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  labs(fill="",y="",title="")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p
#调整P值的位置：
p2<-p+facet_wrap(~ Type, scales="free",nrow = 4)+
  xlab("") + ylab("Abundance")+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)),size=4,method="anova",
                     vjust = 0.5,label.x = 1.2)+
  theme(strip.text.x = element_text(size=10, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")

p2

#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="Tissue.Abundance.addp.boxplot.bytime.split.pdf",plot=p2,
       device='pdf',path=".",width=8,height=6)

###################################################################
#########bar plot cell
library(RColorBrewer)
display.brewer.all()
mycol<-colorRampPalette(brewer.pal(12,'Paired'),alpha = TRUE)(16)
mycol <- colorRampPalette(brewer.pal(7,'Accent'))(7)

b<-exLR_origin.results.cell$est.ab.sum
head(b);dim(b)

a<-as.data.frame(b)

a$sample<-rownames(a)

head(a)


library("reshape")
library("ggplot2")
data_melt<-melt(a,id=c("sample"))
head(data_melt)
colnames(data_melt)<-c("sample","Type","Abundance")
head(data_melt)

c<-data_melt %>% group_by(sample) %>% mutate(RA=Abundance/sum(Abundance)*100)
head(c)


p<-ggplot(data=c,aes(x=sample,y=RA,fill=Type))+
  geom_bar(stat="identity")+theme_bw()
p
p.final<-p+
  labs(y="Relative Abundance(%)",x="",title="")+
  scale_fill_manual(values=c(mycol))+
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1,
                                   vjust =1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=13,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=13,colour="black"),
        axis.title.y = element_text(size=13,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  scale_x_discrete(limit = unique(a$sample))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.final

ggsave(filename="barplot_EV_origin_cell_AB.pdf",plot=p.final,
       device='pdf',path=".",width=8.3,height=5)

###丰度绝对值
head(c)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(c, x = "Type", y = "Abundance",
               palette =c("lightblue"),
               outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  xlab("") + ylab("Abundance")+
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

ggsave(filename="boxplot_EV_origin_cell_Aboundance.pdf",plot=p,
       device='pdf',path=".",width=6,height=3.5)

###by sex
phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
head(c)

exp_melt.final<-merge(c,phenotype,by = "sample",all=T)
dim(exp_melt.final);head(exp_melt.final)

head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")

p <- ggboxplot(exp_melt.final, x = "Type", y = "Abundance",
               color = "Sex", palette = "nejm",
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
  labs(fill="Sex",x="",y="Abundance",title="") 
#p
# Use only p.format as label. Remove method name.
p2<-p + stat_compare_means(aes(group = Sex),label = "p.format",size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="cell.Abundance.addp.boxplot.bysex.pdf",plot=p2,
       device='pdf',path=".",width=6,height=3.5)
####split sex
head(exp_melt.final)
co<-c("DeepSkyBlue","Tomato")
p<-ggplot(data=exp_melt.final, aes(x=Sex,y=Abundance))+geom_boxplot(aes(fill=Sex),outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  labs(fill="",y="",title="")+scale_fill_manual(values=co)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p
#调整P值的位置：
p2<-p+facet_wrap(~ Type, scales="free",nrow = 4)+
  xlab("") + ylab("Abundance")+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)),size=4,method="t.test",
                     vjust = 0.5,label.x = 1.2)+
  theme(strip.text.x = element_text(size=10, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")

p2

#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="Cell.Abundance.addp.boxplot.bysex.split.pdf",plot=p2,
       device='pdf',path=".",width=5,height=6)

####split Time
head(exp_melt.final)
#co<-c("DeepSkyBlue","Tomato") scale_fill_manual(values=co)
p<-ggplot(data=exp_melt.final, aes(x=Time,y=Abundance))+
  geom_boxplot(aes(fill=Time),outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  labs(fill="",y="",title="")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p
#调整P值的位置：
p2<-p+facet_wrap(~ Type, scales="free",nrow = 4)+
  xlab("") + ylab("Abundance")+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)),size=4,method="anova",
                     vjust = 0.5,label.x = 1.2)+
  theme(strip.text.x = element_text(size=10, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")

p2

#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="Cell.Abundance.addp.boxplot.bytime.split.pdf",plot=p2,
       device='pdf',path=".",width=5,height=6)
