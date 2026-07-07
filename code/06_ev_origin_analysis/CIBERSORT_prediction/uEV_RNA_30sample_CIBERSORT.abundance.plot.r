setwd("/public/home/huhui/project/rna_urine_EV/result/STAR/EV_origin/EV-origin-master")
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\自己用CIBERSORT预测-参考矩阵还是EV-origin的")
library(e1071)


####################################################################
load("D:\\项目\\TCGA_PAAD\\TCGA_Gtex\\Download_from_Gtex_database\\gtex.count.exp.Rdata")
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\自己用CIBERSORT预测-参考矩阵还是EV-origin的")
list.files()
results<-read.delim("uEV_RNA_30sample_CIBERSORT-Results.txt",header=T,row.names = 1)
dim(results)
head(results)

results<-results[,c("Adipose.Tissue",	"Bladder",	"Brain",	
                    "Colon",	"Esophagus",	"Heart",	"Kidney",	"Liver",	
                    "Lung",	"Muscle",	"Nerve",	"Pancreas",	"Pituitary",	"Skin",	"Small.Intestine","Stomach")]
rownames(results) <- gsub(".", " ", rownames(results), fixed = TRUE)
head(results)
#write.table(results,file="ExoRbase3.0.EVPs_in_Urine_origins_abs.tissue.txt",
 #     sep='\t',quote=FALSE,row.names=T)
#########bar plot
library(RColorBrewer)
display.brewer.all()
col <- colorRampPalette(brewer.pal(8,'Accent'))(7)
mycol<-colorRampPalette(brewer.pal(12,'Paired'),alpha = TRUE)(16)

mycol<-colorRampPalette(brewer.pal(12,'Set3'),alpha = TRUE)(16)

library(stringi)


b<-results
dim(b)

# 获取行名
row_names <- rownames(b)

# 使用stringi的stri_sort函数进行自然排序
sorted_row_names <- stri_sort(row_names, numeric = TRUE)

# 按照排序后的行名重排数据框（如果只是重新排序行，不改变列）
b_sorted <- b[sorted_row_names, ]

# 查看排序后的行名
rownames(b_sorted)

b<-b_sorted[64:125,]
b<-b_sorted[1:63,]
b<-b_sorted

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
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1,
                                   vjust =1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=13,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=11,colour="black"),
        axis.title.y = element_text(size=13,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  scale_x_discrete(limit = unique(a$sample))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.final

ggsave(filename="barplot_uEV_RNA_30sample_CIBERSORT-Results.tissue.pdf",plot=p.final,
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

b<-results.cell$est.ab.sum
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
