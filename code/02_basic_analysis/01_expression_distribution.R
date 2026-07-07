# =============================================================================
# 02. Expression Distribution Analysis (TPM-based)
# Description: Boxplots, violin plots, density plots of TPM expression
#   stratified by time point, sex, and donor
# Input: all_sample_tpm.txt, sample_info.txt
# Output: Distribution plots (PDF)
# =============================================================================

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
df<-exp.tpm[,c(1:30)]

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
dim(exp_melt.final);dim(exp_melt.final.1)

head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt, x = "sample", y = "tpm",
              palette =c("lightblue"),
              outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  xlab("") + ylab("Exp. (log2(TPM+1))")+
  theme(axis.text.x = 
          element_text(size = 14,angle = 45,
                       hjust = 1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=14,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=12),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=16),legend.key.size = unit(4,'mm'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ylim(0,12)
p
#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="ggplot2.boxplot.all_sample_log2TPM_exp.new.pdf",plot=p,
       device='pdf',path=".",width=10,height=3)

####by sex
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt.final, x = "sample", 
               y = "tpm",
               fill="Sex",
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
ggsave(filename="ggplot2.boxplot.bysex.all_sample_log2TPM_exp.pdf",plot=p,
       device='pdf',path=".",width=10,height=3.5)
################add p
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")

p <- ggboxplot(exp_melt.final, x = "Time", y = "tpm",
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
  labs(fill="type",x="",y="Expression (log2(TPM+1))",title="") 
#p
# Use only p.format as label. Remove method name.
p2<-p + stat_compare_means(aes(group = Sex),label = "p.format",size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="ggplot2.addp.boxplot.bysex.all_sample_log2TPM_exp.pdf",plot=p2,
       device='pdf',path=".",width=6,height=3.5)

####
p<-ggplot(exp_melt.final, aes(x=Sex, y=tpm,color=Sex) )+
  geom_violin(trim=FALSE,position = "dodge",scale="width") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.shape = NA)+ #绘制箱线图
  scale_color_manual(values = c("#E31A1C","#1F78B4"))+
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
#my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p +
  stat_compare_means(aes(group = Sex),label = "p.format",
                     size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2
#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="ggplot2.violinplot.bysex.log2TPM_exp.pdf",plot=p2,
       device='pdf',path=".",width=3,height=3.5)
##paired
head(exp_melt.final)
exp_melt.final %>% group_by(Sex) %>% dplyr::summarise(count = median(tpm))

library(dplyr)
b1<-nocar_ann %>% subset(CD44_type=="CD44+") %>% 
  group_by(sample_type) %>% dplyr::summarise(count = n())

compare_means(len ~ supp, data = ToothGrowth, 
              group.by = "dose", paired = TRUE)
# Box plot facetted by "dose"
p <- ggpaired(ToothGrowth, x = "supp", y = "len",
              color = "supp", palette = "jama", 
              line.color = "gray", line.size = 0.4,
              facet.by = "dose", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", paired = TRUE)

####by time
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt.final, x = "sample", y = "tpm",
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
ggsave(filename="ggplot2.boxplot.bytime.all_sample_log2TPM_exp.pdf",plot=p,
       device='pdf',path=".",width=10,height=3.5)

################add p
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")

library(dplyr)
exp_melt.final.1<-exp_melt.final %>% subset(scaled_count>=-0.3 & scaled_count<=0.2) 
dim(exp_melt.final);dim(exp_melt.final.1)


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
ggsave(filename="ggplot2.violinplot.bytime.log2TPM_exp.pdf",plot=p2,
       device='pdf',path=".",width=4,height=3.5)
#
head(exp_melt.final)
exp_melt.final.1.a <- exp_melt.final %>% subset(Sex=="female")
exp_melt.final.1.b <- exp_melt.final %>% subset(Sex=="male")

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

ggsave(filename="ggplot2.female.boxplot.noP.bytime.log2TPM_exp.pdf",plot=p2,
    device='pdf',path=".",width=5.5,height=3.5)
###################################
#############################################
####同一个人不同时间点比较,TPM
phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
dim(exp_melt.final.1);dim(exp_melt.final)
unique(exp_melt.final$Donor)
exp_melt.final.2 <- exp_melt.final %>% subset(Donor %in% c("4","5","6","11","12","13")) 
#exp_melt.final.2.1<-exp_melt.final.2 %>% subset(scaled_count>=-0.3 & scaled_count<=0.2)
dim(exp_melt.final.2);dim(exp_melt.final.2.1)
#exp_melt.final.1.a <- exp_melt.final.1 %>% subset(Sex=="female")
#exp_melt.final.1.b <- exp_melt.final.1 %>% subset(Sex=="male")

p <- ggboxplot(exp_melt.final.2, x = "Time", y = "tpm",
               color = "Time", palette = "nejm",
               outlier.colour = NA,
               facet.by = "Donor")+
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
                           size=3.5,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="ggplot2.6donor.split.boxplot.bytime.log2TPM_exp.pdf",plot=p2,
       device='pdf',path=".",width=5,height=5.5)
#######################################################################
####同一个人不同时间点比较,raw count log2(X+1)
head(data1);dim(data1)
dat<-log2(data1+1)
library("reshape")
dat$geneid<-rownames(dat)
exp_melt<-melt(dat,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","log2_count")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype,by = "sample",all=T)
#exp_melt.final.1<-exp_melt.final %>% subset(scaled_count>=-0.3 & scaled_count<=0.2) 
head(exp_melt.final)
dim(exp_melt.final);dim(exp_melt.final.1)

#phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
#head(phenotype);dim(phenotype)
#dim(exp_melt.final.1);dim(exp_melt.final)
unique(exp_melt.final$Donor)
exp_melt.final.2 <- exp_melt.final %>% subset(Donor %in% c("4","5","6","11","12","13")) 
#exp_melt.final.2.1<-exp_melt.final.2 %>% subset(scaled_count>=-0.3 & scaled_count<=0.2)
dim(exp_melt.final.2);dim(exp_melt.final.2.1)
#exp_melt.final.1.a <- exp_melt.final.1 %>% subset(Sex=="female")
#exp_melt.final.1.b <- exp_melt.final.1 %>% subset(Sex=="male")

p <- ggboxplot(exp_melt.final.2, x = "Time", y = "log2_count",
               color = "Time", palette = "nejm",
               outlier.colour = NA,
               facet.by = "Donor")+
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
  labs(fill="type",x="",y="Gene expression (log2_count)",title="") 
#p
# Use only p.format as label. Remove method name.
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p + stat_compare_means(comparisons = my_comparisons,label = "p.format",
                           size=3.5,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="ggplot2.6donor.split.boxplot.bytime.log2_count_exp.pdf",plot=p2,
       device='pdf',path=".",width=5,height=5.5)

####violinplot.bytime log2count
p<-ggplot(exp_melt.final, aes(x=Time, y=log2_count,color=Time) )+
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
  labs(fill="",x="",y="Gene expression",title="")
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p +
  stat_compare_means(comparisons = my_comparisons,label = "p.format",
                     size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2
#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="ggplot2.violinplot.bytime.log2_count_exp.pdf",plot=p2,
       device='pdf',path=".",width=4,height=3.5)

library(dplyr)
exp_melt.final %>% group_by(Time) %>% dplyr::summarise(median = median(log2_count))

####log2 count
#
head(exp_melt.final)
exp_melt.final.a <- exp_melt.final %>% subset(Sex=="female")
exp_melt.final.b <- exp_melt.final %>% subset(Sex=="male")

p <- ggboxplot(exp_melt.final.b, x = "Time", y = "log2_count",
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
  labs(fill="type",x="",y="Gene expression (log2 count)",title="") 
#p
# Use only p.format as label. Remove method name.
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p + stat_compare_means(comparisons = my_comparisons,label = "p.format",
                           size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="ggplot2.male.boxplot.bytime.log2_count_exp.pdf",plot=p2,
       device='pdf',path=".",width=5.5,height=3.5)
exp_melt.final.a %>% group_by(Time) %>% 
  dplyr::summarise(median = median(log2_count))
exp_melt.final.b %>% group_by(Time) %>% 
  dplyr::summarise(median = median(log2_count))

###violin sex log2 count 
p<-ggplot(exp_melt.final, aes(x=Sex, y=log2_count,color=Sex) )+
  geom_violin(trim=FALSE,position = "dodge",scale="width") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.shape = NA)+ #绘制箱线图
  scale_color_manual(values = c("#E31A1C","#1F78B4"))+
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
  labs(fill="",x="",y="Gene expression",title="")
#my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p +
  stat_compare_means(aes(group = Sex),label = "p.format",
                     size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2
#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="ggplot2.violinplot.bysex.log2_count_exp.pdf",plot=p2,
       device='pdf',path=".",width=3,height=3.5)
exp_melt.final %>% group_by(Sex) %>% 
  dplyr::summarise(median = median(log2_count))

###
p <- ggboxplot(exp_melt.final, x = "Time", y = "log2_count",
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
  labs(fill="type",x="",y="Gene expression (log2 count)",title="")
#p
# Use only p.format as label. Remove method name.
p2<-p + stat_compare_means(aes(group = Sex),label = "p.format",size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="ggplot2.addp.boxplot.bysex.all_sample_log2_count_exp.pdf",plot=p2,
       device='pdf',path=".",width=6,height=3.5)

exp_melt.final %>% group_by(Time,Sex) %>% 
  dplyr::summarise(median = median(log2_count))

########总图
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt, x = "sample", y = "log2_count",
               palette =c("lightblue"),
               outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  xlab("") + ylab("Gene expression (log2 count)")+
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
ggsave(filename="ggplot2.boxplot.all_sample_log2_count_exp.pdf",plot=p,
       device='pdf',path=".",width=10,height=3.5)

####by sex
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt.final, x = "sample", 
               y = "log2_count",
               fill="Sex",
               palette =c("jco"),
               outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  xlab("") + ylab("Gene expression (log2 count)")+
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
#pdf(file = "ggplot2.boxplot.all_sample_log2_count_exp.pdf",12,5)
ggsave(filename="ggplot2.boxplot.bysex.all_sample_log2_count_exp.pdf",plot=p,
       device='pdf',path=".",width=10,height=3.5)

####by time
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt.final, x = "sample", y = "log2_count",
               fill="Time",
               palette =c("jco"),
               outlier.colour = NA,show.legend = FALSE)+
  theme_bw()+
  xlab("") + ylab("Gene expression (log2 count)")+
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
ggsave(filename="ggplot2.boxplot.bytime.all_sample_log2_count_exp.pdf",plot=p,
       device='pdf',path=".",width=10,height=3.5)

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
p<-ggline(data_median.2, x = "Time", y = "median",
          color = "Donor", palette = "jco")+
  theme_classic()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x="",y="Expression (log2(TPM+1))",title="")
p
p2<-p+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="lineplot.6donor.bytime.log2tpm_exp.pdf",plot=p2,
       device='pdf',path=".",width=5,height=3.5)

####violinplot.bytime median count
p<-ggplot(data_median.2, aes(x=Time, y=median,color=Time) )+
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
  labs(fill="",x="",y="Gene expression (median count)",title="")
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p +
  stat_compare_means(comparisons = my_comparisons,label = "p.format",
                     size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2
#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="ggplot2.violinplot.bytime.median_count_exp.pdf",plot=p2,
       device='pdf',path=".",width=4,height=3.5)

library(dplyr)
data_median.2 %>% group_by(Time) %>% dplyr::summarise(median = median(median))

####median count
#
head(data_median.2)
data_median.2.a <- data_median.2 %>% subset(Sex=="female")
data_median.2.b <- data_median.2 %>% subset(Sex=="male")

p <- ggboxplot(data_median.2, x = "Time", y = "log2_count",
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
  labs(fill="type",x="",y="Gene expression (log2 count)",title="") 
#p
# Use only p.format as label. Remove method name.
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p + stat_compare_means(comparisons = my_comparisons,label = "p.format",
                           size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="ggplot2.male.boxplot.bytime.log2_count_exp.pdf",plot=p2,
       device='pdf',path=".",width=5.5,height=3.5)
exp_melt.final.a %>% group_by(Time) %>% 
  dplyr::summarise(median = median(log2_count))
exp_melt.final.b %>% group_by(Time) %>% 
  dplyr::summarise(median = median(log2_count))

###violin sex log2 count 
p<-ggplot(data_median.2, aes(x=Sex, y=median,color=Sex) )+
  geom_violin(trim=FALSE,position = "dodge",scale="width") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.shape = NA)+ #绘制箱线图
  scale_color_manual(values = c("#E31A1C","#1F78B4"))+
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
  labs(fill="",x="",y="Gene expression (median count)",title="")
#my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p +
  stat_compare_means(aes(group = Sex),label = "p.format",
                     size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2
#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="ggplot2.violinplot.bysex.median_count_exp.pdf",plot=p2,
       device='pdf',path=".",width=3,height=3.5)
data_median.2 %>% group_by(Sex) %>% 
  dplyr::summarise(median = median(median))

###
p <- ggboxplot(data_median.2, x = "Time", y = "median",
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
  labs(fill="type",x="",y="Gene expression (median count)",title="")
#p
# Use only p.format as label. Remove method name.
p2<-p + stat_compare_means(aes(group = Sex),label = "p.format",size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

#pdf(file = "ggplot2.boxplot.all_sample_scale_count_exp.pdf",12,5)
ggsave(filename="ggplot2.addp.boxplot.bysex.all_sample_median_count_exp.pdf",plot=p2,
       device='pdf',path=".",width=6,height=3.5)

#####################################################################
############样本聚类热图
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

library(genefilter)
library(Biobase)
data1<-df[apply(df,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(data1)
#log(TPM+1)
#data1<-log2(data1+1)

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
               main="All samples TPM correlation",
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
file <- paste(".", "/all_sample_cluster.png", sep="")
png(filename=file, height = 2200, width = 2200, res = 300, units = "px")
print(p)
dev.off()

file <- paste(".", "/all_sample_cluster.pdf", sep="")
pdf(file=file, height = 8, width = 8)
print(p)
dev.off()

###############################################################################
#################################
####相关性画图
file <- paste(".", "/all_sample_pearson.cor.xls", sep="")
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


file <- paste(".", "/all_sample_cor.png", sep="")
png(filename=file, height = 3000, width = 3000, res = 300, units = "px")

pairs(log10(data1+1), pch = ".", upper.panel = panel.cor, diag.panel = panel.hist, 
      lower.panel = panel.lm, main = "Pearson Correlation of all samples",gap=0)
dev.off()


file <- paste(".", "/all_sample_cor.pdf", sep="")
pdf(file=file, height = 10, width = 10)
pairs(log10(data1+1), pch = ".", upper.panel = panel.cor, diag.panel = panel.hist, 
      lower.panel = panel.lm, main = "Pearson Correlation of all samples",gap=0)
dev.off()

#################################################################3
####density plot
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
ggsave(filename="density.plot.all_sample.tpm.pdf",plot=test,
       device='pdf',path=".",width=5,height=3.2)
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
ggsave(filename="density.byTime.plot.all_sample.tpm.pdf",plot=test2,
       device='pdf',path=".",width=11,height=3.5)
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
  facet_wrap( ~ exp_melt.final$Donor, drop=TRUE,scale="free",
              ncol=6)+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  
test3
ggsave(filename="density.byDonor.plot.all_sample.tpm.pdf",plot=test3,
       device='pdf',path=".",width=10,height=4)
