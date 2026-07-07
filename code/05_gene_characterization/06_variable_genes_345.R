# =============================================================================
# 06_variable_genes_345
# Description: 345 Highly Variable Genes Analysis
# =============================================================================

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
data_sd<-data_melt.final %>% group_by(geneid,Donor) %>% 
  dplyr::summarise(sd = sd(count))
head(data_sd)

data_mean<-data_melt.final %>% group_by(geneid,Donor) %>% 
  dplyr::summarise(mean = mean(count))
head(data_mean)

data_mean<-cbind(data_mean,expr[data_mean$geneid,c("gene_name","gene_biotype")])
head(data_mean);dim(data_mean)
#########################################################33
m1<-read.delim("data_mean.add_sd.TPM.donors.txt",header=T,row.names = NULL)
head(m1)
Donor4<-subset(m1,Donor=="4")
Donor5<-subset(m1,Donor=="5")
Donor6<-subset(m1,Donor=="6")
Donor11<-subset(m1,Donor=="11")
Donor12<-subset(m1,Donor=="12")
Donor13<-subset(m1,Donor=="13")


summary(Donor4$sd)
summary(Donor5$sd)
summary(Donor6$sd)
summary(Donor11$sd)
summary(Donor12$sd)
summary(Donor13$sd)

#threshold_sd <- 2  # 设置标准差的阈值
selected_genes.Donor4 <- Donor4[Donor4$sd > 10, ]
selected_genes.Donor4 <- selected_genes.Donor4[selected_genes.Donor4$mean > 1, ]

selected_genes.Donor5 <- Donor5[Donor5$sd > 10, ]
selected_genes.Donor5 <- selected_genes.Donor5[selected_genes.Donor5$mean > 1, ]

selected_genes.Donor6 <- Donor6[Donor6$sd > 10, ]
selected_genes.Donor6 <- selected_genes.Donor6[selected_genes.Donor6$mean > 1, ]

selected_genes.Donor11 <- Donor11[Donor11$sd > 10, ]
selected_genes.Donor11 <- selected_genes.Donor11[selected_genes.Donor11$mean > 1, ]

selected_genes.Donor12 <- Donor12[Donor12$sd > 10, ]
selected_genes.Donor12 <- selected_genes.Donor12[selected_genes.Donor12$mean > 1, ]

selected_genes.Donor13 <- Donor13[Donor13$sd > 10, ]
selected_genes.Donor13 <- selected_genes.Donor13[selected_genes.Donor13$mean > 1, ]

#selected_genes.D43 <-  D43[D43$sd < summary(D43$sd)[4], ]
#selected_genes.D43 <- selected_genes.D43[selected_genes.D43$mean > 1, ]

dim(selected_genes.Donor4);dim(selected_genes.Donor5);
dim(selected_genes.Donor6);dim(selected_genes.Donor11)
dim(selected_genes.Donor12);dim(selected_genes.Donor13);

#write.table(selected_genes.Donor4,
 #           "selected_genes.Donor4.txt", 
  #          sep="\t", col.names=T,quote = F,row.names = T)

####venn
library(venn)
library(VennDiagram)

venn_list <- list(selected_genes.Donor4$geneid, selected_genes.Donor5$geneid, 
                  selected_genes.Donor6$geneid, selected_genes.Donor11$geneid, 
                  selected_genes.Donor12$geneid, selected_genes.Donor13$geneid)   # 制作韦恩图搜所需要的列表文件
names(venn_list) <- c("Donor4","Donor5","Donor6","Donor11","Donor12","Donor13")    # 把列名赋值给列表的key值
venn_list = purrr::map(venn_list,na.omit)      # 删除列表中每个向量中的NA

#作图
venn(venn_list,
     zcolor='style', # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
     opacity = 0.3,  # 调整颜色透明度
     box = F,        # 是否添加边框
     ilcs = 1,     # 数字大小
     sncs = 1.2,        # 组名字体大小
     ilabels = "counts"
)

# 更多参数 ?venn查看

# 查看交集详情,并导出结果
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )
write.table(inter, "venn.result.csv", row.names = FALSE, sep = ',', quote = FALSE)

#######表达波动较大基因6个数据集取交集
genes.2888 <- Reduce(intersect,list(
  selected_genes.Donor4$geneid, selected_genes.Donor5$geneid, 
  selected_genes.Donor6$geneid, selected_genes.Donor11$geneid, 
  selected_genes.Donor12$geneid, selected_genes.Donor13$geneid))
length(genes.2888)
#########################热图
####protein coding
gene.345<-read.delim("345gene.common.inuEV.protein.txt",header=T,row.names = NULL)

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library("ggplot2")
data(geneList, package="DOSE")

#keygene=df.c[which(df.c$annotation1=="Promoter (<=1kb)"),]

Gene_list<- bitr(gene.345$Symbol,fromType="SYMBOL",toType="ENSEMBL",OrgDb="org.Hs.eg.db")

phenotype.1
df1<-df.pro[Gene_list$ENSEMBL,phenotype.1$sample]
#df1$geneid<-rownames(df1)
dim(df1)

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
data.m<-na.omit(data.m)
dim(data.m)
library("ComplexHeatmap")
library(circlize)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)

pdf(file = "./heatmap.pro.345gene.expTPM.noscaled.pdf", 5.5, 5)
final_data<-na.omit(final_data)
mycol <- colorRamp2(c(0, 100, 200), c("blue", "white", "red"))
Heatmap(as.matrix(final_data),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()

pdf(file = "./heatmap.pro.345gene.expTPM.scaled.pdf", 5.5, 5)
mycol <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()

####lnc   2888个基因中没有lnc
df1<-df.lnc[genes.1,]
df1$geneid<-rownames(df1)
dim(df)
df2<-df.lnc[genes.1,]
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
        row_dend_width = unit(8, "mm"),column_split = phenotype$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)

pdf(file = "./heatmap.lnc.2888gene.expTPM.scaled.pdf", 6.5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()
##############################################
######box plot
#data1<-df[apply(df,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(data1)

data1<-df.pro[Gene_list$ENSEMBL,]
#log(TPM+1)
data1<-log2(data1+1)
normalized_data<-as.data.frame(data1)
normalized_data$geneid<-rownames(normalized_data)
normalized_data

library("reshape")
library(dplyr)
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","tpm")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype,by = "sample",all=T)
head(exp_melt.final)
exp_melt.final.1<-exp_melt.final %>% subset(Donor %in% c("4","5","6","11","12","13")) 
dim(exp_melt.final);
dim(exp_melt.final.1)

head(exp_melt.final.1)

library("ggpubr")
library("reshape")
library("ggplot2")

####by time
head(exp_melt.final.1)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt.final.1, x = "Donor", y = "tpm",
               fill="Donor",
               palette =c("jco"),
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
ggsave(filename="./boxplot.bydonor.2888gene.log2TPM.pdf",plot=p,
       device='pdf',path=".",width=6,height=4.5)

################add p
p <- ggboxplot(exp_melt.final.1, x = "Donor", y = "tpm",
               color ="Donor",
               palette = "nejm",show.legend = FALSE)+
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
                     ref.group = "4")                    # Pairwise comparison against reference
p1
ggsave(filename="./boxplot.bydonor.345gene.log2TPM.addP.4.pdf",plot=p1,
       device='pdf',path=".",width=6,height=3.6)

p1<-p+stat_compare_means(method = "anova", label.y = 18)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "12")                    # Pairwise comparison against reference
p1
ggsave(filename="./boxplot.bydonor.345gene.log2TPM.addP.12.pdf",plot=p1,
       device='pdf',path=".",width=6,height=3.6)


my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p +
  stat_compare_means(comparisons = my_comparisons,label = "p.format",
                     size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2
ggsave(filename="./boxplot.bytime.2gene.log2TPM.pdf",plot=p2,
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
######################################################333
####density plot protein coding gene
head(exp_melt.final.1)
head(exp_melt.final)
exp_melt.final.1$log10tmp<-log10(exp_melt.final.1$tpm)
test<-ggplot(exp_melt.final.1, aes(x=tpm, color=Donor)) +
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
ggsave(filename="density.plot.2888gene.tpm.6donor.pdf",plot=test,
       device='pdf',path=".",width=4.5,height=2.8)
###by time

head(exp_melt.final.1)
#exp_melt.final$log10tmp<-log10(exp_melt.final$tpm)
head(exp_melt.final.1)
exp_melt.final.1$Donor<-as.factor(exp_melt.final.1$Donor)
test2<-ggplot(exp_melt.final.1, aes(x=tpm, color=Donor)) +
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
  facet_grid( ~ exp_melt.final.1$Time, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title="Protein coding")
test2
ggsave(filename="density.plot.2888gene.byTime.plot.tpm.pdf",plot=test2,
       device='pdf',path=".",width=9,height=2.8)

###by donor

head(exp_melt.final.1)
#exp_melt.final$log10tmp<-log10(exp_melt.final$tpm)
head(exp_melt.final.1)
exp_melt.final.1$Donor<-as.factor(exp_melt.final.1$Donor)
test3<-ggplot(exp_melt.final.1, aes(x=tpm, color=Time)) +
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
  facet_grid( ~ exp_melt.final.1$Donor, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill=c("white")))+
  theme(strip.text = element_text(size = 12,face = 'bold',colour = "gray2"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(title="Protein coding")
test3
ggsave(filename="density.plot.345gene.byDonor.plot.tpm.pdf",plot=test3,
       device='pdf',path=".",width=9,height=2.3)


###VennDiagram最多5个圈
library(VennDiagram)
library(grid)
T<-venn.diagram(list(D1=selected_genes.D1$geneid,D15=selected_genes.D15$geneid,
                     D29=selected_genes.D29$geneid,D43=selected_genes.D43$geneid),
                filename=NULL,scaled = TRUE,cex=1.4,cat.cex=1.5,
                euler.d = TRUE,main="Overlap between DEG EMT gene",
                main.cex=1.8,col = "transparent",
                fill = c("cornflowerblue", "green", "yellow"),
                cat.col = c("darkblue", "darkgreen", "orange"),
                cat.pos = c(0,-2,180))
grid.draw(T)

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

##########heatmap protein_coding 波动较大的基因
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

