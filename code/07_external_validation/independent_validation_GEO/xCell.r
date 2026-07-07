#install.packages("installr")
#library(installr)
#updateR() 

#install.pakages("devtools")
#install.packages('quadprog')
#BiocManager::install(c("GSVA","GSEABase"), version = "3.8")
#??װxCell
#devtools::install_github('dviraran/xCell')

#载入需要的R包
library(tidyverse)
library(ggplot2)
library(reshape2)
library(corrplot)
library(xCell)
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\独立验证集GEO\\GSE72922")
getwd()
list.files()

expr<-read.delim("./geneMatrix.healthy.txt",header=T,row.names = 1)
#expr<-read.delim("./exp.tpm.for_tissue_origin.txt",header=T)
head(expr)



#去掉含有NA的行
exp.tpm.1 =na.omit(exp.tpm)


#读取基因表达矩阵
#exp<-read.table("TCGA-KIRC-expr.txt",header=T,row.names=1,sep="\t")
#xCell计算免疫细胞矩阵
xcell<-xCellAnalysis(expr,          rnaseq=TRUE)
#保存xCell计算结果
write.table(xcell,file="xCell.txt",row.names=T,col.names=T,sep="\t",quote=F)
save(xcell,file="xcell.RData")
#heatmap for reflect the relative proportions of different cell types using xCell
str(as.data.frame(xcell))
xcell<-as.data.frame(xcell)
#xcell$celltype<-rownames(xcell)

library("ComplexHeatmap")
library(circlize)
mycol <- colorRamp2(c(0, 0.15, 0.3), c("blue", "white", "red"))

pdf(file = "./heatmap.xcell.Relative.aboundance.noscaled.pdf", 5.5, 3)
mycol <- colorRamp2(c(0, 0.05, 0.3), c("SteelBlue", "white", "red"))
Heatmap(as.matrix(xcell),cluster_columns = FALSE,cluster_rows = FALSE,
        col=mycol,name = "Relative aboundance", 
        row_dend_width = unit(8, "mm"),
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = TRUE)
dev.off()
######################################################################################3
head(phenotype)
#info$stim<-factor(info$stim,levels=c("NC","D1.5","D7"))
rownames(phenotype)<-phenotype$sample


annot_df <- data.frame(Sex = phenotype.1$Sex, Time = phenotype.1$Time)
rownames(annot_df)<-rownames(phenotype.1)
col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
           Time = c("D1" = "Firebrick1","D15" = "RoyalBlue1",
                    "D29" = "orange","D43" = "LimeGreen") )

ha <- HeatmapAnnotation(Sex = phenotype$Sex, Time = phenotype$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))
#
pdf(file = "./heatmap.xcell.Relative.aboundance.addlegend.pdf", 7, 9)
mycol <- colorRamp2(c(0, 0.02, 0.1), c("SteelBlue", "white", "red"))
Heatmap(as.matrix(xcell),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "Relative aboundance", 
        row_dend_width = unit(8, "mm"),
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = TRUE)
dev.off()

Heatmap(as.matrix(xcell),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "Relative aboundance", 
        row_dend_width = unit(8, "mm"),
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = TRUE)

Heatmap(as.matrix(xcell),cluster_columns = FALSE,cluster_rows = TRUE,
        col=mycol,name = "Relative aboundance", 
        row_dend_width = unit(8, "mm"),column_split = phenotype$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = TRUE,top_annotation = ha)

exp<-test1
#目标基因集
head(gene335)
#genelist<-c("IGF2","CRABP1","PCP4","DEFB1","HOXB8","MMP7","LYPD2","AOC1","RAB4B","PFN2","PIGV","RAD23A")
genelist<-gene335$SYMBOL
#提取基因集的表达矩阵
goal_exp<-filter(exp,rownames(exp) %in% genelist)
#合并目标基因集表达矩阵和免疫细胞矩阵
combine<-rbind(goal_exp,xcell)
#计算相关系数
comcor<-cor(t(combine))
#计算显著性差异
comp<-cor.mtest(comcor,conf.level=0.95)
pval<-comp$p
#获取目标基因相关性矩阵
t<-as.data.frame(comcor)
goalcor<-select(as.data.frame(comcor),rownames(t)[1:336])%>%rownames_to_column(var="celltype")
goalcor<-filter(goalcor,!(celltype %in% genelist))
##长宽数据转换
goalcor<-melt(goalcor,id.vars="celltype")
colnames(goalcor)<-c("celltype","Gene","correlation")
head(goalcor)
#获取目标基因集pvalue矩阵
genelist<-rownames(t)[1:336]
pval<-select(as.data.frame(pval),genelist)%>%rownames_to_column(var="celltype")
pval<-filter(pval,!(celltype %in% genelist))
#长宽数据转换
pval<-melt(pval,id.vars="celltype")
colnames(pval)<-c("celltype","gene","pvalue")
#将pvalue和correlation两个文件合并
final<-left_join(goalcor,pval,by=c("celltype"="celltype","Gene"="gene"))
head(final)
####绘图
#添加一列,来判断pvalue值范围
final$sign<-case_when(final$pvalue<0.05 &final$pvalue>0.01 &abs(final$correlation)>=0.6 ~"*", 
                      final$pvalue<0.01 &final$pvalue>0.001 &abs(final$correlation)>=0.6 ~"**",
                      final$pvalue<0.001 &abs(final$correlation)>=0.6 ~"***",
                      final$pvalue>0.05 ~"")
final.1<-final[which(abs(final$correlation)>=0.6),]
ggplot(data=final.1,aes(x=Gene,y=celltype))+
  geom_tile(aes(fill=correlation),colour="white",size=1)+
  scale_fill_gradient2(low="#2b8cbe",mid="white",high="#e41a1c")+
  geom_text(aes(label=sign),colour="black")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  guides(fill=guide_legend(title="* p<0.05\n\n** p<0.01\n\n*** p<0.001\n\ncorrelation"))
ggsave("xcell.335gene.correlation.pdf",width=20,height=7)
###############################################################################
###12个基因和一个lnc与xcell细胞相对含量的相关性
genelist<-c("BEAN1","CROCC","DLX1","EPN1","GSC2","GUCY2D","HES3","MID1IP1","MPDZ",
            "SOGA1","VEGFA",
            "ZNF623","LINC01814")
#提取基因集的表达矩阵
goal_exp<-filter(exp,rownames(exp) %in% genelist)
#合并目标基因集表达矩阵和免疫细胞矩阵
combine<-rbind(goal_exp,xcell)
#计算相关系数
comcor<-cor(t(combine))
#计算显著性差异
comp<-cor.mtest(comcor,conf.level=0.95)
pval<-comp$p
#获取目标基因相关性矩阵
t<-as.data.frame(comcor)
goalcor<-select(as.data.frame(comcor),genelist)%>%rownames_to_column(var="celltype")
goalcor<-filter(goalcor,!(celltype %in% genelist))
##长宽数据转换
goalcor<-melt(goalcor,id.vars="celltype")
colnames(goalcor)<-c("celltype","Gene","correlation")
head(goalcor)
#获取目标基因集pvalue矩阵
pval<-select(as.data.frame(pval),genelist)%>%rownames_to_column(var="celltype")
pval<-filter(pval,!(celltype %in% genelist))
#长宽数据转换
pval<-melt(pval,id.vars="celltype")
colnames(pval)<-c("celltype","gene","pvalue")
#将pvalue和correlation两个文件合并
final<-left_join(goalcor,pval,by=c("celltype"="celltype","Gene"="gene"))
head(final)
####绘图
#添加一列,来判断pvalue值范围
final$sign<-case_when(final$pvalue<0.05 &final$pvalue>0.01 &abs(final$correlation)>=0.6 ~"*", 
                      final$pvalue<0.01 &final$pvalue>0.001 &abs(final$correlation)>=0.6 ~"**",
                      final$pvalue<0.001 &abs(final$correlation)>=0.6 ~"***",
                      final$pvalue>0.05 ~"")
final$sign<-case_when(final$pvalue<0.05 &final$pvalue>0.01 &abs(final$correlation)>=0.5 ~"*", 
                      final$pvalue<0.01 &final$pvalue>0.001 &abs(final$correlation)>=0.5 ~"**",
                      final$pvalue<0.001 &abs(final$correlation)>=0.5 ~"***",
                      final$pvalue>0.05 ~"")
final.1<-final[which(abs(final$correlation)>=0.5),]
ggplot(data=final.1,aes(x=Gene,y=celltype))+
  geom_tile(aes(fill=correlation),colour="white",size=1)+
  scale_fill_gradient2(low="#2b8cbe",mid="white",high="#e41a1c")+
  geom_text(aes(label=sign),colour="black")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  guides(fill=guide_legend(title="* p<0.05\n\n** p<0.01\n\n*** p<0.001\n\ncorrelation"))
ggsave("xcell.12gene.1lnc.correlation.pdf",width=6.5,height=6)
