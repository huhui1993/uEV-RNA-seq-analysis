# =============================================================================
# 02b_pca_analysis_protein_coding
# Description: PCA Analysis - Protein Coding Genes Only
# =============================================================================

#教程：https://mp.weixin.qq.com/s/ijj2Ww9hEfIXfL-7K2IOWg
list.files()
#expr <- read.delim("exp.for.PCA.xlsx",header = TRUE)  ###行为基因名，列为样本名
library(FactoMineR)
library(factoextra) # 用于提取PCA分析结果信息，其实不用也可以。
library(FactoInvestigate)

expr<-read.delim("all_sample_count.addAnotation.txt",header=T,row.names = 1)
head(expr);dim(expr)
expr<-expr[which(expr$gene_biotype=="protein_coding"),]
expr<-expr[,c(3:32)]
head(expr)
expr[expr==0]<-0.01

type<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",
                 header=T,row.names = 1)

getwd()
#pca分析
expr<-expr[apply(expr,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(expr)
expr<-expr[apply(expr,1,function(x){XXXX<-FALSE;if(any(x>5)){XXXX<-TRUE};return(XXXX)}),];dim(expr)
res.pca <- PCA(t(expr), 
               scale.unit = TRUE, 
               ###如果变量个数高于样本数，则最只能输出样本数-1个主成分(nrow(env)-1)##
               ncp = ncol(expr),  
               graph = FALSE)


#提取分析结果
#factoextra包中有打包函数用于提取FactoMineR包函数分析结果。

# 3.1 提取特征值-每个主成分对样本方差的解释度。
eig.val = get_eigenvalue(res.pca) 
#eig.val = get_eig(res.pca)
eig.val

# 3.2 提取样本结果
ind.res =  get_pca_ind(res.pca) # factoextra提取分析结果，其实可以不用。
ind.res$coord
ind.res$cos2
ind.res$contrib

# 3.3 提取变量结果
var.res <- get_pca_var(res.pca) # 提取响应变量结果矩阵列表
var.res$coord
var.res$cor
var.res$cos2
var.res$contrib

#FactoInvestigate包的Investigate()自动输出FactoMineR包PCA分析的结果解释,输出文件可以是word,html和PDF。

# 3.5 FactoInvestigate包的Investigate()输出结果报告
#?Investigate # 查看函数帮助信息
#Investigate(res.pca,
#            file = "Investigate.Rmd",openFile = FALSE,
 #           document = c("word_document"))


#dimdesc()可用于鉴定与给定主成分具有统计学意义相关性的响应变量。

# 3.6 主成分维度描述
#基于对各主成分的因子分析，找出该主成分最能表征的定性(样本分类信息)和定量变量。
#?dimdesc
dim = dimdesc(res.pca,
              axes = 1:3,
              proba = 0.05) # proba设置显著性阈值
dim$call # 设置信息
dim$Dim.1$quanti #定量响应变量与PC1的相关性r值及p值。
dim$Dim.1$quali #附加定性响应变量与PC1的决定系数R^2及p值。
dim$Dim.2$quali
#4 绘图
#4.1 特征值绘图-碎石图
pdf(file = "fviz_eig.higher5.pdf",6,4)
fviz_eig(res.pca,addlabels = TRUE) # 对每个主成分对方差的解释度绘制柱形图。
#fviz_screeplot(res.pca)# 对每个主成分对方差的解释度绘制柱形图,图中没有解释度值。
nouse = dev.off()

#4.2  样本散点图

#fviz_pca_ind()绘制样本散点图。样本图主要调整内容：1)点的大小、颜色与形状；2)字体(标签与坐标轴等)大小与颜色；3)图主题设置；4)分组置信区间类型。
# 颜色设置
library(ggsci)
library(scales)
mypal = pal_d3("category10")(10)
show_col(mypal)
mypal = mypal[-8]
show_col(mypal)
mypal

library(RColorBrewer)
display.brewer.all()
col <- colorRampPalette(brewer.pal(8,'Accent'))(7)
mycol<-colorRampPalette(brewer.pal(12,'Set3'),alpha = TRUE)(26)

# 4.2.1 绘制特定主成分样本图
ind1 = fviz_pca_ind(res.pca,
                    axes = c(1, 2), #绘制主成分1和2.
                    geom =c("point","text"), # 样本以点+标签形式展示。
                    ##设置样本点属性##
                    pointshape=16, # 点形状
                    pointsize=2,# 点大小
                    col.ind = type$type, # 样本颜色以类别区分
                    palette = mypal, # 样本depth各分组对应颜色。
                    alpha.ind = 1, # 颜色透明度
                    ##设置样本标签属性##
                    repel = TRUE,# 使用ggrepel避免标签重叠，默认为FALSE
                    labelsize=3,# 标签大小，样本标签颜色默认与样本点保持一致。
                    ##不添加置信圈##
                    addEllipses = FALSE,
)
ind1

# 4.2.2 图中保留指定样本
## 可根据样本标签保留样本、或者根据cos2、contrib值大小选择保留样本的数量。
ind2.1 = fviz_pca_ind(res.pca, 
                      geom.ind = c("point","text"), 
                      pointsize =2,
                      habillage=factor(type$type), # 函数点颜色和形状设置，只能使用同一分类因子。
                      #label = "none", #设置不显示样本标签
                      #col.ind = env$depth, # color by groups
                      palette = mypal, # 使用预先设定好的颜色，保证图色系的统一。
                      mean.point = FALSE, #设置不显示分组均点。默认以更大的点显示。
                      select.ind = list(contrib = 20), # 保留contrib最高的20个的样本。
                      legend.position = "bottom",
)
ind2.1
## ggpubr::ggpar()设置图形参数
ind2.2 = ggpubr::ggpar(ind1, # 设置图形参数的对象
                       title = "Principal Component Analysis",
                       subtitle = "group",
                       #caption = "Source: factoextra",
                       xlab = "PC1 (25.1%)", ylab = "PC2 (10.7%)",
                       legend.title = "Type",
                       ggtheme = theme_bw(), 
                       xlim = range(min(ind.res$coord[,1])-1,max(ind.res$coord[,1])+1 ),
                       ylim = range(min(ind.res$coord[,2])-1,max(ind.res$coord[,2])+1 ),
                       #palette = mypal
)
ind2.2

# 4.2.3 设置多边形分组阴影
##没有找到fviz_pca_ind()直接根据分组信息设置点形状的参数，
ind3 = fviz_pca_ind(res.pca, 
                    geom.ind = "point", # 样本点不添加样本标签
                    #pointshape= 16, ##此函数不能根据分组信息设置不同的点形状(反正我没发现，要设置需要改函数)
                    habillage=factor(type$type), #此参数颜色和性状同时设置
                    pointsize=2,# 点大小
                    #col.ind = env$depth, # 此参数只能设置颜色
                    palette = mypal, # 使用预先设定好的颜色，保证图色系的统一。
                    mean.point = TRUE, #默认以更大的点显示分组均值。
                    star.plot = TRUE, # 分组均值点与组内样本使用线段连接。
                    ##绘制多边形分组阴影##
                    addEllipses = TRUE, 
                    ellipse.type = "convex", #可以自行按分组计算凸点，绘制多边形。
                    ##设置主题##
                    ggtheme = theme_bw()
)
ind3

# 将图组合输出到本地
library(ggpubr)
ind = ggarrange(ind2.2,ind3, # 子图，按照想要的顺序排列
                nrow = 1,ncol = 2, # 2行2列的子图矩阵，按行排列
                labels=c("ind2.2","ind3"), #子图标签
                hjust = c(0,0.25,0,0.25),# 默认-0.5，子图标签水平设置，负值为右，正值为左。设置一个值，则每个子图的位置一样，一个向量，则为每个子图单独设置标签位置。
                vjust = c(1.5,1.5,0,0),# 默认1.5,垂直调整子图标签位置，值越大，越向下。
                widths = rep(0.8,4), # 对每个子图的相关宽/高进行设置
                heights = rep(0.8,4)) 
ggsave(filename = "pca.final.higher5.pdf",plot = ind,width = 12,height = 5,device = "pdf",family = "Times")

#####################################################
# 4.2.4 设置分组椭圆区间
## 假设为多元正态分布
ind4 = fviz_pca_ind(res.pca, 
                    geom.ind = "point",  # 设置只显示样本点
                    habillage = factor(type$type),
                    pointsize = 2,
                    palette = mypal, 
                    addEllipses = TRUE, 
                    star.plot = TRUE,
                    ellipse.type = "norm", 
                    ellipse.level = 0.95, # 正态概率下聚集椭圆的大小
                    ellipse.alpha = 0.1,
                    legend.title = "Depth",
                    ggtheme = theme(panel.grid = element_blank(),
                                    panel.background = element_rect(fill = "white"),
                                    axis.line = element_line(color = "black")
                    )
)
ind4

## 多变量t-分布
ind5 = fviz_pca_ind(res.pca, 
                    geom.ind = "point",  # 设置只显示样本点
                    habillage = factor(type$type),
                    pointsize = 2,
                    palette = mypal, 
                    addEllipses = TRUE, 
                    star.plot = TRUE,
                    ellipse.type = "t", # 多变量t-分布
                    ellipse.level = 0.95,
                    ellipse.alpha = 0.1,
                    legend.title = "Depth",
                    ggtheme = theme(panel.grid = element_blank(),
                                    panel.background = element_rect(fill = "white"),
                                    axis.line = element_line(color = "black")
                    )
)
ind5

## 围绕分组均值的置信椭圆
ind6 = fviz_pca_ind(res.pca, 
                    geom.ind = "point",  # 设置只显示样本点
                    habillage = factor(type$type),
                    pointsize = 2,
                    palette = mypal, 
                    addEllipses = TRUE, 
                    star.plot = TRUE,
                    ellipse.type = "confidence", # confidence ellipses arround group mean points
                    ellipse.level = 0.95,
                    ellipse.alpha = 0.1,
                    legend.title = "Depth",
                    ggtheme = theme(panel.grid = element_blank(),
                                    panel.background = element_rect(fill = "white"),
                                    axis.line = element_line(color = "black")
                    )
) 
ind6

## 欧式距离
ind7 = fviz_pca_ind(res.pca, 
                    geom.ind = "point",  # 设置只显示样本点
                    habillage = factor(type$type),
                    pointsize = 2,
                    palette = mypal, 
                    addEllipses = TRUE, 
                    star.plot = TRUE,
                    ellipse.type = "euclid", # 与coord_fixed()联用，绘制一个半径等于分组的圆，代表中心的欧式距离。
                    ellipse.level = 0.95,
                    ellipse.alpha = 0.1,
                    legend.title = "Depth",
                    ggtheme = theme(panel.background = element_rect(fill = "white"),
                                    axis.line = element_line(color = "black")
                    )
)+coord_fixed()
ind7

# 将图组合输出到本地
library(ggpubr)
ind.ellipse = ggarrange(ind4,ind5,ind6,ind7, # 子图，按照想要的顺序排列
                        nrow = 2,ncol = 2, # 2行2列的子图矩阵，按行排列
                        labels=c("ind4","ind5","ind6","ind7"), #子图标签
                        hjust = c(-0.2,0.25,-0.2,0.25),# 默认-0.5，子图标签水平设置，负值为右，正值为左。设置一个值，则每个子图的位置一样，一个向量，则为每个子图单独设置标签位置。
                        vjust = c(1.5,1.5,0,3),# 默认1.5,垂直调整子图标签位置，值越大，越向下。
                        widths = rep(0.8,4), # 对每个子图的相关宽/高进行设置
                        heights = rep(0.8,4)) 
ggsave(filename = "ind.ellipse.higher5.pdf",plot = ind.ellipse,width = 12,height = 10,device = "pdf",family = "Times")

#4.3 绘制特征变量图
#factoextra包中fviz_pca_var()用于绘制变量散点或箭头图。变量可用点、箭头和变量标签显示，以点显示的设置于样本没有差异。这里只描述箭头和标签名绘图过程。
# 4.3.1  默认绘图参数
var1 = fviz_pca_var(res.pca)
var1

# 4.3.2 以变量分类设置变量箭头颜色
var2 = fviz_pca_var(res.pca,
                    #geom.var=c("point", "arrow", "text"),
                    #pointshape=16,
                    #pointsize=2,
                    #geom.var=c("arrow", "text"), # 变量以箭头和文本标签显示
                    labelsize=4,
                    arrowsize=0.75, # 箭头大小
                    repel = FALSE,
                    col.var = type$type, # 变量分类
                    palette = mycol,
                    col.circle = "grey50", # correlation circle颜色
                    circlesize = 0.75, # 设置圆圈线的粗细
                    axes.linetype = "dashed", # 过(0,0)点轴线的线型。
                    legend.title = "Type",
                    ggtheme = theme_pubr(legend = "right"),
)
var2 = ggpubr::ggpar(var2, # 设置图形参数的对象
                     xlab = "PC1 (25.1%)", ylab = "PC2 (10.7%)",
)
var2

# 4.3.3 以cos2值为变量设置梯度颜色
var3 = fviz_pca_var(res.pca,
                    #geom.var=c("point", "arrow", "text"),
                    #pointshape=16,
                    #pointsize=2,
                    #geom.var=c("arrow", "text"), # 变量以箭头和文本标签显示
                    labelsize=2,
                    arrowsize=0.75, # 箭头大小
                    repel = FALSE,
                    select.var=list(contrib = 200),
                    col.var = "cos2", # cos2值越高，箭头颜色越黄。
                    gradient.cols = c(mypal[1],mypal[2]), # 颜色顺序对应从低到高的cos2值。
                    #palette = mypal,
                    col.circle = "grey50", # correlation circle颜色
                    circlesize = 0.75, # 设置圆圈线的粗细
                    axes.linetype = "dashed", # 过(0,0)点轴线的线型。
                    legend.title = "Type",
                    ggtheme = theme_bw()
)
var3

# 4.3.4 以contrib值为变量设置梯度颜色
var4 = fviz_pca_var(res.pca,
                    #geom.var=c("point", "arrow", "text"),
                    #pointshape=16,
                    #pointsize=2,
                    #geom.var=c("arrow", "text"), # 变量以箭头和文本标签显示
                    labelsize=2,
                    arrowsize=0.75, # 箭头大小
                    repel = FALSE,
                    select.var=list(contrib = 200),
                    col.var = "contrib", 
                    gradient.cols = "jco", # 可以直接使用ggsci中的颜色。
                    #palette = mypal,
                    col.circle = "grey50", # correlation circle颜色
                    circlesize = 0.75, # 设置圆圈线的粗细
                    axes.linetype = "dashed", # 过(0,0)点轴线的线型。
                    legend.title = "Type",
                    ggtheme = theme_bw()
)
var4

# 将图组合输出到本地
library(ggpubr)
var = ggarrange(var3,var4, 
                nrow = 1,ncol = 2, 
                labels=c("var3","var4"), #子图标签
                hjust = c(-0.2,0.25,-0.2,0.25),# 默认-0.5，
                vjust = c(1.5,1.5,0,3)) 
ggsave(filename = "var.higher5.pdf",plot = var,width = 14,height = 6,device = "pdf",family = "Times")

#4.4 绘制样本-特征变量双序图

#fviz_pca_biplot()用于绘制PCA双序图。

# 4.4.1  绘制双序图--分类变量设置箭头颜色
p1 = fviz_pca_biplot(res.pca, 
                     repel = TRUE,
                     addEllipses = TRUE, 
                     label = "var", # 只显示标量标签
                     labelsize=4, # 标签大小
                     #col.var = "cos2",
                     col.var=factor(colnames(expr)[1:20]),
                     fill.ind = type$type,
                     col.ind = "white",
                     pointshape=21,
                     pointsize = 2,
                     mean.point=FALSE,
                     ggtheme = theme_bw(),
                     legend.title = list(fill = "Depth", color = "Type"),
)+
  ggpubr::fill_palette("d3")+ # 样本及置信圈填充颜色
  ggpubr::color_palette("npg") # 响应变量分类颜色 
p1

# 4.4.2  绘制双序图-连续变量设置箭头颜色梯度
p2 = fviz_pca_biplot(res.pca,
                     repel = FALSE,
                     addEllipses = TRUE, 
                     label = "var",
                     #habillage = factor(env$depth),
                     fill.ind = type$type,
                     col.ind = "white",
                     pointshape=21,
                     pointsize = 2,
                     palette = "jco",
                     mean.point=FALSE,
                     col.var = "contrib",
                     alpha.var ="contrib",
                     ggtheme = theme_bw(),
                     gradient.cols = mypal[1:3], # 响应变量定量分类
                     #col.var=factor(data[1,4:14]),
                     legend.title = list(fill = "type", color = "Contrib",alpha="Contrib"),
)
p2

# 将图组合输出到本地
library(ggpubr)
biplot = ggarrange(p1,p2,
                   ncol = 2, 
                   labels=c("p1","p2")#子图标签
) 
ggsave(filename = "biplot.pdf",plot = p2,width = 6,height = 7,device = "pdf",family = "Times")

# R环境
sessionInfo()
############################################
#prcomp()函数
expr<-read.delim("all_sample_count.addAnotation.txt",header=T,row.names = 1)
head(expr);dim(expr)
expr<-expr[which(expr$gene_biotype=="protein_coding"),]
dim(expr)
expr<-expr[,c(3:32)]
head(expr)
expr[expr==0]<-0.01

type<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",
                 header=T,row.names = 1)
getwd()
#pca分析
#expr<-expr[apply(expr,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(expr)
expr<-expr[apply(expr,1,function(x){XXXX<-FALSE;if(any(x>5)){XXXX<-TRUE};return(XXXX)}),];dim(expr)
library(ggplot2)
library(ggrepel)
iris_input <- t(expr) # 使用R自带iris数据集（150*5，行为样本，列为特征）
rownames(iris_input) 
#<- paste("sample",1:nrow(iris_input),sep = "")# 设置样本名
#head(iris_input) # 查看数据集前几行
iris_input[1:4,1:5]
# 进行PCA（scale. = TRUE表示分析前对数据进行归一化）
#pca1 <- prcomp(iris_input[,-ncol(iris_input)],center = TRUE,scale. = TRUE)

pca1 <- prcomp(iris_input,center = TRUE,scale. = TRUE)

df1 <- pca1$x # 提取PC score
df1 <- as.data.frame(df1) # 注意：如果不转成数据框形式后续绘图时会报错
head(df1)

summ1 <- summary(pca1)
summ1

# 提取主成分的方差贡献率,生成坐标轴标题
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

# 绘制PCA得分图
library(ggplot2)
p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color = type$type))+
  stat_ellipse(aes(fill = type$type),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 3.5)+geom_text_repel(aes(label=rownames(type)),size = 3,
                                         box.padding = unit(0.1,"lines"),point.padding = unit(0.1,"lines"),
                                         segment.color = 'black',show.legend = FALSE,color = 'black')+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","pink","LightSkyBlue","Turquoise"))+
  scale_colour_manual(values = c("purple","orange","pink","LightSkyBlue","Turquoise"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca1
ggsave(filename = "PCA.higher5.pdf",plot = p.pca1,width = 8,height = 7,device = "pdf",family = "Times")

# 绘制PCA得分图
library(ggplot2)
p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2))+
  geom_point(size = 3.5)+geom_text_repel(aes(label=rownames(type)),size = 3,
                                         box.padding = unit(0.1,"lines"),point.padding = unit(0.1,"lines"),
                                         segment.color = 'black',show.legend = FALSE,color = 'black')+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","pink","LightSkyBlue","Turquoise"))+
  scale_colour_manual(values = c("purple","orange","pink","LightSkyBlue","Turquoise"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca1
ggsave(filename = "PCA.higher5.nogroup.protein_coding.pdf",plot = p.pca1,width = 6,height = 6,device = "pdf",family = "Times")
##########################################
library(pheatmap)
expr[1:3,1:4]
pdf(file = "heatmap.higher5.protein_coding.pdf",8,8)
pheatmap(expr,cellwidth=NA,cellheight=NA,scale="row",show_rownames = FALSE,
         legend_labels = "scale exp", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         border_color=NA,fontsize=13,cluster_rows = TRUE,treeheight_col = 20)
dev.off()
########################################
#protein_coding+lnc TPM
expr<-read.delim("all_sample_tpm.addsymbol.txt",header=T,row.names = 1)
head(expr);dim(expr)
expr<-subset(expr,gene_biotype %in% c("protein_coding","lncRNA"))
#expr<-expr[which(expr$gene_biotype=="protein_coding"),]
dim(expr)
expr<-expr[,c(1:30)]
head(expr);dim(expr)
expr[expr==0]<-0.01

type<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",
                 header=T,row.names = 1)
getwd()
#pca分析
#expr<-expr[apply(expr,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(expr)
expr<-expr[apply(expr,1,function(x){XXXX<-FALSE;if(any(x>5)){XXXX<-TRUE};return(XXXX)}),];dim(expr)
library(ggplot2)
library(ggrepel)
iris_input <- t(expr) # 使用R自带iris数据集（150*5，行为样本，列为特征）
rownames(iris_input) 
#<- paste("sample",1:nrow(iris_input),sep = "")# 设置样本名
#head(iris_input) # 查看数据集前几行
iris_input[1:4,1:5]
# 进行PCA（scale. = TRUE表示分析前对数据进行归一化）
#pca1 <- prcomp(iris_input[,-ncol(iris_input)],center = TRUE,scale. = TRUE)

pca1 <- prcomp(iris_input,center = TRUE,scale. = TRUE)

df1 <- pca1$x # 提取PC score
df1 <- as.data.frame(df1) # 注意：如果不转成数据框形式后续绘图时会报错
head(df1)

summ1 <- summary(pca1)
summ1

# 提取主成分的方差贡献率,生成坐标轴标题
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

# 绘制PCA得分图
library(ggplot2)
p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color = type$type))+
  stat_ellipse(aes(fill = type$type),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 3.5)+geom_text_repel(aes(label=rownames(type)),size = 3,
                                         box.padding = unit(0.1,"lines"),point.padding = unit(0.1,"lines"),
                                         segment.color = 'black',show.legend = FALSE,color = 'black')+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","pink","LightSkyBlue","Turquoise"))+
  scale_colour_manual(values = c("purple","orange","pink","LightSkyBlue","Turquoise"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca1
ggsave(filename = "tmpPCA.higher5.pdf",plot = p.pca1,width = 8,height = 7,device = "pdf",family = "Times")

# 绘制PCA得分图
library(ggplot2)
p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2))+
  geom_point(size = 3.5)+geom_text_repel(aes(label=rownames(type)),size = 3,
                                         box.padding = unit(0.1,"lines"),point.padding = unit(0.1,"lines"),
                                         segment.color = 'black',show.legend = FALSE,color = 'black')+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","pink","LightSkyBlue","Turquoise"))+
  scale_colour_manual(values = c("purple","orange","pink","LightSkyBlue","Turquoise"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca1
ggsave(filename = "PCA.higher5.nogroup.protein_coding.lnc.pdf",
       plot = p.pca1,width = 6,height = 6,device = "pdf",family = "Times")
##########################################
library(pheatmap)
expr[1:3,1:4]
pdf(file = "heatmap.higher5.protein_coding.lnc.pdf",8,8)
pheatmap(expr,cellwidth=NA,cellheight=NA,scale="row",show_rownames = FALSE,
         legend_labels = "scale exp", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         border_color=NA,fontsize=13,cluster_rows = TRUE,treeheight_col = 20)
dev.off()

write.table(expr,file="all_sample_count.protein_lnc.higher5.txt",
            quote = F,sep="\t",row.names=T)

#complexheatmap
head(expr)

#final_data<-apply(exp.100,2,as.numeric)
#rownames(final_data)<-rownames(exp.100)
final_data<-expr
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

head(phenotype)
annot_df<-data.frame(type=phenotype$Time)
#info$stim<-factor(info$stim,levels=c("NC","D1.5","D7"))
rownames(phenotype)<-phenotype$sample
ha <- HeatmapAnnotation(Sex = phenotype$Sex, Time = phenotype$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "./pca.new/heatmap.allgene.expTPM.scaled.addlegend.pdf", 7.8, 7.5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()
###########################################
expr[expr==0]<-0.01
#expr<-expr[apply(expr[,c(2:5)],1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(expr)

options(width = 160)
expr.dim <- dim(expr)
expr.dim
expr[1:3,1:4]
N_samp <- expr.dim[2] - 1
N_gene <- expr.dim[1]
expr.pca<- princomp(expr[,-1],cor = TRUE)


pdf(file = "./rna.scatter.pdf", 5.5, 5.5)
biplot(expr.pca, cex = 0.5, pc.biplot = FALSE)
title(main = "Genes Scatter Diagram")
nouse = dev.off()

library(plotrix)
pc.var <- expr.pca$sdev ** 2
load <- loadings(expr.pca)
pc <- data.frame(load[,1:3])
rnam <- rownames(pc)

pdf(file = "./rna.3D.pdf", 7, 7)
library(scatterplot3d)
pc.var.scaled <- 100 * pc.var/sum(pc.var)
pc.var.scaled = sprintf("%.2f", pc.var.scaled)
axis.names <- paste(names(pc.var), pc.var.scaled, sep = "(")
axis.names <- paste(axis.names, "%)", sep = "")
s3d <- scatterplot3d(pc, type = "h", angle = 40, highlight.3d = FALSE, 
                     scale.y = 0.7, pch = 20, main = "PCA 3D figure", 
                     las = 1, color = rainbow(N_samp), xlab = axis.names[1], 
                     ylab = axis.names[2], zlab = axis.names[3])
text(s3d$xyz.convert(pc), labels = 1:N_samp, cex = 1.0)
legend("right", pch = 20, col = rainbow(N_samp), cex = 0.8, legend = rnam)
nouse = dev.off()

pdf(file = "./rna.pca_figure.pdf", 7.5, 7.5)
library(plotrix)
pc.var <- expr.pca$sdev ** 2

whichgap <- which.max(pc.var[-N_samp] - pc.var[-1])
gap.size <- 0.9 * (pc.var[whichgap] - pc.var[whichgap + 1])
gaps <- c(pc.var[whichgap + 1] + 0.05 * gap.size, pc.var[whichgap] - 0.05 * gap.size)
col <- color.gradient(c(0, 1), c(0, 1, 0), c(1, 0), N_samp)
gap.barplot(pc.var, gap=gaps, col = col, main = "Variances of Principal Components", 
            ylim = c(0,expr.pca$sdev[1] ** 2 - gap.size), 
            xlab = "Components",
            ylab = "Variances", 
            xlim = c(0,N_samp + 1), 
            ytics=as.numeric(sprintf("%.3f", c(min(pc.var),max(pc.var),pc.var[whichgap+1],pc.var[whichgap]))))

box(lty = "solid", col = 'black')

load <- loadings(expr.pca)
pc <- data.frame(load[,1:3])
rnam <- rownames(pc)

plot(pc[,-3], pch = 20,cex=2.0, col = rainbow(N_samp), 
     main = "loadings for PC1 and PC2 in samples",
     xlab = axis.names[1],ylab = axis.names[2])
text(pc[,1],pc[,2], cex = 1, adj = c(-0.2,0.5))
legend("right", pch = 20, col = rainbow(N_samp),cex = 0.8,legend = rnam)

plot(pc[,-2], pch = 20,cex=2.0, col = rainbow(N_samp), 
     main = "loadings for PC1 and PC3 in samples",
     xlab = axis.names[1],ylab = axis.names[3])
text(pc[,1],pc[,3], cex = 1, adj = c(-0.2,0.5))
legend("right", pch = 20, col = rainbow(N_samp),cex = 0.8,legend = rnam)

plot(pc[,-1], pch = 20,cex=2.0, col = rainbow(N_samp), 
     main = "loadings for PC2 and PC3 in samples",
     xlab = axis.names[2],ylab = axis.names[3])
text(pc[,2],pc[,3], cex = 1, adj = c(-0.2,0.5))
legend("right", pch = 20, col = rainbow(N_samp),cex = 0.8,legend = rnam)

nouse = dev.off()
