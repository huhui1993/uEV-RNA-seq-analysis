
###################################################################################################
############################################################lnc
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后")
exp.tpm<-read.delim("all_sample_tpm.lnc.txt",header=T,row.names = 1)
expr<-read.delim("all_sample_count.addAnotation.txt",header=T,row.names = 1)

head(expr)
head(exp.tpm);dim(exp.tpm)
expr_matrix<-exp.tpm[,c(1:30)]

# 多维尺度分析（MDS）图:基于样本间距离矩阵，直观展示样本相似性
# 使用edgeR的plotMDS
library(edgeR)
library(limma)

# 计算logCPM
logCPM <- cpm(expr_matrix, log = TRUE, prior.count = 3)

# MDS图
mds_result <- plotMDS(logCPM, plot = FALSE)

mds_data <- data.frame(
  Dimension1 = mds_result$x,
  Dimension2 = mds_result$y,
  TimePoint = type$Time,
  Subject = type$Donor
)

ggplot(mds_data, aes(x = Dimension1, y = Dimension2, 
                     color = Subject)) +
  geom_point(size = 4) +
  geom_path(aes(group = Subject), alpha = 0.5) +
  labs(x = "Leading logFC dim 1", y = "Leading logFC dim 2") +
  theme_minimal()

#####高变异基因热图
# 选取变异最大的基因（如前500个高变异基因）
library(pheatmap)
library(RColorBrewer)

# 计算基因方差，选取高变异基因
gene_vars <- apply(logCPM, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:500]

# 创建热图
annotation_col <- data.frame(
  TimePoint = type$Time,
  Subject = type$Donor,
  row.names = colnames(expr_matrix)
)

pheatmap(logCPM[top_genes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = FALSE,  # 基因太多，不显示基因名
         show_colnames = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Top 500 variable genes")

###########4. 相关性矩阵热图
# 计算样本间相关性矩阵
cor_matrix <- cor(logCPM, method = "spearman")

# 可视化相关性矩阵
pheatmap(cor_matrix,
         annotation_col = annotation_col,
         annotation_row = annotation_col,
         main = "Sample-sample correlation",
         color = colorRampPalette(c("blue", "white", "red"))(100))

pheatmap(cor_matrix,
         main = "Sample-sample correlation",
         color = colorRampPalette(c("blue", "white", "red"))(100))



#################
phenotype<-read.delim("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\sample.type.txt",
                      header=T,row.names = NULL)
head(phenotype);dim(phenotype)
library(dplyr)
###提取特定供体
phenotype.1<-phenotype %>% subset(Donor %in% c("4","5","6","11","12","13"))
head(phenotype.1);dim(phenotype.1)

phenotype.d1<-phenotype %>% subset(Time %in% c("D1"))
head(phenotype.d1);dim(phenotype.d1)

exp.tpm.1<-exp.tpm[c("ENSG00000225528","ENSG00000293064"),]
exp.tpm.1.pro<-exp.tpm.1[,phenotype.1$sample]
dim(exp.tpm.1.pro)
#[1]  2 24
#############D1 sample exp
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\D1")
exp.tpm.d1<-exp.tpm[,phenotype.d1$sample]
dim(exp.tpm.d1);head(exp.tpm.d1)


exp.tpm.d1$geneid<-rownames(exp.tpm.d1)

library("reshape")
data_melt<-melt(exp.tpm.d1,id=c("geneid"))
head(data_melt)
colnames(data_melt)<-c("geneid","sample","tpm")
head(data_melt);dim(data_melt)

head(phenotype.d1)
data_melt.final<-merge(data_melt,phenotype.d1,by = "sample",all=T)
head(data_melt.final);dim(data_melt.final)

library(dplyr)
data_sd<-data_melt.final %>% group_by(geneid) %>% 
  dplyr::summarise(sd = sd(tpm))
head(data_sd)

data_mean<-data_melt.final %>% group_by(geneid) %>% 
  dplyr::summarise(mean = mean(tpm))
head(data_mean)

data_mean<-cbind(data_mean,expr[data_mean$geneid,c("gene_name","gene_biotype")])
head(data_mean);dim(data_mean)


#write.table(data_sd,
 #           "data_sd.TPM.donors.txt", 
  #          sep="\t", col.names=T,quote = F,row.names = F)

#write.table(data_mean,
 #           "data_mean.TPM.donors.txt", 
  #          sep="\t", col.names=T,quote = F,row.names = F)

m1 <- merge(data_mean, data_sd, by = c("geneid"))
head(m1)
write.table(m1,
            "data_mean.add_sd.TPM.12D1samples.lnc.txt", 
            sep="\t", col.names=T,quote = F,row.names = F)
head(m1)

#threshold_sd <- 2  # 设置标准差的阈值
S1 <- m1[m1$sd < 4, ]
S1<-S1[S1$mean > 1,]

S2 <- m1[m1$sd >= 4 & m1$sd <= 10, ]
S2<-S2[S2$mean > 1,]

S3 <- m1[m1$sd > 10, ]
S3<-S3[S3$mean > 1,]

dim(S1);dim(S2);dim(S3)

####第一部分：计算与可视化基因变异全景
#我们将使用变异系数作为核心指标，它消除了表达量级的影响，能公平地比较不同表达水平基因的个体间变异。
# 加载包
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

dim(exp.tpm.d1);head(exp.tpm.d1)
tpm_matrix<-exp.tpm.d1[,c(1:12)]
# 假设您的TPM矩阵名为 tpm_matrix (行=基因，列=样本)
# 1. 计算每个基因的关键统计量
gene_stats <- data.frame(gene = rownames(tpm_matrix),
  median_expression = apply(tpm_matrix, 1, median), # 中位表达水平
  mean_expression = apply(tpm_matrix, 1, mean),
  sd_expression = apply(tpm_matrix, 1, sd),         # 标准差
  cv_expression = apply(tpm_matrix, 1, function(x) sd(x) / mean(x) * 100)
  ) # 变异系数(%)
gene_stats.D1.lnc<-gene_stats
save(gene_stats.D1.lnc,file="gene_stats.D1.lnc.Rdata")

gene_stats<-gene_stats[gene_stats$mean_expression>1,]
  # 按CV从大到小排序
  gene_stats_sorted <- gene_stats[order(-gene_stats$cv_expression), ]
  save(gene_stats_sorted,file="gene_stats_sorted.lnc.RData")
  # 2. 全景可视化：基因表达水平与变异程度的关系
  p1 <- ggplot(gene_stats, aes(x = log10(median_expression + 1), y = cv_expression)) +
    geom_point(alpha = 0.5, size = 0.8, color = "grey60") +
    geom_smooth(method = "loess", color = "red", se = FALSE) + # 添加趋势线
    labs(x = "Log10(Median TPM + 1)", y = "Coefficient of Variation (%)",
         title = "Gene Expression Variability Across 12 Donors",
         subtitle = "Each point represents one gene") +
    theme_minimal()
  print(p1)
  ggsave(filename="./point.Gene_Expression_Variability_Across_12_Donors.lnc.pdf",plot=p1,
         device='pdf',path=".",width=3.5,height=3.5)
#结果解读：此图揭示关键规律。通常，低表达基因的CV值普遍偏高（技术噪声占比大），
  #而高表达基因的CV值降低并趋于稳定。那些“偏离红色趋势线、处于高表达高CV区域”的基因，
  #是您要寻找的真正具有显著生物学个体差异的候选基因。

####第二部分：对基因进行三级分类:我们将基于CV的分布，客观地将基因分为“稳定”、“中等变异”、“高变异”三类
  # 3. 基于变异系数(CV)的三分位数进行基因分类
  cv_quantiles <- quantile(gene_stats$cv_expression, probs = c(0.33, 0.67), na.rm = TRUE)
  
  gene_stats$variability_class <- cut(gene_stats$cv_expression,
                                      breaks = c(-Inf, cv_quantiles[1], cv_quantiles[2], Inf),
                                      labels = c("Stable", "Moderate", "Highly Variable"))
  
  # 查看各类基因数量
  table(gene_stats$variability_class)
  
  # 4. 绘制分类结果
  p2 <- ggplot(gene_stats, aes(x = variability_class, y = log10(median_expression + 1), fill = variability_class)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Variability Class", y = "Log10(Median TPM + 1)",
         title = "Expression Level Distribution Across Variability Classes",
         fill = "Class") +
    theme_minimal() +
    theme(legend.position = "none")
  print(p2)

  ggsave(filename="./violin.Gene_Expression_Variability_Distribution_Across_Variability_Classes.lnc.pdf",plot=p2,
         device='pdf',path=".",width=3.5,height=3.5)
#####第三部分：探索不同类别基因的特征
  # 5. 提取各类基因列表（用于后续功能分析）
  gene_stats<-gene_stats[gene_stats$mean_expression>1,]
  stable_genes <- gene_stats$gene[gene_stats$variability_class == "Stable"]
  
  moderate_genes <- gene_stats$gene[gene_stats$variability_class == "Moderate"]
  
  high_var_genes <- gene_stats$gene[gene_stats$variability_class == "Highly Variable"]
  
  
  gene_stats.high_var_genes<-gene_stats[c(high_var_genes),]
  gene_stats_sorted.high_var <- gene_stats.high_var_genes[order(-gene_stats.high_var_genes$cv_expression), ]
  
  gene_stats.moderate_genes<-gene_stats[c(moderate_genes),]
  gene_stats_sorted.moderate <- gene_stats.moderate_genes[order(-gene_stats.moderate_genes$cv_expression), ]
  
  gene_stats.stable_genes<-gene_stats[c(stable_genes),]
  gene_stats_sorted.stable <- gene_stats.stable_genes[order(-gene_stats.stable_genes$cv_expression), ]
  
  ####给基因add Symbol
  head(expr)
  tmp<-expr[,c("gene_name","gene_biotype")]
  tmp$gene<-rownames(expr)
  head(tmp)
  v1<-merge(tmp[gene_stats_sorted.high_var$gene,],gene_stats_sorted.high_var,by=c("gene"))
  v2<-merge(tmp[gene_stats_sorted.moderate$gene,],gene_stats_sorted.moderate,by=c("gene"))
  v3<-merge(tmp[gene_stats_sorted.stable$gene,],gene_stats_sorted.stable,by=c("gene"))
  
  write.table(v1,file='gene_stats_sorted.high_var.lnc.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
  write.table(v2,file='gene_stats_sorted.moderate.lnc.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
  write.table(v3,file='gene_stats_sorted.stable.lnc.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
  # 6. 功能特征分析（以高变异基因为例，需安装并加载clusterProfiler）
  # BiocManager::install("clusterProfiler")
  # BiocManager::install("org.Hs.eg.db")
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  # 将基因符号转换为Entrez ID（示例）
  high_var_entrez <- bitr(high_var_genes, fromType = "ENSEMBL", 
                          toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)$ENTREZID
  
  # 进行GO富集分析
  go_enrich <- enrichGO(gene = high_var_entrez,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP", # 生物过程
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.1,
                        readable = TRUE)
  # 可视化
  dotplot(go_enrich, title = "GO Enrichment of Highly Variable Genes") + 
    theme(axis.text.y = element_text(size = 9))
  p<-dotplot(go_enrich, title = "GO Enrichment of Highly Variable Genes") + 
    theme(axis.text.y = element_text(size = 9))
  ggsave(filename = "GO.Highly_Variable_Genes.higher5.pdf", 
         plot =p,width = 13, height = 11, units = 'cm')
  bp <- setReadable(go_enrich, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
  write.table(bp,file="GO.BP.Highly_Variable_Genes.higher5.txt",sep='\t',
              quote=FALSE,row.names=FALSE,col.names=TRUE)
  
  ####进行KEGG富集分析
  library(R.utils)
  R.utils::setOption("clusterProfiler.download.method","auto")
  
  ego_KEGG <- enrichKEGG(gene         = high_var_entrez,
                         organism     = 'hsa',
                         pvalueCutoff = 1,
                         qvalueCutoff  = 1)
  ego_KEGG <- setReadable(ego_KEGG, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
  write.table(ego_KEGG,file="./KEGG.Highly_Variable_Genes.higher5.txt",
              sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
  dotplot(ego_KEGG,showCategory=10,title="KEGG Enrichment of Highly Variable Genes")
  p1<-dotplot(ego_KEGG,showCategory=10,title="KEGG Enrichment of Highly Variable Genes")
  ggsave(filename = "./KEGG.Highly_Variable_Genes.higher5.pdf", plot =p1,
         width = 13, height = 11, units = 'cm')
  
  ####6.2 moderate_genes
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  # 将基因符号转换为Entrez ID（示例）
  moderate_entrez <- bitr(moderate_genes, fromType = "ENSEMBL", 
                          toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)$ENTREZID
  
  # 进行GO富集分析
  go_enrich <- enrichGO(gene = moderate_entrez,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP", # 生物过程
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.1,
                        readable = TRUE)
  # 可视化
  dotplot(go_enrich, title = "GO Enrichment of Moderate Genes") + 
    theme(axis.text.y = element_text(size = 9))
  p<-dotplot(go_enrich, title = "GO Enrichment of Moderate Genes") + 
    theme(axis.text.y = element_text(size = 9))
  ggsave(filename = "GO.moderate_Genes.higher5.pdf", 
         plot =p,width = 13, height = 11, units = 'cm')
  bp <- setReadable(go_enrich, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
  write.table(bp,file="GO.BP.moderate_Genes.higher5.txt",sep='\t',
              quote=FALSE,row.names=FALSE,col.names=TRUE)
  
  ####进行KEGG富集分析
  library(R.utils)
  R.utils::setOption("clusterProfiler.download.method","auto")
  
  ego_KEGG <- enrichKEGG(gene         = moderate_entrez,
                         organism     = 'hsa',
                         pvalueCutoff = 0.1,
                         qvalueCutoff  = 1)
  ego_KEGG <- setReadable(ego_KEGG, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
  write.table(ego_KEGG,file="./KEGG.moderate_Genes.higher5.txt",
              sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
  dotplot(ego_KEGG,showCategory=10,title="KEGG Enrichment of Moderate Genes")
  p1<-dotplot(ego_KEGG,showCategory=10,title="KEGG Enrichment of Moderate Genes")
  ggsave(filename = "./KEGG.moderate_Genes.higher5.pdf", plot =p1,
         width = 13, height = 11, units = 'cm')
  
  ####6.3 stable_genes
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  # 将基因符号转换为Entrez ID（示例）
  stable_entrez <- bitr(stable_genes, fromType = "ENSEMBL", 
                          toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)$ENTREZID
  
  # 进行GO富集分析
  go_enrich <- enrichGO(gene = stable_entrez,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP", # 生物过程
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.1,
                        readable = TRUE)
  # 可视化
  dotplot(go_enrich, title = "GO Enrichment of stable Genes") + 
    theme(axis.text.y = element_text(size = 9))
  p<-dotplot(go_enrich, title = "GO Enrichment of stable Genes") + 
    theme(axis.text.y = element_text(size = 9))
  ggsave(filename = "GO.stable_Genes.higher5.pdf", 
         plot =p,width = 13, height = 11, units = 'cm')
  bp <- setReadable(go_enrich, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
  write.table(bp,file="GO.BP.stable_Genes.higher5.txt",sep='\t',
              quote=FALSE,row.names=FALSE,col.names=TRUE)
  
  ####进行KEGG富集分析
  library(R.utils)
  R.utils::setOption("clusterProfiler.download.method","auto")
  
  ego_KEGG <- enrichKEGG(gene         = stable_entrez,
                         organism     = 'hsa',
                         pvalueCutoff = 0.1,
                         qvalueCutoff  = 1)
  ego_KEGG <- setReadable(ego_KEGG, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
  write.table(ego_KEGG,file="./KEGG.stable_Genes.higher5.txt",
              sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
  dotplot(ego_KEGG,showCategory=10,title="KEGG Enrichment of stable Genes")
  p1<-dotplot(ego_KEGG,showCategory=10,title="KEGG Enrichment of stable Genes")
  ggsave(filename = "./KEGG.stable_Genes.higher5.pdf", plot =p1,
         width = 13, height = 11, units = 'cm')
  
  
  # 7. 表达模式可视化（热图展示前50个高变异基因）
  top_n <- 50
  top_var_genes <- gene_stats_sorted.high_var$gene[1:top_n]
  heatmap_data <- tpm_matrix[top_var_genes, ]
  
  #v1<-merge(tmp[gene_stats_sorted.high_var$gene,],gene_stats_sorted.high_var,by=c("gene"))
  #head(v1)
  #rownames(v1)<-v1$gene
  head(expr)
  m1<-cbind(expr[rownames(heatmap_data),c("gene_name")],heatmap_data)
  
  write.table(m1,file='top50gene_stats_sorted.high_var.lnc.txt',sep='\t',
              quote=FALSE,row.names=TRUE,col.names=TRUE)
  # 对基因和样本进行聚类
  pdf("lnc.Highly_Variable_Genes.top50.pdf",4,4)
  pheatmap(log2(heatmap_data + 1),
           scale = "row", # 按行标准化，突出基因在不同样本间的模式
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           show_rownames = FALSE, # 基因太多，不显示名字
           main = paste("Top", top_n, "Most Variable Genes"),
           fontsize_col = 9)
  dev.off()
  
  data.m<-apply(log2(heatmap_data+1),1,scale)
  data.m<-t(data.m)
  colnames(data.m)<-colnames(heatmap_data)
  library("ComplexHeatmap")
  library(circlize)
  mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  pdf("lnc.Highly_Variable_Genes.top50.new.pdf",3.8,4)
  mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  data.m<-na.omit(data.m)
  Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
          col=mycol,name = "Exp\nZ score", 
          row_dend_width = unit(8, "mm"),
          row_names_gp=gpar(fontsize = 10),
          row_title_gp = gpar(col = "black",fontsize = 12),
          show_row_names = FALSE,
          column_title =paste("Top", top_n, "Most Variable Genes") )
  dev.off()
  
  # 7.2 表达模式可视化（热图展示前50个中变异基因）
  top_n <- 50
  top_moderate_genes <- gene_stats_sorted.moderate$gene[1:top_n]
  heatmap_data <- tpm_matrix[top_moderate_genes, ]
  
  head(expr)
  m1<-cbind(expr[rownames(heatmap_data),c("gene_name")],heatmap_data)
  
  write.table(m1,file='top50gene_stats_sorted.moderate.lnc.txt',sep='\t',
              quote=FALSE,row.names=TRUE,col.names=TRUE)
  
  # 对基因和样本进行聚类
  pdf("lnc.moderate_Variable_Genes.top50.pdf",4,4)
  pheatmap(log2(heatmap_data + 1),
           scale = "row", # 按行标准化，突出基因在不同样本间的模式
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           show_rownames = FALSE, # 基因太多，不显示名字
           main = paste("Top", top_n, "moderate Variable Genes"),
           fontsize_col = 9)
  dev.off()
  
  data.m<-apply(log2(heatmap_data+1),1,scale)
  data.m<-t(data.m)
  colnames(data.m)<-colnames(heatmap_data)
  library("ComplexHeatmap")
  library(circlize)
  mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  pdf("lnc.moderate_Variable_Genes.top50.new.pdf",3.8,4)
  mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  data.m<-na.omit(data.m)
  Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
          col=mycol,name = "Exp\nZ score", 
          row_dend_width = unit(8, "mm"),
          row_names_gp=gpar(fontsize = 10),
          row_title_gp = gpar(col = "black",fontsize = 12),
          show_row_names = FALSE,
          column_title =paste("Top", top_n, "moderate Variable Genes") )
  dev.off()
  
  # 7.3 表达模式可视化（热图展示前50个低变异基因）
  top_n <- 50
  #top_stable_genes <- gene_stats_sorted.stable$gene[1:top_n]
  #heatmap_data <- tpm_matrix[top_stable_genes, ]
  total_genes <- nrow(gene_stats_sorted.stable)
  # 方法：计算尾端的索引范围
  # 从 (总行数 - top_n + 1) 到 (总行数)
  top_stable_genes <- gene_stats_sorted.stable$gene[(total_genes - top_n + 1):total_genes]
  # 从原始矩阵中提取这些稳定基因的表达数据
  heatmap_data <- tpm_matrix[top_stable_genes, ]
  
  head(expr)
  m1<-cbind(expr[rownames(heatmap_data),c("gene_name")],heatmap_data)
  
  write.table(m1,file='top50gene_stats_sorted.stable.lnc.txt',sep='\t',
              quote=FALSE,row.names=TRUE,col.names=TRUE)
  # 对基因和样本进行聚类
  pdf("lnc.stable_Genes.top50.pdf",4,4)
  pheatmap(log2(heatmap_data + 1),
           scale = "row", # 按行标准化，突出基因在不同样本间的模式
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           show_rownames = FALSE, # 基因太多，不显示名字
           main = paste("Top", top_n, "Most stable Genes"),
           fontsize_col = 9)
  dev.off()
  
  data.m<-apply(log2(heatmap_data+1),1,scale)
  data.m<-t(data.m)
  colnames(data.m)<-colnames(heatmap_data)
  library("ComplexHeatmap")
  library(circlize)
  mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  pdf("lnc.stable_Genes.top50.new.pdf",3.8,4)
  mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  data.m<-na.omit(data.m)
  Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
          col=mycol,name = "Exp\nZ score", 
          row_dend_width = unit(8, "mm"),
          row_names_gp=gpar(fontsize = 10),
          row_title_gp = gpar(col = "black",fontsize = 12),
          show_row_names = FALSE,
          column_title =paste("Top", top_n, "Most stable Genes") )
  dev.off()

  #############################################################################
  ####计算高中低变异基因中高中低表达基因数目百分比pie图
  #  计算每个基因的关键统计量
  gene_stats <- data.frame(gene = rownames(tpm_matrix),
                           median_expression = apply(tpm_matrix, 1, median), # 中位表达水平
                           mean_expression = apply(tpm_matrix, 1, mean),
                           sd_expression = apply(tpm_matrix, 1, sd),         # 标准差
                           cv_expression = apply(tpm_matrix, 1, function(x) sd(x) / mean(x) * 100)
  ) # 变异系数(%)
  
  # 3. 基于变异系数(CV)的三分位数进行基因分类
  cv_quantiles <- quantile(gene_stats$cv_expression, probs = c(0.33, 0.67), na.rm = TRUE)
  
  gene_stats$variability_class <- cut(gene_stats$cv_expression,
                                      breaks = c(-Inf, cv_quantiles[1], cv_quantiles[2], Inf),
                                      labels = c("Stable", "Moderate", "Highly Variable"))
  
  #第一步：根据TPM值划分基因表达水平
  # 在已有gene_stats数据框中添加表达水平类别
  # 使用mean_expression（平均TPM）作为表达水平指标
  gene_stats$expression_level <- cut(gene_stats$mean_expression,
                                     breaks = c(0, 1, 100, Inf),
                                     labels = c("Low (0-1 TPM)", "Medium (1-100 TPM)", "High (>100 TPM)"),
                                     include.lowest = TRUE,  # 包含0
                                     right = FALSE)  # 右开区间 [0,1), [1,100), [100,Inf)
  
  # 查看各表达水平基因数量
  table(gene_stats$expression_level)
  
  # 同时查看变异类别和表达水平的交叉分布
  cross_table <- table(gene_stats$variability_class, gene_stats$expression_level)
  print(cross_table)
  
  #                   Low (0-1 TPM) Medium (1-100 TPM) High (>100 TPM)
  #Stable                   3693               1793              28
  #Moderate                 4332               1343               6
  #Highly Variable          4657                857               0
  
  ##第二步：绘制表达水平与变异程度的分布图
  # 1. 表达水平与CV值的关系（箱线图）
  library(ggplot2)
  library(RColorBrewer)
  
  p_expr_vs_cv <- ggplot(gene_stats, aes(x = expression_level, y = cv_expression, fill = expression_level)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    scale_fill_brewer(palette = "Set3") +
    scale_y_log10() +  # 由于CV值可能跨度大，使用对数坐标
    labs(x = "Expression Level (TPM)", 
         y = "CV (%) (Log10 scale)",
         title = "Coefficient of Variation Across Expression Levels",
         subtitle = "Each box represents the CV distribution within an expression category") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 15, hjust = 1))
  
  print(p_expr_vs_cv)
  ggsave("./boxplot_CV_vs_Expression_Level.lnc.pdf", p_expr_vs_cv, 
         width = 4, height = 4, device = "pdf")
  
  # 2. 双变量分布散点图（展示每个基因的分布）
  p_scatter <- ggplot(gene_stats, aes(x = log10(mean_expression + 0.1), y = cv_expression, 
                                      color = variability_class)) +
    geom_point(alpha = 0.5, size = 0.8) +
    facet_wrap(~ expression_level, scales = "free_x") +
    scale_color_manual(values = c("Stable" = "#4daf4a", 
                                  "Moderate" = "#ff7f00", 
                                  "Highly Variable" = "#e41a1c")) +
    labs(x = "Log10(Mean TPM + 0.1)", 
         y = "CV (%)",
         title = "Gene Distribution: Expression Level vs Variability",
         color = "Variability Class") +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(color = "black", size = 10),
          axis.text.y = element_text(color = "black",size = 10))
  
  print(p_scatter)
  ggsave("./scatter_Expression_vs_CV_by_Level.lnc.pdf", p_scatter, 
         width = 5.5, height = 4, device = "pdf")
  ####第三步：制作饼图展示百分比分布
  # 计算每个变异类别中表达水平分布的百分比
  library(dplyr)
  
  percent_data <- gene_stats %>%
    group_by(variability_class, expression_level) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(variability_class) %>%
    mutate(percent = count / sum(count) * 100,
           label = paste0(round(percent, 1), "%\n(n=", count, ")")) %>%
    ungroup()
  
  # 查看计算结果
  print(percent_data)
  
  # 定义颜色方案
  expr_colors <- c("Low (0-1 TPM)" = "#41F0AE", 
                   "Medium (1-100 TPM)" = "#FFC080", 
                   "High (>100 TPM)" = "#FF8080")
  
  #"#41F0AE","#FFC080", "#FF8080"
  
  # 绘制饼图（三个变异类别并列）
  library(patchwork)  # 用于组合多个图形
  
  # 创建三个饼图函数
  create_pie_chart <- function(data, class_name) {
    class_data <- data[data$variability_class == class_name, ]
    
    # 计算标签位置
    class_data <- class_data %>%
      arrange(desc(expression_level)) %>%
      mutate(ypos = cumsum(percent) - 0.5 * percent)
    
    ggplot(class_data, aes(x = "", y = percent, fill = expression_level)) +
      geom_bar(width = 1, stat = "identity", color = "white") +
      coord_polar("y", start = 0) +
      geom_text(aes(y = ypos, label = label), 
                color = "black", size = 3, fontface = "bold") +
      scale_fill_manual(values = expr_colors) +
      labs(title = class_name,
           fill = "Expression Level") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
            legend.position = "none")  # 移除每个子图的图例，后面添加共用图例
  }
  
  # 生成三个饼图
  pie_stable <- create_pie_chart(percent_data, "Stable")
  pie_moderate <- create_pie_chart(percent_data, "Moderate")
  pie_high_var <- create_pie_chart(percent_data, "Highly Variable")
  
  # 创建共用图例
  legend_plot <- ggplot(percent_data, aes(x = expression_level, fill = expression_level)) +
    geom_bar() +
    scale_fill_manual(values = expr_colors, 
                      name = "Expression Level (TPM)") +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 9))
  
  # 提取图例
  library(gridExtra)
  get_legend <- function(plot) {
    tmp <- ggplot_gtable(ggplot_build(plot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    tmp$grobs[[leg]]
  }
  shared_legend <- get_legend(legend_plot)
  
  # 组合图形
  library(grid)
  top_row <- grid.arrange(pie_stable, pie_moderate, pie_high_var, 
                          ncol = 3, 
                          top = textGrob("Expression Level Distribution Within Each Variability Class", 
                                         gp = gpar(fontsize = 14, fontface = "bold")))
  
  final_plot <- grid.arrange(top_row, shared_legend, 
                             nrow = 2, 
                             heights = c(0.85, 0.15))
  
  # 保存组合图形
  ggsave("./pie_Expression_Level_Distribution_by_Variability.lnc.pdf", 
         final_plot, width = 7, height = 2.3, device = "pdf")
  ####第四步：制作堆叠条形图（备选，更易比较）
  # 堆叠百分比条形图（作为饼图的补充，更容易比较）
  percent_data<-na.omit(percent_data)
  p_stacked <- ggplot(percent_data, aes(x = variability_class, y = percent, fill = expression_level)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    geom_text(aes(label = paste0(round(percent, 1), "%")), 
              position = position_stack(vjust = 0.5, reverse = TRUE),
              size = 3, color = "black", fontface = "bold") +
    scale_fill_manual(values = expr_colors) +
    labs(x = "Variability Class", 
         y = "Percentage (%)",
         title = "Expression Level Composition Within Each Variability Class",
         fill = "Expression Level (TPM)") +
    theme_minimal() +
    theme(axis.text.x = element_text(face = "bold", size = 11),
          axis.text.y = element_text( color = "black",size = 11),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "bottom")
  
  print(p_stacked)
  ggsave("./stacked_bar_Expression_Level_Distribution.lnc.pdf", p_stacked, 
         width = 4.5, height = 4, device = "pdf")
  ####第五步：创建汇总表格
  # 创建详细的汇总表格
  summary_table <- gene_stats %>%
    group_by(variability_class, expression_level) %>%
    summarise(
      Gene_Count = n(),
      Percentage = round(n() / nrow(gene_stats) * 100, 2),
      Mean_CV = round(mean(cv_expression, na.rm = TRUE), 2),
      Median_TPM = round(median(mean_expression, na.rm = TRUE), 2),
      .groups = "drop"
    )
  
  # 打印表格
  print(summary_table)
  
  # 保存为CSV文件
  write.csv(summary_table, "./summary_Expression_Variability_Table.lnc.csv", row.names = FALSE)
  
  # 更美观的表格展示（使用kableExtra，可选）
  if(requireNamespace("kableExtra", quietly = TRUE)) {
    library(kableExtra)
    summary_table %>%
      kbl(caption = "Summary of Gene Distribution by Variability and Expression Level") %>%
      kable_paper("hover", full_width = FALSE) %>%
      column_spec(1, bold = TRUE) %>%
      collapse_rows(columns = 1, valign = "top")
  } 

##########heatmap D1 sample
df1<-df[,phenotype.d1$sample]
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
library("ComplexHeatmap")
library(circlize)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf(file = "./D1/heatmap.D1sample.expTPM.scaled.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()

ha <- HeatmapAnnotation(Sex = phenotype.d1$Sex, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1") ))

data.m<-na.omit(data.m)
pdf(file = "./D1/heatmap.D1sample.expTPM.scaled.addlegend.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE,top_annotation = ha)

dev.off()

#pdf(file = "./D1/heatmap.D1sample.expTPM.noscaled.pdf", 6.5, 1.6)
#mycol <- colorRamp2(c(0, 2, 4), c("blue", "white", "red"))
#Heatmap(as.matrix(df1),cluster_columns = TRUE,cluster_rows = TRUE,
#        col=mycol,name = "expression", 
#        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
#        row_names_gp=gpar(fontsize = 10),
#        row_title_gp = gpar(col = "black",fontsize = 12),
#        show_row_names = FALSE)
#dev.off()
#########################D1样本低中高表达基因的表达量热图
#####高表达并集
df1<-df.pro[,phenotype.d1$sample]
df1$geneid<-rownames(df1)
#S0102=as.factor(df1[df1[,1]>=100,13])
#S0104=as.factor(df1[df1[,2]>=100,13])
#S0105=as.factor(df1[df1[,3]>=100,13])
#S0106=as.factor(df1[df1[,4]>=100,13])
#S0107=as.factor(df1[df1[,5]>=100,13])
#S0108=as.factor(df1[df1[,6]>=100,13])
#S0111=as.factor(df1[df1[,7]>=100,13])
#S0112=as.factor(df1[df1[,8]>=100,13])
#S0113=as.factor(df1[df1[,9]>=100,13])
#S0115=as.factor(df1[df1[,10]>=100,13])
#S0116=as.factor(df1[df1[,11]>=100,13])
#S0117=as.factor(df1[df1[,12]>=100,13])
##多数据集求并集
genes <- Reduce(union,list(
df1[df1[,1]>=100,13], df1[df1[,2]>=100,13], 
df1[df1[,3]>=100,13], df1[df1[,4]>=100,13],
df1[df1[,5]>=100,13],df1[df1[,6]>=100,13],
df1[df1[,7]>=100,13],df1[df1[,8]>=100,13],
df1[df1[,9]>=100,13],df1[df1[,10]>=100,13],
df1[df1[,11]>=100,13],df1[df1[,12]>=100,13]))
length(genes)
dim(df)
df2<-df.pro[genes,phenotype.d1$sample]
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
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)

pdf(file = "./D1/heatmap.pro.D1sample.high.expTPM.scaled.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()

####
#####高表达交集
df1<-df.pro[,phenotype.d1$sample]
df1$geneid<-rownames(df1)

##多数据集求交集
genes <- Reduce(intersect,list(
  df1[df1[,1]>=100,13], df1[df1[,2]>=100,13], 
  df1[df1[,3]>=100,13], df1[df1[,4]>=100,13],
  df1[df1[,5]>=100,13],df1[df1[,6]>=100,13],
  df1[df1[,7]>=100,13],df1[df1[,8]>=100,13],
  df1[df1[,9]>=100,13],df1[df1[,10]>=100,13],
  df1[df1[,11]>=100,13],df1[df1[,12]>=100,13]))
length(genes)
dim(df)
df2<-df.pro[genes,phenotype.d1$sample]
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
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)

pdf(file = "./D1/heatmap.pro.D1sample.high.expTPM.intersect.scaled.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()
##
##############################################
######box plot
#exp.tpm.d1<-df[apply(df,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(exp.tpm.d1)

exp.tpm.d1<-df.pro[genes,phenotype.d1$sample]
#log(TPM+1)
exp.tpm.d1<-log2(exp.tpm.d1+1)
normalized_data<-as.data.frame(exp.tpm.d1)
normalized_data$geneid<-rownames(normalized_data)
normalized_data

library("reshape")
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","tpm")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype.d1,by = "sample",all=T)
head(exp_melt.final)
#exp_melt.final.1<-exp_melt.final %>% subset(Donor %in% c("4","5","6","11","12","13")) 
dim(exp_melt.final);
dim(exp_melt.final.1)

head(exp_melt.final)

library("ggpubr")
library("reshape")
library("ggplot2")

####by time
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt.final, x = "Donor", y = "tpm",
               fill="Donor",
               
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
ggsave(filename="./D1/boxplot.bydonor.420gene.log2TPM.pdf",plot=p,
       device='pdf',path=".",width=6,height=4.5)

################add p
p <- ggboxplot(exp_melt.final, x = "Donor", y = "tpm",
               color ="Donor",
              show.legend = FALSE)+
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
                     ref.group = "2")                    # Pairwise comparison against reference
p1
ggsave(filename="./D1/boxplot.bydonor.420gene.log2TPM.addP.d2.pdf",plot=p1,
       device='pdf',path=".",width=6,height=3.6)

#p1<-p+stat_compare_means(method = "anova", label.y = 18)+      # Add global p-value
 # stat_compare_means(label = "p.signif", method = "t.test",
#                     ref.group = "13")                    # Pairwise comparison against reference
#p1
#ggsave(filename="./boxplot.bydonor.2888gene.log2TPM.addP.13.pdf",plot=p1,
#       device='pdf',path=".",width=6,height=3.6)
######共同高表达基因功能分析
###########富集分析
genes
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library("ggplot2")
data(geneList, package="DOSE")

#keygene=df.c[which(df.c$annotation1=="Promoter (<=1kb)"),]

Gene_list<- bitr(genes,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

ego_KEGG <- enrichKEGG(gene         = Gene_list$ENTREZID,
                       organism     = 'hsa',
                       pvalueCutoff = 0.1,
                       qvalueCutoff  = 1)
ego_KEGG <- setReadable(ego_KEGG, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
write.table(ego_KEGG,file="./D1/420gene.enrichKEGG.txt",
            sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
dotplot(ego_KEGG,showCategory=40,title="KEGG Enrichment")
p1<-dotplot(ego_KEGG,showCategory=40,title="KEGG Enrichment")
ggsave(filename = "./D1/420gene.enrichKEGG.pdf", plot =p1,
       width = 16, height = 16, units = 'cm')
#########
#自己的代码画KEGG
kegg<-read.delim('./D1/420gene.enrichKEGG.filter.txt', stringsAsFactors = FALSE)
dim(kegg)
head(kegg)

kegg<-kegg[which(kegg$Count>=5),]
library(ggplot2)
data<-kegg
head(data)
data$LogFDR<--log10(data$pvalue)
data <- data[order(data$Count,data$LogFDR,decreasing = c(TRUE, FALSE),method = "radix"), ]
data$Description <- factor(data$Description, levels = data$Description)

pp = ggplot(data,aes(Count,Description))
p.final<-pp + geom_point(aes(color=(LogFDR),size=Count))+theme_bw()+
  labs(x="Count",y="",title="KEGG enrichment")+
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1,vjust =1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=11,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=12,colour="black"),
        axis.title.y = element_text(size=12,colour="black"),
        plot.title = element_text(size=13,colour="black"),
        legend.key.size = unit(5,'mm'))+
  scale_colour_gradient(low="DodgerBlue",high="NavyBlue",name="-log P")+
  scale_y_discrete(labels=c("Chemical carcinogenesis - reactive oxygen species"="Chemical carcinogenesis-\nreactive oxygen species"))

p.final
ggsave(filename="./D1/420gene.enrichKEGG.new.pdf",plot=p.final,
       device='pdf',path=".",width=5.7,height=5.4)

####
ego_bp <- enrichGO(gene         = Gene_list$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
ego_bp <- setReadable(ego_bp, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

write.table(ego_bp,file="./D1/420gene.enrichGO.txt",
            sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
dotplot(ego_bp,showCategory=40,title="GO Enrichment")
#p1<-dotplot(ego_bp,showCategory=40,title="BP Enrichment")
#ggsave(filename = "enrichBP.DEG.final.pdf", 
#      plot =p1,width = 19, height = 19, units = 'cm')
#GO的三种条目分开画
go<-read.delim('./D1/420gene.enrichGO.filter.txt', stringsAsFactors = FALSE)
dim(go)
#go$term <- paste(go$ID, go$Description, sep = ': ')
go$term<-go$Description

go <- go[order(go$ONTOLOGY, go$p.adjust, decreasing = c(TRUE, TRUE),method = "radix"), ]
go$term <- factor(go$term, levels = go$term)

p<-ggplot(go, aes(term, -log10(p.adjust))) +
  geom_col(aes(fill = ONTOLOGY), width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c('#D06660', '#5AAD36', '#6C85F5')) +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x = '', y = '-Log10 P-Value\n',title="GO enrichment")
p
ggsave(filename="./D1/420gene.enrichGO.barplot.pdf",plot=p,
       device='pdf',path=".",width=8,height=7)



#####中表达
df1<-df.pro[,phenotype.d1$sample]
df1$geneid<-rownames(df1)

##多数据集求并集
genes <- Reduce(union,list(
df1[df1[,1]>1 & df1[,1]<100,13],
df1[df1[,2]>1 & df1[,2]<100,13],
df1[df1[,3]>1 & df1[,3]<100,13],
df1[df1[,4]>1 & df1[,4]<100,13],
df1[df1[,5]>1 & df1[,5]<100,13],
df1[df1[,6]>1 & df1[,6]<100,13],
df1[df1[,7]>1 & df1[,7]<100,13],
df1[df1[,8]>1 & df1[,8]<100,13],
df1[df1[,9]>1 & df1[,9]<100,13],
df1[df1[,10]>1 & df1[,10]<100,13],
df1[df1[,11]>1 & df1[,11]<100,13],
df1[df1[,12]>1 & df1[,12]<100,13]))

length(genes)
dim(df)
df2<-df.pro[genes,phenotype.d1$sample]
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
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)

pdf(file = "./D1/heatmap.pro.D1sample.medium.expTPM.scaled.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()
#########lnc
df1<-df.lnc[,phenotype.d1$sample]
dim(df1)
df1$geneid<-rownames(df1)
dim(df1)
##多数据集求并集
genes <- Reduce(union,list(
  df1[df1[,1]>=100,13], df1[df1[,2]>=100,13], 
  df1[df1[,3]>=100,13], df1[df1[,4]>=100,13],
  df1[df1[,5]>=100,13],df1[df1[,6]>=100,13],
  df1[df1[,7]>=100,13],df1[df1[,8]>=100,13],
  df1[df1[,9]>=100,13],df1[df1[,10]>=100,13],
  df1[df1[,11]>=100,13],df1[df1[,12]>=100,13]))
length(genes)
dim(df)
df2<-df.lnc[genes,phenotype.d1$sample]
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
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)

pdf(file = "./D1/heatmap.lnc.D1sample.high.expTPM.scaled.pdf", 5, 5)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data.m<-na.omit(data.m)
Heatmap(as.matrix(data.m),cluster_columns = TRUE,cluster_rows = TRUE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.d1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = FALSE)
dev.off()

######################################################################################3
head(phenotype.1)
annot_df<-data.frame(type=phenotype.1$Time)
#info$stim<-factor(info$stim,levels=c("NC","D1.5","D7"))
rownames(phenotype.1)<-phenotype.1$sample

ha2 = HeatmapAnnotation(bar = phenotype.1$Sex,
                        col = list(bar = c("female" = "Salmon", "male" = "DodgerBlue1")))
data.m<-na.omit(data.m)
#data<-data[rownames(data.m),]
#dim(data);dim(data.m)
#ha = rowAnnotation(bar = data$Type,
 #                  col = list(bar = c("Cytotoxic" = "Orange1", "Exhausted" = "LightYellow4")))

#ha = HeatmapAnnotation(bar = phenotype$Time,
 #                  col = list(type = c(rainbow(length(unique(phenotype$Time)))))
#)

annot_df <- data.frame(Sex = phenotype.1$Sex, Time = phenotype.1$Time)
rownames(annot_df)<-rownames(phenotype.1)
col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
           Time = c("D1" = "Firebrick1","D15" = "RoyalBlue1",
                    "D29" = "orange","D43" = "LimeGreen") )

ha <- HeatmapAnnotation(Sex = phenotype.1$Sex, Time = phenotype.1$Time, 
                        col = list(Sex = c("female" = "Salmon", "male" = "DodgerBlue1"), 
                                   Time = c("D1" = "Brown1","D15" = "RoyalBlue1",
                                            "D29" = "orange","D43" = "LimeGreen") ))

pdf(file = "./2个共有基因表达/heatmap.2gene.expTPM.scaled.addlegend.pdf", 7, 2.2)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = TRUE,top_annotation = ha)

dev.off()
#
pdf(file = "./2个共有基因表达/heatmap.2gene.expTPM.noscaled.addlegend.pdf", 7, 2.2)
mycol <- colorRamp2(c(0, 2, 4), c("blue", "white", "red"))
Heatmap(as.matrix(exp.tpm.1.pro),cluster_columns = FALSE,cluster_rows = FALSE,
        col=mycol,name = "expression", 
        row_dend_width = unit(8, "mm"),column_split = phenotype.1$Donor,
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = TRUE,top_annotation = ha)
dev.off()
############################################################################
######box plot
#exp.tpm.d1<-df[apply(df,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(exp.tpm.d1)
#log(TPM+1)
#exp.tpm.d1<-log2(exp.tpm.d1+1)
exp.tpm.d1<-exp.tpm.1.pro
normalized_data<-as.data.frame(exp.tpm.d1)
normalized_data$geneid<-rownames(normalized_data)
normalized_data

library("reshape")
exp_melt<-melt(normalized_data,id=c("geneid"))
head(exp_melt)
colnames(exp_melt)<-c("geneid","sample","tpm")
head(exp_melt);dim(exp_melt)
exp_melt.final<-merge(exp_melt,phenotype.1,by = "sample",all=T)
head(exp_melt.final)
#exp_melt.final.1<-exp_melt.final %>% subset(scaled_count>=-0.3 & scaled_count<=0.2) 
dim(exp_melt.final);
#dim(exp_melt.final.1)

head(exp_melt.final)

library("ggpubr")
library("reshape")
library("ggplot2")

####by time
head(exp_melt.final)
library("ggpubr")
library("reshape")
library("ggplot2")
p <- ggboxplot(exp_melt.final, x = "geneid", y = "tpm",
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
ggsave(filename="./2个共有基因表达/boxplot.bytime.2gene.log2TPM.pdf",plot=p,
       device='pdf',path=".",width=6,height=4.5)

################add p
p <- ggboxplot(exp_melt.final, x = "geneid", y = "tpm",
               color ="Time",
               palette = "nejm",show.legend = FALSE)+
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
my_comparisons <- list( c("D1", "D15"), c("D15", "D29"), c("D29", "D43") )
p2<-p +
  stat_compare_means(comparisons = my_comparisons,label = "p.format",
                     size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2
ggsave(filename="./2个共有基因表达/boxplot.bytime.2gene.log2TPM.pdf",plot=p2,
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
#
head(exp_melt.final)
exp_melt.final.1.a <- exp_melt.final %>% subset(geneid=="ENSG00000225528")
exp_melt.final.1.b <- exp_melt.final %>% subset(geneid=="ENSG00000293064")

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

ggsave(filename="./2个共有基因表达/ENSG00000225528.boxplot.bytime.log2TPM.addP.pdf",plot=p2,
       device='pdf',path=".",width=3.2,height=3.2)
#
exp_melt.final.1.b <- exp_melt.final %>% subset(geneid=="ENSG00000293064")

p <- ggboxplot(exp_melt.final.1.b, x = "Time", y = "tpm",
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

ggsave(filename="./2个共有基因表达/ENSG00000293064.boxplot.bytime.log2TPM.addP.pdf",plot=p2,
       device='pdf',path=".",width=3.2,height=3.2)
