#实例：
#setwd("./GSE72922")
#setwd("D://项目//单细胞和乳酸化//息肉结肠癌GEO//GSE72922")
setwd("/public/home/huhui/project/rna_urine_EV/result/STAR/EV_origin/GEO_data/GSE72922")
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\独立验证集GEO\\GSE72922")
#BiocManager::install("GEOquery")#若已安装，此步骤跳过
library(GEOquery)
getwd()#查看当前路径
gse<-getGEO('GSE72922',destdir =".")##根据GSE号来下载数据，下载_series_matrix.txt.gz,同时自动下载平台文件GPL14951.soft
save(gse,file = 'gset.Rdata')
load("gset.Rdata")
Gset <- gse[[1]]
#获取样本一般信息,样本名字、类型、状态等
pdata<-pData(Gset)
class(pdata)
View(pdata)
write.table(pdata,file="GSE.sampleInfo.xls",quote = F,sep="\t",row.names=T)

unique(pdata$`disease status:ch1`)
type=as.data.frame(pdata$`disease status:ch1`)
rownames(type)<-rownames(pdata)
head(type)
colnames(type)<-c("group")
type$sample<-rownames(type)
head(type)
write.table(type,file="sampletype.txt",quote = F,sep="\t",row.names=F)
unique(type$group)
#[1] "clear-cell renal cell carcinoma (ccRCC)" "healthy"  
s1=type[which(type$group=="clear-cell renal cell carcinoma (ccRCC)"),]
s2=type[which(type$group=="healthy"),]

write.table(s1$sample,file="sample.CRC.txt",quote = F,sep="\t",row.names=F,col.names = F)
write.table(s2$sample,file="sample.healthy.txt",quote = F,sep="\t",row.names=F,col.names = F)



####基因注释
GPL14951<-getGEO('GPL14951',destdir =".")  ##根据GPL号下载的是芯片设计的信息, soft文件
colnames(Table(GPL14951))#查看基因symbol在哪一列
genename <- Table(GPL14951)[,c(1,12)]#获得ID和Gene Symbol列（第11列）
head(genename)#ID和基因symbol对应关系

#提取表达量矩阵
exprSet <- as.data.frame(exprs(gse[[1]]))# 得到表达矩阵，行名为ID，需要转换
exprSet[1:3,1:3]

#转换ID为gene name
exprSet$ID <- rownames(exprSet)
express <- merge(x = genename, y = exprSet, by = "ID")
express$ID =NULL
express[1:3,1:3]
dim(express)

#####1个探针存在对应多个基因的情况，///去掉
express1<-express[!grepl(" /// ",express$Symbol),]
dim(express1)
express1[1:3,1:3]

#########多个探针对应一个基因取均值
table(duplicated(express1$Symbol)) ####查看有多少基因Symbol重复
#aggregate函数：将基因重复的取均值
dim(express1)
test1<-aggregate(x=express1[,2:(dim(express1)[2])],by=list(express1$Symbol),FUN=mean,na.rm=T) 
dim(test1)

colnames(test1)[1]<-c("geneNames")
test1[1:3,1:3]
#去掉基因symbol为空的
test2<-subset(test1,geneNames!="")

write.table(test2,file="geneMatrix.txt",quote = F,sep="\t",row.names=F)
e1=test2[,c("geneNames",s1$sample)]
e2=test2[,c("geneNames",s2$sample)]

write.table(e1,file="geneMatrix.CRC.txt",quote = F,sep="\t",row.names=F,col.names = F)
write.table(e2,file="geneMatrix.healthy.txt",quote = F,sep="\t",row.names=F,col.names = F)
############################################
#健康与癌症样本内参基因验证
# 1. 设置工作目录与加载必要的包 
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\独立验证集GEO\\GSE72922")

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggsci) # 提供 Nature 风格配色

# 2. 读取前面已生成的表达矩阵和样本信息
exp_data <- read.delim("geneMatrix.txt", header = TRUE, stringsAsFactors = FALSE)
pheno <- read.delim("sampletype.txt", header = TRUE, stringsAsFactors = FALSE)

# 将矩阵整理为标准格式 (行名为基因，列名为样本)
rownames(exp_data) <- exp_data$geneNames
exp_data <- exp_data[, -1]

# 统一分组名称 (将原名称简化为 Cancer 和 Healthy)
pheno$group <- ifelse(pheno$group == "healthy", "Healthy", "Cancer")
# 确保分组因子顺序，Healthy 在前作为对照
pheno$group <- factor(pheno$group, levels = c("Healthy", "Cancer"))

# 确保表达矩阵的列和 phenotype 的行对应
common_samples <- intersect(colnames(exp_data), pheno$sample)
exp <- exp_data[, common_samples]
pheno <- pheno[match(common_samples, pheno$sample), ]
group <- pheno$group

# 3. 定义基因集
gene_novel <- c("BRWD1", "EIF4A1", "HGS", "NUDT16", "PDCD5", "PSMD2", "PTEN", 
                "RGPD6", "SF1", "SH3BGRL2", "SMN1", "TCERG1", "TJP2", "TMEM123", "ZBED5-AS1")
gene_trad <- c("GAPDH", "HMBS", "SDHA", "YWHAZ", "NOP10", "OST4", "SNRPG", "TOMM7")

# 过滤出芯片矩阵中实际存在的基因 (防止因为某些基因不存在导致循环报错)
gene_novel_exist <- intersect(gene_novel, rownames(exp))
gene_trad_exist <- intersect(gene_trad, rownames(exp))
gene_all <- c(gene_novel_exist, gene_trad_exist)

# 4. 执行 Wilcoxon 秩和检验
res <- data.frame()

for(g in gene_all) {
  # 提取该基因在 Healthy 和 Cancer 中的表达值
  exp_healthy <- as.numeric(exp[g, group == "Healthy"])
  exp_cancer <- as.numeric(exp[g, group == "Cancer"])
  
  # 执行检验
  test_res <- wilcox.test(exp_healthy, exp_cancer)
  
  # 汇总结果
  res <- rbind(
    res,
    data.frame(
      Gene = g,
      Pvalue = test_res$p.value,
      Panel = ifelse(g %in% gene_novel_exist, "Novel Panel (15-gene)", "Traditional Panel"),
      Status = ifelse(test_res$p.value > 0.05, "Stable (P > 0.05)", "Unstable (P < 0.05)")
    )
  )
}

# 计算 -log10(P-value) 用于绘图
res$negLog10P <- -log10(res$Pvalue)
res$Gene <- factor(res$Gene, levels = res$Gene[order(res$Panel, res$negLog10P)])

# 打印显著性结果查看
print(res)

# 5. 数据转换：为箱线图准备长数据格式 (Long format)
exp_subset <- exp[gene_all, ]
exp_subset$Gene <- rownames(exp_subset)
exp_long <- pivot_longer(exp_subset, cols = -Gene, names_to = "sample", values_to = "Expression")
exp_long <- merge(exp_long, pheno, by = "sample")
exp_long$Panel <- ifelse(exp_long$Gene %in% gene_novel_exist, "Novel Panel (15-gene)", "Traditional Panel")

# =========================================================================
# 6. Nature Style 绘图并将所有图像合并进一个 PDF 文件中
# =========================================================================

pdf(file = "Figure_F_Healthy_vs_Cancer_Validation.pdf", 
    width = 6, height = 4, onefile = TRUE)

# 图表 1：P-value 摘要柱状图 (直观展示哪些基因 P > 0.05)
p1 <- ggplot(res, aes(x = Gene, y = negLog10P, fill = Status)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
  facet_grid(. ~ Panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("Stable (P > 0.05)" = "#4DBBD5FF", "Unstable (P < 0.05)" = "#E64B35FF")) +
  theme_classic(base_size = 12, base_family = "sans") +
  labs(x = "", y = "-Log10 (P-value)", title = "Stability of Reference Genes (Healthy vs ccRCC)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.background = element_rect(fill = "#E5E5E5", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  annotate("text", x = 1, y = -log10(0.05) + 0.2, label = "P = 0.05", color = "red", hjust = 0, size = 4)

print(p1)

# 图表 2：Novel Panel (15-gene) 的表达量箱线图
p2 <- ggplot(subset(exp_long, Panel == "Novel Panel (15-gene)"), 
             aes(x = group, y = Expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
  scale_fill_npg() + # 使用 ggsci 的 Nature 系列配色
  theme_classic(base_size = 12, base_family = "sans") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3, hjust = 0.5) +
  labs(x = "", y = "Expression Level", title = "Novel 15-Gene Panel Expression") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold.italic"), # 基因名斜体加粗
    legend.position = "bottom"
  )

print(p2)

# 图表 3：Traditional Panel (传统内参) 的表达量箱线图
p3 <- ggplot(subset(exp_long, Panel == "Traditional Panel"), 
             aes(x = group, y = Expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_fill_npg() +
  theme_classic(base_size = 12, base_family = "sans") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4, label.y.npc = "top") + # 这里用星号表示显著性
  labs(x = "", y = "Expression Level", title = "Traditional Reference Genes Expression") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "#F4F4F4", color = NA),
    strip.text = element_text(face = "bold.italic"),
    legend.position = "bottom"
  )

print(p3)

dev.off()
#############################################################
####GSE72922数据-候选与别人报道内参基因稳定性分析
# =========================================================================
# 准备工作与数据读取
# =========================================================================
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\独立验证集GEO\\GSE72922")

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci) # 提供 Nature 风格配色

# 1. 读取前面已生成的表达矩阵和样本信息
exp_data <- read.delim("geneMatrix.txt", header = TRUE, stringsAsFactors = FALSE)
pheno <- read.delim("sampletype.txt", header = TRUE, stringsAsFactors = FALSE)

# 统一分组名称，与前面的分析保持一致
pheno$group <- ifelse(pheno$group == "healthy", "Healthy", "Cancer")

# 矩阵整理 (行名为基因，列名为样本)
rownames(exp_data) <- exp_data$geneNames
exp <- exp_data[, -1]

# 提取公共样本并对齐
common_samples <- intersect(colnames(exp), pheno$sample)
exp.ref <- exp[, common_samples]
pheno.ref <- pheno[match(common_samples, pheno$sample), ]

# =========================================================================
# Step 1: 定义基因集并提取矩阵
# =========================================================================
gene.new <- c("BRWD1","EIF4A1","HGS","NUDT16","PDCD5","PSMD2",
              "PTEN","RGPD6","SF1","SH3BGRL2","SMN1","TCERG1",
              "TJP2","TMEM123","ZBED5-AS1")

gene.old <- c("GAPDH","HMBS","SDHA","YWHAZ",
              "NOP10","OST4","SNRPG","TOMM7")

gene.all <- c(gene.new, gene.old)

# 过滤出矩阵中实际存在的基因 (防止报错)
gene.all.exist <- intersect(gene.all, rownames(exp.ref))

# 提取目标基因的表达矩阵并转为数值型
exp.ref <- exp.ref[gene.all.exist, ]
exp.ref[] <- lapply(exp.ref, as.numeric)

# =========================================================================
# Figure A. Expression abundance comparison (表达丰度比较)
# =========================================================================
# 计算平均表达量
mean.exp <- apply(exp.ref, 1, mean)

mean.df <- data.frame(
  Gene = names(mean.exp),
  MeanTPM = mean.exp,
  Group = ifelse(names(mean.exp) %in% gene.new, "New candidates", "Published controls")
)

# Nature 风格画图 A
pA <- ggplot(mean.df, aes(x = reorder(Gene, MeanTPM), y = MeanTPM, fill = Group)) +
  geom_bar(stat = "identity", color = "black", width = 0.75) +
  coord_flip() +
  scale_fill_npg() + # 使用 Nature 杂志经典配色 (红蓝)
  theme_classic(base_size = 12) +
  labs(x = "Reference Genes", y = "Mean Expression", title = "Expression Abundance (GSE72922)") +
  theme(
    axis.text.y = element_text(face = "italic", color = "black"), # 基因名斜体
    axis.text.x = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(file = "FigureA.Expression_abundance_comparison_GSE72922.pdf", plot = pA, width = 5, height = 5)


# =========================================================================
# Figure B. Stability comparison (整体稳定性评价，CV越低越稳定)
# =========================================================================
# 定义 CV 计算函数
cv.fun <- function(x){
  sd(x) / mean(x) * 100
}

# 计算所有样本的整体 CV
cv.value <- apply(exp.ref, 1, cv.fun)

cv.df <- data.frame(
  Gene = names(cv.value),
  CV = cv.value,
  Group = ifelse(names(cv.value) %in% gene.new, "New candidates", "Published controls")
)

# Nature 风格画图 B
pB <- ggplot(cv.df, aes(x = reorder(Gene, -CV), y = CV, fill = Group)) + # 注意此处 -CV 让最稳定的在最上面
  geom_bar(stat = "identity", color = "black", width = 0.75) +
  coord_flip() +
  scale_fill_npg() +
  theme_classic(base_size = 12) +
  labs(x = "Reference Genes", y = "Coefficient of Variation (CV %)", title = "Overall Stability (GSE72922)") +
  theme(
    axis.text.y = element_text(face = "italic", color = "black"),
    axis.text.x = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(file = "FigureB.Stability_comparison_GSE72922.pdf", plot = pB, width = 5, height = 5)


# =========================================================================
# Figure C. Intra-group variability (组内变异度，替代原纵向变异)
# =========================================================================
long.df <- data.frame()

for(g in rownames(exp.ref)) {
  # 构建包含表达量和分组信息的临时数据框
  tmp <- data.frame(
    Sample = colnames(exp.ref),
    Expression = as.numeric(exp.ref[g, ]),
    group = pheno.ref$group
  )
  
  # 按 Healthy/Cancer 分组计算组内 CV
  intra.cv <- tmp %>%
    group_by(group) %>%
    summarise(
      CV = sd(Expression) / mean(Expression) * 100,
      .groups = 'drop'
    )
  
  # 取两组的平均 CV 作为该基因的组内平均变异度
  long.df <- rbind(
    long.df,
    data.frame(
      Gene = g,
      MeanIntraCV = mean(intra.cv$CV),
      Group = ifelse(g %in% gene.new, "New candidates", "Published controls")
    )
  )
}

# Nature 风格画图 C
pC <- ggplot(long.df, aes(x = reorder(Gene, -MeanIntraCV), y = MeanIntraCV, fill = Group)) +
  geom_bar(stat = "identity", color = "black", width = 0.75) +
  coord_flip() +
  scale_fill_npg() +
  theme_classic(base_size = 12) +
  labs(x = "Reference Genes", y = "Mean Intra-group CV (%)", title = "Intra-group Variability (GSE72922)") +
  theme(
    axis.text.y = element_text(face = "italic", color = "black"),
    axis.text.x = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(file = "FigureC.Intragroup_variability_GSE72922.pdf", plot = pC, width = 5, height = 5)

print("All figures successfully generated in Nature style!")
