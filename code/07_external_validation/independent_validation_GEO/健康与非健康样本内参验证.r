############################################
#健康与癌症样本内参基因验证
# 1. 设置工作目录与加载必要的包 
setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\EV-origin\\独立验证集GEO\\GSE125442")

# 1. 设置工作目录 (请根据您的实际路径修改)
# setwd("D:\\您的项目路径\\GSE125442")

# 加载所需 R 包
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggsci) # 提供 Nature 风格配色

# =========================================================================
# 2. 数据读取与预处理 (针对 GSE125442 特有格式)
# =========================================================================
# 读取带有注释信息的表达矩阵
raw_data <- read.delim("GSE125442_gene_exp.txt", header = TRUE, stringsAsFactors = FALSE)

# 提取样本列名 (包含 "_RPKM" 的列)
sample_cols <- grep("_RPKM", colnames(raw_data), value = TRUE)

# 构建只包含基因名和表达量的纯净矩阵
exp_data <- raw_data[, c("gene", sample_cols)]

# 去除基因名为 NA 或空的行
exp_data <- subset(exp_data, gene != "" & !is.na(gene))

# 处理可能存在的重复基因名 (取平均值)
exp_data <- aggregate(. ~ gene, data = exp_data, FUN = mean)

# 将 gene 列设置为行名
rownames(exp_data) <- exp_data$gene
exp <- exp_data[, -1]

# 提取并构建表型(分组)信息
samples <- colnames(exp)
pheno <- data.frame(sample = samples, stringsAsFactors = FALSE)

# 根据列名包含 'ctrl' 还是 'Utest' 来打标签
pheno$group <- ifelse(grepl("ctrl", pheno$sample, ignore.case = TRUE), "Healthy", "Cancer")

# 强制转换分组因子，保证 Healthy 在前 (作为基准对照)
pheno$group <- factor(pheno$group, levels = c("Healthy", "Cancer"))
group <- pheno$group

# =========================================================================
# 3. 定义基因集并取交集
# =========================================================================
gene_novel <- c("BRWD1", "EIF4A1", "HGS", "NUDT16", "PDCD5", "PSMD2", "PTEN", 
                "RGPD6", "SF1", "SH3BGRL2", "SMN1", "TCERG1", "TJP2", "TMEM123", "ZBED5-AS1")
gene_trad <- c("GAPDH", "HMBS", "SDHA", "YWHAZ", "NOP10", "OST4", "SNRPG", "TOMM7")

# 过滤出矩阵中实际检测到的基因
gene_novel_exist <- intersect(gene_novel, rownames(exp))
gene_trad_exist <- intersect(gene_trad, rownames(exp))
gene_all <- c(gene_novel_exist, gene_trad_exist)

# =========================================================================
# 4. 执行 Wilcoxon 秩和检验
# =========================================================================
res <- data.frame()

for(g in gene_all) {
  # 提取该基因在 Healthy 和 Cancer 中的表达值
  exp_healthy <- as.numeric(exp[g, group == "Healthy"])
  exp_cancer <- as.numeric(exp[g, group == "Cancer"])
  
  # 执行 Wilcoxon 检验
  test_res <- wilcox.test(exp_healthy, exp_cancer)
  
  # 汇总统计结果
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

# 计算 -log10(P-value) 用于柱状图可视化
res$negLog10P <- -log10(res$Pvalue)

# 按 P 值排序，使图表更加美观
res$Gene <- factor(res$Gene, levels = res$Gene[order(res$Panel, res$negLog10P)])

# 在控制台打印显著性结果供核对
print("=== Wilcoxon Test Results ===")
print(res)

# =========================================================================
# 5. 数据宽转长格式 (为 ggplot2 箱线图准备)
# =========================================================================
exp_subset <- exp[gene_all, ]
exp_subset$Gene <- rownames(exp_subset)

# 转换为长格式
exp_long <- pivot_longer(exp_subset, cols = -Gene, names_to = "sample", values_to = "Expression")
exp_long <- merge(exp_long, pheno, by = "sample")

# 标记属于哪个 Panel
exp_long$Panel <- ifelse(exp_long$Gene %in% gene_novel_exist, "Novel Panel (15-gene)", "Traditional Panel")

# =========================================================================
# 6. 绘制 Nature Style 图表并输出为连打 PDF
# =========================================================================
pdf(file = "Figure_F_Healthy_vs_Cancer_Validation_GSE125442.pdf", 
    width = 6, height = 4, onefile = TRUE)

# ---------------- 图表 1：P-value 摘要柱状图 ----------------
p1 <- ggplot(res, aes(x = Gene, y = negLog10P, fill = Status)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 1) +
  facet_grid(. ~ Panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("Stable (P > 0.05)" = "#4DBBD5FF", "Unstable (P < 0.05)" = "#E64B35FF")) +
  theme_classic(base_size = 12, base_family = "sans") +
  labs(x = "", y = "-Log10 (P-value)", title = "Stability of Reference Genes (GSE125442: Healthy vs ccRCC)") +
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

# ---------------- 图表 2：Novel Panel (15-gene) 箱线图 ----------------
p2 <- ggplot(subset(exp_long, Panel == "Novel Panel (15-gene)"), 
             aes(x = group, y = Expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
  scale_fill_npg() + # 引入 Nature 系列配色 (红蓝)
  theme_classic(base_size = 12, base_family = "sans") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4, hjust = 0.5) +
  labs(x = "", y = "Expression (RPKM)", title = "Novel 15-Gene Panel Expression (GSE125442)") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold.italic"),
    legend.position = "bottom"
  )

print(p2)

# ---------------- 图表 3：Traditional Panel 箱线图 ----------------
p3 <- ggplot(subset(exp_long, Panel == "Traditional Panel"), 
             aes(x = group, y = Expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.6) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_fill_npg() +
  theme_classic(base_size = 12, base_family = "sans") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4, label.y.npc = "top") +
  labs(x = "", y = "Expression (RPKM)", title = "Traditional Reference Genes Expression (GSE125442)") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "#F4F4F4", color = NA),
    strip.text = element_text(face = "bold.italic"),
    legend.position = "bottom"
  )

print(p3)

dev.off()
