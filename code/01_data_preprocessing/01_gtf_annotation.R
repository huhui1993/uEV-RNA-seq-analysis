# =============================================================================
# 01. GTF Gene Annotation
# Description: Parse GTF file to extract gene ID, gene symbol, and biotype
# Input: GRCh38 GTF file, all_sample_count.txt
# Output: Gene annotation table with biotype classification
# =============================================================================

list.files()
library(refGenome)
ens <- ensemblGenome()
read.gtf(ens, "Homo_sapiens.GRCh38.111.gtf")###导入gtf文件 比较耗时
class(ens)
my_gene <- getGenePositions(ens)

#这里面不要用genecode的GTF文件导入会使R崩溃**
  
colnames(my_gene)  
gene_id<- my_gene[,21]
gene_name <- my_gene[,16]
gene_biotype <- my_gene[,18]
a <- cbind(gene_id,gene_name,gene_biotype)
write.table(a,file="GRCh38.111.gtf.geneID2Symbol2biotype.txt",
            quote = F,sep="\t",row.names=F)

expr<-read.delim("all_sample_count.txt",header=T,row.names = 1)
a.1<-a[rownames(expr),]  
expr.2<-cbind(a.1$gene_name,a.1$gene_biotype,expr)
write.table(expr.2,file="all_sample_count.addAnotation.txt",
            quote = F,sep="\t",row.names=T)
#statistic
colSums( expr.2 != 0)
write.table(as.data.frame(colSums( expr.2 != 0)),file="Number_of_expressed_gene_in_each_sample.txt",
            quote = F,sep="\t",row.names=F)
head(expr.2)
newname<-gsub("a.1[^[:alnum:]///' ]","",colnames(expr.2))
colnames(expr.2)<-newname

library(dplyr)  # 调包
df<-c()
tmp<-expr.2[,c(2,32)]
tmp2<-tmp[which(tmp[,2] != 0),]
df0 <- tmp2 %>% group_by(gene_biotype) %>%  # 按Category统计每个类别的数据点个数
  summarise(count=n(),.groups = 'drop') %>% as.data.frame()  # 定义为数据框类型

df<-df0 #随意赋值
numseq=seq(3,32)

for (i in 3:32) {
tmp<-expr.2[,c(2,i)]
tmp2<-tmp[which(tmp[,2] != 0),]
df0 <- tmp2 %>% group_by(gene_biotype) %>%  # 按Category统计每个类别的数据点个数
  summarise(count=n(),.groups = 'drop') %>% as.data.frame()  # 定义为数据框类型
df=merge(df,df0,by = "gene_biotype",all=T)}

df  # 查看统计结果
dim(df)
colnames(df)<-c("gene_biotype","test",colnames(expr.2)[3:32])
write.table(df,file="Number_of_expressed_gene_in_each_sample.in_each_biotype.txt",
            quote = F,sep="\t",row.names=F)
df.1<-df[,-2]
head(df.1)


library("reshape")
library("ggplot2")
data_melt<-melt(a,id=c("cluster"))
head(data_melt)
colnames(data_melt)<-c("cluster","Type","cell_percent")
head(data_melt)
p<-ggplot(data=data_melt,aes(x=Type,y=cell_percent,fill=cluster))+
  geom_bar(stat="identity")+theme_bw()
p.final<-p+
  labs(y="Cell count percentage in each cluster (%)",x="",title="")+
  theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1,
                                   vjust =0.5,colour = "black",lineheight=100),
        axis.text.y = element_text(size=13,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=13,colour="black"),
        axis.title.y = element_text(size=13,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  scale_x_discrete(limit = unique(data_melt$Cluster))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.final
library(RColorBrewer)
display.brewer.all()
col <- colorRampPalette(brewer.pal(8,'Accent'))(7)
mycol<-colorRampPalette(brewer.pal(12,'Set3'),alpha = TRUE)(13)
p.final.new<-p.final+scale_fill_manual(values=mycol)
scale_fill_manual(values=c('CornflowerBlue',
                           'LightSlateBlue','DodgerBlue',
                           'DeepSkyBlue','SteelBlue','LightSteelBlue',
                           'DarkTurquoise','CadetBlue',
                           'Aquamarine',
                           'MediumSeaGreen','Cyan'))
ggsave(filename="allsample_barplot_cell_percentage_in_cluster.pdf",plot=p.final.new,
       device='pdf',path=".",width=4.8,height=6)