# =============================================================================
# Sankey Diagram for Sample Distribution
# Description: Sankey diagram showing sample distribution across time points and donors
# Input: sample_info.txt
# Output: Sankey diagram (PDF/PNG)
# =============================================================================

list.files()
library(ggalluvial)
library(ggplot2)
library(dplyr)
#data
phenotype<-read.delim("sample.type.txt",header=T,row.names = NULL)
head(phenotype);dim(phenotype)
phenotype
#分组计算频数
LIHCData <- group_by(phenotype,Donor,Sex,Time) %>% summarise(., count = n())
#查看宽数据格式
head(LIHCData)
LIHCData$Donor<-as.factor(LIHCData$Donor)
ggplot(as.data.frame(LIHCData),
       aes(axis1 = Donor, axis2 = Sex, axis3 = Time,
           y= count)) +
  scale_x_discrete(limits = c("Donor", "Sex", "Time"), 
                   expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Donor)) +
  geom_stratum()+
  geom_label(stat="stratum",aes(label=after_stat(stratum))) +
  theme_minimal() +
  ggtitle("Urine EV Samples for RNA-Seq")

ggplot(as.data.frame(LIHCData),aes(y=count,axis1 = Donor, axis2 = Sex, axis3 = Time))+
  geom_alluvium(aes(fill=Donor))+
  geom_stratum(width=1/6,fill="black",color="grey")+
  geom_label(stat="stratum",aes(label=after_stat(stratum)))

#to_lodes_form生成alluvium和stratum列，主分组位于key列中
LIHC_long <- to_lodes_form(data.frame(LIHCData),
                           key = "Demographic",
                           axes = 1:3)
head(LIHC_long)

# 绘制桑基图
ggplot(data = LIHC_long,
       aes(x = Demographic, stratum = stratum, alluvium = alluvium,
           y = count, label = stratum)) +
  geom_alluvium(aes(fill = stratum)) +
  geom_stratum() + geom_text(stat = "stratum") +
  theme_minimal() +
  ggtitle("Urine EV Samples for RNA-Seq")

library(RColorBrewer)
display.brewer.all()
col <- colorRampPalette(brewer.pal(8,'Accent'))(8)
col2<-colorRampPalette(brewer.pal(12,'Paired'),alpha = TRUE)(12)
mycol<-c(col2,col)
#scale_fill_manual(values=c(mycol))

p<-ggplot(data = LIHC_long,
       aes(x = Demographic, stratum = stratum, alluvium = alluvium,
           y = count, label = stratum)) +
  geom_alluvium(aes(fill = stratum)) +
  geom_stratum() + geom_text(stat = "stratum") +
  theme_minimal() +
  ggtitle("Urine EV Samples for RNA-Seq")+
  scale_fill_manual(values=c(mycol))
p
ggsave(filename="样本分布桑基图.new2.pdf",plot=p,
       device='pdf',path=".",width=5.5,height=5.5)








library(ggthemes)
library(ggsci)
g=ggplot(LIHC_long, aes(x = Demographic,y=count,fill =stratum, 
                  stratum = stratum, alluvium = alluvium,label = stratum)) +
  geom_col(width = 0.4,color=NA)+
  geom_flow(width = 0.4,alpha = 0.2,knot.pos = 0)   +#knot.pos可以使连线更直
  #geom_alluvium( width = 0.4,alpha = 0.2,knot.pos = 0)+ 与geom_flow效果相似
  scale_fill_manual(values =c(mycol))+
  theme_map()+
  theme(axis.text.x=element_text(size=20,vjust = 5),
        legend.position = 'none')

g=ggplot(LIHC_long, aes(x = Demographic,y=count,fill =stratum, 
                        stratum = stratum, alluvium = alluvium,label = stratum)) +
  geom_col(width = 0.4,color=NA)+
  geom_flow(width = 0.4,alpha = 0.2,knot.pos = 0)   +#knot.pos可以使连线更直
  #geom_alluvium( width = 0.4,alpha = 0.2,knot.pos = 0)+ 与geom_flow效果相似
  scale_fill_manual(values =c(mycol))+
  geom_text(stat = "stratum")+ theme_map()+
  ggtitle("Urine EV Samples for RNA-Seq")

g
ggsave(filename="样本分布桑基图new.pdf",plot=g,
       device='pdf',path=".",width=5.5,height=5.5)
