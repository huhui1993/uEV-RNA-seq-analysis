setwd("D:\\项目\\231021661.报告--尿液外泌体转录组\\自己重新跑后\\D1")
list.files(pattern = "*.txt")
stat.count<-read.delim("all_result.D1.sample.low.medium.high.gene.count.txt",
                       header=T,row.names = NULL)
library("ggpubr")
library("reshape")
library("ggplot2")
library(dplyr)

stat.count

b <- stat.count %>% mutate(low.percent=low/sum(low)*100)
b <- b %>% mutate(medium.percent=medium/sum(medium)*100)
b <- b %>% mutate(high.percent=high/sum(high)*100)
b
#b <-  mutate_if(b, is.numeric, round, digit=3)
b.1<-b[,c("Sample","low.percent","medium.percent","high.percent")] 

data_melt<-melt(b.1,id=c("Sample"))
head(data_melt)
colnames(data_melt)<-c("Sample","type","percentage")
head(data_melt)


colors <- colorRampPalette(c("darkblue","DodgerBlue", "skyblue"))(3)

colors <- colorRampPalette(c("#41F0AE","#FFC080", "#FF8080"))(3)
p6 <- ggplot(data_melt, aes(x=2, y=percentage, fill=type)) +
  geom_bar(position = 'fill', stat = 'identity') +
  facet_wrap(Sample~., strip.position="bottom",ncol = 6) + # ncol= 1, strip.position="right"
  theme_void() +
  scale_fill_manual(values=setNames(colors, levels(data_melt$type))) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 12)) +
  coord_polar(theta = 'y') +
  xlim(0.5, 2.5) +
  labs(x=NULL, y=NULL,fill=NULL)
p6
ggsave(filename="./D1.pie.exp.3type.plot.pdf",plot=p6,
       device='pdf',path=".",width=5,height=3)

