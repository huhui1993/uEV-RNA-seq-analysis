# =============================================================================
# Donut Chart for Sample Composition
# Description: Multi-layer donut chart showing sample distribution
# Input: sample_info.txt
# Output: Donut chart (PDF)
# =============================================================================

list.files()
#1.1 安装相关包webr
library("devtools")
devtools::install_github("cardiomoon/moonBook")
devtools::install_github("cardiomoon/webr")
#加载软件包：
require(ggplot2)
require(moonBook)
require(webr)
#1.2 webr包的使用首先，我们来看基础饼图的绘制。
#这里我们构建一个数据，其中包含5个浏览器及其对应的市场份额。
browser=c("MSIE","Firefox","Chrome","Safari","Opera")
share=c(50,21.9,10.8,6.5,1.8)
df=data.frame(browser,share)
df
#我们只需要这样一行代码即可绘制饼图：
PieDonut(df,aes(browser,count=share))
#其中aes()指定了我们调用哪列数据，browser是分类变量，而count指向我们的share变量，是数值变量。
#我们可以微调一下。r0是核心空白圆的半径，
#start是起绘点的角度（按顺时针），
#而labelpositionThreshold指定了占比超过多少就要把标签外置，这样更美观。
PieDonut(df,aes(browser,count=share),r0=0.5,start=3*pi/2,labelpositionThreshold=0.1)
#如果你想调整饼图或者环图的半径，可以使用 r1和r2参数。
#此外，使用showPieName=FALSE可以去掉最中心的名称
PieDonut(df,aes(browser,count=share),r0=0.5,
         start=3*pi/2,labelpositionThreshold=0.1,
         showPieName=FALSE)

#我们接下来展示饼环图的绘制。
#这里，我们调用R的mtcars包来做演示。
#值得注意的是，这里我们只指定了gear和carb两个分类变量，
#而其占比是根据两个变量中每个元素出现的频数来统计的
#（也就是我们有很多行数据，有一些重复出现的数据行）。
PieDonut(mtcars,aes(gear,carb),start=3*pi/2)
#PieDonut()函数有非常多可以调整的参数，其实开发者是很用心的，
#但美中不足的是他还没有把调整配色的参数开发好，
#我们不修改源代码的话是不能调整配色的

#比较流行一种作图方式是把要突出的那块“馅饼”切下来。
#这里，我们可以使用一个explode参数来实现，同时要指定explodeDonut=TRUE，
#此外，可以通过maxx参数来调整其突出的幅度。
#explode=3意味着把第三种馅饼拿出来。
PieDonut(mtcars,aes(gear,carb),start=3*pi/2,explode=3,explodeDonut=TRUE,maxx=1.7)

#当然，我们也可以只把第3，6，9份甜甜圈拿出来。只需要增加selected参数。
#但可能由于在这个数据中分类变量全是数字代表的，在执行上存在一些问题。
#我们改为使用moonbook包中的acs数据来演示。
#这是 857 名急性冠脉综合征 (ACS) 患者的人口统计和实验室数据。
#如果你想根据诊断显示吸烟状况的分布，请使用以下代码：
PieDonut(acs,aes(Dx,smoking),explode=1,selected=c(3,6,9),explodeDonut=TRUE)

#注意，这里环图的比例是相对于上一级分类（父级分类）的，
#如果我们想展示子级分类的绝对比例，则可以增加参数ratioByGroup=FALSE。
PieDonut(acs,aes(Dx,smoking),explode=1,selected=c(3,6,9),
         explodeDonut=TRUE,ratioByGroup=FALSE)
#####
##如果是已经统计好了的数据，那么我们就需要在前面的代码基础上指定count参数。
#示例数据
mydata <- data.frame(ownership=c(rep("private", 3), rep("public",3),rep("mixed", 3)),  
                     landuse=c(rep(c("residential", "recreation", "commercial"),3)),  
                     acres=c(108,143,102, 300,320,500, 37,58,90))
mydata
PieDonut(mydata,aes(ownership,landuse,count=acres))

#####展示自己的数据
phenotype
PieDonut(phenotype,aes(Time,Sex),start=3*pi/2,ratioByGroup=FALSE)
PieDonut(phenotype,aes(Donor,Time),start=3*pi/2,ratioByGroup=FALSE)
#########################################################
if(!require(plotly)) 
  install.packages("plotly")
if(!require(sunburstR)) 
  install.packages("sunburstR")

library("sunburstR")
mydata <- data.frame(ownership=c(rep("private", 3), rep("public",3),rep("mixed", 3)),  
                     landuse=c(rep(c("residential", "recreation", "commercial"),3)),  
                     acres=c(108,143,102, 300,320,500, 37,58,90))
mydata
plotseq <- data.frame(path = paste(mydata$ownership,mydata$landuse,sep = "-"), 
                      num = mydata$acres)
sunburst(plotseq)
#修改为自定义配色
sunburst(plotseq,colors = list(range = RColorBrewer::brewer.pal(9, "Set3")))
#另一种绘图方案
sund2b(plotseq)
#我们可以使用showLabels = TRUE展示分类的标签，还可以修改根变量的标签。
sund2b(plotseq, showLabels = TRUE, rootLabel = "myplot")
