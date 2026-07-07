library("ComplexHeatmap")

list<-colnames(df1)
all_result=data.frame()
for(i in 1:length(list)){
  t<-df1[,i]
  a<-length(t[t>0 & t<=1])
  b<-length(t[t>1 & t<100])
  c<-length(t[t>=100])
  d<-cbind(a,b,c)
  all_result<-rbind(all_result,d)
}
##########protein coding
df1<-df.pro[,phenotype.d1$sample]
head(df1)
dim(df1)
colnames(df1)
#[1] "S0102" "S0104" "S0105" "S0106" "S0107" "S0108" "S0111" "S0112" "S0113" "S0115" "S0116" "S0117"
##high exp
df1$geneid<-rownames(df1)
lt=list(
            S0102=as.factor(df1[df1[,1]>=100,13]),
            S0104=as.factor(df1[df1[,2]>=100,13]),
            S0105=as.factor(df1[df1[,3]>=100,13]),
            S0106=as.factor(df1[df1[,4]>=100,13]),
            S0107=as.factor(df1[df1[,5]>=100,13]),
            S0108=as.factor(df1[df1[,6]>=100,13]),
            S0111=as.factor(df1[df1[,7]>=100,13]),
            S0112=as.factor(df1[df1[,8]>=100,13]),
            S0113=as.factor(df1[df1[,9]>=100,13]),
           S0115=as.factor(df1[df1[,10]>=100,13]),
           S0116=as.factor(df1[df1[,11]>=100,13]),
           S0117=as.factor(df1[df1[,12]>=100,13])
)

#lt = list(a = sample(letters, 5),
 #         b = sample(letters, 10),
  #        c = sample(letters, 15))
m1 = make_comb_mat(lt)
m1
m = make_comb_mat(lt)
UpSet(m)
m2<-m[comb_degree(m) >= 10]
UpSet(m2)

m3<-m[comb_degree(m) <= 1]
UpSet(m3)

pdf(file="./D1/D1.sample.highexp.upset.pdf",7,5)
UpSet(m2)
dev.off()
###
##medium exp
df1$geneid<-rownames(df1)

lt=list(
  S0102=as.factor(df1[df1[,1]>1 & df1[,1]<100,13]),
  S0104=as.factor(df1[df1[,2]>1 & df1[,2]<100,13]),
  S0105=as.factor(df1[df1[,3]>1 & df1[,3]<100,13]),
  S0106=as.factor(df1[df1[,4]>1 & df1[,4]<100,13]),
  S0107=as.factor(df1[df1[,5]>1 & df1[,5]<100,13]),
  S0108=as.factor(df1[df1[,6]>1 & df1[,6]<100,13]),
  S0111=as.factor(df1[df1[,7]>1 & df1[,7]<100,13]),
  S0112=as.factor(df1[df1[,8]>1 & df1[,8]<100,13]),
  S0113=as.factor(df1[df1[,9]>1 & df1[,9]<100,13]),
  S0115=as.factor(df1[df1[,10]>1 & df1[,10]<100,13]),
  S0116=as.factor(df1[df1[,11]>1 & df1[,11]<100,13]),
  S0117=as.factor(df1[df1[,12]>1 & df1[,12]<100,13])
)
#lt = list(a = sample(letters, 5),
#         b = sample(letters, 10),
#        c = sample(letters, 15))
m1 = make_comb_mat(lt)
m1
m = make_comb_mat(lt)
UpSet(m)
m2<-m[comb_degree(m) >= 10]
UpSet(m2)
pdf(file="./D1/D1.sample.mediumexp.upset.pdf",9,5)
UpSet(m2)
dev.off()
##low exp
df1$geneid<-rownames(df1)
lt=list(
  S0102=as.factor(df1[df1[,1]>0 & df1[,1]<=1,13]),
  S0104=as.factor(df1[df1[,2]>0 & df1[,2]<=1,13]),
  S0105=as.factor(df1[df1[,3]>0 & df1[,3]<=1,13]),
  S0106=as.factor(df1[df1[,4]>0 & df1[,4]<=1,13]),
  S0107=as.factor(df1[df1[,5]>0 & df1[,5]<=1,13]),
  S0108=as.factor(df1[df1[,6]>0 & df1[,6]<=1,13]),
  S0111=as.factor(df1[df1[,7]>0 & df1[,7]<=1,13]),
  S0112=as.factor(df1[df1[,8]>0 & df1[,8]<=1,13]),
  S0113=as.factor(df1[df1[,9]>0 & df1[,9]<=1,13]),
  S0115=as.factor(df1[df1[,10]>0 & df1[,10]<=1,13]),
  S0116=as.factor(df1[df1[,11]>0 & df1[,11]<=1,13]),
  S0117=as.factor(df1[df1[,12]>0 & df1[,12]<=1,13])
)
#lt = list(a = sample(letters, 5),
#         b = sample(letters, 10),
#        c = sample(letters, 15))
m1 = make_comb_mat(lt)
m1
m = make_comb_mat(lt)
UpSet(m)
m2<-m[comb_degree(m) >= 10]
UpSet(m2)
pdf(file="./D1/D1.sample.lowexp.upset.pdf",9,5)
UpSet(m2)
dev.off()
##################################################3
########lncRNA
df1<-df.lnc[,phenotype.d1$sample]
head(df1)
dim(df1)

colnames(df1)
#[1] "S0102" "S0104" "S0105" "S0106" "S0107" "S0108" "S0111" "S0112" "S0113" "S0115" "S0116" "S0117"
##high exp
df1$geneid<-rownames(df1)
lt=list(
  S0102=as.factor(df1[df1[,1]>=100,13]),
  S0104=as.factor(df1[df1[,2]>=100,13]),
  S0105=as.factor(df1[df1[,3]>=100,13]),
  S0106=as.factor(df1[df1[,4]>=100,13]),
  S0107=as.factor(df1[df1[,5]>=100,13]),
  S0108=as.factor(df1[df1[,6]>=100,13]),
  S0111=as.factor(df1[df1[,7]>=100,13]),
  S0112=as.factor(df1[df1[,8]>=100,13]),
  S0113=as.factor(df1[df1[,9]>=100,13]),
  S0115=as.factor(df1[df1[,10]>=100,13]),
  S0116=as.factor(df1[df1[,11]>=100,13]),
  S0117=as.factor(df1[df1[,12]>=100,13])
)

#lt = list(a = sample(letters, 5),
#         b = sample(letters, 10),
#        c = sample(letters, 15))
m1 = make_comb_mat(lt)
m1
m = make_comb_mat(lt)
UpSet(m)
m2<-m[comb_degree(m) >= 10]
UpSet(m2)

m3<-m[comb_degree(m) <= 1]
UpSet(m3)


pdf(file="./D1/D1.lncRNA.sample.highexp.upset.pdf",9,5)
UpSet(m)
dev.off()
