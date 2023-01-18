## read the FTIR spectrum data of handlebar grips.
f1=read.csv("Edata.csv")
f1=f1[,-1]
f1=f1[,67:937] #650-4000 cm-1
data=data.frame(na.omit(f1))
pdim=ncol(data)-1
x_index=rep(0,pdim)
for(i in 1:pdim){
  x_index[i]=as.numeric(substring(colnames(data)[i],2,))
}
all_y=c(2:7,9:30,32,34:37,40,41,44:49) ## sample selection
x=filter(data,data$y %in% all_y)
colnames(x)=x_index

## run HC for all sample pairs
r=NULL
for(j in all_y){
  for(k in all_y[all_y<j]){
    r0=HC(j,k,x) 
    r0=c(j,k,r0)
    r=append(r,r0)
  }
}
rm=matrix(r,ncol=pdim+2,byrow=TRUE)
colnames(rm)=c("S1","S2",x_index)
df=data.frame(rm)
#write.csv(df,"HC features.csv") 

##run JDINAC
y1=13;y2=16
y_ind=c(y1,y2)
vs=NULL 
for(i in c(1:length(y_ind))){
  num=y_ind[i]
  sam=df[df$S1==num|df$S2==num,]
  c=colSums(sam[,3:(pdim+2)])  
  c=data.frame(index,c)
  all=c[c$c>0,1] 
  vs=c(vs,all)
}
vs=unique(vs)
vs=sort(vs,decreasing = FALSE)
vs=vs[vs!=4000]
x_ind=ifelse(x_index %in% vs,1,0)
x_ind[pdim+1]=1
EDGE=edge_construct(vs)

sam_train=NULL
for(i in y_ind){
  sam=data[data$y==i,]
  train=sam[,x_ind==1]
  sam_train=rbind(sam_train,train)
}
dim(sam_train)  #30*80；30*7
sam_test=sam_train

## classlabel
datafit=sam_train  #训练数据,含y
datapre=sam_test    #测试数据,含y
classlabel_train=ifelse(datafit$y %in% y_ind1,0,1) 
classlabel_test=ifelse(datapre$y %in% y_ind1,0,1)#
datafit=datafit[,-ncol(sam_train)]   
datapre=datapre[,-ncol(sam_test)]   

#bandwidth h0,h1
sam_0=data[data$y %in% y1,x_ind==1]  #
bws_0=NULL
for(i in 1:ncol(sam_0[,-1])){
  bw=bw.nrd0(sam_0[,i])
  bws_0=c(bws_0,bw)
}
sam_1=data[data$y %in% y,x_ind==1]
bws_1=NULL
for(i in 1:ncol(sam_1[,-1])){
  bw=bw.nrd0(sam_1[,i])
  bws_1=c(bws_1,bw)
}

source('D:\\红外光谱数据\\代码\\固定带宽JDINAC.R')
EDGE=edge_construct(vs)
IA=nrow(EDGE)  #number of candidate interactions
difnet=jdinac.h(h0=bws_0,h1=bws_1,EDGE=EDGE,classLabel=classlabel_train,DataFit=datafit,DataPre=datapre,nsplit=10,nfolds=5) 
ypre=ifelse(difnet$yPre>0.5,1,0)
sum(ypre!=classlabel_test)/nrow(datapre)   #error
ypre=data.frame(ypre)
ypre['true_class']=classlabel_test
ypre['label']=sam_test$y
ypre=cbind(ypre,difnet$yPre,difnet$preY)
Vars=data.frame(difnet$Vars)
Vars['v1']=vs[Vars[,3]]  
Vars['v2']=vs[Vars[,4]]  
eset=difnet$Eset
eset=data.frame(eset)
eset['v1']=vs[eset[,1]]
eset['v2']=vs[eset[,2]] 
egs=nrow(eset) # number of interactions selected by JDINAC
re=egs/IA  #screening efficiency
vf=sort(unique(c(eset[,"v1"],eset[,"v2"])))
vn=length(vf) #number of nodes

#JDINAC选中数据
dt=data[data$y %in% c(y1,y2),c(which(x_index %in% vf),ncol(data))]
med1=colMeans(dt[dt$y %in% y1,-ncol(dt)])
med1=colMeans(dt[dt$y %in% y2,-ncol(dt)])
dm=cbind(med1,med2,vf)
colnames(dm)=c(paste0("E",y1),paste0("E",y2),"VF")
cor=c()
for(e in 1:egs){
  v1=eset[e,"v1"];v2=eset[e,"v2"]
  r1=which(dm[,3]==v1)
  r2=which(dm[,3]==v2)
  prod=(dm[r1,1]-dm[r2,1])*(dm[r1,2]-dm[r2,2])
  cd=sign(prod)
  cor=c(cor,cd)
}  


## differential interaction network
library(igraph)
library(ggraph)
library(tidygraph)
#file="D:\\红外光谱数据\\JDINAC补充实验\\个体差异网络\\30&32.xlsx"
#eset=read.xlsx(file,sheetName = 'eset')
ed=eset[,c("v1","v2","numb")]  #构建网络的dataframe
nodes=as_data_frame(nw,what="vertices")
degs=degree(nodes)
ed=filter(ed,  !is.na(v2))
nw=graph_from_data_frame(ed,directed=FALSE)
graph_gt <- as_tbl_graph(nw)
layout=create_layout(graph_gt, layout = 'stress')
p1=ggraph(layout)+geom_edge_link(color="gray") + geom_node_point(size=6.3,shape=21,alpha=0.2)+ 
  geom_node_text(aes(label=name),size=2)+scale_color_discrete()+scale_edge_width(range=c(0.2,3))+theme_graph()+
  theme(legend.position = "none")
p1

## bandiwdth sensitivity
## table 2
vs=c(1269,1273,1728,1732,1736,1740)
EDGE=edge_construct(vs)
source('D:\\红外光谱数据\\代码\\固定带宽JDINAC.R')
difnet=jdinac.h(h0=bws_0,h1=bws_1,EDGE=EDGE,classLabel=classlabel_train,DataFit=datafit,DataPre=datapre,nsplit=30,nfolds=5) 
ypre=ifelse(difnet$yPre>0.5,1,0)
sum(ypre!=classlabel_test)/nrow(datapre)   #error
eset=difnet$Eset
eset=data.frame(eset)
eset['v1']=vs[eset[,1]]
eset['v2']=vs[eset[,2]] 
eset['freq']=eset['numb']/60  #frequency

Vars=data.frame(difnet$Vars) #coefficients
Vars['v1']=vs[Vars[,3]]  
Vars['v2']=vs[Vars[,4]]  
CI=eset[,c("numb","v1","v2")]
CI["mean"]="mean";CI["SD"]="SD"
for(i in 1:nrow(CI)){
  v1=CI[i,"v1"];v2=CI[i,"v2"]
  coes=Vars[Vars$v1==v1 & Vars$v2==v2,"coef"]
  pvalue=t.test(coes,0,var.equal=TRUE,alternative = "two.sided")$p.value
  CI[i,"mean"]=mean(coes)
  CI[i,"SD"]=sd(coes)
  CI[i,"p-value"]=pvalue
}

CI$mean=as.numeric(CI$mean)
CI$SD=as.numeric(CI$SD)
CI["LCI"]=CI$mean-1.96*CI$SD
CI["UCI"]=CI[,"mean"]+1.96*CI[,"SD"]
CI["interaction"]=paste0("(",CI$v1,",",CI$v2,")")
CI$interaction=factor(CI$interaction,levels=CI$interaction,ordered = TRUE)

##figure 5
p5=ggplot(aes(x=LCI,xend=UCI,y=interaction),data=CI)+
  geom_dumbbell(colour_x = "#FFB6C1",colour_xend = "#4169E1",
                size_x = 2,size_xend = 2,size=1,color="#66CDAA")+
  geom_point(aes(x=mean,y=interaction),color="#EE7621",size=2)+
  theme_bw()+
  xlab("Coefficients")+
  ylab("Interactions")
p5


## figure 7(b)

