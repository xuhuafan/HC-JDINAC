# read the FTIR spectrum data of handlebar grips.
f1=read.csv("Edata.csv")
data=data.frame(na.omit(f1))
pdim=ncol(data)-1
x_index=rep(0,pdim)
for(i in 1:pdim){
  x_index[i]=as.numeric(substring(colnames(data)[i],2,)) ##wavenumbers
}
all_y=unique(data$y) ##sample labels


# run HC for all sample pairs：
r=NULL
for(j in all_y){
  for(k in all_y[all_y>j]){
    r0=HC(j,k,data) 
    r0=c(j,k,r0)
    r=append(r,r0)
  }
}
rm=matrix(r,ncol=pdim+2,byrow=TRUE)
colnames(rm)=c("S1","S2",colnames(data)[-(pdim+1)])
df=data.frame(rm)
#write.csv(df,"HC feature.csv") 


# compute FDR of feature selection
df=read.csv("HC feature.csv")
y1=13;y2=16
fsi=df[df$S1==y1 & df$S2==y2,-c(1,2)]
fs=x_index[fsi==1]
rate=length(fs)/pdim
##range of the valid feature according to expert experience
range=c(700:800,855:895,949:989,1000:1040,1053:1093,1103:1143,1220:1344,1400:1500,
        1565:1605,1610:1650,1702:1742,1774:1814,2493:2533,2835:2875,2905:2945,2938:2978)
FDR=length(setdiff(fs,range))/length(fs)


# run JDINAC
y1=13;y2=16
vs=diff_fea(y1,y2,df) # differential features set
EDGE=edge_construct(vs) # candidate interaction set
## input data
x_ind=ifelse(x_index %in% vs,1,0)
x_ind[pdim+1]=1
sam_train=data[data$y %in% y_ind,x_ind==1]
sam_test=sam_train
datafit=sam_train  
datapre=sam_test  
classlabel_train=ifelse(datafit$y %in% y1,0,1) 
classlabel_test=ifelse(datapre$y %in% y1,0,1) 
datafit=datafit[,-ncol(sam_train)]   
datapre=datapre[,-ncol(sam_test)]   
## bandwidth
sam_0=data[data$y %in% y1,x_ind==1]  
bws_0=apply(sam_0[,-ncol(sam_0)],2,bw.nrd0)
sam_1=data[data$y %in% y2,x_ind==1]
bws_1=apply(sam_1[,-ncol(sam_1)],2,bw.nrd0)

source('HC-JDINAC.R')
difnet=jdinac.h(h0=bws_0,h1=bws_1,EDGE=EDGE,classLabel=classlabel_train,DataFit=datafit,DataPre=datapre,nsplit=10,nfolds=5) 
ypre=ifelse(difnet$yPre>0.5,1,0) ##predicted class
sum(ypre!=classlabel_test)/nrow(datapre)   #classification error
ypre=data.frame(ypre)
ypre['true_class']=classlabel_test
ypre['label']=sam_test$y
Vars=data.frame(difnet$Vars)  #coefficients
Vars['v1']=vs[Vars[,3]]  
Vars['v2']=vs[Vars[,4]]  
eset=difnet$Eset  #differential interaction
eset=data.frame(eset)
eset['v1']=vs[eset[,1]]
eset['v2']=vs[eset[,2]] 

IA=nrow(EDGE)  #number of candidate interactions
egs=nrow(eset) # number of interactions/edges selected by JDINAC
re=egs/IA  #screening efficiency
vf=sort(unique(c(eset[,"v1"],eset[,"v2"])))
vn=length(vf) #number of nodes
fdr=length(unique(c(which((eset$v1 %in% range)==0),which((eset$v2 %in% range)==0))))/egs


## compute the sign of edges
dt=data[data$y %in% y_ind,c(which(x_index %in% vf),ncol(data))]
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
  cd[cd==0]=1
  cor=c(cor,cd)
}  
## plot a graph to dispaly the sign of edges
vf=c(1126,1740,1281,1728,1724,2920,1277,1736,1743,2912)
dd=dm[dm[,3] %in% vf,]
Sample=c(rep("E13",length(vf)),rep("E16",length(vf)))
df=cbind(Sample,rep(dd[,3],2),rbind(dd[,1],dd[,2]))
df=as.data.frame(df)
colnames(df)=c("Sample","VF","Abs")
df$VF=as.character(df$VF)
df$Abs=as.numeric(df$Abs)
df$VF<- factor(df$VF,levels = c("1126","1740","1281","1728","1724","2920","1277","1736","1743","2912"),ordered = TRUE )
d1=df[df$VF %in% c("1126","1740"),]
d2=df[df$VF %in% c("1281","1728"),]
d3=df[df$VF %in% c("1724","2920"),]
d4=df[df$VF %in% c("1277","1736"),]
d5=df[df$VF %in% c("1743","2912"),]
p=ggplot()+ 
  geom_point(aes(x = VF,y = Abs,color=Sample,shape=Sample),data = df,size=2)+
  xlab(expression(Wavenumbers(cm^-1)))+ylab("Absorbance")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),legend.position = "right")+
  theme(text=element_text(size=12, family="serif"))+
  geom_line(aes(x = VF,y = Abs,group=Sample,color=Sample),size=1.2,data = d1)+
  geom_line(aes(x = VF,y = Abs,group=Sample,color=Sample),size=1.2, data= d2)+
  geom_line(aes(x = VF,y = Abs,group=Sample,color=Sample),size=1.2,data = d3)+
  geom_line(aes(x = VF,y = Abs,group=Sample,color=Sample),size=1.2,data = d4)+
  geom_line(aes(x = VF,y = Abs,group=Sample,color=Sample),size=1.2,data = d5)+
  geom_vline(xintercept = seq(2.5,8.5,2),linetype="dashed",color="grey")+
  geom_vline(xintercept = 6.5,color="DarkGrey")+
  annotate("text",x = 3.5, y = 0.3, label = 'positive',size=3.5,color="LightSteelBlue4")+
  annotate("text",x = 8.5, y = 0.3, label = 'negative',size=3.5,color="LightSteelBlue4")


# plot differential interaction network
library(igraph)
library(ggraph)
library(tidygraph)
ed=eset[,c("v1","v2")]  #dataframe for construct network
nodes=as_data_frame(nw,what="vertices")
degs=degree(nodes)
nw=graph_from_data_frame(ed,directed=FALSE)
graph_gt <- as_tbl_graph(nw)
layout=create_layout(graph_gt, layout = 'stress')
p1=ggraph(layout)+geom_edge_link(color="gray") + geom_node_point(size=6.3,
            shape=21,aes(fill=ifelse(degs>=5,"red","steelblue")),alpha=0.2)+ 
  geom_node_text(aes(label=name),size=2)+scale_color_discrete()+
  scale_edge_width(range=c(0.2,3))+theme_graph()+
  theme(legend.position = "none")
p1


# stability of coefficients
y1=13;y2=16
y_ind=c(y1,y2)
vs=c(1269,1273,1728,1732,1736,1740)
EDGE=edge_construct(vs)
x_ind=ifelse(x_index %in% vs,1,0)
x_ind[pdim+1]=1
sam_train=data[data$y %in% y_ind,x_ind==1]
sam_test=sam_train
datafit=sam_train  
datapre=sam_test  
classlabel_train=ifelse(datafit$y %in% y1,0,1) 
classlabel_test=ifelse(datapre$y %in% y1,0,1) 
datafit=datafit[,-ncol(sam_train)]   
datapre=datapre[,-ncol(sam_test)]   
## bandwidth
sam_0=data[data$y %in% y1,x_ind==1]  
bws_0=apply(sam_0[,-ncol(sam_0)],2,bw.nrd0)
sam_1=data[data$y %in% y2,x_ind==1]
bws_1=apply(sam_1[,-ncol(sam_1)],2,bw.nrd0)
source('HC-JDINAC.R')
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
## plot dumbbell of coefficients
library(ggplot2)
dbc=ggplot(aes(x=LCI,xend=UCI,y=interaction),data=CI)+
  geom_dumbbell(colour_x = "#FFB6C1",colour_xend = "#4169E1",
                size_x = 2,size_xend = 2,size=1,color="#66CDAA")+
  geom_point(aes(x=mean,y=interaction),color="#EE7621",size=2)+
  theme_bw()+
  xlab("Coefficients")+
  ylab("Interactions")
dbc


# density plot
y_ind=c(30,32)
fs=c(2846,2927,1365,1462)
x_ind=ifelse(x_index %in% fs,1,0)
x_ind[pdim+1]=1
dt=data[data$y %in% y_ind,x_ind==1]
Sample=c(rep("E30",15),rep("E32",15))
dt["Sample"]=Sample
dtp=gather(dt, X, abs, X1365:X2927) 
dtp$X=factor(dtp$X,levels=c("X2846","X2927","X1365","X1462"),ordered = T,labels=c("2846 (S)","2927 (S)","1365 (US)","1462 (US)"))
p1<-ggplot(data=dtp,aes(x=abs))+  
  geom_density( aes(fill=Sample),alpha=0.4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=12,  family="serif"),legend.position = "right")+
  xlab("Absorbance")+
  ylab("Density")+
  facet_grid(X~.)

p2=ggplot(data=dt,aes(x=X2846, y=X2927,group=Sample)) +
  geom_point(size=0.2,aes(color=Sample)) +
  stat_density_2d(aes(alpha = ..level.., fill = Sample), geom = "polygon",)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),legend.position = "right",
        legend.key.size = unit(10, "pt"))+
  annotate("text",x = 0.07, y = 0.28, label = 'Selected\n interaction',size=3,color="LightSteelBlue4")+
  theme(text=element_text(size=12, family="serif"))


# community analysis
c0=c(3,13,26,29);c1=c(16,21,30,44)
y_ind=c(c0,c1)
## eset13_16: differential interaction of E13,E16
vs=sort(unique(c(eset13_16[,"v1"],eset13_16[,"v2"]))) 
rs=c();cs=c()
for(i in 1:nrow(eset13_16)){
  r=which(vs==eset13_16[i,"v1"])
  rs=c(rs,r)
  cl=which(vs==eset13_16[i,"v2"])
  cs=c(cs,cl)
  
}
EDGE=cbind(rs,cs)
x_ind=ifelse(x_index %in% vs,1,0)
x_ind[pdim+1]=1
sam_train=data[data$y %in% y_ind,x_ind==1]
sam_test=sam_train
datafit=sam_train  
datapre=sam_test  
classlabel_train=ifelse(datafit$y %in% c0,0,1) 
classlabel_test=ifelse(datapre$y %in% c0,0,1) 
datafit=datafit[,-ncol(sam_train)]   
datapre=datapre[,-ncol(sam_test)]   
## bandwidth
sam_0=data[data$y %in% c0,x_ind==1]  
bws_0=apply(sam_0[,-ncol(sam_0)],2,bw.nrd0)
sam_1=data[data$y %in% c1,x_ind==1]
bws_1=apply(sam_1[,-ncol(sam_1)],2,bw.nrd0)
source('HC-JDINAC.R')
difnet=jdinac.h(h0=bws_0,h1=bws_1,EDGE=EDGE,classLabel=classlabel_train,DataFit=datafit,DataPre=datapre,nsplit=10,nfolds=5) 
ypre=ifelse(difnet$yPre>0.5,1,0)
sum(ypre!=classlabel_test)/nrow(datapre)   #error
eset=difnet$Eset
eset=data.frame(eset)
eset['v1']=vs[eset[,1]]
eset['v2']=vs[eset[,2]] 
Vars=data.frame(difnet$Vars) #coefficients
Vars['v1']=vs[Vars[,3]]  
Vars['v2']=vs[Vars[,4]]  
ypre=data.frame(ypre)
ypre['true_class']=classlabel_test
ypre['label']=sam_test$y
ypre=cbind(ypre,difnet$yPre,difnet$preY)
colnames(ypre)=c("pred_class","true_class","label","pred",paste0("prob",1:20))
## compute confusing rate
prob=ypre[,-c(1:4)]
nums=rep(0,nrow(prob))
for(i in 1:ncol(prob)){
  nums=nums+(abs(prob[,i]-0.5)<0.1)
}
CR=sum(nums)/(nrow(ypre)*20)