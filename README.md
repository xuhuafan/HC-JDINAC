# HC-JDINAC
This package provides the R function to implement the HC-JDINAC method for micro trace evidence examination, as used in the manuscript, "HC-JDINAC: Identifying Important Individual Features and Interaction Differences from Highly Similar FTIR Spectra".
  
## HC-LDA
### Description
Perfrom HC-LDA,FAIR-LDA,SFLDA and FANS. Compute classification error (CE) and FDR like Table 3 in paper.
### Examples
```
##Given y1,y2, perform HC-LDA, FAIR-LDA, FANS and SFLDA. 
data=read.csv("Edata.csv")
y1=13;y2=16  ## perform binary classification between any two individauls by changing y1,y2
result_hc=HC_LDA(y1,y2,data,alpha=0.1)  ## HC-LDA result
result_hc$res  ##testing classification error(CE) and FDR
result_fair=FAIR_LDA(y1,y2,data,alpha=0.1)  ## FAIR-LDA result
result_fair$res  ##testing classification error(CE) and FDR
result_fans=FANS(y1,y2,data,nsplit=10,nfolds=5)  ## FANS result
result_fans$res  ##testing classification error(CE) and FDR
result_sflda=SFLDA(y1,y2,data)  ## SFLDA result
result_sflda$res  ##testing classification error(CE) and FDR
```

## JDINAC with binned bivariate kernel density estimation.
### Description
Perform JDINAC with binned bivariate kernel density estimation using kernel $K_\textbf{h}$ applied to the rescaling of the product kernel $K(\textbf{x}) = K(x_1)K(x_2)$ by the vector of bandwidths $\textbf{h} = (h_1, h_2)$:
$$
{\color{orange}{
K_\textbf{h}(\textbf{x})=K_{h_1}(x_1)K_{h_2}(x_2),}}
$$
where $K$ is a symmetric probability density function and $K_h(x) = K(x/h)/h$ is a rescaling of $K$ by the bandwidth $h>0$, $\textbf{x}=(x_1,x_2)$, $x_1$ and $x_2$ are different features.
### Examples
```
##perform HC-JDINAC on given two individuals y1,y2
y1=13;y2=16
data=read.csv("individual.csv")
dfh=read.csv("HC feature.csv")  ##HC feature selection result
result=jdinac_bd(y1,y2,data,dfh,k0=1,k1=1)
eset=result$eset ## differential edge/interaction selected by JDINAC
err=result$err ## classification error
p=diffnet_plot(eset) ## differential interaction network

##perform HC-JDINAC on given two communities c0,c1
c0=c(3,13,26,29);c1=c(16,21,30,44)
data=read.csv("community.csv")
dfh=read.csv("HC feature.csv")  ##HC feature selection result
result=jdinac_bd(y1,y2,data,dfh,k0=1,k1=1)
eset=result$eset ## differential edge/interaction selected by JDINAC
err=result$err ## classification error
p=diffnet_plot(eset) ## differential interaction network

##perform communitities c0,c1 identification using differential interactions of sample pair y1,y2. (see Table 6 in paper)
y1=13;y2=16
data=read.csv("individual.csv")
dfh=read.csv("HC feature.csv")  ##HC feature selection result
result1=jdinac_bd(y1,y2,data,dfh,k0=1,k1=1)
eset1=result$eset ## individual differential edge/interaction selected by JDINAC
c0=c(3,13,26,29);c1=c(16,21,30,44)
dt=read.csv("community-individual.csv")
res=com_ind(eset,c0,c1,dt)  ##testing classification error and confusing rate
```

## JDINAC_m.R
### Description
Perform JDINAC with bivariate kernel density estimation using bandwidth matrix.
### Usage
`jdinac.m(h0,h1,EDGE,classLabel,DataFit,DataPre,nsplit=10,nfolds=5) `
### Arguments
* h0: bandwidth matrix list for class 0, consisting of bandwidth matrixs for each edge (interaction) in EDGE.
* h1: bandwidth matrix list for class 1, consisting of bandwidth matrixs for each edge (interaction) in EDGE.
* EDGE: array indices, which edge will be tested.
* classLabel: must be 0 or 1, e.g. 1 for cases and 0 for controls.
* DataFit: training data matrices containing one sample per row, one variable per column.
* DataPre: testing data matrices containing one sample per row, one variable per column.
* nsplit: randomly split the dataset 2*nsplit times.
* nfolds: number of folds for Cross-validation in penalized logistic regression.
### Values
The output will be a list containing the following components:
* yPre: vector of average predicted probablity for each sample in testing data DataPre.
* Eset: matrix of differential edge set. The first two columns are the array indices; the edge between row-th variable and col-th variable. The 3rd column (numb) is the differential dependency weight; not normalized.
* Vars: matrix of nonzero coefficients in all data split. The first column (t) is the t-th data split. The second column (coef) is the nonzero coefficients (excluding intercept). The last two columns are the array indices, the edge between row-th variable and col-th variable.
* preY: matrix of predicted probablity for each sample testing data DataPre for each data split. Total 2\*nsplit columns.
### Examples
```
## perform HC-JDINAC on sample pair E13,E16 or communites C0,C1
y1=13;y2=16
data=read.csv("individual.csv")
pdim=ncol(data)-2
train=data[data$y %in% c(y1,y2) & data$type=="train",-ncol(data)]
test=data[data$y %in% c(y1,y2) & data$type=="test",-ncol(data)]
###communities data
#y1=c(3,13,26,29);y2=c(16,21,30,44)
#data=read.csv("community.csv")
#pdim=ncol(data)-3
#train=data[data$type=="train",1:(pdim+1)]
#test=data[data$type=="test",1:(pdim+1)]

index=rep(0,pdim)
for(i in 1:pdim){index[i]=as.numeric(substring(colnames(data)[i],2,))}
dfh=read.csv("HC feature.csv")  ##HC feature selection result
dh=colSums(dfh[dfh[,1] %in% c(y1,y2) | dfh[,2] %in% c(y1,y2),-c(1,2)])
vs=index[which(dh>=1)]
vs=vs[vs<4000]
p=length(vs)
x_ind=ifelse(index %in% vs,1,0)  ## feature index selected by HC 
x_ind[pdim+1]=1

E=diag(0,p)
colnames(E)=vs
rownames(E)=vs
f_same=function(vs,x,num){row=ifelse(vs %in% x,num,0)}
pvc=c(600:750,1076:1116,1183:1223,1232:1272,1310:1350,1409:1449,2830:2870,2904:2944,2948:2988)  
caco3=c(698:732,833:896,1064:1104,1420:1460,1776:1816,2493:2533,2855:2895,2962:3002)
PAEs=c(685:790,1020:1093,1104:1144,1266:1306,1361:1401,1442:1482,1560:1620,1710:1750,2843:2894,2912:2981,3017:3057,3417:3457)
E[f_same(vs,pvc,1)==1,f_same(vs,pvc,1)==1]=1
E[f_same(vs,caco3,2)==2,f_same(vs,caco3,2)==2]=1
E[f_same(vs,PAEs,3)==3,f_same(vs,PAEs,3)==3]=1
E[lower.tri(E,diag=T)]=0  
EDGE=which(E!=0, arr.ind=T)  

## input data
sam_train=train[,x_ind==1]
sam_test=test[,x_ind==1]
classlabel_train=ifelse(sam_train$y %in% y1,0,1) 
classlabel_test=ifelse(sam_test$y %in% y1,0,1) 
datafit=sam_train[,-ncol(sam_train)]   
datapre=sam_test[,-ncol(sam_test)]   
## bandwidth matrix
joint_bw0=list()
joint_bw1=list()
for(e in 1:nrow(EDGE)){
  rnum=EDGE[e,1];cnum=EDGE[e,2]
  edt0=sam_0[,c(rnum,cnum)]
  bw0=Hns(edt0)
  joint_bw0[[e]]=bw0
  edt1=sam_1[,c(rnum,cnum)]
  bw1=Hns(edt1)
  joint_bw1[[e]]=bw1
}

difnet=jdinac.m(h0=joint_bw0,h1=joint_bw1,EDGE=EDGE,classLabel=classlabel_train,DataFit=datafit,DataPre=datapre,nsplit=10,nfolds=5) 
ypre=ifelse(difnet$yPre>0.5,1,0) ##predicted class
mean(ypre!=classlabel_test)   #classification error
ypre=data.frame(ypre)
ypre['true_class']=classlabel_test
ypre['label']=sam_test$y
ypre=cbind(ypre,difnet$yPre,difnet$preY)
Vars=data.frame(difnet$Vars)  #coefficients
Vars['v1']=vs[Vars[,3]]  
Vars['v2']=vs[Vars[,4]]  
eset=difnet$Eset  #differential interaction
eset=data.frame(eset)
eset['v1']=vs[eset[,1]]
eset['v2']=vs[eset[,2]] 
```




