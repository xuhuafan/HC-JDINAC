#' @description function jdinac.m: perform JDINAC with bivariate kernel density estimation using bandwidth matrix.
#'
#' @param h0: bandwidth matrix list for class 0, consisting of bandwidth matrixs for each edge (interaction) in EDGE.
#' @param h1: bandwidth matrix list for class 1, consisting of bandwidth matrixs for each edge (interaction) in EDGE.
#' @param EDGE: array indices, which edge will be tested.
#' @param classLabel: must be 0 or 1, e.g. 1 for cases and 0 for controls.
#' @param DataFit: training data matrices containing one sample per row, one variable per column.
#' @param DataPre: testing data matrices containing one sample per row, one variable per column.
#' @param nsplit: randomly split the dataset 2*nsplit times.
#' @param nfolds: number of folds for Cross-validation in penalized logistic regression.
#'
#' @return The output will be a list containing the following components:
#' @export yPre: vector of average predicted probablity for each sample in testing data DataPre.
#' @export Eset: matrix of differential edge set. The first two columns are the array indices; the edge between row-th variable and col-th variable. The 3rd column (numb) is the differential dependency weight; not normalized.
#' @export Vars: matrix of nonzero coefficients in all data split. The first column (t) is the t-th data split. The second column (coef) is the nonzero coefficients (excluding intercept). The last two columns are the array indices, the edge between row-th variable and col-th variable.
#' @export preY: matrix of predicted probablity for each sample testing data DataPre for each data split. Total 2*nsplit columns.

library(ks)
library(akima)
library(glmnet)
denPre2D.m <- function(edge,H0,H1,classLabel,DataFit,newData){
  Index0 <- which(classLabel==0)
  Index1 <- which(classLabel==1)
  x <- DataFit[c(Index0,Index1),edge[1]]   
  y <- DataFit[c(Index0,Index1),edge[2]]    

  x0 <- newData[,edge[1]]
  y0 <- newData[,edge[2]]
  N <- max(c(x,y))
  
  coI <- 1:length(Index0)
  caI <- length(Index0)+(1:length(Index1))
  caSmooth <- try(kde(x=cbind(x[caI], y[caI]), H=H1),silent=TRUE)
  if('try-error' %in% class(caSmooth))            
  {
    caSmooth <- kde(x=cbind(x[caI], y[caI]),H=H1,gridsize = 100)                              
  }
  caP <- caSmooth$estimate
  coSmooth <- try(kde(x=cbind(x[coI], y[coI]), H=H0),silent=TRUE)
  if('try-error' %in% class(coSmooth))            
  {
    coSmooth <- kde(x=cbind(x[coI], y[coI]),H=H0,gridsize = 100)                             
  }
  coP <- coSmooth$estimate
  
  # make sure there are no zeros in the smooth function (since we will take a log of that)
  # caP[caP==0] <- min(caP[caP>0])/100
  caP <- pmax(caP, 1e-10)
  # coP[coP==0] <- min(coP[coP>0])/100
  coP <- pmax(coP, 1e-10)
  
  cafit <- bicubic(x=caSmooth$eval.points[[1]], y=caSmooth$eval.points[[2]], z=caP, x0=x0,y0=y0)
  cofit <- bicubic(x=coSmooth$eval.points[[1]], y=coSmooth$eval.points[[2]], z=coP, x0=x0,y0=y0)
  cafit$z <- pmax(cafit$z, 1e-10)
  cofit$z <- pmax(cofit$z, 1e-10)
  denPre <- log(cafit$z/cofit$z)
  denPre
}

jdinac.m <- function(h0,h1,EDGE,classLabel,DataFit,DataPre,nsplit=10,nfolds=5){
  if(missing(DataPre)) {
    DataPre <- DataFit
    zPre <- zFit
  } 
 if(nsplit<=0){       
    stop("nsplit must be positive integer")    
  } 
  else {
    preY <- NULL
    vset <- NULL
    vars=NULL
    for(i in 1:nsplit){
      size0 <- sum(classLabel==0)
      size1 <- sum(classLabel==1)
      index0=which(classLabel==0)
      index1=which(classLabel==1)
      sn0 <- round(size0/2)
      sn1 <- round(size1/2)
      set.seed(10*i+456)
      splitid = c(sample(index0,sn0),sample(index1,sn1))
      cLabel <- classLabel[splitid]   
      denX=matrix(nrow=nrow(rbind(DataFit[-splitid,],DataPre)),ncol=0)
      for(e in 1:nrow(EDGE)){
        H0=h0[[e]]
        H1=h1[[e]]
        edge=EDGE[e,]
        denx=denPre2D.m(edge,H0,H1,classLabel=cLabel,DataFit=DataFit[splitid,],newData=rbind(DataFit[-splitid,],DataPre))
        denX=cbind(denX,denx)
      }
      y <- classLabel[-splitid]  
      cv.fit <- cv.glmnet(x=data.matrix(denX[1:length(y),]), y=as.matrix(y), family = "binomial", nfolds=nfolds) 
      yp <- predict(cv.fit,newx=data.matrix(denX[-(1:length(y)),]), s="lambda.min",type="response")
      preY <- cbind(preY,yp)  
      
      coes<-coef(cv.fit,s="lambda.min")[-1]  
      coefs <- which(coes!=0)
      vset <- c(vset,coefs)   
      
      var=cbind(coefs,rep(2*i-1,length(coefs)),coes[coefs])
      vars=rbind(vars,var)
      
      cLabel <- classLabel[-splitid]  
      denX=matrix(nrow=nrow(rbind(DataFit[splitid,],DataPre)),ncol=0)
      for(e in 1:nrow(EDGE)){
        H0=h0[[e]]
        H1=h1[[e]]
        edge=EDGE[e,]
        denx=denPre2D.m(edge,H0,H1,classLabel=cLabel,DataFit=DataFit[-splitid,],newData=rbind(DataFit[splitid,],DataPre))
        denX=cbind(denX,denx)
      }
      y <- classLabel[splitid]
      cv.fit <- cv.glmnet(x=data.matrix(denX[1:length(y),]), y=as.matrix(y), family = "binomial", nfolds=nfolds) 
      yp <- predict(cv.fit,newx=data.matrix(denX[-(1:length(y)),]), s="lambda.min",type="response")
      preY <- cbind(preY,yp)  
      
      coes<-coef(cv.fit,s="lambda.min")[-1]  
      coefs <- which(coes!=0)
      vset <- c(vset,coefs)   
      var=cbind(coefs,rep(2*i,length(coefs)),coes[coefs]) 
      vars=rbind(vars,var) 
    } 
    yPre <- rowMeans(preY)  
    numb <- table(vset)   
    Vid <- as.numeric(rownames(numb))  
    Eset <- cbind(EDGE[Vid,],as.numeric(numb))  
    Eset <- Eset[order(Eset[,3],decreasing=T),]  
    colnames(Eset) <- c("row","col","numb")    
    Vars=cbind(vars,EDGE[vars[,1],])
    Vars=Vars[,-1]
    colnames(Vars)=c("t","coef","row","col")    
  } 
  list(yPre=yPre,Eset=Eset,Vars=Vars,preY=preY) 
}

#' @description function jdinac_Md: perform JDINAC with binned bivariate kernel density estimation on given two classes y1,y2;
#' @param y1: the first class to be classified.
#' @param y2: the second class to be classified.
#' @param data: data matrices containing one sample per row, column "y" is class label, column "type" indicates "train" or "test", and the other columns are variables.
#' @param dfh: data matrices containing one pairwise HC test per row, the first two columns are class labels tested by HC, and the other columns are 0,1 indicator, where 1 indicates the corresponding feature is selected by HC.
#' @return The output will be a list containing the following components:
#' @export err: testing classification error.
#' @export Eset: matrix of differential edge set. The first two columns are the array indices; the edge between row-th variable and col-th variable. The 3rd column (numb) is the differential dependency weight; not normalized. The last two columns are features(wavenumbers) of each edge set.
jdinac_md=function(y1,y2,data,dfh){
  pdim=870
  train=data[data$y %in% c(y1,y2) & data$type=="train",1:(pdim+1)]
  test=data[data$y %in% c(y1,y2) & data$type=="test",1:(pdim+1)]
  index=rep(0,pdim)
  for(i in 1:pdim){index[i]=as.numeric(substring(colnames(data)[i],2,))}
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
  sam_train=train[,x_ind==1]
  sam_test=test[,x_ind==1]
  classlabel_train=ifelse(sam_train$y %in% y1,0,1)
  classlabel_test=ifelse(sam_test$y %in% y1,0,1)
  datafit=sam_train[,-ncol(sam_train)]
  datapre=sam_test[,-ncol(sam_test)]
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
  err=mean(ypre!=classlabel_test)   #classification error
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
  list(err=err,eset=eset)
}

#' @description function group_plot: Wavenumbers are grouped and dumbbell graphs and bar graphs are drawn to analyze the distribution of differential interactions.
#' @param d: dataframe containing one differential edge per row, and columns "v1","v2" repersent the two nodes connected by the differntial edge.
#' @return The output will be a network diagram.
#' @export p1: the dumbbell graph of differential interactions.
#' @export p2: the bar graph of the frequency of the differential interactions.
library(dplyr)
library(ggplot2)
library(ggalt)
group_plot=function(d){
  d["var1"]="var1"
  for(i in 1:nrow(d)){
    v1=d[i,'v1']
    if(v1 %in% c(1072:1099)){d[i,"var1"]=1084}
    else if(v1 %in% c(1103:1144)){d[i,"var1"]=1124}
    else if(v1 %in% c(1238:1246)){d[i,"var1"]=1242}
    else if(v1 %in% c(1269:1306)){d[i,"var1"]=1286}
    else if(v1 %in% c(1460:1482)){d[i,"var1"]=1462}
    else if(v1 %in% c(1420:1459)){d[i,"var1"]=1440}
    else if(v1 %in% c(1710:1750)){d[i,"var1"]=1730}
    else if(v1 %in% c(2904:2944)){d[i,"var1"]=2924}
    else if(v1 %in% c(1361:1401)){d[i,"var1"]=1381}
    else {d[i,"var1"]=v1}
  }
  d["var2"]="var2"
  for(i in 1:nrow(d)){
    v1=d[i,'v2']
    if(v1 %in% c(1072:1099)){d[i,"var2"]=1084}
    else if(v1 %in% c(1103:1144)){d[i,"var2"]=1124}
    else if(v1 %in% c(1238:1246)){d[i,"var2"]=1242}
    else if(v1 %in% c(1269:1306)){d[i,"var2"]=1286}
    else if(v1 %in% c(1460:1482)){d[i,"var2"]=1462}
    else if(v1 %in% c(1420:1459)){d[i,"var2"]=1440}
    else if(v1 %in% c(1710:1750)){d[i,"var2"]=1730}
    else if(v1 %in% c(2904:2944)){d[i,"var2"]=2924}
    else if(v1 %in% c(1361:1401)){d[i,"var2"]=1381}
    else {d[i,"var2"]=v1}
  }
  group_d = d %>%
    group_by(var1,var2) %>%
    summarise(freq=sum(numb),.groups="keep")
  group_d=group_d[order(group_d$freq,decreasing = T),]
  group_d["db"]=1:nrow(group_d)
  group_d["fq"]=group_d$freq/sum(group_d$freq)
  group_d$var1=as.numeric(group_d$var1)
  group_d$var2=as.numeric(group_d$var2)

  p1=ggplot(aes(x=var1,xend=var2,y=db),data=group_d)+
    geom_dumbbell(colour_x = "#FFB6C1",colour_xend = "#4169E1",
                  size_x = 1.5,size_xend = 1.5,size=1,color="gray")+
    xlab(expression(Wavenumbers(cm^-1)))+ylab("")+
    scale_x_continuous(breaks = c(1124,1286,1462,1730,2000,2500,2924))+
    theme_bw()+coord_flip()+
    theme(text=element_text(size=14,  family="serif"),
          axis.title.y = element_text(margin =margin(0,0.2,0,0,'cm')))
  p2=ggplot(aes(x=db,y=fq),data=group_d)+
    geom_bar(stat="identity",fill="#FC4E07")+
    theme_bw()+xlab("Interactions")+ylab("Frequency") +
    theme(text=element_text(size=14,  family="serif"),
          axis.title.y = element_text(margin =margin(0,0.8,0,0,'cm')))
  return(list(p1=p1,p2=p2))
}
