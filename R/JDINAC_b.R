#' @title Perform JDINAC with binned bivariate kernel density estimation 
#'
#' @param h0: bandwidth vector for class 0, consisting of bandwidths for each variable in DataFit.
#' @param h1: bandwidth vector for class 1, consisting of bandwidths for each variable in DataFit.
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
#' @examples 
#' ## perform HC-JDINAC on sample pair E13,E16 or communites C0,C1
#' y1=13;y2=16
#' data=read.csv("individual.csv")
#' pdim=ncol(data)-2
#' train=data[data$y %in% c(y1,y2) & data$type=="train",-ncol(data)]
#' test=data[data$y %in% c(y1,y2) & data$type=="test",-ncol(data)]
#' 
#' ## communities data
#' #y1=c(3,13,26,29);y2=c(16,21,30,44)
#' #data=read.csv("community.csv")
#' #pdim=ncol(data)-3
#' #train=data[data$type=="train",1:(pdim+1)]
#' #test=data[data$type=="test",1:(pdim+1)]
#' 
#' index=rep(0,pdim)
#' for(i in 1:pdim){index[i]=as.numeric(substring(colnames(data)[i],2,))}
#' 
#' dfh=read.csv("HC feature.csv")  ##HC feature selection result
#' dh=colSums(dfh[dfh[,1] %in% c(y1,y2) | dfh[,2] %in% c(y1,y2),-c(1,2)])
#' vs=index[which(dh>=1)]
#' vs=vs[vs<4000]
#' p=length(vs)
#' x_ind=ifelse(index %in% vs,1,0)  ## feature index selected by HC 
#' x_ind[pdim+1]=1
#' 
#' E=diag(0,p)
#' colnames(E)=vs
#' rownames(E)=vs
#' f_same=function(vs,x,num){row=ifelse(vs %in% x,num,0)}
#' pvc=c(600:750,1076:1116,1183:1223,1232:1272,1310:1350,1409:1449,2830:2870,2904:2944,2948:2988)  
#' caco3=c(698:732,833:896,1064:1104,1420:1460,1776:1816,2493:2533,2855:2895,2962:3002)
#' PAEs=c(685:790,1020:1093,1104:1144,1266:1306,1361:1401,1442:1482,1560:1620,1710:1750,2843:2894,2912:2981,3017:3057,3417:3457)
#' E[f_same(vs,pvc,1)==1,f_same(vs,pvc,1)==1]=1
#' E[f_same(vs,caco3,2)==2,f_same(vs,caco3,2)==2]=1
#' E[f_same(vs,PAEs,3)==3,f_same(vs,PAEs,3)==3]=1
#' E[lower.tri(E,diag=T)]=0  
#' EDGE=which(E!=0, arr.ind=T)  
#' 
#' ## input data
#' sam_train=train[,x_ind==1]
#' sam_test=test[,x_ind==1]
#' classlabel_train=ifelse(sam_train$y %in% y1,0,1) 
#' classlabel_test=ifelse(sam_test$y %in% y1,0,1) 
#' datafit=sam_train[,-ncol(sam_train)]   
#' datapre=sam_test[,-ncol(sam_test)]   
#' ## bandwidth
#' sam_0=train[train$y %in% y1,x_ind==1]  
#' bws_0=apply(sam_0[,-ncol(sam_0)],2,bw.nrd0)
#' sam_1=train[train$y %in% y2,x_ind==1]
#' bws_1=apply(sam_1[,-ncol(sam_1)],2,bw.nrd0)
#' 
#' difnet=jdinac(h0=bws_0,h1=bws_1,EDGE=EDGE,classLabel=classlabel_train,DataFit=datafit,DataPre=datapre,nsplit=10,nfolds=5) 
#' ypre=ifelse(difnet$yPre>0.5,1,0) ##predicted class
#' mean(ypre!=classlabel_test)   #classification error
#' ypre=data.frame(ypre)
#' ypre['true_class']=classlabel_test
#' ypre['label']=sam_test$y
#' Vars=data.frame(difnet$Vars)  #coefficients
#' Vars['v1']=vs[Vars[,3]]  
#' Vars['v2']=vs[Vars[,4]]  
#' eset=difnet$Eset  #differential interaction
#' eset=data.frame(eset)
#' eset['v1']=vs[eset[,1]]
#' eset['v2']=vs[eset[,2]] 


#' ##Example2: perform community identification using differential interactions of sample pair E13,E16
#' y1=c(3,13,26,29);y2=c(16,21,30,44)
#' data=read.csv("community-individual.csv")
#' pdim=ncol(data)-3
#' train=data[data$type=="train",1:(pdim+1)]
#' test=data[data$type=="test",1:(pdim+1)]
#' 
#' index=rep(0,pdim)
#' for(i in 1:pdim){index[i]=as.numeric(substring(colnames(data)[i],2,))}
#' 
#' ##Construct interaction edge using differential interactions "eset" of sample pair E13,E16
#' vs=sort(unique(c(eset$v1,eset$v2)))  
#' rs=c();cs=c()
#' for(i in 1:nrow(eset)){
#'   v1=eset[i,"v1"];v2=eset[i,"v2"]
#'   r=which(vs==v1);c=which(vs==v2)
#'   rs=c(rs,r);cs=c(cs,c)
#' }
#' eset["r"]=rs;eset["c"]=cs
#' EDGE=eset[,c("r","c")]
#' 
#' p=length(vs)
#' x_ind=ifelse(index %in% vs,1,0) 
#' x_ind[pdim+1]=1
#' 
#' E=diag(0,p)
#' colnames(E)=vs
#' rownames(E)=vs
#' f_same=function(vs,x,num){row=ifelse(vs %in% x,num,0)}
#' pvc=c(600:750,1076:1116,1183:1223,1232:1272,1310:1350,1409:1449,2830:2870,2904:2944,2948:2988)  
#' caco3=c(698:732,833:896,1064:1104,1420:1460,1776:1816,2493:2533,2855:2895,2962:3002)
#' PAEs=c(685:790,1020:1093,1104:1144,1266:1306,1361:1401,1442:1482,1560:1620,1710:1750,2843:2894,2912:2981,3017:3057,3417:3457)
#' E[f_same(vs,pvc,1)==1,f_same(vs,pvc,1)==1]=1
#' E[f_same(vs,caco3,2)==2,f_same(vs,caco3,2)==2]=1
#' E[f_same(vs,PAEs,3)==3,f_same(vs,PAEs,3)==3]=1
#' E[lower.tri(E,diag=T)]=0  
#' EDGE=which(E!=0, arr.ind=T)  
#' 
#' ## input data
#' sam_train=train[,x_ind==1]
#' sam_test=test[,x_ind==1]
#' classlabel_train=ifelse(sam_train$y %in% y1,0,1) 
#' classlabel_test=ifelse(sam_test$y %in% y1,0,1) 
#' datafit=sam_train[,-ncol(sam_train)]   
#' datapre=sam_test[,-ncol(sam_test)]   
#' ## bandwidth
#' sam_0=train[train$y %in% y1,x_ind==1]  
#' bws_0=apply(sam_0[,-ncol(sam_0)],2,bw.nrd0)
#' sam_1=train[train$y %in% y2,x_ind==1]
#' bws_1=apply(sam_1[,-ncol(sam_1)],2,bw.nrd0)
#' 
#' difnet=jdinac(h0=bws_0,h1=bws_1,EDGE=EDGE,classLabel=classlabel_train,DataFit=datafit,DataPre=datapre,nsplit=10,nfolds=5) 
#' ypre=ifelse(difnet$yPre>0.5,1,0) ##predicted class
#' mean(ypre!=classlabel_test)   #classification error
#' ypre=data.frame(ypre)
#' ypre['true_class']=classlabel_test
#' ypre['label']=sam_test$y
#' Vars=data.frame(difnet$Vars)  #coefficients
#' Vars['v1']=vs[Vars[,3]]  
#' Vars['v2']=vs[Vars[,4]]  
#' eset=difnet$Eset  #differential interaction
#' eset=data.frame(eset)
#' eset['v1']=vs[eset[,1]]
#' eset['v2']=vs[eset[,2]] 
library(KernSmooth)
library(akima)
library(glmnet)
denPre2D <- function(edge,h0,h1,classLabel,DataFit,newData){
  Index0 <- which(classLabel==0)
  Index1 <- which(classLabel==1)

  x <- DataFit[c(Index0,Index1),edge[1]]  
  y <- DataFit[c(Index0,Index1),edge[2]]    

  x0 <- newData[,edge[1]]
  y0 <- newData[,edge[2]]
  N <- max(c(x,y))
  
  coI <- 1:length(Index0)
  caI <- length(Index0)+(1:length(Index1))
  
  bw1=h1[c(edge[1],edge[2])]
  bw0=h0[c(edge[1],edge[2])]
  caWidth <- c(bw.nrd0(x[caI]),bw.nrd0(y[caI]))  
  coWidth <- c(bw.nrd0(x[coI]),bw.nrd0(y[coI]))   
  
  gridSize <- ceiling(N / c(min(caWidth[1], coWidth[1]), min(caWidth[2], coWidth[2])))
  gridSize <- pmax(gridSize,10) 
  #make sure there are at least 100 points in total(10*10)

  caSmooth <- bkde2D(x=cbind(x[caI], y[caI]), bandwidth=bw1,  gridsize=gridSize)
  caP <- caSmooth$fhat
  coSmooth <- bkde2D(x=cbind(x[coI], y[coI]), bandwidth=bw0, gridsize=gridSize)
  coP <- coSmooth$fhat
  
  # make sure there are no zeros in the smooth function (since we will take a log of that)
  # caP[caP==0] <- min(caP[caP>0])/100
  caP <- pmax(caP, 1e-10)
  # coP[coP==0] <- min(coP[coP>0])/100
  coP <- pmax(coP, 1e-10)
  
  cafit <- bicubic(x=caSmooth$x1, y=caSmooth$x2, z=caP, x0=x0,y0=y0)
  cofit <- bicubic(x=coSmooth$x1, y=coSmooth$x2, z=coP, x0=x0,y0=y0)
  cafit$z <- pmax(cafit$z, 1e-10)
  cofit$z <- pmax(cofit$z, 1e-10)
  denPre <- log(cafit$z/cofit$z)
  denPre
}

jdinac <- function(h0,h1,EDGE,classLabel,DataFit,DataPre,nsplit=10,nfolds=5){
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
      
      #set.seed(100*i+456)
      splitid = c(sample(index0,sn0),sample(index1,sn1))
      cLabel <- classLabel[splitid]     
      denX <- apply(EDGE,1,denPre2D,h0=h0,h1=h1,classLabel=cLabel,DataFit=DataFit[splitid,],newData=rbind(DataFit[-splitid,],DataPre))
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
      denX <- apply(EDGE,1,denPre2D,h0=h0,h1=h1,classLabel=cLabel,DataFit=DataFit[-splitid,],
                    newData=rbind(DataFit[splitid,],DataPre))
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
    Eset <- cbind(EDGE[Vid,],numb)  
    Eset <- Eset[order(Eset[,3],decreasing=T),]  
    colnames(Eset) <- c("row","col","numb")
    
    Vars=cbind(vars,EDGE[vars[,1],])
    Vars=Vars[,-1]
    colnames(Vars)=c("t","coef","row","col")
  } 
  list(yPre=yPre,Eset=Eset,Vars=Vars,preY=preY) 
}
