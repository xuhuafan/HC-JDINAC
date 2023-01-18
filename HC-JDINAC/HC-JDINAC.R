#' Perform HC feature selection
#'
#' \code{HC} returns the HC feature selection result.
#'
#' @param y1: One sample label to be tested by HC
#' @param y2: The other sample label to be tested by HC
#' @param x: Data matrix
#' @param alpha: Threshold value of t test in HC
#'
#' @return: The output will be a vector of 0,1 with a length equal to the number of variables, where 1 indicates that the variable is valid.
HC=function(y1,y2,x,alpha=0.1){
  p_t=NULL
  m=ncol(x)-1  #number of variables
  x_index=colnames(x)[1:m]
  for(i in 1:m){  
    z1=x[x$y %in% y1,i]
    z2=x[x$y %in% y2,i]
    if(length(unique(z1))==1 & length(unique(z2))==1){p=1}
    else{
      p=t.test(z1,z2,var.equal=FALSE,alternative = "two.sided")$p.value}
    p_t=append(p_t,p)}
  p_t=data.frame(x_index,p_t)
  colnames(p_t)=c("wave","pvalue")
  p_t=p_t[order(p_t$pvalue,decreasing=FALSE),]
  ## compute HC value
  p=sort(p_t$pvalue)
  HC=NULL
  for(i in 1:m){
    if(p[i]*(1-p[i])==0){
      HC[i]=sqrt(m)*(i/m-p[i])/sqrt(p[i]*(1-p[i])+0.01)
    }
    else
      HC[i]=sqrt(m)*(i/m-p[i])/sqrt(p[i]*(1-p[i]))
  }
  p_t["HC"]=HC
  ## alpha=0.1
  HC1=HC[1:ceiling(alpha*m)]
  t=which.max(HC1) ##reject first t hypothesis
  t_index=p_t$wave[1:t]  #valid features
  v=ifelse(x_index %in% t_index,1,0)  
 return(v)
}

#' Interaction/edge construction according to synergies between wavenumbers/features
#'
#' \code{edge_construct} returns interactive edges between features.
#'
#' @param vs: Valid single features selected by HC.
#'
#' @return: The output will be a matrix with each row corresponding to an interaction.
edge_construct(vs){
  q=length(vs)
  E=diag(0,q)
  colnames(E)=vs
  rownames(E)=vs
  f=function(x,num){
    row=ifelse(vs %in% x,num,0)
  }
  ## characteristic absorption peaks of pvc:
  pvc=c(600:750,1076:1116,1183:1223,1232:1272,1310:1350,1409:1449,2830:2870,2904:2944,2948:2988)  
  E[f(pvc,1)==1,f(pvc,1)==1]=1
  ## characteristic absorption peaks of caco3:
  caco3=c(698:732,833:896,1064:1104,1420:1460,1776:1816,2493:2533,2855:2895,2962:3002)
  E[f(caco3,2)==2,f(caco3,2)==2]=1
  ## characteristic absorption peaks of PAEs:
  PAEs=c(685:790,1020:1093,1104:1144,1266:1306,1361:1401,1442:1482,1560:1620,1710:1750,2843:2894,2912:2981,3017:3057,3417:3457)
  E[f(PAEs,3)==3,f(PAEs,3)==3]=1
  E[lower.tri(E,diag=T)]=0  
  EDGE=which(E!=0, arr.ind=T) 
  return(EDGE)
}

library(KernSmooth)
library(akima)
library(glmnet)
denPre2D.h <- function(edge,h0,h1,classLabel,DataFit,newData,method=c("integers","bandwidth")[2]){
  Index0 <- which(classLabel==0)
  Index1 <- which(classLabel==1)
  x <- DataFit[c(Index0,Index1),edge[1]]   
  y <- DataFit[c(Index0,Index1),edge[2]]   
  x0 <- newData[,edge[1]]
  y0 <- newData[,edge[2]]
  N <- max(c(x,y))
  coI <- 1:length(Index0)
  caI <- length(Index0)+(1:length(Index1))
  caWidth <- c(bw.nrd0(x[caI]),bw.nrd0(y[caI]))  
  coWidth <- c(bw.nrd0(x[coI]),bw.nrd0(y[coI]))   
  gridSize <- switch(method,integers  = c(N, N),bandwidth = ceiling(N / c(min(caWidth[1], coWidth[1]), min(caWidth[2], coWidth[2]))))
  gridSize <- pmax(gridSize,10) 
  #make sure there are at least 100 points in total(10*10)
  caSmooth <- bkde2D(x=cbind(x[caI], y[caI]), bandwidth=h1,  gridsize=gridSize)
  caP <- caSmooth$fhat
  coSmooth <- bkde2D(x=cbind(x[coI], y[coI]), bandwidth=h0, gridsize=gridSize)
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
  ## log(f_ij/g_ij)
  denPre <- log(cafit$z/cofit$z)
  denPre
}

#' Perfrom JDINAC
#'
#' \code{JDINAC.h} returns the results of JDINAC method, including predicted probablity, selected interactions and their weights, and non-zero coefficients.
#'
#' @param h0: The bandwidth of class 0 samples.
#' @param h1: The bandwidth of class 1 samples.
#' @param EDGE: The interaction terms to consider in logistic regression.
#' @param classLabel: The class labels of DataFit.
#' @param DataFit: The data used to train JDINAC model.
#' @param DataPre: The data to be predicted.
#' @param nsplit: The number of data split.
#' @param nfolds: The folds of cross validation.
#'
#' @return: The output will be a named list containing the predicted probablity of DataPre, the selected interactions and their weight, and non-zero coefficients.
jdinac.h <- function(h0,h1,EDGE,classLabel,DataFit,DataPre,nsplit=10,nfolds=5){
  
  if(missing(DataPre)) {
    DataPre <- DataFit
  } 
  else if(nsplit<=0){       
    stop("nsplit must be positive integer")    
  } 
  else {
    preY <- NULL
    vset <- NULL
    vars=NULL
    for(i in 1:nsplit){
      size0 <- sum(classLabel==0)
      size1 <- sum(classLabel==1)
      sn0 <- round(size0/2)
      sn1 <- round(size1/2)
      index0 <- which(classLabel==0)
      index1 <- which(classLabel==1)
      set.seed(10*i+123)
      splitid <- c(sample(index0,sn0),sample(index1,sn1))
      
      cLabel <- classLabel[splitid] 
      denX <- apply(EDGE,1,denPre2D.h,h0=h0,h1=h1,classLabel=cLabel,DataFit=DataFit[splitid,],newData=rbind(DataFit[-splitid,],DataPre))
      y <- classLabel[-splitid]  
      cv.fit <- cv.glmnet(x=denX[1:length(y),], y=y, family = "binomial", nfolds=nfolds) 
      yp <- predict(cv.fit,newx=denX[-(1:length(y)),], s="lambda.min",type="response")
      preY <- cbind(preY,yp)
      
      coes<-coef(cv.fit,s="lambda.min")[-1]  
      coefs <- which(coes!=0)
      vset <- c(vset,coefs)  
      
      var=cbind(coefs,rep(2*i-1,length(coefs)),coes[coefs]) 
      vars=rbind(vars,var)
      
      cLabel <- classLabel[-splitid]        
      denX <- apply(EDGE,1,denPre2D.h,h0=h0,h1=h1,classLabel=cLabel,DataFit=DataFit[-splitid,],
                    newData=rbind(DataFit[splitid,],DataPre))
      y <- classLabel[splitid]
      cv.fit <- cv.glmnet(x=denX[1:length(y),], y=y, family = "binomial", nfolds=nfolds) 
     
      yp <- predict(cv.fit,newx=denX[-(1:length(y)),], s="lambda.min",type="response")
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
  list(yPre=yPre,Eset=Eset,Vars=Vars) 
}
