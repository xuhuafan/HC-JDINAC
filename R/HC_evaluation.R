#' @description function HC_LDA: perform HC-LDA on given two classes y1,y2, where dividing the data into training set and test set according to 8:7.
#'
#' @param y1: the first class to be tested by HCT.
#' @param y2: the second class to be tested by HCT.
#' @param data: data matrices containing one sample per row, column "y" is class label, and the other columns are variables.
#' @param alpha: significance level of two-sample test.
#'
#' @return The output will be a list containing the following components:
#' @export sel_wave: features selected by HC on the training set.
#' @export res: vector of length two containing testing classification error and FDR.
library(MASS)
HCT=function(y1,y2,data,alpha=0.1){
  Zs=NULL
  p_t=NULL
  for(i in 1:(ncol(data)-1)){
    z1=data[data$y %in% y1,i]
    z2=data[data$y %in% y2,i]
    if(sum(z1-z2)==0 |fivenum(z1)[5]-fivenum(z1)[1]==0
       |fivenum(z2)[5]-fivenum(z2)[1]==0 )
    {pt=1;zs=0}
    else
    {t_test=t.test(z1,z2,var.equal=TRUE,alternative = "two.sided")
    zs=t_test$statistic
    pt=t_test$p.value}
    Zs=append(Zs,zs)
    p_t=append(p_t,pt)
  }
  res=data.frame(1:(ncol(data)-1),p_t,Zs)
  p=sort(p_t)
  HC=NULL
  N=length(p)
  for(i in 1:N){
    if(p[i]*(1-p[i])==0){
      HC[i]=sqrt(N)*(i/N-p[i])/sqrt(p[i]*(1-p[i])+0.01)
    }
    else
      HC[i]=sqrt(N)*(i/N-p[i])/sqrt(p[i]*(1-p[i]))
  }
  HC1=HC[1:ceiling(alpha*N)]
  t=which.max(HC1)
  th_z=(sort(abs(Zs),decreasing = T))[t]
  res["abs_Z"]=abs(Zs)
  res["HC"]=HC[rank(p_t)]
  v=ifelse(res$abs_Z >= th_z,1,0)
  res["v"]=v
  colnames(res)=c("col","p-value","Z-stat","abs_Z","HC","sel")
  return(list(res=res,th_z=th_z))
}

HC_LDA=function(y1,y2,data,alpha=0.1){
  dt=data[data$y %in% c(y1,y2),]
  x_index=rep(0,ncol(dt)-1)
  for(i in 1:(ncol(dt)-1)){
    x_index[i]=as.numeric(substring(colnames(dt)[i],2,))
  }
  #ind1=which(dt$y == y1)
  #ind2=which(dt$y == y2)
  ##8:7 split train set and test set
  #set.seed(y1+111)
  #sam_ind1=sample(ind1,length(ind1)*8/15)
  #set.seed(y2+111)
  #sam_ind2=sample(ind2,length(ind2)*8/15)
  #dat=dt[c(sam_ind1,sam_ind2),]
  sam1=data[data$y %in% y1,]
  sam2=data[data$y %in% y2,]
  set.seed(y1+111)
  s1=sample(nrow(sam1),ceiling(nrow(sam1)*0.5))
  train1=data.frame(sam1[s1,])
  set.seed(y2+111)
  s2=sample(nrow(sam2),ceiling(nrow(sam2)*0.5))
  train2=data.frame(sam2[s2,])
  train=rbind(train1,train2)
  test=data.frame(rbind(sam1[-s1,],sam2[-s2,]))

  HCT_res=HCT(y1,y2,train,alpha=0.1)$res
  sel=HCT_res$sel;sel_col=which(sel==1);sel_wave=x_index[sel_col]
  train=train[,c(sel_col,ncol(dt))]
  test=test[,c(sel_col,ncol(dt))]
  lda.fit=lda(y~., data=train)
  lda.pred=predict(lda.fit,test)
  lda.class=lda.pred$class
  err=mean(lda.class!=test$y)
  ##compute FDR
  ##vector "reg" is the true positivie region annotated by experts
  if(y1==13 & y2==16){
    reg=c(700:800,855:895,949:989,1000:1040,1053:1093,1103:1143,
          1220:1320,1304:1344,1400:1500,1565:1605,1610:1650,1702:1743,
          1774:1814,2493:2533,2835:2875,2905:2945,2938:2978)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)}
  else if(y1==7 & y2==21){
    reg=c(700:800,855:895,949:989,1000:1040,1055:1095,1100:1140,
      1220:1320,1304:1344,1400:1500,1565:1605,1610:1650,1702:1742,
      1774:1814,2493:2533,2840:2880,2905:2945,2938:2978,3714:3754)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)
    }
  else if(y1==30 & y2==32){
    reg=c(700:800,855:895,949:989,1000:1040,1055:1095,1100:1140,
          1220:1320,1312:1352,1400:1500,1555:1595,1702:1742,
          1774:1814,2493:2533,2832:2872,2905:2945,2938:2978)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)
    }
  else{FDR="unknown"}
  return(list(sel_wave=sel_wave,res=c(err,FDR)))
}


#' @description function FAIR_LDA: perform FAIR-LDA on given two classes y1,y2, where dividing the data into training set and test set according to 8:7.
#'
#' @param y1: the first class to be tested by FAIR.
#' @param y2: the second class to be tested by FAIR.
#' @param data: data matrices containing one sample per row, column "y" is class label, and the other columns are variables.
#' @param alpha: significance level of two-sample test.
#'
#' @return The output will be a list containing the following components:
#' @export sel_wave: features selected by HC on the training set.
#' @export res: vector of length two containing testing classification error and FDR. 
library(MASS)
FAIR=function(y1,y2,data,alpha=0.1){
  Zs=NULL
  p_t=NULL
  for(i in 1:(ncol(data)-1)){
    z1=data[data$y %in% y1,i]
    z2=data[data$y %in% y2,i]
    if(sum(z1-z2)==0 |fivenum(z1)[5]-fivenum(z1)[1]==0
       |fivenum(z2)[5]-fivenum(z2)[1]==0 )
    {pt=1;zs=0}
    else
    {t_test=t.test(z1,z2,var.equal=TRUE,alternative = "two.sided")
    zs=t_test$statistic
    pt=t_test$p.value}
    Zs=append(Zs,zs)
    p_t=append(p_t,pt)
  }
  res=data.frame(1:(ncol(data)-1),p_t,Zs)
  v=ifelse(res$p_t <= alpha,1,0)
  res["v"]=v
  colnames(res)=c("col","p-value","Z-stat","sel")
  return(res)
}

FAIR_LDA=function(y1,y2,data,alpha=0.1){
  dt=data[data$y %in% c(y1,y2),]
  x_index=rep(0,ncol(dt)-1)
  for(i in 1:(ncol(dt)-1)){
    x_index[i]=as.numeric(substring(colnames(dt)[i],2,))
  }
  sam1=data[data$y %in% y1,]
  sam2=data[data$y %in% y2,]
  set.seed(y1+111)
  s1=sample(nrow(sam1),ceiling(nrow(sam1)*0.5))
  train1=data.frame(sam1[s1,])
  set.seed(y2+111)
  s2=sample(nrow(sam2),ceiling(nrow(sam2)*0.5))
  train2=data.frame(sam2[s2,])
  train=rbind(train1,train2)
  test=data.frame(rbind(sam1[-s1,],sam2[-s2,]))
  
  FAIR_res=FAIR(y1,y2,train,alpha=0.1)
  sel=FAIR_res$sel;sel_col=which(sel==1);sel_wave=x_index[sel_col]
  train=train[,c(sel_col,ncol(dt))]
  test=test[,c(sel_col,ncol(dt))]
  lda.fit=lda(y~., data=train)
  lda.pred=predict(lda.fit,test)
  lda.class=lda.pred$class
  err=mean(lda.class!=test$y)
  ##compute FDR
  ##vector "reg" is the true positivie region annotated by experts
  if(y1==13 & y2==16){
    reg=c(700:800,855:895,949:989,1000:1040,1053:1093,1103:1143,
          1220:1320,1304:1344,1400:1500,1565:1605,1610:1650,1702:1743,
          1774:1814,2493:2533,2835:2875,2905:2945,2938:2978)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)}
  else if(y1==7 & y2==21){
    reg=c(700:800,855:895,949:989,1000:1040,1055:1095,1100:1140,
          1220:1320,1304:1344,1400:1500,1565:1605,1610:1650,1702:1742,
          1774:1814,2493:2533,2840:2880,2905:2945,2938:2978,3714:3754)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)
  }
  else if(y1==30 & y2==32){
    reg=c(700:800,855:895,949:989,1000:1040,1055:1095,1100:1140,
          1220:1320,1312:1352,1400:1500,1555:1595,1702:1742,
          1774:1814,2493:2533,2832:2872,2905:2945,2938:2978)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)
  }
  else{FDR="unknown"}
  return(list(sel_wave=sel_wave,res=c(err,FDR)))
}


#' @description function FANS：perform FANS between two classes y1,y2, where dividing the data into training set and test set according to 8:7.
#'
#' @param y1: the first class to be classified.
#' @param y2: the second class to be classified.
#' @param data: data matrices containing one sample per row, column "y" is class label, and the other columns are variables.
#' @param nsplit: randomly split the dataset 2*nsplit times.
#' @param nfolds: number of folds for Cross-validation in penalized logistic regression.
#'
#' @return The output will be a list containing the following components:
#' @export sel_wave: features selected by HC on the training set.
#' @export res: vector of length two containing testing classification error and FDR. 

library(glmnet)
library(ks)
FANS=function(y1,y2,data,nsplit=10,nfolds=5){
  preY <- NULL
  vset <- NULL
  vars=NULL
  x_index=rep(0,870)
  for(i in 1:870){
    x_index[i]=as.numeric(substring(colnames(data)[i],2,))
  }
  sam1=data[data$y %in% y1,]
  sam2=data[data$y %in% y2,]
  set.seed(y1+111)
  s1=sample(nrow(sam1),ceiling(nrow(sam1)*0.5))
  train1=data.frame(sam1[s1,])
  set.seed(y2+111)
  s2=sample(nrow(sam2),ceiling(nrow(sam2)*0.5))
  train2=data.frame(sam2[s2,])
  train=rbind(train1,train2)
  test=data.frame(rbind(sam1[-s1,],sam2[-s2,]))
  train$y=ifelse(train$y == y1,0,1)
  test$y=ifelse(test$y == y1,0,1)
  DataFit=train[,-ncol(train)]
  DataPre=test[,-ncol(test)]
  h0=apply(train[train$y==0,-ncol(train)],2,bw.nrd0)
  h1=apply(train[train$y==1,-ncol(train)],2,bw.nrd0)
  classLabel=train$y
  size0 <- sum(classLabel==0)
  size1 <- sum(classLabel==1)
  index0=which(classLabel==0)
  index1=which(classLabel==1)
  sn0 <- round(size0/2)
  sn1 <- round(size1/2)
  for(i in 1:nsplit){
    set.seed(i+123)
    splitid0=sample(index0,sn0)
    splitid1=sample(index1,sn1)
    splitid=c(splitid0,splitid1)
    D10=DataFit[splitid0,]
    D11=DataFit[splitid1,]
    D2=DataFit[-splitid,]
    y=classLabel[-splitid]
    Z=matrix(nrow = nrow(D2),ncol=ncol(DataFit))
    newZ=matrix(nrow = nrow(DataPre),ncol=ncol(DataPre))
    for(j in 1:ncol(DataFit)){
      f=kde(D10[,j],h=h0[j])
      g=kde(D11[,j],h=h1[j])
      fj=dkde(D2[,j],f)
      fj=pmax(fj, 1e-8)
      gj=dkde(D2[,j],g)
      gj=pmax(gj, 1e-8)
      dr=fj/gj
      Z[,j]=log(fj/gj)
      newf=dkde(DataPre[,j],f)
      newf=pmax(newf, 1e-8)
      newg=dkde(DataPre[,j],g)
      newg=pmax(newg, 1e-8)
      newZ[,j]=log(newf/newg)
    }
    cv.fit=cv.glmnet(Z,y,family="binomial",nfolds=nfolds)
    yp <- predict(cv.fit,newx=newZ,s="lambda.min",type="response")
    preY <- cbind(preY,yp)
    coes<-coef(cv.fit,s="lambda.min")[-1]
    coefs <- which(coes!=0)
    vset <- c(vset,coefs)
    D10=DataFit[setdiff(index0,splitid0),]
    D11=DataFit[setdiff(index1,splitid1),]
    D2=DataFit[splitid,]
    y=classLabel[splitid]
    Z=matrix(nrow = nrow(D2),ncol=ncol(DataFit))
    newZ=matrix(nrow = nrow(DataPre),ncol=ncol(DataPre))
    for(j in 1:ncol(DataFit)){
      f=kde(D10[,j],h=h0[j])
      g=kde(D11[,j],h=h1[j])
      fj=dkde(D2[,j],f)
      fj=pmax(fj, 1e-8)
      gj=dkde(D2[,j],g)
      gj=pmax(gj, 1e-8)
      dr=fj/gj
      Z[,j]=log(fj/gj)
      newf=dkde(DataPre[,j],f)
      newf=pmax(newf, 1e-8)
      newg=dkde(DataPre[,j],g)
      newg=pmax(newg, 1e-8)
      newZ[,j]=log(newf/newg)
    }
    cv.fit=cv.glmnet(Z,y,family="binomial",nfolds=nfolds)
    yp <- predict(cv.fit,newx=newZ,s="lambda.min",type="response")
    preY <- cbind(preY,yp)
    coes<-coef(cv.fit,s="lambda.min")[-1]
    coefs <- which(coes!=0)
    vset <- c(vset,coefs)
  }
  yPre <- rowMeans(preY)
  predclass=ifelse(yPre>0.5,1,0)
  err=mean(predclass!=test$y)
  vs <- sort(unique(vset))
  sel_wave=x_index[vs]
  ##compute FDR
  ##vector "reg" is the true positivie region annotated by experts
  if(y1==13 & y2==16){
    reg=c(700:800,855:895,949:989,1000:1040,1053:1093,1103:1143,
          1220:1320,1304:1344,1400:1500,1565:1605,1610:1650,1702:1743,
          1774:1814,2493:2533,2835:2875,2905:2945,2938:2978)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)}
  else if(y1==7 & y2==21){
    reg=c(700:800,855:895,949:989,1000:1040,1055:1095,1100:1140,
          1220:1320,1304:1344,1400:1500,1565:1605,1610:1650,1702:1742,
          1774:1814,2493:2533,2840:2880,2905:2945,2938:2978,3714:3754)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)
  }
  else if(y1==30 & y2==32){
    reg=c(700:800,855:895,949:989,1000:1040,1055:1095,1100:1140,
          1220:1320,1312:1352,1400:1500,1555:1595,1702:1742,
          1774:1814,2493:2533,2832:2872,2905:2945,2938:2978)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)
  }
  else{FDR="unknown"}
  list(sel_wave=sel_wave,res=c(err,FDR))
}

#' @description function SFLDA：perform SFLDA between two classess y1,y2, where dividing the data into training set and test set according to 8:7.
#'
#' @param y1: the first class to be classified.
#' @param y2: the second class to be classified.
#' @param data: data matrices containing one sample per row, column "y" is class label, and the other columns are variables.
#'
#' @return The output will be a list containing the following components:
#' @export sel_wave: features selected by HC on the training set.
#' @export res: vector of length two containing testing classification error and FDR. 

library(CVXR)
library(foreach)
library(doParallel)
library(caret)
library(splines)
penalty = function(beta, lambda = 0, eta = 0) {
  p = length(beta)
  smooth =  eta/2 * p_norm(diff(x = beta, differences = 1), 2)
  lasso = lambda * p_norm(beta, 1)
  smooth + lasso
}

loss_mean=function(data,beta,X){
  Y=data.matrix(data[,-ncol(data)])
  ind1=which(data$y == 1)
  ind0=which(data$y == 0)
  n=nrow(Y)
  N=ncol(Y)
  p=4
  Nm=150
  knots1=seq(0,1,length=Nm+2) 
  XB1<-bs(X,knots=knots1[-c(1,(Nm+2))],degree=p-1,intercept=TRUE,Boundary.knots=knots1[c(1,(Nm+2))])
  Beta=solve(t(XB1)%*%(XB1))%*%t(XB1) 
  betahat=Y%*%t(Beta)
  etahat=betahat%*%t(XB1)
  delta=colMeans(etahat[ind1,])-colMeans(etahat[ind0,])
  gamma=cov(etahat)
  list(delta=delta,loss=quad_form(beta, gamma)/2
       - t(delta) %*% beta)
}

CV5=function(dt,lambda,eta,beta,X,k=5){
  set.seed(123)
  folds <- createFolds(y=dt$y,k)
  ERR=c()
  for(f in 1:k){
    fold_test <-dt[folds[[f]],]
    fold_train <- dt[-folds[[f]],]
    Loss=loss_mean(fold_train,beta,X)
    obj = Loss$loss + penalty(beta, lambda, eta)
    prob = Problem(Minimize(obj))
    result = solve(prob)
    beta_hat = result$getValue(beta)
    predx = as.matrix(fold_test[,-ncol(dt)])
    Tx = apply(predx,1,function(x){(x %*% beta_hat - t(Loss$delta) %*% beta_hat)^2 - (x %*% beta_hat)^2})
    predclass = ifelse(Tx<0,1,0)
    err=mean(as.numeric(fold_test[,ncol(dt)]!=predclass))
    ERR=c(ERR,err)
  }
  return(mean(ERR))
}

SFLDA=function(y1,y2,data){
  x_index=rep(0,ncol(data)-1)
  for(i in 1:(ncol(data)-1)){
    x_index[i]=as.numeric(substring(colnames(data)[i],2,))
  }
  X=as.numeric(x_index-x_index[1])/(x_index[length(x_index)]-x_index[1])
  beta = Variable(ncol(data)-1)
  sam1=data[data$y %in% y1,]
  sam2=data[data$y %in% y2,]
  set.seed(y1+111)
  s1=sample(nrow(sam1),ceiling(nrow(sam1)*0.5))
  train1=data.frame(sam1[s1,])
  set.seed(y2+111)
  s2=sample(nrow(sam2),ceiling(nrow(sam2)*0.5))
  train2=data.frame(sam2[s2,])
  train=rbind(train1,train2)
  train$y=ifelse(train$y %in% y1,1,0)
  test=data.frame(rbind(sam1[-s1,],sam2[-s2,]))
  test$y=ifelse(test$y %in% y1,1,0)  
  mdata=rbind(train,test)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  res <- foreach(lam = rep(c(2^c(-1:-10),0),11), eta = c( rep(c(2^c(-1:-10),0),rep(11,11))), 
                 .packages=c('CVXR',"caret","splines"),.export=c("CV5","loss_mean","penalty"),.combine = "rbind") %dopar% CV5(train,lam,eta,beta,X,k=5)   
  stopImplicitCluster()
  RES=data.frame(lam = rep(c(2^c(-1:-10),0),11), eta = c( rep(c(2^c(-1:-10),0),rep(11,11))))
  RES['err']=res[,1]
  opt_para=which(RES$err==min(RES$err))
  num=length(opt_para)
  Beta_hat=matrix(0,nrow=ncol(train)-1,ncol=num)
  ERR=c();trainERR=c()
  for(i in 1:num){
    l=RES[opt_para[i],1]
    e=RES[opt_para[i],2]
    Loss = loss_mean(train,beta,X)
    obj = Loss$loss + penalty(beta, l, e)
    prob = Problem(Minimize(obj))
    result = solve(prob)
    beta_hat = result$getValue(beta)
    Beta_hat[,i]=beta_hat
    predx = as.matrix(mdata[,-ncol(test)])
    Tx = apply(predx,1,function(x){(x %*% beta_hat - t(Loss$delta) %*% beta_hat)^2 - (x %*% beta_hat)^2})
    predclass = ifelse(Tx<0,1,0)
    err=as.numeric(mdata[,ncol(test)]!=predclass)
    trainerr=mean(err[1:nrow(train)])
    testerr=mean(err[-(1:nrow(train))])
    trainERR=c(trainERR,trainerr)
    ERR=c(ERR,testerr)
    if(trainerr==0){break}}
  opt_k=which.min(trainERR)
  opt_l=RES[opt_para[opt_k],1]
  opt_e=RES[opt_para[opt_k],2]
  fres=c(y1,y2,opt_l,opt_e,ERR[opt_k],trainERR[opt_k])
  BETA=Beta_hat[,opt_k]
  BETA2=abs(BETA)
  if(max(BETA2)<1e-4){fs=ifelse(BETA2>quantile(BETA2,0.95),1,0)}
  else{fs=ifelse(abs(BETA) > 0.0001,1,0)}
  sel_col=which(fs==1)
  sel_wave=x_index[sel_col]
  ##compute FDR
  ##vector "reg" is the true positivie region annotated by experts
  if(y1==13 & y2==16){
    reg=c(700:800,855:895,949:989,1000:1040,1053:1093,1103:1143,
          1220:1320,1304:1344,1400:1500,1565:1605,1610:1650,1702:1743,
          1774:1814,2493:2533,2835:2875,2905:2945,2938:2978)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)}
  else if(y1==7 & y2==21){
    reg=c(700:800,855:895,949:989,1000:1040,1055:1095,1100:1140,
          1220:1320,1304:1344,1400:1500,1565:1605,1610:1650,1702:1742,
          1774:1814,2493:2533,2840:2880,2905:2945,2938:2978,3714:3754)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)
  }
  else if(y1==30 & y2==32){
    reg=c(700:800,855:895,949:989,1000:1040,1055:1095,1100:1140,
          1220:1320,1312:1352,1400:1500,1555:1595,1702:1742,
          1774:1814,2493:2533,2832:2872,2905:2945,2938:2978)
    FDR=length(setdiff(sel_wave,reg))/length(sel_wave)
  }
  else{FDR="unknown"}
  
  list(sel_wave=sel_wave,res=c(ERR[opt_k],FDR))
}

