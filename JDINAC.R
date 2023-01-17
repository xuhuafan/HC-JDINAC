library(KernSmooth)
library(akima)
library(glmnet)

# D1密度估计，求解D2对数似然比：
########################################
denPre2D <- function(edge,classLabel,DataFit,newData,method=c("integers","bandwidth")[2]){
  
  ## 提出不同类别的样本：
  Index0 <- which(classLabel==0)
  Index1 <- which(classLabel==1)
  
  ##？DataFit：用于密度估计的部分数据D1：
  x <- DataFit[c(Index0,Index1),edge[1]]
  y <- DataFit[c(Index0,Index1),edge[2]]
  
  ##？newData:用于L1-logistics回归的部分数据D2:
  x0 <- newData[,edge[1]]
  y0 <- newData[,edge[2]]
  N <- max(c(x,y))
  
  ## 0类样本编码：（1、2、n0)
  coI <- 1:length(Index0)
  ## 1类样本编码:(n0+1、n1+2、n0+n1)
  caI <- length(Index0)+(1:length(Index1))
  
  ## 分别对两类求带宽：
  caWidth <- c(bw.nrd0(x[caI]),bw.nrd0(y[caI]))
  coWidth <- c(bw.nrd0(x[coI]),bw.nrd0(y[coI]))
  
  gridSize <- switch(method,integers  = c(N, N),bandwidth = ceiling(N / c(min(caWidth[1], coWidth[1]), min(caWidth[2], coWidth[2]))))
  ## 以10向下截断，10以下的取10：
  gridSize <- pmax(gridSize,10) 
  #make sure there are at least 100 points in total(10*10)
  
  ## 2维的密度估计矩阵:  #带宽增加0.001
  caSmooth <- bkde2D(x=cbind(x[caI], y[caI]), bandwidth=caWidth,  gridsize=gridSize)
  ## f_ij估计：
  caP <- caSmooth$fhat
  coSmooth <- bkde2D(x=cbind(x[coI], y[coI]), bandwidth=coWidth, gridsize=gridSize)
  ## g_ij估计：
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
  ## D2上的log(f_ij/g_ij)
  denPre <- log(cafit$z/cofit$z)
  denPre
}

# JDINAC流程：
########################################
jdinac.z <- function(EDGE,classLabel,DataFit,DataPre,zFit,zPre,nsplit=10,nfolds=5){
#### nsplit=T,即重复次数：
  
  if(class(zFit)!="matrix" | class(zPre)!="matrix")
    stop("zFit and zPre must be 'matrix'")
  if(nsplit<=0){       
    stop("nsplit must be positive integer")  
  } 
  else {
    preY <- NULL
    vset <- NULL
    ## repeat step1-3:
    for(i in 1:nsplit){
      size0 <- sum(classLabel==0)
      size1 <- sum(classLabel==1)
      sn0 <- round(size0/2)
      sn1 <- round(size1/2)
      ## 测试集、训练集对半分：
      splitid <- c(sample(1:size0,sn0),sample((1:size1)+size0,sn1))
      
      ## 训练集D1：
      cLabel <- classLabel[splitid]
      ## 利用denPre2D函数估计联合密度：
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[splitid,],newData=rbind(DataFit[-splitid,],DataPre))
      y <- classLabel[-splitid]
      
      ## 测试集D2:
      zfit <- zFit[-splitid, drop=F]
      ## cv求解最优惩罚参数—lambda:
      cv.fit <- cv.glmnet(x=cbind(zfit,denX[1:length(y),]), y=y, family = "binomial", nfolds=nfolds) 
      yp <- predict(cv.fit,newx=cbind(zPre,denX[-(1:length(y)),]), s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      coefs <- which(coef(cv.fit,s="lambda.min")[-(1:(1+ncol(zFit)))] !=0)
      vset <- c(vset,coefs)
      
      cLabel <- classLabel[-splitid]        
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[-splitid,],newData=rbind(DataFit[splitid,],DataPre))
      y <- classLabel[splitid]
      zfit <- zFit[splitid, ,drop=F]
      cv.fit <- cv.glmnet(x=cbind(zfit,denX[1:length(y),]), y=y, family = "binomial", nfolds=nfolds) 
      yp <- predict(cv.fit,newx=cbind(zPre,denX[-(1:length(y)),]), s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      coefs <- which(coef(cv.fit,s="lambda.min")[-(1:(1+ncol(zFit)))] !=0)
      vset <- c(vset,coefs)     
    } 
    
    ## y的平均预测(step 4):
    yPre <- rowMeans(preY) 
    numb <- table(vset)
    Vid <- as.numeric(rownames(numb))
    
    ## 统计权重(step 5)：
    Eset <- cbind(EDGE[Vid,],numb)
    Eset <- Eset[order(Eset[,3],decreasing=T),]
    colnames(Eset) <- c("row","col","numb")
  } 
  list(yPre=yPre,Eset=Eset) 
}

# 包含更多情况的JDINAC：
#############################################
jdinac <- function(EDGE,classLabel,DataFit,DataPre,zFit=NULL,zPre=NULL,nsplit=10,nfolds=5){
  
  ## 如果缺少测试集，那么测试集=训练集：
  if(missing(DataPre)) {
    DataPre <- DataFit
    zPre <- zFit
  } 
  ## 如果zFit和zPre非空：
  if(!is.null(zFit) & !is.null(zPre)){
    Pre.Z <- jdinac.z(EDGE=EDGE,classLabel=classLabel,DataFit=DataFit,DataPre=DataPre, zFit=zFit,zPre=zPre,nsplit=nsplit,nfolds=nfolds)
    ## y的类别预测值
    yPre <- Pre.Z$yPre 
    ## 具有依赖协同作用的成对基因的差异化权重
    Eset <- Pre.Z$Eset 
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
      index0=which(classLabel==0)
      index1=which(classLabel==1)
      sn0 <- round(size0/2)
      sn1 <- round(size1/2)
      splitid <- c(sample(index0,sn0),sample(index1,sn1))
      
      cLabel <- classLabel[splitid]        
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[splitid,],newData=rbind(DataFit[-splitid,],DataPre))
      y <- classLabel[-splitid]
      cv.fit <- cv.glmnet(x=denX[1:length(y),], y=y, family = "binomial", nfolds=nfolds) 
      
      yp <- predict(cv.fit,newx=denX[-(1:length(y)),], s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      coefs <- which(coef(cv.fit,s="lambda.min")[-1] !=0)
      vset <- c(vset,coefs)
      ##
      coes<-coef(cv.fit,s="lambda.min")[-1]  #除截距项外的系数
      var=cbind(coefs,rep(2*i-1,length(coefs)),coes[coefs]) #特征索引、次数、系数
      vars=rbind(vars,var)
      ##
    
      cLabel <- classLabel[-splitid]        
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[-splitid,],
                    newData=rbind(DataFit[splitid,],DataPre))
      y <- classLabel[splitid]
      cv.fit <- cv.glmnet(x=denX[1:length(y),], y=y, family = "binomial", nfolds=nfolds) 
     
      yp <- predict(cv.fit,newx=denX[-(1:length(y)),], s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      coefs <- which(coef(cv.fit,s="lambda.min")[-1] !=0)
      vset <- c(vset,coefs)   
      ##
      coes<-coef(cv.fit,s="lambda.min")[-1]  #除截距项外的系数
      var=cbind(coefs,rep(2*i,length(coefs)),coes[coefs]) #特征索引、次数、系数
      vars=rbind(vars,var)
      ##
    } 
    yPre <- rowMeans(preY) 
    numb <- table(vset)
    Vid <- as.numeric(rownames(numb))  
    Eset <- cbind(EDGE[Vid,],numb)
    Eset <- Eset[order(Eset[,3],decreasing=T),]
    colnames(Eset) <- c("row","col","numb")
    ##
    Vars=cbind(vars,EDGE[vars[,1],])
    Vars=Vars[,-1]
    colnames(Vars)=c("t","coef","row","col")
    ##
  } 
  list(yPre=yPre,Eset=Eset,Vars=Vars,preY=preY) 
}



jdinac.p <- function(EDGE,classLabel,DataFit,DataPre,zFit=NULL,zPre=NULL,nsplit,nfolds=5){
  
  ## 如果缺少测试集，那么测试集=训练集：
  if(missing(DataPre)) {
    DataPre <- DataFit
    zPre <- zFit
  } 
  ## 如果zFit和zPre非空：
  if(!is.null(zFit) & !is.null(zPre)){
    Pre.Z <- jdinac.z(EDGE=EDGE,classLabel=classLabel,DataFit=DataFit,DataPre=DataPre, zFit=zFit,zPre=zPre,nsplit=nsplit,nfolds=nfolds)
    ## y的类别预测值
    yPre <- Pre.Z$yPre 
    ## 具有依赖协同作用的成对基因的差异化权重
    Eset <- Pre.Z$Eset 
  } 
  else if(nsplit<=0){       
    stop("nsplit must be positive integer")    
  } 
  else {
    preY <- NULL
    vars=NULL
    vs=colnames(DataFit[,-ncol(DataFit)])
    vx1=vs[EDGE[,1]];vx2=vs[EDGE[,2]]
    cnames=paste0("(",vx1,",",vx2,")")
      size0 <- sum(classLabel==0)
      size1 <- sum(classLabel==1)
      index0=which(classLabel==0)
      index1=which(classLabel==1)
      sn0 <- round(size0/2)
      sn1 <- round(size1/2)
      set.seed(nsplit+123)
      splitid <- c(sample(index0,sn0),sample(index1,sn1))
      
      cLabel <- classLabel[splitid]        
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[splitid,],newData=rbind(DataFit[-splitid,],DataPre))
      colnames(denX)=cnames
      y <- classLabel[-splitid]
      cv.fit <- cv.glmnet(x=denX[1:length(y),], y=y, family = "binomial", nfolds=nfolds) 
      
      yp <- predict(cv.fit,newx=denX[-(1:length(y)),], s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      #coefs <- which(coef(cv.fit,s="lambda.min")[-1] !=0)
      #vset <- c(vset,coefs)
      ##
      coes<-coef(cv.fit,s="lambda.min")  #除截距项外的系数
      vars=cbind(vars,coes)
      ##
      
      cLabel <- classLabel[-splitid]        
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[-splitid,],
                    newData=rbind(DataFit[splitid,],DataPre))
      colnames(denX)=cnames
      y <- classLabel[splitid]
      cv.fit <- cv.glmnet(x=denX[1:length(y),], y=y, family = "binomial", nfolds=nfolds) 
      
      yp <- predict(cv.fit,newx=denX[-(1:length(y)),], s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      #coefs <- which(coef(cv.fit,s="lambda.min")[-1] !=0)
      #vset <- c(vset,coefs)   
      ##
      coes<-coef(cv.fit,s="lambda.min")  #除截距项外的系数
      vars=cbind(vars,coes)
      ##
    
  } 
  list(Vars=vars,preY=preY) 
}






