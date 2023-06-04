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
