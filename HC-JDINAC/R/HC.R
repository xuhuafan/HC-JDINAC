#' @title Perfrom Higher Criticism Thresholding feature selection
#' 
#' @param x1 data matrices for first class containing one sample per row, one variable per column. 
#' @param x2 data matrices for second class containing one sample per row, one variable per column. 
#' @param alpha: significance level of two-sample test.
#'
#' @return: The output will be a list containing the following components:
#' @export res: data frame, the first column is the variable index; the second column through fourth columns are p-value, Z-statistic, the absolute value of the Z-statistic and HC objective score, respectively; the last column is feature selection indicator of 0,1, where 1 indicates the corresponding feature is selected by HCT and 0 indicates it is not selected.
#' @export th_z: the higher criticism threshold t^HC.
#' @examples ### Perform HCT on sample pair E13,E16.
#' data=read.csv("Edata.csv")
#' data=data.frame(na.omit(data))
#' y1=13;y2=16
#' x1=data[data$y == y1,-ncol(data)]
#' x2=data[data$y == y2,-ncol(data)]
#' HCT_res=HCT(x1,x2,alpha=0.1)$res
#' sel=HC_res$sel
#' vs=colnames(data)[which(sel==1)]  #selected features


HCT=function(x1,x2,alpha=0.1){
  Zs=NULL
  p_t=NULL
  for(i in 1:(ncol(x1)-1)){
    z1=x1[,i]
    z2=x2[,i]
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
  res=data.frame(1:ncol(x1),p_t,Zs)
  ## 计算p值：
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
  v=ifelse(res$abs_z >= th_z,1,0)
  res["v"]=v
  colnames(res)=c("col","p-value","Z-stat","abs_Z","HC","sel")
  return(list(res=res,th_z=th_z))
}


