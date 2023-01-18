#' Perform HC feature selection
#'
#' \code{HC} returns the HC feature selection result.
#'
#' @param y1 One sample label to be tested by HC
#' @param y2 The other sample label to be tested by HC
#' @param x Data matrix
#' @param alpha Threshold value of t test in HC
#'
#' @return The output will be a vector of 0,1 with a length equal to the number of variables, where 1 indicates that the variable is valid.
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
