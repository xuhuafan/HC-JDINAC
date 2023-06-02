# HC-JDINAC

This package provides the R function to implement the HC-JDINAC method for micro trace evidence examination.  

## HC
### Description
Perfrom Higher Criticism Thresholding feature selection
### Usage
HCT(x1,x2,alpha=0.1)
### Arguments
* x1: data matrices for first class containing one sample per row, one variable per column. 
* x2: data matrices for second class containing one sample per row, one variable per column. 
* alpha: significance level of two-sample test.

### Value
The output will be a list containing the following components: \
* res: data frame, the first column is the variable index; the second column through fourth columns are p-value, Z-statistic, the absolute value of the Z-statistic and HC objective score, respectively; the last column is feature selection indicator of 0,1, where 1 indicates the corresponding feature is selected by HCT and 0 indicates it is not selected. \
* th_z: the higher criticism threshold t^HC. \
### References

### Examples
'''
##Perform HCT on sample pair E13,E16.
data=read.csv("Edata.csv")
data=data.frame(na.omit(data))
y1=13;y2=16
x1=data[data$y == y1,-ncol(data)]
x2=data[data$y == y2,-ncol(data)]
HCT_res=HCT(x1,x2,alpha=0.1)$res
sel=HC_res$sel
vs=colnames(data)[which(sel==1)]  #selected features

##Perform HC-LDA on sample pair E13,E16.
data=read.csv("Edata.csv")
data=data.frame(na.omit(data))
y1=13;y2=16
dt=data[data$y %in% c(y1,y2),]
ind1=which(dt$y == y1)
ind2=which(dt$y == y2)
##8:7 split train set and test set
sam_ind1=sample(ind1,ceiling(length(ind1)*0.5))
sam_ind2=sample(ind2,ceiling(length(ind2)*0.5))
x1=dt[sam_ind1,-ncol(dt)]
x2=dt[sam_ind2,-ncol(dt)]
HCT_res=HCT(x1,x2,alpha=0.1)$res
sel=HC_res$sel;sel_col=which(sel==1)
train=dt[c(sam_ind1,sam_ind2),c(sel_col,ncol(dt))]
test=dt[-c(sam_ind1,sam_ind2),c(sel_col,ncol(dt))]
lda.fit=lda(y~., data=train)
lda.pred=predict(lda.fit,test) 
lda.class=lda.pred$class
err=mean(lda.class!=test$y)  ##classification error
'''


HC.R - the function to implement Higher Criticism Thresholding feature selection.
JDINAC_b.R - the function to implement JDINAC with binned bivariate kernel density estimation. 
JDINAC_m.R - the function to implement JDINAC with with bivariate kernel density estimation using bandwidth matrix.

Edata.csv - FTIR spectral dataset of 41 handlebar grip individuals; column "y" represents individual labels, the other columns represent varibles; 15 observations for each individual. 
HC feature.csv - feature selection result of HC pairwise test on Edata.csv. 
individual.csv - dataset for individual HC-JDINAC analysis; containing training data and testing data.
community.csv - dataset for community HC-JDINAC analysis; containing training data and testing data.
community-individual.csv - dataset for community HC-JDINAC analysis using sample pair E13,E16's differential interactions.


