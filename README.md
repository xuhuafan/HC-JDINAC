# HC-JDINAC
This package provides the R function to implement the HC-JDINAC method for micro trace evidence examination, as used in the manuscript, "HC-JDINAC: Identifying Important Individual Features and Interaction Differences from Highly Similar FTIR Spectra".
  
## HC_evaluation.R
### Description
Perfrom HC-LDA, FAIR-LDA, SFLDA and FANS. Compute classification error (CE) and FDR like Table 3 in paper.
### Examples
```
##Given individual y1,y2, perform HC-LDA, FAIR-LDA, SFLDA and FANS.
data=read.csv("Edata.csv")
y1=13;y2=16  ## perform binary classification between any two individauls by changing y1,y2
result_hc=HC_LDA(y1,y2,data,alpha=0.1)  ## HC-LDA result
result_hc$res  ##testing classification error(CE) and FDR
result_fair=FAIR_LDA(y1,y2,data,alpha=0.1)  ## FAIR-LDA result
result_fair$res  ##testing classification error(CE) and FDR
result_sflda=SFLDA(y1,y2,data)  ## SFLDA result
result_sflda$res  ##testing classification error(CE) and FDR
result_fans=FANS(y1,y2,data,nsplit=10,nfolds=5)  ## FANS result
result_fans$res  ##testing classification error(CE) and FDR
```

## JDINAC_b.R
### Description
Perform JDINAC with binned bivariate kernel density estimation using kernel $K_\textbf{h}$ applied to the rescaling of the product kernel $K(\textbf{x}) = K(x_1)K(x_2)$ by the vector of bandwidths $\textbf{h} = (h_1, h_2)$.
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
gp=group_plot(eset)
p1=gp$p1  #dumbbell diagram like Fig.5 in paper
p2=gp$p2  #bar chart like Fig.5 in paper

##perform HC-JDINAC on given two communities C0,C1
C0=c(3,13,26,29);C1=c(16,21,30,44)
data=read.csv("community.csv")
dfh=read.csv("HC feature.csv")  ##HC feature selection result
result=jdinac_bd(y1,y2,data,dfh,k0=1,k1=1)
eset=result$eset ## differential edge/interaction selected by JDINAC
err=result$err ## classification error
p=diffnet_plot(eset) ## differential interaction network

##perform communitities C0,C1 identification using differential interactions of sample pair y1,y2. (see Table 6 in paper)
y1=13;y2=16
data=read.csv("individual.csv")
dfh=read.csv("HC feature.csv")  ##HC feature selection result
eset1=jdinac_bd(y1,y2,data,dfh,k0=1,k1=1)$eset  ## individual differential edges/interactions
C0=c(3,13,26,29);C1=c(16,21,30,44)
dt=read.csv("community-individual.csv")
res=com_ind(eset1,C0,C1,dt)  ##testing classification error and confusing rate
```

## JDINAC_m.R
### Description
Perform JDINAC with bivariate kernel density estimation using bandwidth matrix with kernel $K$ and bandwidth matrix $\textbf{H}$, where $K_{\textbf{H}}(\textbf{x}) = |\textbf{H}|^{-1/2}K(\textbf{H}^{-1/2}\textbf{x})$ and $\textbf{x}=(x_1,x_2)$.
### Examples
```
##perform HC-JDINAC on given two individuals y1,y2
y1=13;y2=16
data=read.csv("individual.csv")
dfh=read.csv("HC feature.csv")  ##HC feature selection result
result=jdinac_md(y1,y2,data,dfh)
eset=result$eset ## differential edge/interaction selected by JDINAC
err=result$err ## classification error
p=diffnet_plot(eset) ## differential interaction network
gp=group_plot(eset)
p1=gp$p1  #dumbbell diagram like Fig.5 in paper
p2=gp$p2  #bar chart like Fig.5 in paper
```




