# HC-JDINAC

This package provides the R function to implement the HC-JDINAC method for micro trace evidence examination.  

HC.R - the function to implement Higher Criticism Thresholding feature selection.
JDINAC_b.R - the function to implement JDINAC with binned bivariate kernel density estimation. 
JDINAC_m.R - the function to implement JDINAC with with bivariate kernel density estimation using bandwidth matrix.

Edata.csv - FTIR spectral dataset of 41 handlebar grip individuals; column "y" represents individual labels, the other columns represent varibles; 15 observations for each individual. 
HC feature.csv - feature selection result of HC pairwise test on Edata.csv. 
individual.csv - dataset for individual HC-JDINAC analysis; containing training data and testing data.
community.csv - dataset for community HC-JDINAC analysis; containing training data and testing data.
community-individual.csv - dataset for community HC-JDINAC analysis using sample pair E13,E16's differential interactions.


