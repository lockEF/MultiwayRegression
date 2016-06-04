\name{MultiwayClassification-package}
\alias{MultiwayClassification}
\docType{package}
\title{
Linear Classification of Multi-Way Data
}
\description{
Description: Performs linear classification for data with multi-way structure.  The distance-weighted discrimination (DWD) or support vector machine (SVM) classification objectives are optimized under the assumption that the multi-way coefficients have low rank. 
}
\details{
\tabular{ll}{
Package: \tab MultiwayClassification-package\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-06-03\cr
License: \tab GPL-3\cr
}

}
\author{
Eric F. Lock and Tianmeng Lyu

Maintainer: Eric F. Lock <elock@umn.edu>
}

\examples{
	data(IFNB_Data) ##Load gene expression time course data (?IFNB_Data for more info)
	
##multi-way DWD	
results.mw <- mul.dwd(DataArray,y=Class, rank=1) #estimate rank 1 model
names(results.mw)
##Compute projection onto the classification direction for each individual:
DWD_scores <- c()
for(i in 1:length(Class)) DWD_scores[i] = sum(DataArray[i,,]*results.mw$beta)+results.mw$int

##multi-way SVM
results.mw <- mul.svm(DataArray,y=Class,R=20, rank=1) #estimate rank 1 model
##Compute projection onto the classification direction for each individual:
SVM_scores <- c()
for(i in 1:length(Class)) SVM_scores[i] = sum(DataArray[i,,]*results.mw$beta)+results.mw$int
}


\keyword{ package }