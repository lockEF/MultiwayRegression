\name{MultiwayRegression-package}
\alias{MultiwayRegression}
\docType{package}
\title{
Perform tensor-on-tensor regression
}
\description{
Functions to predict one multi-way array (i.e., a tensor) from another multi-way array. 
}
\details{
\tabular{ll}{
Package: \tab MultiwayRegression-package\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2017-01-09\cr
License: \tab GPL-3\cr
}

}
\author{
Eric F. Lock

Maintainer: Eric F. Lock <elock@umn.edu>
}

\examples{
data(SimData) ##loads simulated X: 100 x 15 x 20 and Y: 100 x 5 x 10 
Results <- rrr(X,Y,R=2)  ##Fit rank 2 model with no regularization
Y_pred <- ctprod(X,Results$B,2)  ##Array of fitted values
}


\keyword{ package }