\name{rrr}
\alias{rrr}
\title{
Penalized reduced rank regression for tensors
}
\description{
Fits a linear model to estimate one multi-way array from another, under the restriction that the coefficient array has given PARAFAC rank.  By default, estimates are chosen to minimize a least-squares objective;  an optional penalty term allows for $L_2$ regularization of the coefficient array.  
}
\usage{
rrr(X,Y,R=1,lambda=0,annealIter=0,convThresh=10^(-5), seed=0)
}
\arguments{
  \item{X}{
  	A predictor array of dimension N x P_1 x ... x P_L.  
}
  \item{Y}{
  	An outcome array of dimension N x Q_1 x ... X Q_M.  
}
  \item{R}{
  	Assumed rank of the P_1 x ... x P_L x Q_1 x ... x Q_M coefficient array. 
}
  \item{lambda}{
  	Ridge ($L_2$) penalty parameter for the coefficient array.
}
 \item{annealIter}{
   Number of tempering iterations to improve initialization 
}
 \item{convThresh}{
  	Converge threshold for the absolute difference in the objective function between two iterations
}
 \item{seed}{
  	Random seed for generation of initial values.
}
}
\value{
%%  ~Describe the value returned 
\item{U}{
 List of length L. U[[l]]: P_l x R gives the coefficient basis for the l'th mode of X.}
\item{V}{
	List of length M. V[[m]]: Q_m x R gives the coefficient basis for the m'th mode of Y.}
\item{sse}{
	Vector giving the sum of squared residuals at each iteration.}
\item{sseR}{
	Vector giving the value of the objective (sse+penalty) at each iteration.}
}

\references{
}
\author{
Eric F. Lock 
}

\examples{
data(SimData) ##loads simulated X: 100 x 15 x 20 and Y: 100 x 5 x 10 
Results <- rrr(X,Y,R=2)
}


