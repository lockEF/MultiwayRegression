\name{rrrBayes}
\alias{rrrBayes}
\title{
Bayesian inference for reduced rank regression
}
\description{
Performs Bayesian inference for a linear model to estimate one multi-way array from another, under the restriction that the coefficient array has given PARAFAC rank.  By default, estimates are chosen to minimize a least-squares objective;  an optional penalty term allows for $L_2$ regularization of the coefficient array.  
}
\usage{
rrrBayes(X,Y,Inits,Xmat.new,R=1,lambda=0,Samples=1000, thin=1,seed=0)
}
\arguments{
  \item{X}{
  	A predictor array of dimension N x P_1 x ... x P_L for the training data.  
}
  \item{Y}{
  	An outcome array of dimension N x Q_1 x ... x Q_M for the training data.  
}
  \item{Inits}{
  	Initial values.  Inits$U gives a list of length L where Inits$U[[l]]: P_l x R gives the coefficient basis for the l'th mode of X.  Inits$V gives a list of length M where Inits$V[[m]]: Q_m x R gives the coefficient basis for the m'th mode of Y. Can be the output of rrr(...).      
}
 \item{Xmat.new}{
  	Predictor array of dimension M x P_1 x ... x P_L.  Each row gives the entries for a new  P_1 x ... x P_L predictor observation in vectorized form.        
}
  \item{R}{
  	Assumed rank of the P_1 x ... x P_L x Q_1 x ... x Q_M coefficient array. 
}
  \item{lambda}{
  	Ridge ($L_2$) penalty parameter for the coefficient array, inversely proportional to the variance of the coefficients under a Gaussian prior.  
}
 \item{Samples}{
   Length of the MCMC sampling chain.
}
 \item{thin}{
  	Thinning value, for thin=j, only every j'th observation in the MCMC chain is saved.
}
 \item{seed}{
  	Random seed for generation of initial values.
}
}
\value{
%%  ~Describe the value returned 
An array of dimension (Samples/thin) x M x Q_1 x ... x Q_M, giving (Samples/thin) samples from the posterior predictive of the outcome array predicted by Xmat.new.  
}

\references{
}
\author{
Eric F. Lock 
}




