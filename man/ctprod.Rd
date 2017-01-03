\name{ctprod}
\alias{ctprod}
\title{
Compute the contracted tensor product between two multiway arrays.  
}
\description{
Computes the contracted tensor product between two multiway arrays.
}
\usage{
ctprod(A,B,K)
}
\arguments{
  \item{A}{
An array of dimension P_1 x ... x P_L x R_1 x ... x R_K.
}
 \item{B}{
An array of dimension R_1 x ... x R_K x Q_1 x ... x Q_M.
}
\item{K}{
A positive integer, giving the number of modes to collapse.  
}
}
\value{
An array C of dimension P_1 x ... x P_L x Q_1 x ... x Q_M, given by the contracted tensor product of A and B.  
}

\author{
Eric F. Lock
}



