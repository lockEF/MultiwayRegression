\name{mulstan}
\alias{mulstan}
\title{
Standardize a Multi-way Array
}
\description{
Standardizes the first dimension of a three-way array.
}
\usage{
mulstan(multi)
}
\arguments{
  \item{xmul}{
An array of dimensions N x P x M. 
}
}
\value{
A standardized array (multi[,i,j] is scaled to have mean 0 and standard deviation 1 for all i,j).
}

\author{
Eric F. Lock, Tianmeng Lyu, and Lynn E. Eberly
}



