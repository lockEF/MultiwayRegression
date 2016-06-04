\name{mul.dwd}
\alias{mul.dwd}
\title{
Multi-way Distance Weighted Discrimination
}
\description{
Optimizes the DWD objective for a multi-way dataset and a vector of class labels, under the assumption that the coefficient array is of low-rank.
}
\usage{
mul.dwd(xmul,y,rank=1,C=100)
}
\arguments{
  \item{xmul}{
  	An array of dimensions N x P x M.  The initial dimension (N) gives the cases to be classified.
}
  \item{y}{
  	A vector of length N giving the class label ('-1' or '1') for each case.  
}
  \item{rank}{
  	Assumed rank of the P x M coefficient matrix.  
}
  \item{C}{
  	Penalty parameter for the DWD objective (adjusted by median pairwise distance)
}
}
\value{
%%  ~Describe the value returned 
\item{beta}{
	P x M matrix of coefficients}
\item{int}{
	Intercept}
If r=1:
\item{w}{
	Vector of length P, giving weights in the second dimension}
\item{v}{
	Vector of length M, giving weights in the third dimension}
}

\references{
Hanwen, H., Lu, X., Liu, Y.,  Haaland, P., & Marron, J. S. (2012). R/DWD: distance-weighted discrimination for classification, visualization and batch adjustment. Bioinformatics, 28(8):1182-3.
}
\author{
Eric F. Lock, Tianmeng Lyu, and Lynn E. Eberly
}

\seealso{
\code{\link{kdwd}}
}
\examples{
data(IFNB_Data) ##Load gene expression time course data (?IFNB_Data for more info)
results.mw <- mul.dwd(DataArray,y=Class, rank=1) #estimate rank 1 model
names(results.mw)
##Compute projection onto the classification direction for each individual:
DWD_scores <- c()
for(i in 1:length(Class)) DWD_scores[i] = sum(DataArray[i,,]*results.mw$beta)+results.mw$int
}


