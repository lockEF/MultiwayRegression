\name{mul.svm}
\alias{mul.svm}
\title{
Multi-way Support Vector Machines
}
\description{
Optimizes the SVM objective for a multi-way dataset and a vector of class labels, under the assumption that the coefficient array is of low-rank.
}
\usage{
mul.svm=function(xmul,y,R=50,rank=1,C=1)
}
\arguments{
  \item{xmul}{
	An array of dimensions N x P x M.  The initial dimension (N) gives the cases to be classified.
}
  \item{y}{
	A vector of length N giving the class label ('-1' or '1') for each case.  
}
\item{R}{
	Number of intitializations.
}
  \item{rank}{
	Assumed rank of the P x M coefficient matrix.  
}
  \item{C}{
	Penalty parameter for the SVM objective.
}
}
\value{
%%  ~Describe the value returned 
\item{beta}{
	P x M matrix of coefficients}
\item{int}{
	Intercept}
\item{obj}{
	Value of the SVM objective}
If r=1:
\item{w}{
	Vector of length P, giving weights in the second dimension}
\item{v}{
	Vector of length M, giving weights in the third dimension}
}

\references{
Karatzoglou, A., Meyer, D., & Hornik, K. (2006).   Support Vector Machines in R. Journal of Statistical Software, 15(9):1-28.
}
\author{
Eric F. Lock, Tianmeng Lyu, and Lynn E. Eberly
}

\seealso{
\code{\link{ksvm}}
}
\examples{
data(IFNB_Data) ##Load gene expression time course data (?IFNB_Data for more info)
results.mw <- mul.svm(DataArray,y=Class,R=20, rank=1) #estimate rank 1 model
##Compute projection onto the classification direction for each individual:
SVM_scores <- c()
for(i in 1:length(Class)) SVM_scores[i] = sum(DataArray[i,,]*results.mw$beta)+results.mw$int
}


