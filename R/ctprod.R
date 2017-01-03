ctprod <- function(A,B,K){
  da <- dim(A)
  la <- length(da)
  db <- dim(B)
  lb <- length(db)
  Amat <- array(A, dim=c(prod(da[1:(la-K)]),prod(da[(la-K+1):la])))
  Bmat <- array(B, dim=c(prod(db[1:K]),prod(db[(K+1):lb])))
  Cmat <- Amat %*% Bmat
  C <- array(Cmat,dim=c(da[1:(la-K)],db[(K+1):lb]))
  return(C)
} 