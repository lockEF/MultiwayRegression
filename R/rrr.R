
rrr <- function(X,Y,R=1,lambda=0,annealIter=0,convThresh=10^(-5), seed=0){
  if(seed>0) set.seed(seed)
  lambdaFin = lambda
  L = length(dim(X))-1
  N=dim(X)[1]
  P = dim(X)[2:(L+1)]
  M = length(dim(Y))-1
  Q = dim(Y)[2:(M+1)]
  Xmat = array(X,dim=c(N,prod(P)))
  Ymat = array(Y,dim=c(N,prod(Q)))
  Yvec = as.vector(Y)
  
  U = list()
  for(l in 1:L) U[[l]] = matrix(rnorm(P[l]*R),ncol=R)
  V = list()
  for(l in 1:M) V[[l]] = matrix(rnorm(Q[l]*R),ncol=R)
  Vmat = matrix(nrow=prod(Q),ncol=R)
  for(r in 1:R){
    Vr = lapply(V, function(x) x[,r])
    Vmat[,r] = as.vector(array(apply(expand.grid(Vr), 1, prod), dim=Q))
  }
#  B = array(apply(expand.grid(U), 1, prod), dim=P)
  Bmat = matrix(nrow=R,ncol=prod(P))
  for(r in 1:R){
    Ur = lapply(U, function(x) x[,r])
    Br = array(apply(expand.grid(Ur), 1, prod), dim=P)  
    Bmat[r,] = array(Br,dim=c(1,prod(P)))
  }
  Y_pred = as.vector(Xmat%*%t(Bmat)%*%t(Vmat))
  Xarrays = list()
  for(l in 1:L){
    Xarrays[[l]] <- array(dim=c(N,P[l],prod(P[-l])))
    perm = c(l,c(1:L)[-l])
    for(i in 1:N){
      X_slice = array(Xmat[i,],dim=c(P))
      X_slice_perm = aperm(X_slice,perm)
      Xarrays[[l]][i,,] = array(X_slice_perm,dim=c(P[l],prod(P[-l])))
    }
  }  
  Ymats = list()
  for(m in 1:M){
    perm = c(m+1,c(1:(M+1))[-(m+1)])
    Y_perm = aperm(Y,perm) 
    Ymats[[m]] = array(Y_perm,dim=c(Q[m],N*prod(Q[-m])))
  }
  sse=c()
  sseSig = c()
  sseR = c()
  j=0
  conv=convThresh+1
  while(conv>convThresh){
    j=j+1
    if(j<=annealIter){
      lambda = 100*(annealIter-j)/annealIter+lambdaFin
    }
    for(l in 1:L){
      ###Matrix B
      B_red = matrix(nrow=prod(P[-l]),ncol=R)
      for(r in 1:R){
        Ured_r <- lapply(U[-l], function(x) x[,r]) 
        B_red[,r] = as.vector(array(apply(expand.grid(Ured_r), 1, prod), dim=P[-l])) 
      }
      X_red = matrix(nrow=N,ncol=P[l]*r)
      perm = c(l,c(1:L)[-l])
      for(i in 1:N){X_red[i,] = as.vector(Xarrays[[l]][i,,] %*% B_red)} 
      C = matrix(nrow=dim(X_red)[1]*prod(Q),ncol=dim(X_red)[2])
      for(r in 1:R){
        index = ((r-1)*P[l]+1):(r*P[l])
        C[,index] = kronecker(Vmat[,r],X_red[,index])
      }
      lambdaMat = kronecker(t(Vmat)%*%Vmat * t(B_red)%*%B_red,lambda*diag(P[l]))
      CP <- crossprod(X_red)*kronecker(t(Vmat)%*%Vmat, matrix(rep(1,P[l]^2),nrow=P[l]))
      regMat = solve(CP+lambdaMat) 
      Ulvec <- regMat%*%(t(C)%*%Yvec)
      U[[l]] = matrix(Ulvec,nrow=P[l])
      Bmat = matrix(nrow=R,ncol=prod(P))
      for(r in 1:R){
        Ur = lapply(U, function(x) x[,r])
        Br = array(apply(expand.grid(Ur), 1, prod), dim=P)  
        Bmat[r,] = array(Br,dim=c(1,prod(P)))
      }
      Y_pred = as.vector(Xmat%*%t(Bmat)%*%t(Vmat))
      sse[(M+L)*(j-1)+l] = sum((Yvec-Y_pred)^2)
      sseR[(M+L)*(j-1)+l] = sum((Yvec-Y_pred)^2)+lambda*sum((t(Bmat)%*%t(Vmat))^2)
    }
    for(m in 1:M){
      BRmat = matrix(nrow = prod(c(P,Q[-m])),ncol=R)
      D = matrix(nrow=N*prod(Q[-m]),ncol=R)
      for(r in 1:R){
        Vecs <- lapply(c(U,V[-m]), function(x) x[,r])
        Br = array(apply(expand.grid(Vecs), 1, prod), dim=c(prod(P),prod(Q[-m])))
        BRmat[,r] = as.vector(Br)
        D[,r] = as.vector(Xmat %*% Br)
      }
      V[[m]] = Ymats[[m]]%*%D%*%solve(t(D)%*%D+lambda*t(BRmat)%*%BRmat) ##now for X, use Xorig
      Vmat = matrix(nrow=prod(Q),ncol=R)
      for(r in 1:R){
        Vr = lapply(V, function(x) x[,r])
        Vmat[,r] = as.vector(array(apply(expand.grid(Vr), 1, prod), dim=Q))
      }
      Y_pred = as.vector(Xmat%*%t(Bmat)%*%t(Vmat))
      sse[(M+L)*(j-1)+L+m] = sum((Yvec-Y_pred)^2)
      sseR[(M+L)*(j-1)+L+m] = sum((Yvec-Y_pred)^2)+lambda*sum((t(Bmat)%*%t(Vmat))^2)
    }
  if(j>1) conv = abs(sseR[length(sseR)]-sseR[length(sseR)-L-M])/sseR[length(sseR)] 
  }
  B = array(t(Bmat)%*%t(Vmat),dim = c(P,Q))
  return(list(U=U,V=V,B=B,sse=sse,sseR=sseR))
}



