
rrrBayes <- function(X,Y,Inits,X.new,R=1,lambda=0,Samples=1000, thin=1,seed=0){
   if(seed>0) set.seed(seed)
#  lambda = sum(Y^2)/sum(X^2)*lambda
  U = Inits$U
  V = Inits$V
  L = length(dim(X))-1
  N=dim(X)[1]
  P = dim(X)[2:(L+1)]
  M = length(dim(Y))-1
  Q = dim(Y)[2:(M+1)]
  Xmat = array(X,dim=c(N,prod(P)))
  Ymat = array(Y,dim=c(N,prod(Q)))
  Xmat.new = array(X.new,dim=c(dim(X.new)[1],prod(P)))
  Yvec = as.vector(Y)
  Vmat = matrix(nrow=prod(Q),ncol=R)
  for(r in 1:R){
    Vr = lapply(V, function(x) x[,r])
    Vmat[,r] = as.vector(array(apply(expand.grid(Vr), 1, prod), dim=Q))
  }
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
  sig2 = c()
  tau2 = c()
  Ypred_samps = array(dim=(c(floor(Samples/thin), dim(Xmat.new)[1]*prod(Q))))
  for(j in 1:Samples){
  sig2[j] = 1/rgamma(1,length(Yvec)/2,.5*sum((Yvec-Y_pred)^2))
    for(l in 1:L){
      ###Matrix B
      B_red = matrix(nrow=prod(P[-l]),ncol=R)
      for(r in 1:R){
        Ured_r <- lapply(U[-l], function(x) x[,r]) 
        B_red[,r] = as.vector(array(apply(expand.grid(Ured_r), 1, prod), dim=P[-l])) 
      }
      X_red = matrix(nrow=N,ncol=P[l]*r)
      for(i in 1:N){X_red[i,] = as.vector(Xarrays[[l]][i,,] %*% B_red)} 
      C = matrix(nrow=dim(X_red)[1]*prod(Q),ncol=dim(X_red)[2])
      for(r in 1:R){
        index = ((r-1)*P[l]+1):(r*P[l])
        C[,index] = kronecker(Vmat[,r],X_red[,index])
      }
      lambdaMat = kronecker(t(Vmat)%*%Vmat * t(B_red)%*%B_red,lambda*diag(P[l]))
      CP <- crossprod(X_red)*kronecker(t(Vmat)%*%Vmat, matrix(rep(1,P[l]^2),nrow=P[l]))
      regMat = solve(CP+lambdaMat) 
      UlvecMean = regMat%*%(t(C)%*%Yvec)
      UlvecVar = sig2[j]*regMat
      Ulvec = mvrnorm(1,UlvecMean,UlvecVar)
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
      regMat = solve(t(D)%*%D+lambda*t(BRmat)%*%BRmat)
      V[[m]] = Ymats[[m]]%*%D%*%regMat 
      V[[m]] = V[[m]] + mvrnorm(Q[m],rep(0,R),sig2[j]*regMat)
      for(r in 1:R){
        Vr = lapply(V, function(x) x[,r])
        Vmat[,r] = as.vector(array(apply(expand.grid(Vr), 1, prod), dim=Q))
      }
      Y_pred = as.vector(Xmat%*%t(Bmat)%*%t(Vmat))
      sse[(M+L)*(j-1)+L+m] = sum((Yvec-Y_pred)^2)
      sseR[(M+L)*(j-1)+L+m] = sum((Yvec-Y_pred)^2)+lambda*sum((t(Bmat)%*%t(Vmat))^2)
    }    
     ##standardization ###scale to have each vector the same norm.
      if(j%%thin==0){  
        pred.mean = as.vector(Xmat.new%*%t(Bmat)%*%t(Vmat))
        Ypred_samps[floor(j/thin),] = as.vector(Xmat.new%*%t(Bmat)%*%t(Vmat))+rnorm(length(pred.mean),0,sqrt(sig2[j]))
      }
      ###scale to avoid potential singularity issues  
      colNorms = matrix(nrow = L+M, ncol = R)
      for(l in 1:L){ colNorms[l,] = sqrt(colSums(U[[l]]^2))}
      for(m in 1:M){ colNorms[L+m,] = sqrt(colSums(V[[m]]^2))}
      scaleFactor = apply(colNorms, 2, prod)^(1/(L+M))
      for(l in 1:L){ U[[l]] = scale(U[[l]],center=FALSE,scale=colNorms[l,]/scaleFactor)}
      for(m in 1:M){ V[[m]] = scale(V[[m]],center=FALSE,scale=colNorms[L+m,]/scaleFactor)}
      #update Vmat
      for(r in 1:R){
        Vr = lapply(V, function(x) x[,r])
        Vmat[,r] = as.vector(array(apply(expand.grid(Vr), 1, prod), dim=Q))
      }
  }
  Ypred_samps <- array(Ypred_samps,dim=c(floor(Samples/thin),dim(Xmat.new)[1],Q) )
  return(Ypred_samps)
}



