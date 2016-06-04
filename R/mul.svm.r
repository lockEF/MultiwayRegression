mul.svm=function(xmul,y,R=50,rank=1,C=1)
{ 
  if (rank==1)
  {
    res=matrix(rep(NA,R),ncol=1)
    res[,1]=apply(as.matrix(c(1:R)),1,FUN=function(t){mul.svmstep(xmul,y,t,rank=1,C=C)[[4]]})
    sel=which(res[,1]==max(res[,1]))
    out=mul.svmstep(xmul,y,sel,rank=1,C=C)
    return(list('beta'=out$beta,'w'=out$w,'v'=out$v,'int'=out$int,'obj'=out$obj))
  }
  if (rank>1)
  {
    res=matrix(rep(NA,R),ncol=1)
    res[,1]=apply(as.matrix(c(1:R)),1,FUN=function(t){mul.svmstep(xmul,y,t,rank=rank,C=C)[[3]]})
    sel=which(res[,1]==max(res[,1]))
    out=mul.svmstep(xmul,y,sel,rank=rank,C=C)
    return(list('beta'=out$beta,'int'=out$int,'obj'=out$obj))
  }
}

mul.svmstep=function(xmul,y,seed,rank,C)
{
  if (rank==1)
  {
    set.seed(seed)
    n1=dim(xmul)[1]
    n2=dim(xmul)[2]
    n3=dim(xmul)[3]
    v=runif(n3,0,1)
    w=runif(n2,0,1)
    dataw=matrix(0,nrow=n1,ncol=n2)
    datav=matrix(0,nrow=n1,ncol=n3)
    ite=0
    diff=1
    y=as.factor(y)
    while (diff>0.001)
    {
      ite=ite+1
      vpre=v
      wpre=w
      betapre=kronecker(v,w)
      ###fix v, update w
      v=v/sqrt(sum(v^2))
      for (i in 1:n1)
      {
        dataw[i,]=xmul[i,,]%*%v
      }
      traw=cbind.data.frame(y,dataw)
      fitw=ksvm(y ~., data=traw,scaled=FALSE,kernel="vanilladot",C=C)
      w=as.vector(t(fitw@coef[[1]]) %*% as.matrix(traw[fitw@alphaindex[[1]],-1]))
      
      intw= - fitw@b #rho is the negative intercept!!
      ###fix w, update v
      w=w/sqrt(sum(w^2))
      for (i in 1:n1)
      {
        datav[i,]=t(t(xmul[i,,])%*%w)
      }
      trav=cbind.data.frame(y,datav)
      fitv=ksvm(y ~., data=trav,scaled=FALSE,kernel="vanilladot",C=C)
      v = as.vector(t(fitv@coef[[1]]) %*% as.matrix(trav[fitv@alphaindex[[1]],-1]))
      intv= - fitv@b #rho is the negative intercept!!
      objv=fitv@obj
      
      beta=kronecker(v,w)    
      diff=dist(rbind(as.vector(beta),as.vector(betapre)))
    }
    result=list('beta'=w%*%t(v),'w'=w,'v'=v,'int'=intv,'obj'=objv)
    return(result)
  }
  
  if (rank>1)
  {
    set.seed(seed)
    r=rank
    n1=dim(xmul)[1]
    n2=dim(xmul)[2]
    n3=dim(xmul)[3]
    w0=matrix(NA,nrow=r,ncol=n2)
    v0=matrix(NA,nrow=r,ncol=n3)
    for (i in 1:r)
    {
      w0[i,]=runif(n2,0,1)
    }
    for (i in 1:r)
    {
      v0[i,]=runif(n3,0,1)
    }
    dataw=matrix(0,nrow=n1,ncol=r*n2)
    datav=matrix(0,nrow=n1,ncol=r*n3)
    ite=0
    diffbeta=1
    #obj=c(0,0)
    y=as.factor(y)
    beta=matrix(0,ncol=n3,nrow=n2)
    for (i in 1:r)
    {
      beta=beta+w0[i,]%*%t(v0[i,])
    }
    while(diffbeta>0.001)
    {
      ite=ite+1
      betapre=beta
      
      decomp=svd(beta,r,r)
      vs=matrix(NA,nrow=r,ncol=n3)
      for (i in 1:r)
      {
        vs0=decomp$v[,i]*sqrt(decomp$d[i])
        vs[i,]=vs0/sqrt(sum(vs0^2))
      }
      
      ###fix v, update w
      for (i in 1:n1)
      {
        datawi=c()
        for (j in 1:r)
        {
          datawi=c(datawi,xmul[i,,]%*%vs[j,])
        }
        dataw[i,]=datawi 
      }
      traw=cbind.data.frame(y,dataw)
      fitw=ksvm(y ~., data=traw,scaled=FALSE,kernel="vanilladot",C=C)
      wall=as.vector(t(fitw@coef[[1]]) %*% as.matrix(traw[fitw@alphaindex[[1]],-1]))
      w=matrix(NA,nrow=r,ncol=n2)
      for (i in 1:r)
      {
        w[i,]=wall[((i-1)*n2+1):(i*n2)]
      }
      intw= - fitw@b #rho is the negative intercept!!
      #obj=c(obj,fitw@obj)
      
      beta=matrix(0,ncol=n3,nrow=n2)
      for (i in 1:r)
      {
        beta=beta+w[i,]%*%t(vs[i,])
      }
      
      decomp=svd(beta,r,r)
      ws=matrix(NA,nrow=r,ncol=n2)
      for (i in 1:r)
      {
        ws0=decomp$u[,i]*sqrt(decomp$d[i])
        ws[i,]=ws0/sqrt(sum(ws0^2))
      }
      
      ###fix w, update v
      for (i in 1:n1)
      {
        datavi=c()
        for (j in 1:r)
        {
          datavi=c(datavi,t(xmul[i,,])%*%ws[j,])
        }
        datav[i,]=datavi 
      }
      trav=cbind.data.frame(y,datav)
      fitv=ksvm(y ~., data=trav,scaled=FALSE,kernel="vanilladot",C=C)
      vall = as.vector(t(fitv@coef[[1]]) %*% as.matrix(trav[fitv@alphaindex[[1]],-1]))
      v=matrix(NA,nrow=r,ncol=n3)
      for (i in 1:r)
      {
        v[i,]=vall[((i-1)*n3+1):(i*n3)]
      }
      intv= - fitv@b #rho is the negative intercept!!
      
      beta=matrix(0,ncol=n3,nrow=n2)
      for (i in 1:r)
      {
        beta=beta+ws[i,]%*%t(v[i,])
      }
      #obj=c(obj,fitv@obj)
      objv=fitv@obj
      
      diffbeta=dist(rbind(as.vector(beta),as.vector(betapre)))
    }
    mulrankr=rep(0,n2*n3)
    for (i in 1:r)
    {
      mulrankr=mulrankr+kronecker(v[i,],ws[i,])
    }
    result=list('beta'=beta,'int'=intv,'obj'=objv)
    return(result)
  }
}