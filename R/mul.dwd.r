mul.dwd=function(xmul,y,rank=1,C=100)
{
  n1=dim(xmul)[1]
  n2=dim(xmul)[2]
  n3=dim(xmul)[3]
  xfull <- matrix(nrow=n1,ncol=n2*n3)
  for(i in 1:n1){
    xfull[i,] =  as.vector(xmul[i,,])
  }
  if (rank==1)
  {
    k=sum(y==-1)
    D=median(as.vector(as.matrix(dist(xfull))[(k+1):nrow(xfull),1:k]))
    v0=runif(n3,0,1)
    w0=runif(n2,0,1)
    v=v0/sqrt(sum(v0^2))
    w=w0/sqrt(sum(w0^2)) 
    dataw=matrix(0,nrow=n1,ncol=n2)
    datav=matrix(0,nrow=n1,ncol=n3)
    ite=0
    diff=1
    y=as.factor(y)
    while(diff>0.001)
    {
      ite=ite+1
      vpre=v
      wpre=w
      betapre=kronecker(vpre, wpre)
      ###fix v, update w
      for (i in 1:n1)
      {
        dataw[i,]=xmul[i,,]%*%v
      }
      traw=cbind.data.frame(y,dataw)
      d=median(as.vector(as.matrix(dist(dataw))[(k+1):nrow(dataw),1:k]))
      c=C*d^2/D^2
      fitw=kdwd(y~.,data=traw,type="bdwd",C=c,scaled=FALSE)
      #fitw@w would be a matrix if n<p, so need to transform it into a vector.
      w=as.vector(fitw@w[[1]])
      intw=fitw@b0[[1]]
      #dist():compute the distances between the rows of a data matrix
      ###fix w, update v
      for (i in 1:n1)
      {
        datav[i,]=t(t(xmul[i,,])%*%w)
      }
      trav=cbind.data.frame(y,datav)
      d=median(as.vector(as.matrix(dist(datav))[(k+1):nrow(datav),1:k]))
      c=C*d^2/D^2
      fitv=kdwd(y~.,data=trav,type="bdwd",C=c,scaled=FALSE)
      v=as.vector(fitv@w[[1]])
      intv=fitv@b0[[1]]
      diffint=dist(rbind(intw,intv))
      beta=kronecker(v, w)
      diff=dist(rbind(as.vector(beta),as.vector(betapre)))
    }
    result=list('beta'=w%*%t(v),'w'=w,'v'=v,'int'=intv)
    return(result)
  }
  
  if (rank>1)
  {
    r=rank
    x=xfull
    k=sum(y==-1)
    D=median(as.vector(as.matrix(dist(x))[(k+1):nrow(x),1:k]))
    w0=matrix(NA,nrow=r,ncol=n2)
    v0=matrix(NA,nrow=r,ncol=n3)
    for (i in 1:r)
    {
      w0[i,]=runif(n2,0,1)
      v0[i,]=runif(n3,0,1)
    }
    dataw=matrix(0,nrow=n1,ncol=r*n2)
    datav=matrix(0,nrow=n1,ncol=r*n3)
    ite=0
    diffbeta=1
    y=as.factor(y)
    svdall=c(0,0)
    beta=matrix(0,ncol=n3,nrow=n2)
    for (i in 1:r)
    {
      beta=beta+w0[i,]%*%t(v0[i,])
    }
    
    beta=beta/sqrt(sum(beta^2))
    
    while(diffbeta>0.001)
    {
      ite=ite+1
      betapre=beta
      
      decomp=svd(beta,r,r)
      vs=t(decomp$v[,1:r])
#      vs=matrix(NA,nrow=r,ncol=n3)
#      for (i in 1:r)
#      {
#        vs0=decomp$v[,i]*sqrt(decomp$d[i])
#        vs[i,]=vs0/sqrt(sum(vs0^2))
#      }
      
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
      d=median(as.vector(as.matrix(dist(dataw))[(k+1):nrow(dataw),1:k]))
      c=C*d^2/D^2
      fitw=kdwd(y~.,data=traw,type="bdwd",C=c,scaled=FALSE)
      w=matrix(NA,nrow=r,ncol=n2)
      for (i in 1:r)
      {
        w[i,]=as.vector(fitw@w[[1]])[((i-1)*n2+1):(i*n2)]
      }
      intw=fitw@b0[[1]]
      
      beta=matrix(0,ncol=n3,nrow=n2)
      for (i in 1:r)
      {
        beta=beta+w[i,]%*%t(vs[i,])
      }
      
      decomp=svd(beta,r,r)
      ws=t(decomp$u[,1:r])
  #    ws=matrix(NA,nrow=r,ncol=n2)
  #    for (i in 1:r)
  #    {
  #      ws0=decomp$u[,i]*sqrt(decomp$d[i])
  #      ws[i,]=ws0/sqrt(sum(ws0^2))
  #    }
      
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
      d=median(as.vector(as.matrix(dist(datav))[(k+1):nrow(datav),1:k]))
      c=C*d^2/D^2
      fitv=kdwd(y~.,data=trav,type="bdwd",C=c,scaled=FALSE)
      v=matrix(NA,nrow=r,ncol=n3)
      for (i in 1:r)
      {
        v[i,]=as.vector(fitv@w[[1]])[((i-1)*n3+1):(i*n3)]
      }
      intv=fitv@b0[[1]]
      
      beta=matrix(0,ncol=n3,nrow=n2)
      for (i in 1:r)
      {
        beta=beta+ws[i,]%*%t(v[i,])
      }
      
      diffbeta=dist(rbind(as.vector(beta),as.vector(betapre)))
    }
    
    mulrankr=rep(0,n2*n3)
    for (i in 1:r)
    {
      mulrankr=mulrankr+kronecker(v[i,],ws[i,])
    }
    result=list('beta'=beta,'int'=intv)
    return(result)
  }
}