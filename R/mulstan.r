mulstan=function(multi)
{
  multis=multi
  n1=dim(multi)[1]
  n2=dim(multi)[2]
  n3=dim(multi)[3]
  for (i in 1:n3)
  {
    for (j in 1:n2)
    {
      multis[,j,i]=(multi[,j,i]-mean(multi[,j,i]))/sd(multi[,j,i])
    }
  }
  return(multis)
}