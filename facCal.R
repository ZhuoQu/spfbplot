facCal=function(n,dim,iter=100)
{
temp=cm(n,dim,iter,100)
fac=temp[1]*(temp[2]-dim+1)/(dim*temp[2])
ptm <- proc.time()
M=temp[2]
cutoff=qf(0.993,dim,M-dim+1)
return(c(fac,cutoff))
}
