Bootstrap<-function(bsTims,sparse_set1,sett)
{
  if (class(sparse_set1)=="matrix")
  {
    newsparse_set1<-list()
    newsparse_set1[[1]]<-sparse_set1
    sparse_set1<-newsparse_set1
  }
  p<-length(sparse_set1)
  n<-nrow(sparse_set1[[1]])/2
  train_bootstrap<-list()
  mean_fit<-matrix(NA,nrow=n,ncol=p*length(sett))
  mean_variance<-matrix(NA,nrow=n,ncol=p*length(sett))
  var_fit_difference<-matrix(NA,nrow=n,ncol=p*length(sett))
  
  
  for (bt in 1:bsTims)
  {
    cat(bt,"\n")
   mfpca_b<-mFPCA(sparse_set1,sett,Bt=TRUE)
    train_bootstrap[[bt]]<-multi_FPCA(mfpca_b,sett) 
  }
  
  
  for (rw in 1:n)
  {
    for (cl in 1: (p*length(sett)) )
    {
      use_fitdif<-c()
      use_fit<-c()
      use_var<-c()
      for (bt in 1:bsTims)
      {
        use_fitdif<-c(use_fitdif,train_bootstrap[[bt]]$fit_difference[rw,cl])
        use_fit<-c(use_fit,train_bootstrap[[bt]]$fit[rw,cl])
        use_var<-c(use_var,train_bootstrap[[bt]]$variance[rw,cl])
      }
      var_fit_difference[rw,cl]<-var(use_fitdif,na.rm=TRUE)
      mean_fit[rw,cl]<-median(use_fit,na.rm=TRUE)
      mean_variance[rw,cl]<-median(use_var,na.rm=TRUE)
    }
  }
  
  ######################## results of overall bootstrap theta
  fit_bs<-mean_fit
  variance_bs<-mean_variance+var_fit_difference
  
  func<-list(fit=fit_bs,variance=variance_bs)
  return (func)
}


########################################## trimmed_mean is made to remove the extremely large variance part in bootstrap
trimmed_mean<-function(x,tr,na.rm)
{
  newx<-sort(x,decreasing = FALSE,na.last=NA)
  if (tr[2]=="symmetric")
    {func<-mean(newx,trim=as.numeric(tr[1]))
} else if (tr[2]=="upper")
 {
  end<-quantile(newx,probs=1-as.numeric(tr[1]),na.rm)
  cand_out<-boxplot(newx,plot=FALSE)$out
  out<-max(end,min(cand_out,na.rm=TRUE))
  func<-mean(newx[which(newx<out)])
  } else if (tr[2]=="lower")
  {
    st<-quantile(newx,probs=as.numeric(tr[1]),na.rm)
    cand_out<-boxplot(newx,plot=FALSE)$out
    out<-max(st,min(cand_out))
    func<-mean(newx[which(newx>out)])
    
  }
  return (func)
}

