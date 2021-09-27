mFPCA<-function(data,sett,Bt)
{
  if (class(data)=="matrix")
  {
  datanew<-list()
  datanew[[1]]<-data
  data<-datanew
  } 
    p=length(data)
  nr<-nrow(data[[1]])/2
  nc<-ncol(data[[1]])
  
  
  org_data<-lapply(1:p,function(i){data[[i]][1:nr,]})
  
  sparse_data<-lapply(1:p,function(i){data[[i]][-(1:nr),]})
  
  new_org_data<-org_data
  newsparse_data<-sparse_data
  
  time_index<-lapply(1:p, function(i){matrix(TRUE,nrow=nr,ncol=nc)})
  
  for (j in 1:p)
  {for (i in 1:nr)
  {
    time_index[[j]][i,which(is.na(sparse_data[[j]][i,]))]<-FALSE
  }
  }
  
  obj<-lapply(1:p, function(i){
    Ly<-sparse_data[[i]]
    
    Lt<-list(sett)
    res<-funData(argvals = Lt,X=Ly)
    return(res) }  )
  expres<-list(type="uFPCA")
  mFData<-multiFunData(obj)
  newmFData<-mFData
  if (Bt==TRUE)
  {
  x<-1:nr
  sample_nb<-sample(x,nr,replace=TRUE)
  sample_index<-sort(sample_nb)
  
     for (i in 1:p)
    {newmFData[[i]]@X<-mFData[[i]]@X[sample_index,]
  newsparse_data[[i]]<-sparse_data[[i]][sample_index,]
  new_org_data[[i]]<-org_data[[i]][sample_index,]
    }
  
  } else {
    sample_index<-1:nr
  }

  estimated_observed_p<-sapply(1:p,function(i){apply(newsparse_data[[i]],1,function(j){1-sum(is.na(j))/length(sett)})})
  
  min_observed_points<-length(sett)*min(estimated_observed_p)
  
  M<-ifelse(min_observed_points<=3,4,9)
  
  fit_face_sparsecase<-MFPCA(newmFData,M,uniExpansions =lapply(1:p,function(i){expres}),fit=TRUE)
  
  eigenfun<-fit_face_sparsecase$functions
  
  gamma<-diag(fit_face_sparsecase$values)
  meanfun<-fit_face_sparsecase$meanFunction
  
  mean<-lapply(1:p, function(i){as.vector(meanfun[[i]]@X)})
  
  fitY<-lapply(1:p,function(i){fit_face_sparsecase$fit[[i]]@X})
  residual<-lapply(1:p,function(i){org_data[[i]]-fitY[[i]]})
  
  score<-fit_face_sparsecase$scores
  meas_error<-unlist(lapply(residual, function(i){var(as.vector(i))}))
  
  Phi<-c()
  for (nb in 1:p)
  {
    Phi<-cbind(Phi,eigenfun[[nb]]@X) 
  }
  result<-list(eigenfun=Phi,eigenvalue=gamma,sigma_sqerror=meas_error,mean=mean,fit=fitY,score=score,time=time_index,sample_index=sample_index,org=new_org_data,sparse=newsparse_data)
  
  return (result)
}

