library(funData)
library(MFPCA)
library(base)
multi_FPCA<-function(mfpca,sett){
  #### here mfpca is the result of function mFPCA
  ## mfpca is a list containing 5 elements: eigenfun,eigenvalue,sigma_sqerror,mean,fit,sample_index,eigenscore,org,sparse
  gamma<-mfpca$eigenvalue #9*9 matrix
  mean<-mfpca$mean
  sample_index<-mfpca$sample_index
  eigenfun<-mfpca$eigenfun
  #ufpcapredict<-lapply(1:p, function(i){mean[[i]]+eigenscore%*%eigenfun[[i]]@X})
  
  fitY<-mfpca$fit ### after permutation
  residual<-lapply(1:p,function(i){mfpca$org[[i]]-fitY[[i]]})  ### after permutation 3 lists
  
  Phi<-mfpca$eigenfun  ### 9*150 matrix
  
  len<-function(l,i)
  {
    func<-length(which(mfpca$time[[l]][sample_index[i],]==TRUE))
    return (func)
  }
  
  phi<-function(i)  ### 9*135 matrix for example
  {
    t_ind<-c()
    for (nb in 1:p)
    {t_ind<-c(t_ind,which(mfpca$time[[nb]][sample_index[i],]==TRUE))
    }
    func<-Phi[,t_ind]
    return (func)
  }
  
  H<-function(i)  ### 9*135 matrix for example
  {
    gamma%*%phi(i)
  }
  
   # observed_y_minus_mu<-function(i)
   # {
   #   t_ind<-c()
   #   val<-c()
   #   for (nb in 1:p)
   #   {
   #   observe<-mfpca$sparse[[nb]][i,which(mfpca$time[[nb]][sample_index[i],]==TRUE)]
   #   mean<-mfpca$mean[[nb]][which(mfpca$time[[nb]][sample_index[i],]==TRUE)]
   #   val<-c(val,observe-mean)
   #   }
   #   return (val)
   # }
  
    SigmaY<-function(i) 
  {
      diag_vect<-c()
      for (nb in 1:p)
      {diag_vect<-c(diag_vect,rep(mfpca$sigma_sqerror[nb],len(nb,i)))}
      Sigma_m<-diag(diag_vect)
    t(phi(i))%*%gamma%*%phi(i)+Sigma_m
  }
  
   # update_eigenscore<-matrix(NA,nrow=n,ncol=M)
   # for (NR in 1:n)
   # {
   #   #update_eigenscore[NR,]<-H(NR)%*%solve(SigmaY(NR))%*%observed_y_minus_mu(NR)
   #   update_eigenscore[NR,]<-solve(phi(i)%*%t(phi(i))+mfpca$sigma_sqerror[1]*solve(gamma))%*%phi(i)%*%observed_y_minus_mu(NR)
   # }
   # 
   score<-mfpca$score
   
   # update_fit<-list()
   # for (i in 1:p)
   # {
   #  update_fit[[i]]<-mean[[i]]+update_eigenscore%*%Phi[,((i-1)*length(sett)+1):(i*length(sett))]
   # }
   # fitY
  Omega<-function(i)
  {
    func<-gamma-H(i)%*%solve(SigmaY(i))%*%t(H(i))
    return (func)
  }
  
  cond_variance<-function(i)
  {
    cov_matrix<-t(Phi)%*%Omega(i)%*%Phi
    func<-diag(cov_matrix)
    return (func)
  } #### 50*3 vector
  
  cond_variance_fulltheta<-t(sapply(1:n,function(i){cond_variance(i)}))
  
  fit_sim<-c()
  for (i in 1:p)
  {
    fit_sim<-cbind(fit_sim,fitY[[i]])
  }
  
  fit_est<-matrix(NA,ncol=p*length(sett),nrow=n)
  
  fit_est[sample_index,]<-fit_sim
  
  var_est<- matrix(NA,nrow=nrow(cond_variance_fulltheta),ncol=ncol(cond_variance_fulltheta))
  
  var_est[sample_index,]<-cond_variance_fulltheta
  
  residual_sim<-c()
  for (i in 1:p)
  {residual_sim<-cbind(residual_sim,residual[[i]])
  }
  residual_est<-matrix(NA,nrow=n,ncol=p*length(sett))
  residual_est[sample_index,]<-residual_sim
  
sumtable<-list(fit=fit_est,variance=var_est,fit_difference=residual_est)
  ### now the results transform to the ascending order
  return (sumtable)
}

