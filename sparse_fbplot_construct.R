### fit is a n*(t*p) matrix of the fitted data such as cbind(n*length(time),n*length(time),...,n*length(time)).
#### sparse is a n*(t*p) matrix of the sparse data such as cbind(n*length(time),n*length(time),...,n*length(time)).
##p is the number of variables
## x is the xoordinate of time or index. If not given, 1,...length(x) is provided.
##med is a logical vector with length n such as c(False, False,...,True,False,...False), where TRUE means the median. 
## med can be obtained from the msplot_modified function.
##medlabel chooses whether to show the label of the median.
##remove is a logical vector with length n such as c(False, False,...,True,...,True,False,...False), where TRUE means the outlier.
## remove can be obtained from the msplot_modified function.
## outlabel chooses whether to show the label of the outlier.
## plot chooses whether to show the plot.
## fullout chooses whether to show all the outliers or just outliers outside the 50% central region.

sparse_fbplot<-function (fit, sparse,p=3,x = NULL,depth = NULL, med=NULL,medlabel=TRUE,remove=NULL,outlabel=TRUE,plot = TRUE, 
                         prob = 0.5, color = 6, outliercolor.fb = 2, barcol = 4, outliercolor.dir=3, fullout = FALSE, 
                         factor = 1.5,smoothtp="3RS3R") ## if remove is given, then we use MS plot to detect outliers already
{
  tp = dim(fit)[2]/p
  n = dim(fit)[1]
  if (length(x) == 0) {
    x = 1:tp
  }
  
  index<-order(depth, decreasing = TRUE)
  
  if (length(med)==0)
  {
    med = depth == max(depth)
  }
  sparse_ds<-function(i)
  {
    subfit<-t(fit[,((i-1)*tp+1):(i*tp)])
    subsparse<-t(sparse[,((i-1)*tp+1):(i*tp)])
    medavg = matrix(subfit[,med==1 ], ncol = sum(med), nrow = tp)
    y = apply(medavg, 1, mean)
    
    outpoint.dir = which(remove == 1)
    if (length(outpoint.dir)>0)
    {index<-setdiff(index,outpoint.dir)
    } 
    
    m = ceiling(length(index) * prob)
    center = subfit[, index[1:m]]
    sp_center<-subsparse[, index[1:m]]
    
    inf = apply(center, 1, min) ### shown in figure
    sup = apply(center, 1, max) ### shown in figure
    if (prob == 0.5) {
      dist = factor * (sup - inf)
      upper = sup + dist
      lower = inf - dist
      outly = (subfit <= lower) + (subfit >= upper)
      outcol = colSums(outly)
      removefb = (outcol > 0)
      colum = 1:n
      outpoint.fb = colum[removefb == 1]
      outfb<-setdiff(outpoint.fb,outpoint.dir)
    }
    
    sparse_density_center<-t(sapply(1:nrow(sp_center),function(nnr)
    {
      resf<-data.frame(obs=length(which(!is.na(sp_center[nnr,])))/m,spa=length(which(is.na(sp_center[nnr,])))/m)
      return (resf)
    }))
    return(list(outpointfb = outfb,outpointdir=outpoint.dir,sparse_density_ct=sparse_density_center))
  }
  #####################################################################################
  sparse_density_ct<-lapply(1:p,function(i){sparse_ds(i)$sparse_density_ct})
  outpointfb<-lapply(1:p,function(i){sparse_ds(i)$outpointfb})
  outpointdir<-sparse_ds(1)$outpointdir
  
  ########################################################################################
  subfb<-function(i)
  { 
    subfit<-t(fit[,((i-1)*tp+1):(i*tp)])
    subsparse<-t(sparse[,((i-1)*tp+1):(i*tp)])
    
    medavg = matrix(subfit[,med ], ncol = sum(med), nrow = tp)
    y = apply(medavg, 1, mean)
    if (length(outpointdir)>0)
    {index<-setdiff(index,outpointdir)
    } 
    m = ceiling(length(index) * prob)
    center = subfit[, index[1:m]]
    sp_center<-subsparse[, index[1:m]]
    sp_ct<- sparse_density_ct[((i-1)*length(time)+1):(i*length(time))]
    outcenter = subfit[, -index[1:m]]
    inf = apply(center, 1, min) ### shown in figure
    sup = apply(center, 1, max) ### shown in figure
    
    if (length(outpointfb[[i]])==0)
    {good = subfit[, index]
    } else 
    {
      good = subfit[, setdiff(index,outpointfb[[i]])]
    }
    maxcurve = apply(good, 1, max) ### shown in figure
    mincurve = apply(good, 1, min) ## shown in figure
    
    yrange<-c(min(fit), max(fit))
    
    plotout<-function(OP,outlabel)
    { if (length(OP)>0)
    {  out<-subfit[,OP] ##
    sparse_out<-subsparse[, OP] ##
    co=ifelse(sum(OP)==sum(outpointdir),outliercolor.dir,outliercolor.fb)
    if (length(OP)==1)
    {
      sparse_outindex<-which(is.na(sparse_out)) ## index where there are observationums in the midcurves
      nsparse_outindex<-which(!is.na(sparse_out))
      if (length(sparse_outindex)==0)
      {
        lines(x, out[nsparse_outindex],pch=16,lty = 2, lwd = 1,col=co)
      }  else {
        for (ck in 1:length(nsparse_outindex))
        {
          if ((nsparse_outindex[ck]+1)<=length(x))
          {
            segments(x[nsparse_outindex[ck]], out[nsparse_outindex[ck]],x1=x[nsparse_outindex[ck]+1], y1=out[nsparse_outindex[ck]+1],pch=16,lty = 2, lwd = 1,col=co)
          } else {
            segments(x[nsparse_outindex[ck]], out[nsparse_outindex[ck]],x1=x[nsparse_outindex[ck]], y1=out[nsparse_outindex[ck]],pch=16,lty = 2, lwd = 1,col=co)
          }
        }
        for (cl in 1:length(sparse_outindex))
        {
          if ((sparse_outindex[cl]+1)<=length(x))
          {
            segments(x[sparse_outindex[cl]],out[(sparse_outindex[cl])],x1=x[sparse_outindex[cl]+1],y1=out[sparse_outindex[cl]+1],lwd=1,lty=2,col=grey(0.5))
          } else
          {segments(x[sparse_outindex[cl]],out[(sparse_outindex[cl])],x1=x[sparse_outindex[cl]],y1=out[sparse_outindex[cl]],lwd=1,lty=2,col=grey(0.5))}   
        } 
      }
      if (outlabel==TRUE)
      {
        labindex<-ifelse(i==1,max(nsparse_outindex[length(nsparse_outindex)],sparse_outindex[length(sparse_outindex)],na.rm=TRUE),
                         min(nsparse_outindex[1],sparse_outindex[1],na.rm=TRUE))
        text(x[labindex],out[labindex],OP,pos=4,cex=0.5,col=co)
      }
    } else {
      for (nums in 1:length(OP))
      {
        sparse_outindex<-which(is.na(sparse_out[,nums])) ## index where there are observationums in the midcurves
        nsparse_outindex<-which(!is.na(sparse_out[,nums]))
        if (length(sparse_outindex)==0)
        {
          lines(x, out[nsparse_outindex,nums],pch=16,lty = 2, lwd = 1,col=co)
        }  else {
          for (ck in 1:(length(nsparse_outindex)))
          {
            if ((nsparse_outindex[ck]+1)<=nrow(sparse_out))
            {
              segments(x[nsparse_outindex[ck]], out[nsparse_outindex[ck],nums],x1=x[nsparse_outindex[ck]+1], y1=out[nsparse_outindex[ck]+1,nums],pch=16,lty = 2, lwd = 1,col=co)
            } else {
              segments(x[nsparse_outindex[ck]], out[nsparse_outindex[ck],nums],x1=x[nsparse_outindex[ck]], y1=out[nsparse_outindex[ck],nums],pch=16,lty = 2, lwd = 1,col=co)
            }
          }
          for (cl in 1:length(sparse_outindex))
          {
            if ((sparse_outindex[cl]+1)<=nrow(out))
            {
              segments(x[sparse_outindex[cl]],out[(sparse_outindex[cl]),nums],x1=x[sparse_outindex[cl]+1],y1=out[sparse_outindex[cl]+1, nums],lwd=1,lty=2,col=grey(0.5))
            } else
            {segments(x[sparse_outindex[cl]],out[(sparse_outindex[cl]),nums],x1=x[sparse_outindex[cl]],y1=out[sparse_outindex[cl], nums],lwd=1,lty=2,col=grey(0.5))}   
          } 
        }
        if (outlabel==TRUE)
        {
          labindex<-ifelse(i==1,max(nsparse_outindex[length(nsparse_outindex)],sparse_outindex[length(sparse_outindex)],na.rm=TRUE),
                           min(nsparse_outindex[1],sparse_outindex[1],na.rm=TRUE))
          if (OP[nums]!=36)
          {
            xxlab=x[labindex]          
          } else
          {
            xxlab=ifelse(i==1,x[labindex]+0.2,x[labindex]-0.57)
          }
          text(xxlab,out[labindex,nums],OP[nums],pos=ifelse(i==1,4,2),cex=0.5,col=co)
        } 
      }
    }
    }
    }
    par(mai=c(0.6,0.6,0.4,0.1),mar=c(4.5, 3.5, 3, 1))
    
    if (i==p)
    {
      par(mai=c(0.6,0.6,0.4,0.1),mar=c(4, 3.5, 3, 0.5))
    }
    if (plot) {
      if (p==1)
      {
      xlab="Time"
      ylab="Values"}
      else {
        varname<-c("variable ",i)
        xlab="Time"
        ylab="Value"
      }
      #par(mai=c(0.55,0.6,0.4,0.1))
      plot(x, y, lty = 1, lwd = 2, col = "white", type = "l", xlim=c(range(x)[1],range(x)[2]+1/40*(range(x)[2]-range(x)[1])), 
           ylim=yrange,ylab=ylab,
           xlab=xlab,main=paste(ifelse(sum(is.na(sparse))>0,"Sparse ",""),ifelse(length(remove)==0,"","Two-Stage "),"Functional Boxplot",sep=""))
      if (fullout==FALSE) {
        plotout(outpointdir,outlabel)
        plotout(outpointfb[[i]],outlabel)
      }
      
      barval = (x[1] + x[tp])/2
      bar = which(sort(c(x, barval)) == barval)[1]
      lines(c(x[bar], x[bar]), c(maxcurve[bar], sup[bar]), 
            col = barcol, lwd = 2)
      lines(c(x[bar], x[bar]), c(mincurve[bar], inf[bar]), 
            col = barcol, lwd = 2)
    }
    
    x_adj<-x
    x_adj[1]<-x[1]-1/80*(range(x)[2]-range(x)[1])
    x_adj[length(x)]<-x[length(x)]+1/80*(range(x)[2]-range(x)[1])
    xx = c(x_adj, x_adj[order(x, decreasing = TRUE)])
    supinv = sup[order(x, decreasing = TRUE)]
    yy<-c(inf,supinv)
    
    if (plot) {
      med_center<-inf+0.5*(sup-inf)
      sep<-inf+unlist(sparse_density_ct[[i]][,1])*(sup-inf)
      if (abs(sum(sep)-sum(sup))>0.2)
      {sep<-smooth(sep,kind=smoothtp)
      }
      sep_inv<-sep[order(x, decreasing = TRUE)]
      yy_lower = c(inf, sep_inv)
      yy_upper<-c(sep,supinv)
      polygon(xx, yy_lower, col = "magenta", border = "magenta", 
              lwd = 2)
      polygon(xx, yy_upper, col = "grey", border = "grey", 
              lwd = 2)
      polygon(xx, yy, border = barcol, lwd = 2)
      lines(x,maxcurve,type="l",lwd=2,col=barcol)
      lines(x,mincurve,type="l",lwd=2,col=barcol)
      
      
      sparse_midindex<-which(is.na(subsparse[,index[1]])) ## index where there are observations in the midcurves
      nsparse_midindex<-which(!is.na(subsparse[,index[1]]))
      #points(x[nsparse_midindex], subfit[nsparse_midindex, index[1]], type="p",pch=16,lty = 1)
      lines(x, subfit[, index[1]],pch=16,lty = 1, lwd = 2,col="white")
      if (length(sparse_midindex)==0)
      {
        lines(x, subfit[, index[1]],pch=16,lty = 1, lwd = 2,col=gray(0))
      } else 
      {
        for (cl in 1:length(sparse_midindex))
        {
          if ((sparse_midindex[cl]+1)<=nrow(subfit))
          {
            segments(x[sparse_midindex[cl]],subfit[(sparse_midindex[cl]),index[1]],x1=x[sparse_midindex[cl]+1],y1=subfit[sparse_midindex[cl]+1, index[1]],lwd=1.5,lty=1,col="gray")
          }
          else
          {segments(x[sparse_midindex[cl]],subfit[(sparse_midindex[cl]),index[1]],x1=x[sparse_midindex[cl]],y1=subfit[sparse_midindex[cl], index[1]],lwd=1.5,lty=1,col="gray")
          }   
        }
        
        for (ck in 1:(length(nsparse_midindex)))
        {
          if ((nsparse_midindex[ck]+1)<=nrow(subfit))
          {
            segments(x[nsparse_midindex[ck]], subfit[nsparse_midindex[ck], index[1]],x1=x[nsparse_midindex[ck]+1], y1=subfit[nsparse_midindex[ck]+1,index[1]],pch=16,lty = 1, lwd = 2,col=gray(0))
          } else {
            segments(x[nsparse_midindex[ck]], subfit[nsparse_midindex[ck], index[1]],x1=x[nsparse_midindex[ck]], y1=subfit[nsparse_midindex[ck],index[1]],pch=16,lty = 1, lwd = 2,col=gray(0))
          }
        }
        lines(x,med_center,lty=3,lwd=1.5,col="cyan")
      }
      if (medlabel==TRUE) 
      {text(max(x)+1/25*(range(x)[2]-range(x)[1]),subfit[length(x), index[1]],labels=which(med),cex=0.7)}
      if (fullout==TRUE) {
        
        plotout(outpointdir,outlabel)
        plotout(outpointfb[[i]],outlabel)
      }
    }  
  }
  #par(mfrow=c(1,p),mai=c(0.8,0.8,0.4,0.1))
  for (kk in 1:p)
  {
    func<-subfb(kk)
  }
  return(list(depth = depth, outpointfb = outpointfb, outpointdir=outpointdir, sparse_density_ct=sparse_density_ct,medcurve = which(med)))
  
}
