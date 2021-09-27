library("plotfunctions")
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
## when colorrange=NULL, intensity_sparse_fbplot function will obtain the range in all variables and normalize the intensity by dividing by the respective maximum in each variable.
## If input the obtained colorrange and run the intensity_sparse_fbplot again, the normalization is obtained by dividing by the maximum of all variables.

intensity_sparse_fbplot<-function (fit, sparse,p=2,x = NULL,depth = NULL, med=NULL,medlabel=TRUE,remove=NULL,plot = TRUE, 
                                   prob = 0.5, color = 6, outliercolor.fb = 2, barcol = 4, outliercolor.dir=3,fullout = FALSE, 
                                   factor = 1.5,colorrange=NULL,colorsplit=40,contour=TRUE,legend=TRUE) ## if remove is given, then we use MS plot to detect outliers already
{
  index<-order(depth, decreasing = TRUE)
  tp<-length(x)
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
      colum = 1:nrow(fit)
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
  select_col<-c("magenta","tomato","gold","yellow","white")
  grRd <- colorRampPalette(select_col,space="rgb")  ####################################################################
  if (legend==TRUE&&p>1&&sum(is.na(sparse))>0)
  {
    m<-matrix(1:(p+1),nrow=1,ncol=p+1,byrow=TRUE)
    
    if (p==2)
    {layout(mat=m,widths = c(rep((1-0.07)/p,p),0.07))
    } else if (p==3)
    {
      layout(mat=m,widths = c(rep((1-0.035)/p,p),0.035))
    }
    
  }
  
  yrange<-c(min(fit,na.rm=TRUE), max(fit,na.rm=TRUE))
  
  subfb<-function(i,plot)
  { 
    subfit<-t(fit[,((i-1)*tp+1):(i*tp)])
    subsparse<-t(sparse[,((i-1)*tp+1):(i*tp)])
    
    # subfit<-10*t(fit[,((i-1)*tp+1):(i*tp)])
    # subsparse<-10*t(sparse[,((i-1)*tp+1):(i*tp)])
    
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
    
    plotout<-function(OP)
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
            segments(x[sparse_outindex[cl]],out[(sparse_outindex[cl])],x1=x[sparse_outindex[cl]+1],y1=out[sparse_outindex[cl]+1],lwd=1,lty=2,col="gray")
          } else
          {segments(x[sparse_outindex[cl]],out[(sparse_outindex[cl])],x1=x[sparse_outindex[cl]],y1=out[sparse_outindex[cl]],lwd=1,lty=2,col="gray")}   
        } 
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
              segments(x[sparse_outindex[cl]],out[(sparse_outindex[cl]),nums],x1=x[sparse_outindex[cl]+1],y1=out[sparse_outindex[cl]+1, nums],lwd=1,lty=2,col="gray")
            } else
            {segments(x[sparse_outindex[cl]],out[(sparse_outindex[cl]),nums],x1=x[sparse_outindex[cl]],y1=out[sparse_outindex[cl], nums],lwd=1,lty=2,col="gray")}   
          } 
        }
      }
    }
    }
    }
    par(mai=c(0.6,0.6,0.4,0.1),mar=c(4.5, 3.5, 3, 1))
    
    if (i==p)
    {
      par(mai=c(0.6,0.6,0.4,0.1),mar=c(4, 4, 3, 2))
    }
    
    if (i==p&&legend==TRUE)
    {
      par(mai=c(0.52,0.5,0.4,0.02),mgp=c(2.2,1,0))
    }
    
    rangeupper<-ifelse(medlabel==TRUE,range(x)[2]+1/40*(range(x)[2]-range(x)[1]),range(x)[2])
    if (plot) {
      if (p==1)
      {
      xlab="Time"
      ylab="Values"}
      else {
        varname<-c("variable ",i)
        xlab="Time"
        ylab="Values"
      }
      plot(x, y, lty = 1, lwd = 2, col = 1, type = "n", xlim=c(range(x)[1],rangeupper), 
           ylim=yrange,ylab=ylab,
           #yaxt="n",
           xlab=xlab,main=paste("Intensity Sparse ",ifelse(length(remove)==0,"","Two-Stage "),"Functional Boxplot",sep=""))
      # if (i==1 |i ==2)
      #  {ylabel<-c(-15,-10,-5,0,5,10)
      # } else if (i==3)
      # {
      #   ylabel<-c(-10,-5,0,5,10)
      # }
      # ticks<-10*ylabel
      # axis(2,at=ticks,label=ylabel)
      #axis(2, seq(yrange[1],yrange[2],length.out=4),las=2)
      if (fullout==FALSE) {
        plotout(outpointdir)
        plotout(outpointfb[[i]])
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
    
    solo_y<-unique(c(center))
    x_y_match<-function(jk)
    {
      xtime<-rep(x[jk],length(unique(center[jk,])))
      yvalue<-unique(center[jk,])
      zpoint<-c()
      zna_point<-c()
      for (tpp in 1:length(xtime))
      {
        mat_index<-which(center[jk,]==yvalue[tpp])
        obs_index<-which(!is.na(sp_center[jk,]))
        sp_ind<-setdiff(mat_index,obs_index)
        obs_ind<-intersect(mat_index,obs_index)
        zpoint[tpp]<-ifelse(length(obs_ind)>0,length(obs_ind),0)
        zna_point[tpp]<-ifelse(length(sp_ind)>0,length(sp_ind),0)
      }
      wl<-data.frame(x=xtime,y=yvalue,zobs=zpoint,zna=zna_point)
      return (wl)
    }
    point<-lapply(1:length(x),x_y_match)
    xpoint<-unlist(lapply(point, function(ls){ls["x"]}))
    ypoint<-unlist(lapply(point, function(ls){ls["y"]}))
    nobs_point<-unlist(lapply(point, function(ls){ls["zna"]}))
    pp_na<-ppp(rep(xpoint,nobs_point),rep(ypoint,nobs_point),poly=list(x=xx,y=yy))
    #plot(pp_na, pch=20, cols="grey70", main="sparse value point process")  # Plot points
    #Q <- quadratcount(pp_na, nx= 41, ny=30)
    #plot(Q, add=TRUE)  # Add quadrat grid
    Q.d <- density(pp_na,ajdust=1,dimyx=c(200,200))
    #sq<-summary(Q.d)
    norm_Q.d<-Q.d
    #norm_Q.d$v<-Q.d$v/max(Q.d$v,na.rm=TRUE)
    #Q.d.normalize<-Q.d
    #Q.d.normalize$v<-Q.d$v/sq$integral
    if (length(colorrange)==0)
    {
      #rg<-range(Q.d.normalize$v,na.rm=TRUE)
      rg<-range(Q.d[['v']],na.rm=TRUE)
      norm_Q.d$v<-Q.d$v/max(Q.d$v,na.rm=TRUE)
    } else {
      rg<-range(unlist(colorrange))
      norm_Q.d<-Q.d
      norm_Q.d$v<-Q.d$v/max(rg,na.rm=TRUE)
    }
    #pp_na_mark<-ppp(rep(xpoint,nobs_point),rep(ypoint,nobs_point),poly=list(x=xx,y=yy),marks=rg)
    CO<- colourmap(grRd(colorsplit), range = c(0,1))
    
    if (plot) {
      
      if (sum(nobs_point)>0)
      {
        plot(norm_Q.d, las=1,add=TRUE,col=CO)  # Plot density raster
        if (contour==TRUE)
        {contour(norm_Q.d, add=TRUE,drawlabels = FALSE,levels = c(0,0.2,0.4,0.6,0.8,1))
        }
        polygon(xx, yy, border = barcol, lwd = 2)
      } else
      {
        polygon(xx, yy, col = "magenta", border = barcol, lwd = 2)
      }
      #}
      #plot(pp_obs, pch=20, cols="grey70", main=NULL)  # Plot points
      lines(x,maxcurve,type="l",lwd=2,col=barcol)
      lines(x,mincurve,type="l",lwd=2,col=barcol)
      
      
      sparse_midindex<-which(is.na(subsparse[,index[1]])) ## index where there are observations in the midcurves
      nsparse_midindex<-which(!is.na(subsparse[,index[1]]))
      #points(x[nsparse_midindex], subfit[nsparse_midindex, index[1]], type="p",pch=16,lty = 1)
      if (length(sparse_midindex)==0)
      {
        lines(x, subfit[, index[1]],pch=16,lty = 1, lwd = 2,col=gray(0))
      } else 
      {
        for (cl in 1:length(sparse_midindex))
        {
          if ((sparse_midindex[cl]+1)<=nrow(subfit))
          {
            segments(x[sparse_midindex[cl]],subfit[(sparse_midindex[cl]),index[1]],x1=x[sparse_midindex[cl]+1],y1=subfit[sparse_midindex[cl]+1, index[1]],lwd=2,lty=1,col="gray")
          }
          else
          {segments(x[sparse_midindex[cl]],subfit[(sparse_midindex[cl]),index[1]],x1=x[sparse_midindex[cl]],y1=subfit[sparse_midindex[cl], index[1]],lwd=2,lty=1,col="gray")
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
        #lines(x,med_center,lty=2,lwd=1,col="cyan")
      }
      if (medlabel==TRUE) 
      {text(max(x)+1/25*(range(x)[2]-range(x)[1]),subfit[length(x), index[1]],labels=which(med),cex=0.9)}
      if (fullout==TRUE) {
        plotout(outpointdir)
        plotout(outpointfb[[i]])
      }
    }
    return (rg)
  }
  
  rg<-list()
  for (kk in 1:p)
  {
    func<-subfb(kk,plot)
    rg[[kk]]<-c(func[1],func[2])
    
  }
  if (legend==TRUE&&sum(is.na(sparse))>0)
  {
    par(mai=c(0.52,0.2,0.95,0.1))
    plot(1, type = "n", axes=FALSE, xlab="", ylab="",main="%")
    # plot_colors<-c("magenta","red","orange","yellow","green","cyan","blue","slateblue","black","gray")
    
    gradientLegend(valRange=c(0,100),color=select_col,length=1,depth=0.45, side=4,dec=0,inside=TRUE,n.seg=4,pos=c(0.3,0,0.7,1.07),coords=FALSE,cex=1.03)
  }
  return(list(depth = depth, outpointfb = outpointfb, outpointdir=outpointdir, colorrange=rg,sparse_density_ct=sparse_density_ct,medcurve = which(med)))
  
}

