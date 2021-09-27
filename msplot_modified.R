## data is an array where dim(data)=c(n,length(time),p)
##n is the number of samples
##p is the number of variables
msplot_modified<-function(data,depth.dir="RP",plot=FALSE,col.out="green",col.med="purple",
                col.normal="black",col="white",dirout=TRUE,plot.type="scatter",Met="multi")
{
  ###pairwise plots of variation of outlyingness (VO) against mean outlyingness (MO)###
  
  
  temp=dim(data)
  r1=NULL
  ##Univariate cases##
  if (length(temp)==2)
  { 
    n=temp[1]
    result=DirOut(data,depth.dir=depth.dir)
    if (dirout)
    {
      D=result$D
      factor=facCal_num(temp[1],2) 
      fac1=factor$fac1
      cutoff1=factor$fac2 #cut off value for testing/outlier detection#
      num=sum(fac1*D>cutoff1) #number of outliers#
      cutoff=cutoff1/fac1
      
      mo <-result$out_avr
      vo <-result$out_var
      out.dir=which(result$D>cutoff)
      medcurve=which.min(result$D)
      
      pch=rep(20,n)
      pch[out.dir]=19
      pch[medcurve]=17
      
      outlabel=rep('',nrow(data))
      outlabel[out.dir]=out.dir
      outlabel[medcurve]=medcurve
      
      outcol=rep("white",nrow(data))
      outcol[out.dir]=as.numeric(1:length(out.dir))
      outcol[medcurve]="black"
      
      if (plot)
      {
        M=cbind(mo,vo)
        ans=cov.rob(M,method="mcd")
        L=solve(chol.default(solve(ans$cov)))
        theta=1:dim(data)[2]/(dim(data)[2]-1)
        circle=abind(sin(theta*2*pi),cos(theta*2*pi),along=2)
        x=ans$center[1]+(cutoff^(1/2)*L%*%t(circle))[1,]
        y=ans$center[2]+(cutoff^(1/2)*L%*%t(circle))[2,]
        elip.data=data.frame(x=x,y=y,a=rep(1,length(x)))
        
        col.point=rep(col.normal,n)
        col.point[out.dir]=col.out
       
        ms.data=data.frame(x=mo,y=vo,out=col.point,pch=pch)
        
        # p<-ggplot(data=ms.data,aes(x=x,y=y))+geom_point(col=col.point)+
        #   geom_path(data=elip.data,aes(x=x,y=y),show.legend=FALSE,colour="lightblue")+
        #   xlab("MO")+ylab("VO")+labs(title="MS-Plot")+theme(plot.title = element_text(hjust = 0.5))
        # 
        
        plot(mo,vo,type="n",ylab="VO",xlab="MO",ylim=range(vo),xlim=range(mo),main=paste("MS-Plot on","Estimated Data", ifelse(Met=="uni","(uni-method)","(multi-method)"),sep=" "))
        for (i in setdiff(1:length(outlabel),c(out.dir,medcurve)))
        {
          points(mo[i],vo[i],pch=pch[i],col=col.normal,cex=0.7)}
         
        }
        for (i in out.dir)
        {
          points(mo[i],vo[i],pch=pch[i],col=col.out,cex=0.7)
        }
        points(mo[medcurve],vo[medcurve],col=col.med,pch=pch[medcurve])
        text(mo+0.08,vo+0.05*10^{-5},labels = outlabel,col=outcol,cex=0.5)
        legend("top",col=c("black","red","purple"),legend=c("normal","outlier","median"),pch=c(20,19,pch[medcurve]),cex=0.7)
      }
      
      return(r1=list(mo=mo,vo=vo,out.dir=out.dir,medcurve=medcurve
                     #,p=p
      ))
        # if (length(out.dir)==0)
        # {
        #   p<-ggplot(data=ms.data,aes(x=x,y=y,colour=factor(out,labels="F")))+geom_point(size=2)+
        #     guides(colour=guide_legend(title="Outlying"))+scale_color_manual(values=3)+scale_shape_manual(values=20)+
        #     geom_path(data=elip.data,aes(x=x,y=y),show.legend=FALSE,colour="lightblue")+
        #     xlab("MO")+ylab("VO")+labs(title="MS-Plot")+theme(plot.title = element_text(hjust = 0.5))
        # }
        
        
        #plot(mo,vo,type="p",col=col.point,pch=pch,xlim=c(min(x,M[,1]),max(x,M[,1])),
        #     ylim=c(min(y,M[,2]),max(y,M[,2])),...)
        #lines(x,y,type="l",lty=1,col="lightblue")
      
      
        
        return(r1=list(mo=mo,vo=vo,out.dir=out.dir,medcurve=medcurve))
    
    if (!dirout)
    {
      mo <-result$out_avr
      xlim=c(min(mo)-0.1*(max(mo)-min(mo)),max(mo)+0.1*(max(mo)-min(mo)))
      vo <-result$out_var
      fo <-mo^2+vo
      ms.data=data.frame(x=mo,y=vo,col=col)
      p<-ggplot(data=ms.data,aes(x,y))+geom_point(colour=col,show.legend=FALSE)+
        xlab("MO")+ylab("VO")+labs(title="MS-Plot")+theme(plot.title = element_text(hjust = 0.5))
      if (plot)
      {return(list(mo=mo,vo=vo,p=p))}
      else 
      {return(list(mo=mo,vo=vo))}
    }
  } 
   else if (length(temp)==3)
  {
    
    factor=facCal_num(temp[1],dim=temp[3]+1)
    fac2=factor$fac1
    cutoff2=factor$fac2   #cut off value for testing/outlier detection#
    
    d=temp[3]
    n=temp[1]
    

    
    if (d>=2)
    {
      result=DirOut(data,depth.dir=depth.dir)
      
      if (dirout)
      {
        cutoff=cutoff2/fac2
        mo <-result$out_avr
        vo <-result$out_var
        out.dir=which(result$D>cutoff)
        medcurve=which.min(result$D)
        ###
        #hist(result$D[which(result$D<=2000)])
        if (plot)
        {
          M=cbind(mo,vo)
          # ans=cov.rob(M,method="mcd")
          # L=solve(chol.default(solve(ans$cov)))
          # theta=0:200/200
          # v=(0:200-100)/100
          # #circle=abind((1-v^2)^(1/2)*sin(theta*2*pi),(1-v^2)^(1/2)*cos(theta*2*pi),v,along=2)
          # x=as.vector((1-v^2)^(1/2)%*%t(sin(theta*2*pi)))
          # y=as.vector((1-v^2)^(1/2)%*%t(cos(theta*2*pi)))
          # z=as.vector(v%*%t(rep(1,201)))
          # circle=abind(x,y,z,along=2)
          # 
          # x1=ans$center[1]+(cutoff^(1/2)*L%*%t(circle))[1,]
          # y1=ans$center[2]+(cutoff^(1/2)*L%*%t(circle))[2,]
          # z1=ans$center[3]+(cutoff^(1/2)*L%*%t(circle))[3,]
          # 
          col.point=rep(col.normal,n)
          col.point[medcurve]=col.med
          col.point[out.dir]=col.out
          pch=rep(20,n)
          pch[out.dir]=19
          pch[medcurve]=17
          # xlim=c(min(x1,M[,1]),max(x1,M[,1]))
          # ylim=c(min(y1,M[,2]),max(y1,M[,2]))
          # zlim=c(min(z1,M[,3]),max(z1,M[,3]))
          outcol=rep("white",nrow(data))
          outcol[out.dir]=1:length(out.dir)
          outcol[medcurve]="black"
          
          if (plot.type=="parallel")
          {
            data.ms=data.frame(MO=mo,VO=vo,col=col.point) 
            p<-ggparcoord(data.ms,columns=columns,groupColumn =(d+2),showPoints = TRUE)
            return(r1=list(mo=mo,vo=vo,out.dir=out.dir,medcurve=medcurve,p=p))
            
          }
          
          if (plot.type=="scatter")
          {
            MO<-(apply(mo^2,1,sum))^(1/2)
            outlabel=rep('',nrow(data[,,1]))
            outlabel[out.dir]=out.dir
            outlabel[medcurve]=medcurve
            ms.data=data.frame(x=MO,y=vo,col=col.point,outlabel=outlabel)
            #p<-ggplot(data=ms.data,aes(x,y))+geom_point(pch=pch,aes(colour=col),show.legend=TRUE)+
              ##scale_colour_manual(values=c(2,3))+
             # geom_text(aes(label=as.character(outlabel)),hjust=0,vjust=0,size=2)+ylim(0,64)+
              #xlab("||MO||")+ylab("VO")+labs(title="MS-Plot")+theme(plot.title = element_text(hjust = 0.5))
            plot(MO,vo,type="n",ylab="VO",xlab="||MO||",main=paste("MS-Plot on",ifelse(Met=="multi","fitted data (multi-variables)","original data")))
            for (i in setdiff(1:length(outlabel),c(out.dir,medcurve)))
            {
              
              points(MO[i],vo[i],pch=pch[i],col=col.normal,cex=0.71)
            }
            for (i in out.dir)
            {
              points(MO[i],vo[i],pch=pch[i],col=col.out,cex=0.7)
            }
            points(MO[medcurve],vo[medcurve],col=col.med,pch=pch[medcurve])
            text(MO+(range(MO)[2]-range(MO)[1])/70,vo+runif(length(vo),(range(vo)[2]-range(vo)[1])/60,(range(vo)[2]-range(vo)[1])/45),labels = outlabel,col="green",cex=0.5)
           }
          
          return(r1=list(mo=mo,vo=vo,out.dir=out.dir,medcurve=medcurve
                         #,p=p
                         ))
        }
        
        #plot(MO,vo,xlab="MO",ylab="VO",type="n",pch=pch,col=col,cex=0.1,...)
        #points(MO,vo,pch=19,cex=1,col=col.point,...)
      }
      return(r1=list(mo=mo,vo=vo,out.dir=out.dir,medcurve=medcurve))
    }
    if (!dirout)
    {
      mo<-result$out_avr
      MO <-apply(mo^2,1,sum)^(1/2)
      xlim=c(0,max(mo)*1.05)
      vo <-result$out_var
      
      if (plot==TRUE)
      {
        if (plot.type=="parallel")
        {
          data.ms=data.frame(MO=mo,VO=vo) 
          p<-ggparcoord(data.ms,columns=1:(d+1),showPoints = TRUE)
        }
        
        if (plot.type=="scatter")
        {
          ms.data=data.frame(x=MO,y=vo,col=col)
          p<-ggplot(data=ms.data,aes(x,y))+geom_point(colour=col,show.legend=FALSE)+
            xlab("||MO||")+ylab("VO")+labs(title="MS-Plot")+theme(plot.title = element_text(hjust = 0.5))
        }
        
        return(list(mo=mo,vo=vo,p=p))
      }
      return(list(mo=mo,vo=vo))
      
      #plot(mo,vo,type="p",pch=19,xlim=xlim,col=col,...)
    }
  }
  
}
  




