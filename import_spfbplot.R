setwd("~/Desktop/KAUST/project/multivariate_sparse_functional_data/code")
library(compiler)
library(fds)
library(dplyr)
library(ggplot2)
library(gdata)
library(mfaces)
library(face)
library(ggpubr)
library(spatstat)
library(roahd)
library(refund)
library(fdapace)
library(readr)
library(abind)
library(doParallel)
library(fda.usc)
library(plotfunctions)
library(RandomFields)
source("calcSigma.R")
source("facCal.R")
source("facCal_num.R")
source("Dirout.R")
source("msplot_modified.R")

source("sparse_fbplot_construct.R")
source("intensity_sparse_fbplot_construct.R")

fit_for_depth<-function(what,sett)
{
  func<-list()
  for (i in 1:length(sett))
  {
    select_index<-c()
    for (nb in 1:p)
    {select_index<-c(select_index,i+(nb-1)*length(sett))
    }
    func[[i]]<-what[,select_index]
  }
  return (func)
}
