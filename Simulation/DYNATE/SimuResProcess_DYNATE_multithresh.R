rm(list=ls())
#setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/SATET/Results")
library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr) 
library(foreach)
library(doParallel)
source("../R_Func/DYNATE_Func.R")
source("Simu_COMMON_PARS.R")
savedatafile="/hpc/home/xl110/SATET/Results/"
loadpop="/hpc/group/chsi/xl110/SATET/"
loadpop2="/hpc/group/xielab/xl110/SATET/"

pis=0.2

load(paste0(loadpop2,"Output_mutithresh",pis*10,".RData"))
output_cov <- output_cov%>%filter(thresh_val%in% multi_thresh)
output <- output%>%filter(thresh_val%in% multi_thresh)

cl <- makeCluster(40)
registerDoParallel(cl)
load(paste0(loadpop,"sim_map",pis*10,".RData"))

t1 <- Sys.time()
out=foreach(thv=multi_thresh,.combine=rbind)%:%
  foreach(ss.factor=seq(4),.combine=rbind)%:%
  foreach(covar=c(FALSE,TRUE),.combine=rbind)%:%
  foreach(test=c("FET","score"),.combine=rbind)%:%
  foreach(seeds=seq(100),.combine=rbind,
          .packages=c("data.table","stringi","tidyverse","dplyr"))%dopar%{
            if(covar){
              struct_map <- output_cov%>%filter(Test==test,
                                                SS.Factor==ss.factor,Seed==seeds,
                                                thresh_val==thv)
            }else{
              struct_map <- output%>%filter(Test==test,
                                            SS.Factor==ss.factor,Seed==seeds,
                                            thresh_val==thv)
            }
            SATET_multithresh(struct_map,L=4,sim_map=sim_map,thresh_val=thv,seed=seeds)
          }
save(out,file=paste0(savedatafile,"Output_SATET_2multilayer",pis*10,".RData"))# for thresh 13,15,17,19,21


print(Sys.time()-t1)

stopCluster(cl)


