rm(list=ls())

####### File Specific Parameters 
pi1=0.2# proportion of pathogenic RVs. 

library("MASS")
library("Matrix")
library(SCANG)
library(tidyverse)
library(data.table)
library(dplyr)
library(doParallel)
library(foreach)
source("../R_Func/SCANG_Func.R")
source("Simu_COMMON_PARS_SCANG.R")
loadpop="/hpc/group/chsi/xl110/SATET/"
savedatafile="/hpc/home/xl110/SATET/Results/"

pop_name <- paste0(c("pop","pop_cov"),pi1*10,".RData")
output_name <- paste0(savedatafile,"LOC_SCANG_ROC",pi1*10,".RData")

comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

cl <- makeCluster(50)
registerDoParallel(cl)
ptm <- proc.time()
load(paste0(loadpop,pop_name[1]))
Snp_loc = foreach(r=1:R, .combine=rbind,
                         .packages=c("data.table","dplyr","Matrix",
                                     "SCANG","MASS","tidyverse")) %dopar%
  Get_Alt_info(N1=N1,
                      pop=pop,
                      sampsize.factors=sampsize.factors,
                      sim_map=sim_map,
                      params=params,
                      covars=FALSE,
                      filter=filter,
                      f=f,
                      Rocp=Rocp,
                      seed=seed[r],
                SampleID=SampleID[[r]])

load(paste0(loadpop,pop_name[2]))
Snp_loc_cov = foreach(r=1:R, .combine=rbind,
                    .packages=c("data.table","dplyr","Matrix",
                                "SCANG","MASS","tidyverse")) %dopar%
  Get_Alt_info(N1=N1,
                pop=pop,
                sampsize.factors=sampsize.factors,
                sim_map=sim_map,
                params=params,
                covars=TRUE,
                filter=filter,
                f=f,
                Rocp=Rocp,
                seed=seed[r],
                SampleID=SampleID[[r]])


proc.time()-ptm
stopCluster(cl)


save(Snp_loc,Snp_loc_cov,file=output_name)


q(save="no")
