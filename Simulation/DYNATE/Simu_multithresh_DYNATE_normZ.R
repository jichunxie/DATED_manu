rm(list=ls())

####### File Specific Parameters 
pi1=0.2# proportion of pathogenic RVs. 

library(foreach)
library(doParallel)
library(abind)
library(data.table)
library(dplyr)
library(tidyverse)
library(reshape2)
library(pracma)
set.seed(123)
source("../R_Func/DYNATE_Func.R")
source("../R_Func/SPA_functions.R")
source("Simu_COMMON_PARS.R")

loadpop="/hpc/group/chsi/xl110/SATET/"
savedatafile="/hpc/group/xielab/xl110/SATET/"


pop_name <- paste0(c("pop","pop_cov"),pi1*10,".RData")
output_name <- paste0(savedatafile,"Output_mutithresh",pi1*10,".RData")

comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

############## Main Simulation #############

cl <- makeCluster(40)
registerDoParallel(cl)

load(paste0(loadpop,pop_name[1]))
output=foreach(thv=multi_thresh,.combine=rbind)%:%foreach(r=1:R,.combine=rbind,
               .packages=c("purrr","data.table","dplyr",
                           "pracma","Matrix","reshape2","tidyverse"))%dopar%{
                             out=SATET_simu(N1=N1,pop=pop,
                                            sampsize.factors=sampsize.factors,
                                            covars=FALSE,
                                            SampleID=SampleID[[r]],
                                            thresh_val=thv,seed=r)
                             out <- out$DGS%>%mutate(thresh_val=thv)
                             out
                           }


load(paste0(loadpop,pop_name[2]))
output_cov=foreach(thv=multi_thresh,.combine=rbind)%:%foreach(r=1:R,.combine=rbind,
                                                              .packages=c("purrr","data.table","dplyr",
                                                                          "pracma","Matrix","reshape2","tidyverse"))%dopar%{
                                                                            out=SATET_simu(N1=N1,pop=pop,
                                                                                           sampsize.factors=sampsize.factors,
                                                                                           covars=TRUE,
                                                                                           SampleID=SampleID[[r]],
                                                                                           thresh_val=thv,seed=r)
                                                                            out <- out$DGS%>%mutate(thresh_val=thv)
                                                                            out
                                                                          }



stopCluster(cl)

save(output,output_cov,file=output_name)



q(save="no")
