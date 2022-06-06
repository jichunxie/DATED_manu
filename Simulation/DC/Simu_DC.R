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
source("../R_Func/DATED_Func.R")
source("Simu_COMMON_PARS.R")
source("../R_Func/SPA_functions.R")
#loadpop="/hpc/group/biostat/xl110/SATET/"
loadpop="/hpc/group/chsi/xl110/SATET/"
savedatafile="/hpc/home/xl110/SATET/Results/"


pop_name <- paste0(c("pop","pop_cov"),pi1*10,".RData")
output_name <- paste0(loadpop,"Output_DC",pi1*10,".RData")

comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

cl <- makeCluster(30)
registerDoParallel(cl)
load(paste0(loadpop,pop_name[1]))
output=foreach(r=1:R,.combine=rbind,
          .packages=c("purrr","data.table","dplyr",
                    "pracma","Matrix","reshape2","tidyverse"))%dopar%
  SATET_simu_dc(N1=N1,pop=pop,
             sampsize.factors=sampsize.factors,
             covars=FALSE,
             SampleID=SampleID[[r]],
             thresh_val=thresh_val,seed=r)

load(paste0(loadpop,pop_name[2]))
output_cov=foreach(r=1:R,.combine=rbind,
                .packages=c("purrr","data.table","dplyr",
                            "pracma","Matrix","reshape2","tidyverse"))%dopar%
  SATET_simu_dc(N1=N1,pop=pop,
             sampsize.factors=sampsize.factors,
             covars=TRUE,
             SampleID=SampleID[[r]],
             thresh_val=thresh_val,seed=r)
stopCluster(cl)

save(output,output_cov,file=output_name)


q(save="no")
