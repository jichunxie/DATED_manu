rm(list=ls())

library(foreach)
library(doParallel)
library(abind)
library(data.table)
library(tidyverse)
library(reshape2)


source("../R_Func/DYNATE_Func.R")
source("../R_Func/SPA_functions.R")
loadpop="/hpc/group/chsi/xl110/SATET/"
savedatafile="/hpc/home/xl110/SATET/Results/"

load(paste0(loadpop,"ALS_Preprocessed_Data_Euro_Sample.RData"))
# 2634 cases, 7504 controls, 289503 snps,34553 unique domains

str(geno_case_matrix)
snp_dat <- ALS_data_SNV_no_LoF_euro%>% 
  dplyr::select(c("Sample.Type","ID","subRVIS.Domain.Name",
                  "chrno","loc_adj"))%>%
  dplyr::rename(Sample.Name=ID)
snp_dat <- analysis_Euro_dt%>%dplyr::rename(Sample.Name=ID)

Gmat_case=geno_case_matrix
Gmat_ctrl=geno_ctrl_matrix
rm(All_euro_sub,ALS_data_SNV_no_LoF_euro,
   geno_matrix_all)

snp_dat <- snp_dat%>%group_by(domainID)%>%
  mutate(len=n())%>%
  ungroup()%>%
  filter(len>=5)%>%
  arrange(loc_adj)%>%
  mutate(snpID0=snpID,snpID=rleid(loc_adj),Sample.Type=as.character(Sample.Type))%>%
  arrange(-desc(Sample.Type),Sample.Name)%>%
  mutate(Sample.ID=rleid(Sample.Name)-1)



thresh_val=seq(7,35,2)

# # debug
# tvi <- thresh_val <- 7
# teststat="FET";glm_input=NULL;
# Gmat_case=NULL;Gmat_ctrl=NULL;
# struct_map=NULL
# midp=TRUE;score.cutoff=2;use.SPA.score=TRUE;
# L=5;alpha=0.05;alpha1=NULL



cl <- makeCluster(15)
registerDoParallel(cl)

t1 <- Sys.time()
out1 <- foreach(tvi=thresh_val,.combine=rbind,
                .packages=c("purrr","data.table","dplyr",
                            "pracma","Matrix","reshape2",
                            "tidyverse"))%dopar%{
                              out <- Test_DGS(snp_dat=snp_dat,thresh_val=tvi,      
                                              teststat="FET")# should use SATET in 1128
                              out <- out%>%mutate(thresh_val=tvi)
                              SATET_data(out,L=5,alpha=0.05,thresh_val=tvi)
                            }

t2 <- Sys.time()

out2 <- foreach(tvi=thresh_val,.combine=rbind,
                .packages=c("purrr","data.table","dplyr",
                            "pracma","Matrix","reshape2",
                            "tidyverse"))%dopar%{
                              out <- Test_DGS(snp_dat=snp_dat,thresh_val=tvi,
                                              teststat="score")# should use SATET in 1128
                              out <- out%>%mutate(thresh_val=tvi)
                              SATET_data(out,L=5,alpha=0.05,thresh_val=tvi)
                            }
t3 <- Sys.time()

stopCluster(cl)

t2-t1
t3-t2

# old function using SATET is SATET_func


save(out1,out2,file=paste0(loadpop,"Dataanalysis_SATET_filter_multithresh_new.RData"))


q(save="no")
