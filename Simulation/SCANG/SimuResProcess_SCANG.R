rm(list=ls())
library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr) 
library(foreach)
library(doParallel)
source("../R_Func/SCANG_Func.R")
savedatafile="/hpc/home/xl110/SATET/Results/"
loadpop="/hpc/group/chsi/xl110/SATET/"


Get_SCANG_res <- function(start,end,ss.factor,seeds,Alt_dcn=alt_dcn,covar){
  if(covar){
    snp_loc <- Snp_loc_cov%>%filter(seed==seeds,SS_Factor==ss.factor)
  }else{
    snp_loc <- Snp_loc%>%filter(seed==seeds,SS_Factor==ss.factor)
  }

  alt_snps <- snp_loc$snp_name2[which(snp_loc$snp_name%in%sim_map$alt_snps)]
  null_snps <- snp_loc$snp_name2[which(!(snp_loc$snp_name%in%sim_map$alt_snps))]
  
  snpo <- outer(start,alt_snps,"<=")*outer(end,alt_snps,">=")
  rejo <- alt_snps[colSums(snpo)>0]
  Domain.o <- uniqueN(snp_loc[rejo,"DomainID"])/Alt_dcn
  CV.power.o <- sum(colSums(snpo)>0)/length(sim_map$alt_snps)
  FDP.node.o <- mean(1-rowMeans(snpo))
  
  rejo <- NULL
  for(j in seq_along(start)){
    rejo <- union(rejo,start[j]:end[j])
  }
  alt.ds <- unique(unlist(sim_map$domain_w_sa_alt_snps))
  
  rejdsID.o <- as.numeric(unique(snp_loc$DomainID[rejo]))
  
  FDP.dc.o=uniqueN(intersect(rejdsID.o,alt.ds))/uniqueN(rejdsID.o)
  FPR.dc.p=uniqueN(setdiff(rejdsID.o,alt.ds))/length(sim_map$domain_snp_list)
  Prop.snp.o=length(intersect(alt_snps,rejo))/length(rejo)
  
  return(c("SNP_power"=CV.power.o,
           "SNP_FDP"=1-Prop.snp.o,
           "SNP_FPR"=uniqueN(setdiff(rejo,alt_snps))/uniqueN(sim_map$null_snps),
           "SNP_F1"=2/(1/CV.power.o+1/Prop.snp.o),
           "DC_power"=Domain.o,
           "DC_FDP"=FDP.dc.o,
           "DC_FPR"=FPR.dc.p,
           "DC_F1"=2/(1/Domain.o+1/(1-FDP.dc.o)),
           "length"=median(end-start+1)
  ))
}

# seeds=1;ss.factor=4;a=0.95;test="SCANGO"
# t1 <- Sys.time()
# outsc <- SCANG_sim$SCANG%>%filter(#SS_Factor==ss.factor,
#   #alpha==a,
#   #Test==test,
#   seed==seeds
# )%>%
#   group_by(alpha,SS_Factor,Test)%>%
#   summarize(out=Get_SCANG_res(Start,End,unique(SS_Factor),seeds))%>%
#   mutate(index=c("SNP_power","Node_FDP","DC_power",
#                  "SNP_FDP", "DC_FDP"))%>%ungroup()
# t2 <- Sys.time()

# see=SCANG_sim$SCANG%>%filter(seed==seeds,alpha==0.05)%>%
#   group_by(alpha,SS_Factor,Test,covar=FALSE)%>%
#   summarize(out=Get_SCANG_res(Start,End,unique(SS_Factor),seeds,covar=FALSE))%>%
#   mutate(index=c("SNP_power","SNP_FDP","SNP_FPR",
#                  "DC_power","DC_FDP","DC_FPR","length"))%>%ungroup()

nclust=50
cl <- makeCluster(nclust)
registerDoParallel(cl)
for(pis in c(0.2,0.4,0.6,0.8)){
  load(paste0(loadpop,"sim_map",pis*10,".RData"))
  load(paste0(savedatafile,"Output_SCANG_ROC",pis*10,".RData"))
  load(paste0(savedatafile,"LOC_SCANG_ROC",pis*10,".RData"))
  alt_dcn <- uniqueN(unlist(sim_map$domain_w_sa_alt_snps))
  
  t1=Sys.time()
  outdc=foreach(seeds=seq(100),.combine=rbind,
                .packages=c("data.table","stringi","tidyverse","dplyr"))%dopar%{
                  SCANG_sim$SCANG%>%filter(seed==seeds,alpha==0.05)%>%
                    group_by(alpha,SS_Factor,Test,covar=FALSE)%>%
                    summarize(out=Get_SCANG_res(Start,End,unique(SS_Factor),seeds,covar=FALSE))%>%
                    mutate(index=c("SNP_power","SNP_FDP","SNP_FPR","SNP_F1",
                                   "DC_power","DC_FDP","DC_FPR","DC_F1","length"))%>%ungroup()
                }
  outdc_cov=foreach(seeds=seq(100),.combine=rbind,
                .packages=c("data.table","stringi","tidyverse","dplyr"))%dopar%{
                  SCANG_sim_cov$SCANG%>%filter(seed==seeds,alpha==0.05)%>%
                    group_by(alpha,SS_Factor,Test,covar=TRUE)%>%
                    summarize(out=Get_SCANG_res(Start,End,unique(SS_Factor),seeds,covar=TRUE))%>%
                    mutate(index=c("SNP_power","SNP_FDP","SNP_FPR","SNP_F1",
                                   "DC_power","DC_FDP","DC_FPR","DC_F1","length"))%>%ungroup()
                }
  save(outdc,outdc_cov,file=paste0(savedatafile,"Output_SCANG_ROC2_",pis*10,".RData"))
  print(Sys.time()-t1)
}
stopCluster(cl)

