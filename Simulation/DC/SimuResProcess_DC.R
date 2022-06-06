rm(list=ls())

library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr) 
library(foreach)
library(doParallel)
source("../R_Func/DATED_Func.R")
savedatafile="/hpc/home/xl110/SATET/Results/"
loadpop="/hpc/group/chsi/xl110/SATET/"

Out <- Out_cov <- NULL
for(pis in c(0.2,0.4,0.6,0.8)){
  t1 <- Sys.time()
  load(paste0(loadpop,"Output_DC",pis*10,".RData"))
  load(paste0(loadpop,"sim_map",pis*10,".RData"))
  nonnulls = sim_map$alt_snps
  nulls = sim_map$null_snps
  m_dc = length(sim_map$domain_snp_list)
  nonnulls_dc = unlist(sim_map$domain_w_sa_alt_snps)
  nulls_dc = setdiff(seq(m_dc),nonnulls_dc)
  
  out <- out_cov <- NULL
  
  Padj=output%>%dplyr::select(c("pvals","domainID","Test","SS.Factor","Seed","Covar"))%>%
    distinct_all()%>%group_by(Test,SS.Factor,Seed,Covar)%>%
    mutate(padj=p.adjust(pvals,method="BH"))%>%ungroup()%>%
    left_join(output,by=c("domainID","Test","SS.Factor","Seed","Covar"))
  Padj_cov=output_cov%>%dplyr::select(c("pvals","domainID","Test","SS.Factor","Seed","Covar"))%>%
    distinct_all()%>%group_by(Test,SS.Factor,Seed,Covar)%>%
    mutate(padj=p.adjust(pvals,method="BH"))%>%ungroup()%>%
    left_join(output_cov,by=c("domainID","Test","SS.Factor","Seed","Covar"))
  t2 <- Sys.time()
  for(a in seq(0.05,0.95,0.05)){
    out=Padj%>%filter(padj<=a)%>%group_by(Test,SS.Factor,Covar,Seed)%>%
      summarize("DC_FDP"=uniqueN(intersect(domainID,nulls_dc))/uniqueN(domainID),
                "DC_power"=uniqueN(intersect(domainID,nonnulls_dc))/length(nonnulls_dc),
                "DC_FPR"=uniqueN(intersect(domainID,nulls_dc))/length(nulls_dc),
                "SNP_FDP"=uniqueN(intersect(snpID,nulls))/uniqueN(snpID),
                "SNP_power"=uniqueN(intersect(snpID,nonnulls))/length(nonnulls), #true positive rate
                "SNP_FPR"=uniqueN(intersect(snpID,nulls))/length(nulls),.groups = 'drop')%>%
      ungroup()%>%mutate(alpha=a)%>%bind_rows(out)
    out_cov=Padj_cov%>%filter(padj<=a)%>%group_by(Test,SS.Factor,Covar,Seed)%>%
      summarize("DC_FDP"=uniqueN(intersect(domainID,nulls_dc))/uniqueN(domainID),
                "DC_power"=uniqueN(intersect(domainID,nonnulls_dc))/length(nonnulls_dc),
                "DC_FPR"=uniqueN(intersect(domainID,nulls_dc))/length(nulls_dc),
                "SNP_FDP"=uniqueN(intersect(snpID,nulls))/uniqueN(snpID),
                "SNP_power"=uniqueN(intersect(snpID,nonnulls))/length(nonnulls), #true positive rate
                "SNP_FPR"=uniqueN(intersect(snpID,nulls))/length(nulls),.groups = 'drop')%>%
      ungroup()%>%mutate(alpha=a)%>%bind_rows(out_cov)
    print(a)
  }
  Out <- out%>%group_by(Test,SS.Factor,Covar,alpha) %>%
    summarise(DC_power=mean(DC_power),DC_FPR=mean(DC_FPR),DC_FDP=mean(DC_FDP),
              SNP_power=mean(SNP_power),SNP_FPR=mean(SNP_FPR),SNP_FDP=mean(SNP_FDP))%>%
    ungroup()%>%mutate(p1s=pis)%>%bind_rows(Out)
  Out_cov <- out_cov%>%group_by(Test,SS.Factor,Covar,alpha) %>%
    summarise(DC_power=mean(DC_power),DC_FPR=mean(DC_FPR),DC_FDP=mean(DC_FDP),
              SNP_power=mean(SNP_power),SNP_FPR=mean(SNP_FPR),SNP_FDP=mean(SNP_FDP))%>%
    ungroup()%>%mutate(p1s=pis)%>%bind_rows(Out_cov)
  t3=Sys.time()
  print(t3-t1);print(t3-t2)
}

out_dc=Out%>%bind_rows(Out_cov)
rm(Out,Out_cov)
save(out_dc,file=paste0(savedatafile,"Output_DC_ROC_processed.RData"))

str(out_dc)
dim(out_dc)

q(save="no")


