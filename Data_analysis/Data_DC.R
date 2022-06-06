rm(list=ls())

library(foreach)
library(doParallel)
library(abind)
library(data.table)
library(tidyverse)
library(reshape2)

source("../R_Func/DATED_Func.R")
source("../R_Func/SPA_functions.R")

loadpop="/hpc/group/chsi/xl110/SATET/"
savedatafile="/hpc/home/xl110/SATET/Results/"

load(paste0(loadpop,"ALS_Preprocessed_Data_Euro_Sample.RData"))
snp_dat <- analysis_Euro_dt%>%dplyr::rename(Sample.Name=ID)
snp_dat2 <- ALS_data_SNV_no_LoF_euro%>%dplyr::rename(Sample.Name=ID) 
snp_dat <- snp_dat%>%
   left_join(snp_dat2,by=c("Sample.Name","Variant.ID","Sample.Type","chr","loc_adj"))%>%
   mutate(Polyphen.Humdiv.Prediction=as.character(Polyphen.Humdiv.Prediction))%>%
   mutate(Gene.Name=as.character(Gene.Name))

snp_dat <- snp_dat%>%group_by(domainID)%>%
  mutate(len=n())%>%
  ungroup()%>%
  filter(len>=5)%>%
  #arrange(loc_adj)%>%
  mutate(snpID0=snpID,snpID=rleid(loc_adj),Sample.Type=as.character(Sample.Type))%>%
  arrange(-desc(Sample.Type),Sample.Name)%>%
  mutate(Sample.ID=rleid(Sample.Name)-1)

glm_input <- snp_dat%>%
  dplyr::select(c("Sample.ID","Sample.Type"))%>%
  distinct_all()%>%
  mutate(types=ifelse(Sample.Type=="case",1,0))%>%
  dplyr::select(c("types"))

N <- uniqueN(snp_dat$Sample.ID)
nsnp <- uniqueN(snp_dat$snpID)
# Get Gmat case and control
mat_all <- new("dgTMatrix",
               i = as.integer(snp_dat$Sample.ID),
               j = as.integer(snp_dat$snpID-1), x=rep(1,nrow(snp_dat)), 
               Dim = c(N, nsnp))
N1=sum(glm_input)
N0=N-N1

# snp_dat2 <- ALS_data_SNV_no_LoF_euro%>%dplyr::rename(Sample.Name=ID) 
# snp_dat3 <- snp_dat%>%
#   left_join(snp_dat2,by=c("Sample.Name","Variant.ID","Sample.Type","chr","loc_adj"))%>%
#   mutate(Polyphen.Humdiv.Prediction=as.character(Polyphen.Humdiv.Prediction))%>%
#   mutate(Gene.Name=as.character(Gene.Name))
# 
# snp_dat3 <- count_domain_size(snp_dat3)
# snp_dat0 <- snp_dat3%>%filter(counts>=5)%>%
#   mutate_at("Sample.Name",as.character) %>% 
#   arrange(domainID,snpID)%>%
#   rename(snpID0=snpID,domainID0=domainID)%>%
#   mutate(snpID=rleid(snpID0),domainID=rleid(domainID0))%>%
#   arrange(Sample.Type,Sample.Name)%>%
#   mutate(sampleID=rleid(Sample.Name),typeID=as.integer(Sample.Type=="case"))

# sd <- snp_dat0[,c("sampleID","typeID")]%>%distinct_all()

# mat_all <- new("dgTMatrix",
#           i = as.integer(snp_dat0$sampleID-1),
#           j = as.integer(snp_dat0$snpID-1), x=rep(1,nrow(snp_dat0)),
#           Dim=c(max(snp_dat0$sampleID),max(snp_dat0$snpID)))
# glm_input <- matrix(sd$typeID,ncol=1)

# Transform the data so that we only consider domains with size>thresh_val

### for domain collapsing method
t1 <- Sys.time()
out_dc1<-SATET_dc(snp_dat=snp_dat,
                  mat_all=mat_all,
                  glm_input=glm_input,
                  teststat="FET",
                  score.cutoff=2,
                  alpha=0.05)
t2 <- Sys.time()
out_dc2<-SATET_dc(snp_dat=snp_dat,
                  mat_all=mat_all,
                  glm_input=glm_input,
                  teststat="score",
                  score.cutoff=2,
                  alpha=0.05)
t3 <- Sys.time()

t2-t1
t3-t2


q(save="no")
