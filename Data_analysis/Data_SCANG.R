rm(list=ls())

library(foreach)
library(doParallel)
library(SCANG)
library(data.table)
library(tidyverse)
library(reshape2)
library(Matrix)
source("../R_Func/SCANG_Func.R")

loadpop="/hpc/group/chsi/xl110/SATET/"
savedatafile="/hpc/home/xl110/SATET/Results/"

load(paste0(loadpop,"ALS_Preprocessed_Data_Euro_Sample.RData"))

# [1] "All_euro_sub"             "ALS_data_SNV_no_LoF_euro"
# [3] "analysis_Euro_dt"         "geno_case_matrix"        
# [5] "geno_ctrl_matrix"         "geno_matrix_all"         
# [7]  "Ncases"                  "Nctrls"  
# 2634 cases, 7504 controls, 289503 snps,34553 unique domains


snp_dat <- analysis_Euro_dt%>%dplyr::rename(Sample.Name=ID)

#filter begin
snp_dat <- snp_dat%>%group_by(domainID)%>%
  mutate(len=n())%>%
  ungroup()%>%
  filter(len>=5)%>%
  arrange(loc_adj)%>%
  mutate(snpID0=snpID,snpID=rleid(loc_adj),Sample.Type=as.character(Sample.Type))%>%
  arrange(-desc(Sample.Type),Sample.Name)%>%
  mutate(Sample.ID=rleid(Sample.Name)-1)
#filter end

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


snp_loc <- snp_dat%>%
  dplyr::select(c("domain","chr","loc_adj","snpID"))%>%
  distinct_all()%>%
  arrange(loc_adj)%>%
  dplyr::rename(snp_loc=loc_adj)%>%
  mutate(snp_name2=seq(n()))

#### Find Lmin and Lmax
t1 <- Sys.time()
Ls <- Find_L(snp_loc)
t2 <- Sys.time()
t2-t1
#Ls,4min
Lmin <- 1
Lmax <- 103


phenotypedata <- data.frame(phenotype=glm_input$types)

genotype <- as(mat_all,"dgCMatrix")


# snp_loc <- snp_dat%>%
#   dplyr::select(c("chr","loc_adj","snpID"))%>%
#   distinct_all()%>%
#   dplyr::rename(snp_loc=loc_adj)%>%
#   mutate("loc_order"=order(snp_loc))



#genotype <- genotype0[,snp_loc$loc_order]

# snp_loc <- snp_loc%>%
#   arrange(snp_loc)%>%
#   mutate(snp_name2=seq(n()))


t1 <- Sys.time()

obj_nullmodel <- fit_null_glm_SCANG(phenotype~1,data=phenotypedata,family=binomial)
res_lm <- SCANG(genotype,obj_nullmodel,Lmin,Lmax,filter=2e-5,f=0.5)

ptime <- Sys.time()-t1
ptime

save(snp_loc,res_lm,file="Dataanalysis_filter_SCANG.RData")


q(save="no")
