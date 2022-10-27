
rm(list=ls())

#SATET: avg discovered layer, median length

library(foreach)
library(doParallel)
library(abind)
library(data.table)
library(tidyverse)
library(reshape2)

source("R_Func/DYNATE_Func.R")
load("Results/Dataanalysis_SATET_filter_multithresh_new.RData")
#out1n <- out1;out2n <- out2

#load("Results/Dataanalysis_SATET_multithresh.RData")
load("Data/ALS_Preprocessed_Data_Euro_Sample.RData")

#check <- analysis_Euro_dt%>%dplyr::select(c("Sample.Type",""))

snp_dat <- analysis_Euro_dt%>%dplyr::rename(Sample.Name=ID)

snp_dat <- snp_dat%>%group_by(domainID)%>%
  mutate(len=n())%>%
  ungroup()%>%
  filter(len>=5)%>%
  arrange(loc_adj)%>%
  mutate(snpID0=snpID,snpID=rleid(loc_adj),
         Sample.Type=as.character(Sample.Type))%>%
  arrange(-desc(Sample.Type),Sample.Name)%>%
  mutate(Sample.ID=rleid(Sample.Name)-1)

snp_dat2 <- ALS_data_SNV_no_LoF_euro%>%dplyr::rename(Sample.Name=ID) 
snp_dat3 <- snp_dat%>%
  left_join(snp_dat2,by=c("Sample.Name","Variant.ID","Sample.Type","chr","loc_adj"))%>%
  mutate(Polyphen.Humdiv.Prediction=as.character(Polyphen.Humdiv.Prediction))%>%
  mutate(Gene.Name=as.character(Gene.Name))

# average case/ctrl maf
see=snp_dat3%>%dplyr::select(c("snpID","Case.Maf","Ctrl.Maf"))%>%distinct_all()

  
  
  

snp_dat4 <- snp_dat3%>%mutate(Sample.Type=as.character(Sample.Type))

snp_score <- snp_dat3 %>% dplyr::select(c("domainID","snpID","domain","chrno","Gene.Name",
                                          "Polyphen.Humdiv.Prediction"))%>%
  distinct_all()



#region info:  
snp_loc <- snp_dat3%>%dplyr::select(c("snpID","loc_adj","chr"))%>%
  distinct_all()%>%
  mutate(chr=factor(chr,levels=paste0("chr",1:23)))%>%
  arrange(loc_adj)%>%mutate(snp_name2=seq(n()))

method1="DYNATE-FL";method2="DYNATE-SS"
tvi_max=7;L_max=5
out_satet <- out1%>%mutate("Method"=method1)%>%
  bind_rows(out2%>%mutate("Method"=method2))%>%
  filter(thresh_val<=tvi_max) %>% 
  mutate(L1=ifelse(Layer==1,L1,NA)) %>%
  mutate(L2=ifelse(Layer==2,L2,NA))%>%
  mutate(L3=ifelse(Layer==3,L3,NA))%>%
  mutate(L4=ifelse(Layer==4,L4,NA))
out_satet <- out_satet%>%
  mutate(DGS=rowSums(dplyr::select(out_satet,paste0("L",1:L_max)),
                     na.rm=TRUE))%>%
  dplyr::select(-paste0("L",0:L_max))%>%
  filter(Layer<=L_max)%>%
  mutate(Region =paste(Method,thresh_val,Layer,DGS,sep="-"))%>%
  dplyr::select(-c("pvals1","wt","Layer","DGS",
                   "Seed","Test","thresh_val"))
satet_se <-out_satet%>%left_join(snp_loc,by="snpID")%>%
  group_by(Region,Method,chr)%>%
  mutate(group=cumsum(c(1L, diff(snp_name2)) > 1L))%>%
  group_by(group,.add=TRUE)%>%
  summarise(start=min(snp_name2),end=max(snp_name2))
  

#Odds_ratio_region(regions,"snpID","domain",Rej)

Odds_ratio_region <- function(var,var.name,reg,Rej){
  Vars <- unique((Rej%>%filter_at(reg,all_vars(.==var)))$snpID)
  sd <- data.frame(snp_dat4%>%dplyr::select(var.name))[,1]
  pts <- snp_dat4%>%
    mutate(Var=as.integer(sd%in%Vars))%>%
    group_by(Sample.Name,Sample.Type)%>%
    summarize(Var=ifelse(sum(Var)>0,1,0))
  if(sum(pts$Var==1)==0){
    return(NA)
  }else{
    ft <- table(pts$Var,pts$Sample.Type)
    or <- round(ft[2,1]*ft[1,2]/(ft[1,1]*ft[2,2]),1)
    upper <- round(exp(log(or)+1.96*sqrt(sum(1/c(ft)))),1)
    lower <- round(exp(log(or)-1.96*sqrt(sum(1/c(ft)))),1)
    return(paste0(or,"(",lower,",",upper,")"))
  }
}


tvi_max=14
L_max=5
if(tvi_max>7){method1="DYNATE_MT (FL)";method2="DYNATE_MT (SS)"}else{
  method1="DYNATE-FL";method2="DYNATE-SS"
}
# for output cvs leaf info
satet_leafinfo <- out1%>%mutate("Method"=method1)%>%
  bind_rows(out2%>%mutate("Method"=method2))%>%
  dplyr::select(c("domainID","L0","Method","pvals1"))%>%
  distinct_all()
  
#### DYNATE function
DYNATE_analysis <- function(tvi_max=14,L_max=5){
  if(tvi_max>7){method1="DYNATE_MT (FL)";method2="DYNATE_MT (SS)"}else{
    method1="DYNATE-FL";method2="DYNATE-SS"
  }
  out_satet <- out1%>%mutate("Method"=method1)%>%
    bind_rows(out2%>%mutate("Method"=method2))%>%
    filter(thresh_val<=tvi_max) %>% 
    mutate(L1=ifelse(Layer==1,L1,NA)) %>%
    mutate(L2=ifelse(Layer==2,L2,NA))%>%
    mutate(L3=ifelse(Layer==3,L3,NA))%>%
    mutate(L4=ifelse(Layer==4,L4,NA))
  out_satet <- out_satet%>%
    mutate(DGS=rowSums(dplyr::select(out_satet,paste0("L",1:L_max)),na.rm=TRUE))%>%
    dplyr::select(-paste0("L",0:L_max))%>%
    filter(Layer<=L_max)%>%
    mutate(Region =paste(Method,thresh_val,Layer,DGS,sep="-"))
  
  Rej_satet <- out_satet%>%left_join(snp_score,by=c("snpID"))%>%group_by(Region)%>%
    mutate(Polyphen.Humdiv.Prediction=ifelse(Polyphen.Humdiv.Prediction=="NA","unknown",
                                             Polyphen.Humdiv.Prediction))%>%
    mutate(Gene.Name=gsub('^.|.$', '', Gene.Name))
  
  Rej_satet_flag <- Rej_satet%>%
    dplyr::select(-c("snpID","DGS","Polyphen.Humdiv.Prediction"))%>%
    distinct_all()%>%arrange(Method,thresh_val,Layer)%>%
    group_by(Method,thresh_val)%>%
    mutate(first_dc=!duplicated(domain),
           first_gene=!duplicated(Gene.Name))%>%ungroup()
  
  
  Rej_satet.region2 <- Rej_satet%>%dplyr::select(c("snpID","Region","Method"))%>%
    distinct_all()%>%group_by(Region,Method)%>%
    summarise(Region_len=n())%>%ungroup()%>%
    group_by(Method)%>%
    mutate(Med_len=median(Region_len))%>%ungroup()
  Rej_satet.region <- Rej_satet%>%#dplyr::select(-Delt)%>%
    distinct_all()%>%
    dcast(Region+Method~Polyphen.Humdiv.Prediction)%>%
    mutate(Polyphen.div=paste(benign,possibly,probably,sep="/"))%>%
    left_join(Rej_satet.region2,by=c("Region","Method"))
  
  Rej_satet.dc <- Rej_satet %>% # Not correct for polyphen score
    dcast(domain+Method~Polyphen.Humdiv.Prediction)%>%
    mutate(Polyphen.div=paste(benign,possibly,probably,sep="/"))
  Rej_satet.dc1 <- Rej_satet_flag%>%filter(first_dc)%>%
    group_by(Method,domain)%>%
    summarise(Layer_avg=mean(Layer))
  Rej_satet.dc <- Rej_satet.dc %>%left_join(Rej_satet.dc1,by=c("domain","Method"))
  
  Rej_satet.gene <- Rej_satet %>% 
    dcast(Gene.Name+Method~Polyphen.Humdiv.Prediction)%>%
    mutate(Polyphen.div=paste(benign,possibly,probably,sep="/"))
  Rej_satet.gene1 <- Rej_satet_flag%>%filter(first_gene)%>%
    group_by(Method,Gene.Name)%>%
    summarise(Layer_avg=mean(Layer))
  Rej_satet.gene <- Rej_satet.gene %>%left_join(Rej_satet.gene1,by=c("Gene.Name","Method"))
  
  
  ORs <- NULL 
  for(regions in Rej_satet.region$Region){
    ORs <- c(ORs,Odds_ratio_region(regions,"snpID","Region",Rej=Rej_satet))
  }
  Rej_satet.region <- Rej_satet.region%>%
    mutate(ORCI=ORs,OR=str_split(ORs,"[(]",simplify=TRUE)[,1])%>%
    mutate(OR=as.numeric(OR))
  
  ORs <- NULL 
  for(regions in Rej_satet.gene$Gene.Name){
    ORs <- c(ORs,Odds_ratio_region(regions,"snpID","Gene.Name",Rej=Rej_satet))
  }
  Rej_satet.gene <- Rej_satet.gene%>%
    mutate(ORCI=ORs,OR=str_split(ORs,"[(]",simplify=TRUE)[,1])%>%
    mutate(OR=as.numeric(OR))
  
  ORs <- NULL 
  for(regions in Rej_satet.dc$domain){
    ORs <- c(ORs,Odds_ratio_region(regions,"snpID","domain",Rej=Rej_satet))
  }
  Rej_satet.dc <- Rej_satet.dc%>%
    mutate(ORCI=ORs,OR=str_split(ORs,"[(]",simplify=TRUE)[,1])%>%
    mutate(OR=as.numeric(OR))
  return(list("region"=Rej_satet.region,"dc"=Rej_satet.dc,"gene"=Rej_satet.gene))
}


satet7 <- DYNATE_analysis(7,L_max=5)
Rej_satet.region1=satet7$region
Rej_satet.dc1=satet7$dc
Rej_satet.gene1=satet7$gene



# Simu:domain level, Region level, SNP level AUC ROC curve, by pi
# Data analysis: REgion level, table, OR, supp


#### SCANG
load("Results/Dataanalysis_filter_SCANG.RData")


snp_dat <- ALS_data_SNV_no_LoF_euro%>%
  dplyr::select(c("subRVIS.Domain.Name","Gene.Name","chr","loc_adj"))%>%
  dplyr::rename(snp_loc=loc_adj)%>%distinct_all()

snp_loc_info <- snp_loc%>%left_join(snp_dat,by=c("chr","snp_loc"))

reso <- res_lm$SCANG_O_res
ress <- res_lm$SCANG_S_res
resb <- res_lm$SCANG_B_res

rejb <- rejs <- rejo <- NULL
Rejo <- Rejs <- Rejb <- NULL
for(i in seq(nrow(reso))){
  Rejo <- rbind(Rejo,
                data.frame("snp_name2"=reso[i,2]:reso[i,3],
                           "Method"="SCANG-O",
                           "Region"=paste0("SCANG-O-",i)))
  rejo <- union(rejo,reso[i,2]:reso[i,3])
}
for(i in seq(nrow(ress))){
  Rejs <- rbind(Rejs,
                data.frame("snp_name2"=ress[i,2]:ress[i,3],
                           "Method"="SCANG-S","Region"=paste0("SCANG-S-",i)))
  rejs <- union(rejs,ress[i,2]:ress[i,3])
}
for(i in seq(nrow(resb))){
  Rejb <- rbind(Rejb,
                data.frame("snp_name2"=resb[i,2]:resb[i,3],
                           "Method"="SCANG-B","Region"=paste0("SCANG-B-",i)))
  rejb <- union(rejb,resb[i,2]:resb[i,3])
}

Rej <- rbind(Rejo,Rejs,Rejb)%>%left_join(snp_loc,by="snp_name2")%>%
  left_join(snp_score,by=c("snpID","domain"))%>%group_by(Region)%>%
  mutate(Polyphen.Humdiv.Prediction=ifelse(Polyphen.Humdiv.Prediction=="NA","unknown",
                                           Polyphen.Humdiv.Prediction))

scang_se <- Rej%>%
  group_by(Region,Method,chr)%>%
  mutate(group=cumsum(c(1L, diff(snp_name2)) > 1L))%>%
  group_by(group,.add=TRUE)%>%
  summarise(start=min(snp_name2),end=max(snp_name2))


Rej.region2 <- Rej%>%dplyr::select(c("snp_name2","Region","Method"))%>%
  distinct_all()%>%group_by(Region,Method)%>%
  summarise(Region_len=n())%>%ungroup()%>%
  group_by(Method)%>%
  mutate(Med_len=median(Region_len))
Rej.region <- Rej%>%
  distinct_all()%>%
  dcast(Region+Method~Polyphen.Humdiv.Prediction)%>%
  mutate(Polyphen.div=paste(benign,possibly,probably,sep="/"))%>%
  left_join(Rej.region2,by=c("Region","Method"))

Rej.dc <- Rej %>% # Not correct Polyphen score
  dcast(domain+Method~Polyphen.Humdiv.Prediction)%>%
  mutate(Polyphen.div=paste(benign,possibly,probably,sep="/"))
Rej.gene <- Rej %>% 
  dcast(Gene.Name+Method~Polyphen.Humdiv.Prediction)%>%
  mutate(Polyphen.div=paste(benign,possibly,probably,sep="/"))

scang_dc <- Rej%>%dplyr::select(c("Method","Region","domain"))%>%
  distinct_all()%>%
  group_by(Region)%>%
  summarize(Domain=paste(domain,collapse=", "))



#Odds_ratio_region("SCANG-O-1","snpID","Region")

###need change!
ORs <- NULL 
for(regions in Rej.region$Region){
  ORs <- c(ORs,Odds_ratio_region(regions,"snpID","Region",Rej))
}
Rej.region <- Rej.region%>%
  mutate(ORCI=ORs,OR=str_split(ORs,"[(]",simplify=TRUE)[,1])%>%
  mutate(OR=as.numeric(OR))


ORs <- NULL 
for(regions in Rej.gene$Gene.Name){
  ORs <- c(ORs,Odds_ratio_region(regions,"snpID","Gene.Name",Rej))
}
Rej.gene <- Rej.gene%>%
  mutate(ORCI=ORs,OR=str_split(ORs,"[(]",simplify=TRUE)[,1])%>%
  mutate(OR=as.numeric(OR))


ORs<- NULL 
for(regions in Rej.dc$domain){
  ORs <- c(ORs,Odds_ratio_region(regions,"snpID","domain",Rej))
}
Rej.dc <- Rej.dc%>%
  mutate(ORCI=ORs,OR=str_split(ORs,"[(]",simplify=TRUE)[,1])%>%
  mutate(OR=as.numeric(OR))

nam <- c("DYNATE-FL","DYNATE-SS","SCANG-O","SCANG-S","SCANG-B")
Region <- Rej_satet.region1%>%bind_rows(Rej.region)%>%mutate(Measure="Region")
Domain <- Rej_satet.dc1%>%bind_rows(Rej.dc)%>%mutate(Measure="Domain")
Gene <- Rej_satet.gene1%>%bind_rows(Rej.gene)%>%mutate(Measure="Gene")

ALL <-Region%>%bind_rows(Domain,Gene)%>%mutate(logOR=log(OR))%>%
  mutate(Method=factor(Method,levels=nam))

library(ggpubr)
library(ggplot2)
Q25 <- function(x){quantile(x,0.25)}
Q75 <- function(x){quantile(x,0.75)}
ALL2 <- ALL%>%filter(Measure=="Region")%>%
  #pivot_longer(cols=c("OR","Region_len"))%>%
  pivot_longer(cols=c("logOR","Region_len"))%>%
  mutate(Measure=name)

#ALL2%>%group_by(name,Method)%>%summarize(medor=median(value))

#data.frame(ALL2%>%filter(Measure=="logOR",Method=="SCANG-B"))

data_violin1=ggplot(ALL2%>%filter(Measure=="logOR"),
                   aes(Method,value))+geom_violin()+
  stat_summary(fun = median,
               geom = "pointrange",
               fun.min = Q25,
               fun.max = Q75,color="red")+
  ylab("logOR")+theme_classic()+ggtitle("A")+
  geom_hline(yintercept=0,linetype="dashed")+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1,size=9))

data_violin2=ggplot(ALL2%>%filter(Measure!="logOR"),
                    aes(Method,value))+
  geom_violin()+ggtitle("B")+
  stat_summary(fun = median,
               geom = "pointrange",
               fun.min = Q25,
               fun.max = Q75,color="red")+
  ylab("Region Length")+theme_classic()+
  geom_hline(yintercept=0,linetype="dashed")+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1,size=9))

pdf.options(reset = TRUE, onefile = FALSE)
pdf(paste0("Results/Data_Violin",".pdf"),height=3.5,width=7)
ggarrange(data_violin1,data_violin2,ncol=2)
dev.off()

sum0 <- ALL%>%group_by(Method,Measure)%>%summarise(medOR=median(OR))%>%
  ungroup()%>%
  pivot_wider(names_from=Measure,values_from = medOR)

sum1 <- ALL%>%group_by(Method,Measure)%>%summarise(Op_rate=mean(OR<1))%>%
  ungroup()%>%
  pivot_wider(names_from=Measure,values_from = Op_rate)

sum2 <- ALL%>%group_by(Method,Measure)%>%summarise(Count=n())%>%
  ungroup()%>%
  pivot_wider(names_from=Measure,values_from = Count)

sum3 <- ALL%>%filter(Measure=="Region")%>%group_by(Method)%>%
  summarise(Med_len=median(Region_len))

Sum <- sum1%>%left_join(sum2,by="Method")%>%left_join(sum3,by="Method")%>%
  setnames(c("Method","Domain_opr","Gene_opr","Region_opr",
             "Domain_N","Gene_N","Region_N","Med_len"))%>%
  dplyr::select(c("Method","Region_opr",
                  "Med_len",
                  "Domain_opr","Domain_N","Gene_opr","Gene_N"))


library(xtable)
print(xtable(Sum,digits=3),include.rownames=FALSE)



#### Domain information


see=Domain%>%dplyr::select(-c("benign","possibly",
                          "probably",
                          "Polyphen.div","Measure","OR"))%>%
  filter(substr(Method,6,8)!="_MT")%>%
  pivot_wider(id_cols=c(domain,ORCI),names_from=Method,
              values_from=Method)%>%
  mutate("DC-SS"=ifelse(domain%in% c("SOD1:238186:238186_0",
                                       "TARDBP:-:-_2"),"DC-SS",NA))%>%
  mutate("DC-FL"=ifelse(domain%in% c("SOD1:238186:238186_0",
                                       "TARDBP:-:-_2"),"DC-FL",NA))%>%
  unite(col="Methods","DYNATE-SS","DYNATE-FL",
        "SCANG-S","SCANG-B","SCANG-O","DC-SS","DC-FL",sep=",",na.rm =TRUE)
write.csv(see,"Results/Data-Domain.csv")


see2=Domain%>%dplyr::select(-c("benign","possibly",
                              "probably","unknown",
                              "Polyphen.div","Measure","OR"))%>%
  pivot_wider(id_cols=c(domain,ORCI),names_from=Method,
              values_from=Method)%>%
  mutate("DC-SS"=ifelse(domain%in% c("SOD1:238186:238186_0",
                                       "TARDBP:-:-_2"),"DC-SS",NA))%>%
  mutate("DC-FL"=ifelse(domain%in% c("SOD1:238186:238186_0",
                                       "TARDBP:-:-_2"),"DC-FL",NA))%>%
  unite(col="Methods","DYNATE-SS","DYNATE-FL",
        "DC-SS","DC-FL",sep=",",na.rm =TRUE)%>%
  filter(Methods!="")%>%
  dplyr::select(c("domain","ORCI","Methods"))%>%
  rename(Domain=domain,`OR (CI)`=ORCI)

print(xtable(see2),include.rownames=FALSE)


print(xtable(see),include.rownames=FALSE)


### Plot Domains
domain_info <- snp_dat3%>%dplyr::select(c("snpID","loc","domain","chr"))%>%distinct_all()%>%
  mutate(chr=factor(chr,levels=paste0("chr",1:23)))%>%
  arrange(chr,loc)%>%mutate(snp_name2=seq(n()),loc_adj=loc)%>%
  group_by(domain,chr)%>%
  mutate(group=cumsum(c(1L, diff(snp_name2)) > 1L))%>%
  group_by(group,.add=TRUE)%>%
  summarise(start=min(snp_name2),end=max(snp_name2))%>%
  right_join(see,by="domain")%>%
  mutate(OR=strsplit(ORCI,split="[()]")[[1]][1])%>%
  mutate(OR=log(as.numeric(OR)))%>%ungroup()

sedc=domain_info%>%
  filter(domain%in%c("SOD1:238186:238186_0","TARDBP:-:-_2"))%>%
  mutate(Region=domain)%>%
  dplyr::select(c("Region","group","start","end","chr"))
sedc <- rbind(sedc,sedc)
sedc$Method=rep(c("DC-FL","DC-SS"),each=2)

Se <- satet_se%>%bind_rows(scang_se,sedc)%>%ungroup()%>%
  mutate(#range=ifelse(start==2120,"1",ifelse(start<=270000,"2","3")),
         yaxis=rleid(Method))
domain_info2 <- domain_info %>% dplyr::select(c("domain","start","end","OR","chr"))%>%
  distinct_all()%>%
  mutate(#range=ifelse(start<=3000,"1",ifelse(start<=270000,"2","3")),
         Data="Domain log(OR)",
         yaxis=OR)
  
Se <- Se %>%mutate(Data="Region") %>%bind_rows(domain_info2)%>%
  mutate(Data=factor(Data,levels=c("Region","Domain log(OR)")),
         domain2=ifelse(OR>0,domain,"Protective Domains"))
namdc <- names(table(Se$domain2))
namdc <- c(namdc[namdc!="Protective Domains"],"Protective Domains")
Se <- Se %>% mutate(domain2=factor(domain2,levels=namdc))%>%filter(!is.na(start))%>%
  mutate(chr=factor(chr,levels=paste0("chr",c(1,5,18,21))))
  

Se0 <- Se%>%filter(Data=="Region")

col.value <- c("springgreen2",
               #"slateblue",
               "sandybrown",
               #"olivedrab2",
               "lightpink1","brown2","grey69")

# addlim <- data.frame(chr=c("chr21","chr1","chr18","chr5"),start=c(248420,1940,22200,73850),
#                      end=c(248470,1990,22290,74090),Data=rep("Domain log(OR)",2),
#                      OR=NA,domain2=NA,yaxis=1)
# Se <- Se %>% bind_rows(addlim)
# Se0 <- Se%>%filter(Data=="Region")



p=ggplot(Se,aes(x=start,y=yaxis))+
  facet_grid(cols = vars(chr),scales = "free",space="free_x")+
  #facet_grid_sc(cols = vars(chr),scales = list(x = scales_x),space="free_x")+
  geom_rect(data=Se%>%filter(Data=="Domain log(OR)"),aes(ymin=0,xmin=start,xmax=end,
                                                         ymax=yaxis,color=domain2,fill=domain2),
            ,alpha=0.8)+
  scale_fill_manual(values=col.value)+
  scale_color_manual(values=col.value)+
  geom_linerange(data=Se%>%filter(Data=="Region"),aes(xmax=end,xmin=start,group=Region,x=start,y=yaxis),
                 position = position_dodge(0.5),color="dodgerblue2")+
  geom_hline(yintercept=0,linetype="dashed", color = "black")+
  scale_y_continuous(name="Method",
                     breaks = c(1:7),labels =  unique(Se0$Method),
                     sec.axis = sec_axis(~.*1, name="log(OR)"))+
  scale_x_continuous(breaks = seq(25, 300000, 50))+
  xlab("")+
  #xlab("RV Location ID")+
  theme_classic()+labs(fill="Domain",color="Domain")+
  theme(axis.text.x = element_blank(),
       strip.text.x = element_blank())
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  #      strip.text.x = element_text(angle = 45))



pdf("Results/Landscape.pdf",height=3,width=9)
print(p)
dev.off()



######## EPG5

see <- out_satet%>%left_join(snp_dat3,by="snpID")%>%
  filter(domain=="EPG5:-:-_0")%>%
  filter(!is.na(Rs.Number))

unique(see$Rs.Number)


see2 <- out_satet%>%dplyr::select(c("snpID"))%>%
  distinct_all()%>%left_join(snp_dat3,by="snpID")%>%
  filter(domain=="EPG5:-:-_0")%>%
  dplyr::select(c("snpID","domainID","Variant.ID",
                  "Sample.Name","loc_adj","Variant.Type","Genotype","Sample.Type"))%>%
  ungroup()%>%
  group_by(snpID,domainID,Variant.ID,
           Sample.Name,loc_adj,Variant.Type,Genotype,Sample.Type)%>%
  mutate(count=n())

table(see2$count)

see=see2%>%ungroup()%>%group_by(Variant.ID,Sample.Type)%>%
  summarise(number=n())

case=see%>%filter(Sample.Type=="case")
ctrl=see%>%filter(Sample.Type=="ctrl")
hist(case$number,probability=TRUE)
hist(ctrl$number,probability=TRUE)
# freq plot for pts with EGP5 (bar plot)
# freq plot for pts without EPG5  (0 vs 1,2,...), exclude SOD1 TRDP5?


samp <- snp_dat3%>%
  filter(domain%in%c("SOD1:238186:238186_0","TARDBP:-:-_2"))%>%
  dplyr::select(c("Sample.Name"))%>%distinct_all()
samp <- as.character(samp$Sample.Name)

data2 <- out_satet%>%
  filter(Method=="DYNATE-FL")%>%
  dplyr::select(c("snpID","Region"))%>%
  distinct_all()%>%full_join(snp_dat3,by="snpID")%>%
  filter(!Sample.Name%in%samp)%>%
  dplyr::select(c("snpID","domain","Variant.ID",
                  "Sample.Name","loc_adj","Variant.Type",
                  "Genotype","Sample.Type","Region"))%>%
  ungroup()%>%
  mutate(epg5=ifelse(domain=="EPG5:-:-_0"&!is.na(Region),1,0))%>%
  group_by(Sample.Name,Sample.Type)%>%
  summarize(count=sum(epg5))%>%ungroup()%>%
  group_by(Sample.Type,count)%>%
  summarize(Sample.count=n())%>%
  ungroup()%>%group_by(Sample.Type)%>%
  mutate(Total.count=sum(Sample.count))%>%
  mutate(Sample.prop=Sample.count/Total.count)%>%
  ungroup()%>%add_row(Sample.Type="case",count=2,Sample.count=0,Sample.prop=0)



#### only consider samples with mutations in detected region

data1 <- out_satet%>%dplyr::select(c("snpID"))%>%
  distinct_all()%>%left_join(snp_dat3,by="snpID")%>%
  #filter(!domain%in%c("SOD1:238186:238186_0","TARDBP:−:−_2"))%>%
  dplyr::select(c("snpID","domain","Variant.ID",
                  "Sample.Name","loc_adj","Variant.Type",
                  "Genotype","Sample.Type"))%>%
  ungroup()%>%
  mutate(epg5=ifelse(domain=="EPG5:-:-_0",1,0))%>%
  group_by(Sample.Name,Sample.Type)%>%
  summarize(count=sum(epg5))%>%ungroup()%>%
  group_by(Sample.Type,count)%>%
  summarize(Sample.count=n())%>%
  ungroup()%>%group_by(Sample.Type)%>%
  mutate(Total.count=sum(Sample.count))%>%
  mutate(Sample.prop=Sample.count/Total.count)%>%
  ungroup()%>%add_row(Sample.Type="case",count=2,Sample.count=0,Sample.prop=0)


dtab1 <- data1 %>% pivot_wider(names_from=count,
                               values_from =Sample.count,
                               values_fill = 0)%>%
  dplyr::select(-c("Total.count","Sample.prop"))%>%
  group_by(Sample.Type)%>%
  summarize(`0`=sum(`0`),`1`=sum(`1`),`2`=sum(`2`))%>%
  ungroup()%>%
  column_to_rownames("Sample.Type")


data2 <- out_satet%>%dplyr::select(c("snpID"))%>%
  distinct_all()%>%left_join(snp_dat3,by="snpID")%>%
  filter(!Sample.Name%in%samp)%>%
  dplyr::select(c("snpID","domain","Variant.ID",
                  "Sample.Name","loc_adj","Variant.Type","Genotype","Sample.Type"))%>%
  ungroup()%>%
  mutate(epg5=ifelse(domain=="EPG5:-:-_0",1,0))%>%
  group_by(Sample.Name,Sample.Type)%>%
  summarize(count=sum(epg5))%>%ungroup()%>%
  group_by(Sample.Type,count)%>%
  summarize(Sample.count=n())%>%
  ungroup()%>%group_by(Sample.Type)%>%
  mutate(Total.count=sum(Sample.count))%>%
  mutate(Sample.prop=Sample.count/Total.count)%>%
  ungroup()%>%add_row(Sample.Type="case",count=2,Sample.count=0,Sample.prop=0)

dtab2 <- data2 %>% pivot_wider(names_from=count,
                               values_from =Sample.count,
                               values_fill = 0)%>%
  dplyr::select(-c("Total.count","Sample.prop"))%>%
  group_by(Sample.Type)%>%
  summarize(`0`=sum(`0`),`1`=sum(`1`),`2`=sum(`2`))%>%
  ungroup()%>%
  column_to_rownames("Sample.Type")

chisq.test(dtab2,simulate.p.value=TRUE)#0.0010
chisq.test(dtab1,simulate.p.value = TRUE)#0.0005

library(ggpubr)
p2<-ggplot(data=data2, aes(x=count, y=Sample.prop,fill=Sample.Type)) +
  geom_bar(stat="identity", position=position_dodge())+
  ggtitle("No SOD1/TARDBP (p-value=0.0015)")+
  xlab("Counts of EPG5 mutations")+ylab("Sample Proportions")
p1<-ggplot(data=data1, aes(x=count, y=Sample.prop,fill=Sample.Type)) +
  geom_bar(stat="identity", position=position_dodge())+
  ggtitle("With SOD1/TARDBP (p-value=0.0005)")+
  xlab("Counts of EPG5 mutations")+ylab("Sample Proportions")

pdf("Results/EGP5_investigate.pdf",height=4,width=9)
print(ggarrange(p1,p2))
dev.off()

#Those names are in BED format as described in \cite{gussow2016intolerance}
see.domain <- out_satet%>%dplyr::select(c("snpID"))%>%
  distinct_all()%>%left_join(snp_dat3,by="snpID")%>%
  filter(domain=="EPG5:-:-_0")%>%
  dplyr::select(c("snpID","domainID","Variant.ID",
                  "Sample.Name","loc_adj","Variant.Type","Genotype","Sample.Type"))%>%
  ungroup()%>%
  group_by(snpID,domainID,Variant.ID,
           Sample.Name,loc_adj,Variant.Type,Genotype,Sample.Type)%>%
  mutate(count=n())





