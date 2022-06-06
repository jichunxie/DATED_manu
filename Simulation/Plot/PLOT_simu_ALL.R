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


# pis=0.6;upq=0.95;lowq=0.05;tvi=9;c0=39
Simu_analysis <- function(pis,upq,lowq,L_max=3){
  load(paste0("sim_map",pis*10,".RData"))
  load(paste0("Output_SATET_ROC",pis*10,".RData"))
  load(paste0("Output_SCANG_ROC2_",pis*10,".RData"))
  
  #### SATET
  Out_SATET <- out %>% mutate(alpha=as.integer(alpha*100)/100,
                              DC_F1=2/(1/(1-DC_FDP)+1/DC_power),
                              SNP_F1=2/(1/(1-SNP_FDP)+1/SNP_power))%>%
    filter(alpha==0.05) %>% rename(length=length_med)%>%
    pivot_longer(cols=c("DC_FDP","DC_power","SNP_FDP","SNP_power",
                        "Node_FDP_new","length","DC_F1","SNP_F1"))%>%
    mutate("Method"=ifelse(Test=="FET",
                           ifelse(Layer==1,"LC (FL)","DATED-FL"),
                           ifelse(Layer==1,"LC (SS)","DATED-SS")))%>%
    group_by(Method,SS.Factor,Covar,Layer,name) %>%
    summarise(Upper=quantile(value,upq,na.rm=TRUE),
              Lower=quantile(value,lowq,na.rm=TRUE),
              value=mean(value))%>%
    ungroup()%>%mutate(Pi=pis)

  Out_SATET0 <- Out_SATET%>%filter(Layer%in% c(1,L_max),name!="length")
  Out_SATET_len <- Out_SATET%>%
    filter(Layer%in% c(L_max),name=="length")%>%
    dplyr::select(-("Layer"))
    
  Out_SATET <- Out_SATET%>%filter(!(name%in%c("length")))
  #### SCANG
  Out_sc <- outdc%>%bind_rows(outdc_cov)%>%
    filter(!(index%in%c("Node_FDP","length")),alpha==0.05)%>%
    rename(value=out,name=index,Method=Test,
           Covar=covar,SS.Factor=SS_Factor)%>%
    #distinct_all()
    #group_by(alpha,SS.Factor,Method,Covar) %>%
    #mutate()
    #pivot_wider(names_from=name,values_from=value)%>%
    group_by(SS.Factor,Method,Covar,name)%>%
    summarise(Upper=quantile(value,upq,na.rm=TRUE),
              Lower=quantile(value,lowq,na.rm=TRUE),
              value=mean(value))%>%mutate(Pi=pis)%>%ungroup()%>%
    mutate(Method=ifelse(Method=="SCANGB","SCANG-B",
                         ifelse(Method=="SCANGO","SCANG-O","SCANG-S")))
  Out_sc_len <- outdc%>%bind_rows(outdc_cov)%>%
    filter(index=="length",alpha==0.05)%>%
    rename(value=out,name=index,Method=Test,
           Covar=covar,SS.Factor=SS_Factor)%>%
    group_by(SS.Factor,Method,Covar,name)%>%
    summarise(Upper=quantile(value,upq,na.rm=TRUE),
              Lower=quantile(value,lowq,na.rm=TRUE),
              value=mean(value))%>%mutate(Pi=pis)%>%ungroup()%>%
    mutate(Method=ifelse(Method=="SCANGB","SCANG-B",
                         ifelse(Method=="SCANGO","SCANG-O","SCANG-S")))
  
  
  #### DC
  load(paste0("Output_DC",pis*10,".RData"))
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
  Padj <- Padj%>%bind_rows(Padj_cov)
  out=Padj%>%filter(padj<=0.05)%>%group_by(Test,SS.Factor,Covar,Seed)%>%
      summarize("DC_FDP"=uniqueN(intersect(domainID,nulls_dc))/uniqueN(domainID),
                "DC_power"=uniqueN(intersect(domainID,nonnulls_dc))/length(nonnulls_dc),
                "SNP_FDP"=uniqueN(intersect(snpID,nulls))/uniqueN(snpID),
                "SNP_power"=uniqueN(intersect(snpID,nonnulls))/length(nonnulls),
                .groups = 'drop')%>%ungroup()
  
  Out_DC <- out%>%mutate(DC_F1=2/(1/(1-DC_FDP)+1/DC_power),
                         SNP_F1=2/(1/(1-SNP_FDP)+1/SNP_power))%>%
    pivot_longer(cols=c("DC_FDP","DC_power","DC_F1","SNP_FDP","SNP_power","SNP_F1"))%>%
    mutate(Method=ifelse(Test=="FET","DC-FL","DC-SS"))%>%
    group_by(SS.Factor,Method,Covar,name)%>%
    summarise(Upper=quantile(value,upq,na.rm=TRUE),
              Lower=quantile(value,lowq,na.rm=TRUE),
              value=mean(value))%>%
    ungroup()%>%mutate(Pi=pis)
  
  #### Method comparison
  Out=Out_SATET0%>%bind_rows(Out_DC)%>%bind_rows(Out_sc)%>%
    filter(!(name%in%c("Node_FDP_new","DC_FPR","SNP_FPR")))%>%
    mutate(Measure=ifelse(name=="DC_FDP","Avg Domain-FDP",
                          ifelse(name=="DC_power","Avg Domain-Sensitivity",
                                 ifelse(name=="SNP_FDP","Avg RV-FDP",
                                        ifelse(name=="DC_F1","Avg Domain-F1 score",
                                               ifelse(name=="SNP_F1","Avg RV-F1 score","Avg RV-Sensitivity"))))))
  
  #### Layer-wise comparison
  Out_layer=Out_SATET%>%filter(Layer<=L_max)%>%
    mutate(Layer=as.character(Layer),
           Method=ifelse(Method%in%c("DATED-FL","LC (FL)"),"FL",
                         "SS"),
           Measure=ifelse(name=="Node_FDP_new","Avg Node-FDP",
                          ifelse(name=="DC_FDP","Avg Domain-FDP",
                                 ifelse(name=="DC_power","Avg Domain-Sensitivity",
                                        ifelse(name=="SNP_FDP","Avg RV-FDP",
                                               ifelse(name=="SNP_F1","Avg RV-F1 score",
                                                      ifelse(name=="DC_F1","Avg Domain-F1 score","Avg RV-Sensitivity")))))))
  
  
  return(list("out1"=Out,"out2"=Out_layer,"len"=rbind(Out_sc_len,Out_SATET_len)))
}


Out2 <- Simu_analysis(pis=0.2,upq=0.95,lowq=0.05)
Out4 <- Simu_analysis(pis=0.4,upq=0.95,lowq=0.05)
Out6 <- Simu_analysis(pis=0.6,upq=0.95,lowq=0.05)
Out8 <- Simu_analysis(pis=0.8,upq=0.95,lowq=0.05)

Length <- Out2$len%>%bind_rows(Out4$len,Out6$len,Out8$len)%>%
  mutate(covar=ifelse(Covar,"With Covariates","Without Covariates"),
         Pi=paste0("pi=",Pi),
         Measure="Median Region-Length",
         SS.Factor=as.factor(SS.Factor*1000))%>%
  filter(SS.Factor!=4000)%>%
  mutate(Method=factor(Method,levels=c("DATED-FL",
                                       "DATED-SS",
                                       "SCANG-B",     
                                       "SCANG-O",     
                                       "SCANG-S")))

col.value <- c("red","purple","goldenrod2","darkorange2","khaki3")

Plot_length=ggplot(Length%>%filter(SS.Factor!=4000), 
                   aes(x = SS.Factor, y = value,fill=Method))+
  facet_nested(covar+Pi~Measure)+
  scale_fill_manual(values=col.value)+
  scale_y_continuous(breaks=c(5,25))+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  xlab("N0")+ylab("")+
  theme(strip.text.y = element_text(size=8, angle=270),
        strip.text.x = element_text(size=8),
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))



  

Out <-  Out2$out1 %>% bind_rows(Out4$out1,Out6$out1,Out8$out1)%>%
  filter(SS.Factor!=4000)%>%
  mutate(covar=ifelse(Covar,"With Covariates","Without Covariates"),
         Pi=paste0("pi=",Pi),
         SS.Factor=as.factor(SS.Factor*1000))%>%
  mutate(Measure=factor(Measure,levels=c("Avg RV-FDP","Avg Domain-FDP",
                                         "Avg RV-Sensitivity",
                                         "Avg Domain-Sensitivity",
                                         "Avg RV-F1 score","Avg Domain-F1 score")))%>%
  #mutate(Method=ifelse(Method%in%c("SATET-FL","SATET-SS"),paste0("3-layer ",Method),
  #                     Method))%>%
  mutate(Method=factor(Method,levels=c("LC (FL)",
                                        "LC (SS)",
                                        "DATED-FL",
                                        "DATED-SS",
                                         "DC-FL",
                                         "DC-SS",
                                         "SCANG-B",     
                                         "SCANG-O",     
                                         "SCANG-S")))
Out_layer <- Out2$out2 %>% bind_rows(Out4$out2,Out6$out2,Out8$out2)%>%
  mutate(covar=ifelse(Covar,"With Covariates","Without Covariates"),
         Pi=paste0("pi=",Pi),
         SS.Factor=as.factor(SS.Factor*1000))%>%
  filter(SS.Factor!=4000)

#%>% filter(Layer!="4")

Out_layer2 <- Out_layer%>%filter(Measure%in%c("Avg Node-FDP","Avg RV-FDP",
                                              "Avg RV-Sensitivity","Avg RV-F1 score"))%>%
  mutate(Method=paste0(Layer,"-layer"," DATED","-",Method),
         Measure=factor(Measure,levels=c("Avg Node-FDP","Avg RV-FDP",
                                         "Avg RV-Sensitivity","Avg RV-F1 score")))%>%
  filter(Measure %in% c("Avg Node-FDP","Avg RV-Sensitivity"))


Outl2 <- Out_layer2%>%filter(Measure%in%c("Avg Node-FDP","Avg RV-FDP"))%>%
  dplyr::select(-c("Method","value","Upper","Lower","Layer"))%>%
  distinct_all()%>%
  mutate(FDP=ifelse(Measure=="Avg Node-FDP",0.05,1-as.numeric(substr(Pi,4,6))))

col.value <- c("seagreen4","seagreen3","slategray2","slategray3",
               "grey69","grey56","dodgerblue1","dodgerblue4",
               "cornsilk","cornsilk3","red","purple")
satet.plot <- ggplot(Out_layer2,
                     aes(x=SS.Factor,y=value,fill=Method))+
  facet_nested(covar+Pi~Measure)+
  scale_fill_manual(values=col.value)+
  scale_y_continuous(breaks=c(0.25,0.75))+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  geom_hline(data=Outl2,aes(yintercept=FDP),linetype="dashed")+
  xlab("N0")+ylab("")+
  theme(strip.text.y = element_text(size=8, angle=270),
        strip.text.x = element_text(size=8),
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))

  
satet.plot_S <- ggplot(Out_layer2%>%filter(Pi%in%c("pi=0.4","pi=0.8")),
                     aes(x=SS.Factor,y=value,fill=Method))+
  facet_nested(covar+Pi~Measure)+
  scale_fill_manual(values=col.value)+
  scale_y_continuous(breaks=c(0.25,0.75))+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  geom_hline(data=Outl2%>%filter(Pi%in%c("pi=0.4","pi=0.8")),aes(yintercept=FDP),linetype="dashed")+
  xlab("N0")+ylab("")+
  theme(strip.text.y = element_text(size=8, angle=270),
        strip.text.x = element_text(size=8),
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))




shape.value <- c(4,8,10,12,19,17,1,2,5)
col.value <- c("seagreen4","seagreen3","skyblue3","skyblue2","red","purple",
               "goldenrod2","darkorange2","khaki3")


Out <- Out %>% filter(!(Method%in%c("LC (FL)","LC (SS)")),
                      SS.Factor!=4000)
col.value <- c("red","purple","skyblue3","skyblue2",
               "goldenrod2","darkorange2","khaki3")


Outl <- Out%>%filter(Measure%in%c("Avg GS-FDP","Avg RV-FDP"))%>%
  dplyr::select(-c("Method","value","Upper","Lower"))%>%
  distinct_all()%>%
  mutate(FDP=1-as.numeric(substr(Pi,4,6)))

Plot_compare_RV=ggplot(Out%>%
                         filter(Measure%in%c("Avg RV-FDP",
                                             "Avg RV-Sensitivity",
                                             "Avg RV-F1 score"),
                                substr(Method,1,2)!="DC"), 
                       aes(x = SS.Factor, y = value,fill=Method))+
  facet_nested(covar+Pi~Measure)+
  scale_fill_manual(values=col.value[-c(3,4)])+
  scale_y_continuous(breaks=c(0.25,0.75))+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  geom_hline(data=Outl,aes(yintercept=FDP),linetype="dashed")+
  xlab("N0")+ylab("")+
  theme(strip.text.y = element_text(size=8, angle=270),
        strip.text.x = element_text(size=8),
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))

Plot_compare_RV_S=ggplot(Out%>%
                         filter(Measure%in%c("Avg RV-FDP",
                                             "Avg RV-Sensitivity",
                                             "Avg RV-F1 score"),
                                substr(Method,1,2)!="DC",
                                Pi%in%c("pi=0.4","pi=0.8")), 
                       aes(x = SS.Factor, y = value,fill=Method))+
  facet_nested(covar+Pi~Measure)+
  scale_fill_manual(values=col.value[-c(3,4)])+
  scale_y_continuous(breaks=c(0.25,0.75))+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  geom_hline(data=Outl%>%filter(Pi%in%c("pi=0.4","pi=0.8")),aes(yintercept=FDP),linetype="dashed")+
  xlab("N0")+ylab("")+
  theme(strip.text.y = element_text(size=8, angle=270),
        strip.text.x = element_text(size=8),
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))


Plot_compare_dc=ggplot(Out%>%filter(!(Measure%in%c("Avg RV-FDP","Avg RV-Sensitivity",
                                                   "Avg RV-F1 score"))), 
                       aes(x = SS.Factor, y = value,fill=Method))+
  facet_nested(covar+Pi~Measure)+
  scale_fill_manual(values=col.value)+
  scale_y_continuous(breaks=c(0.25,0.75))+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  xlab("N0")+ylab("")+
  theme(strip.text.y = element_text(size=8, angle=270),
        strip.text.x = element_text(size=8),
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))

Plot_compare_dc_S=ggplot(Out%>%filter(!(Measure%in%c("Avg RV-FDP","Avg RV-Sensitivity","Avg RV-F1 score")),
                                    Pi%in%c("pi=0.4","pi=0.8")), 
                       aes(x = SS.Factor, y = value,fill=Method))+
  facet_nested(covar+Pi~Measure)+
  scale_fill_manual(values=col.value)+
  scale_y_continuous(breaks=c(0.25,0.75))+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  xlab("N0")+ylab("")+
  theme(strip.text.y = element_text(size=8, angle=180),
        strip.text.x = element_text(size=8),
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))


pdf.options(reset = TRUE, onefile = FALSE)
pdf("Comparison_ALL_dc.pdf",height=4,width=7.5)
print(Plot_compare_dc)
dev.off()

pdf.options(reset = TRUE, onefile = FALSE)
pdf("Comparison_ALL_RV.pdf",height=4,width=7.5)
print(Plot_compare_RV)
dev.off()


pdf.options(reset = TRUE, onefile = FALSE)
pdf("Comparison_ALL_dc_S.pdf",height=4,width=7.5)
print(Plot_compare_dc_S)
dev.off()

pdf.options(reset = TRUE, onefile = FALSE)
pdf("Comparison_ALL_RV_S.pdf",height=4,width=7.5)
print(Plot_compare_RV_S)
dev.off()

pdf.options(reset = TRUE, onefile = FALSE)
pdf("Comparison_ALL_length.pdf",height=4,width=5)
print(Plot_length)
dev.off()



pdf.options(reset = TRUE, onefile = FALSE)
pdf("SATET_layerwise_ALL.pdf",height=4,width=7)
print(satet.plot)
dev.off()

pdf.options(reset = TRUE, onefile = FALSE)
pdf("SATET_layerwise_ALL_S.pdf",height=4,width=7)
print(satet.plot_S)
dev.off()


