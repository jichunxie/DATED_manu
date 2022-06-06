rm(list=ls())
#setwd("/Users/xuechan/Documents/Study/Rotation-Jichun/SATET/Results")
library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr) 
library(foreach)
library(doParallel)
library(ggh4x)
source("../R_Func/DATED_Func.R")
savedatafile="/hpc/home/xl110/SATET/Results/"
loadpop="/hpc/group/chsi/xl110/SATET/"

# pis=0.6;upq=0.95;lowq=0.05;L=3;c0=39
#Others
Simu_analysis <- function(pis,upq,lowq,c0=39){
  # c0=39 is from the information of the simulated genotype information
  load(paste0("sim_map", pis * 10, ".RData"))
  
  #### DATED multithresh
  load(paste0("Output_SATET_multilayer", pis * 10, ".RData"))
  out1 <- out
  load(paste0("Output_SATET_2multilayer", pis * 10, ".RData"))
  out <- out%>%bind_rows(out1)#%>%dplyr::select(-c("L5","L6"))
  
  
  load(paste0("Output_SCANG_ROC2_",pis*10,".RData"))
  
  m = sim_map$num_snps
  nonnulls = sim_map$alt_snps
  nonnulls_omni = sim_map$alt_sparse_snps
  nulls = sim_map$null_snps
  m_dc = length(sim_map$domain_snp_list)
  nonnulls_dc = unlist(sim_map$domain_w_sa_alt_snps)
  nulls_dc = setdiff(seq(m_dc), nonnulls_dc)
  
  smap_res <- out %>% mutate(Lmax=floor(log(c0/thresh_val,2))+1) %>%filter(thresh_val>=5,
                                                                           thresh_val<=11)
  Lres <- max(smap_res$Lmax)
  # Multilayer analysis
  Res0 <- smap_res %>% dplyr::select(-c("snpID", "alt_snps")) %>%
    filter(Layer<=Lmax)%>%distinct_all() %>%
    pivot_longer(cols = paste0("L", seq(Lres)), names_to = "name") %>%
    filter(name == paste0("L", Layer)) %>%
    group_by(Test, Layer, value, seed,thresh_val) %>% # value, the DGS index at specific layer
    mutate(FD_node = as.integer(sum(alt_leaf) == 0),
           FD_node_new = 1 - sum(alt_leaf) / n(),
           TP_node = as.integer(sum(alt_leaf) != 0),
           TP_node_new = sum(alt_leaf) / n() ) %>% ungroup() %>%
    dplyr::select(
      c("Test","Layer","value","FD_node","FD_node_new","TP_node",
        "TP_node_new","thresh_val","Covar","SS.Factor","pi1",
        "alt_leaf","Covar","seed")
    ) %>%
    distinct_all()
  
  Res1 <-
    smap_res %>% filter(Layer <= Lmax) %>% group_by(Test, SS.Factor, 
                                                    Covar, seed,thresh_val) %>%
    summarize(
      "DC_FDP" = uniqueN(intersect(domainID, nulls_dc)) / uniqueN(domainID),
      "DC_power" = uniqueN(intersect(domainID, nonnulls_dc)) /
        length(nonnulls_dc),
      "DC_FPR" = uniqueN(intersect(domainID, nulls_dc)) / length(nulls_dc),
      "SNP_FDP" = uniqueN(intersect(snpID, nulls)) / uniqueN(snpID),
      "SNP_power" = uniqueN(intersect(snpID, nonnulls)) / length(nonnulls),
      #true positive rate
      "SNP_FPR" = uniqueN(intersect(snpID, nulls)) / length(nulls),
      .groups = 'drop'
    ) %>% ungroup() %>% mutate(DC_F1=2/(1/(1-DC_FDP)+1/DC_power),
                               SNP_F1=2/(1/(1-SNP_FDP)+1/SNP_power))
  Res2 <-
    Res0 %>% group_by(Test, SS.Factor, Covar, seed,thresh_val) %>%
    summarize(
      "Node_FDP" = mean(FD_node),
      "Node_FDP_new" = mean(FD_node_new),
      "Leaf_FDP" = 1 - mean(alt_leaf),
      .groups = 'drop'
    ) %>%
    ungroup() %>%
    left_join(Res1, by = c("Test", "SS.Factor", "Covar", "seed","thresh_val"))
  
  out_DATED_mt <-
    Res2 %>% mutate(alpha = 0.05) %>% 
    pivot_longer(cols = c("DC_FDP", "DC_power", "DC_F1",
                          "SNP_FDP", "SNP_power","SNP_F1","Node_FDP_new")) %>%
    mutate("Method" = ifelse(
      Test == "FET",
      paste0("DATED-FL (",thresh_val,")"),
      paste0("DATED-SS (",thresh_val,")")
    )) %>%
    group_by(Method, SS.Factor, Covar, name,thresh_val) %>%
    summarise(
      Upper = quantile(value, upq, na.rm = TRUE),
      Lower = quantile(value, lowq, na.rm = TRUE),
      value = mean(value)
    ) %>%
    ungroup() %>% mutate(Pi = pis)
  
  #### SCANG
  # Out_sc <- outdc%>%bind_rows(outdc_cov)%>%
  #   filter(index!="Node_FDP",alpha==0.05)%>%
  #   rename(value=out,name=index,Method=Test,
  #          Covar=covar,SS.Factor=SS_Factor)%>%
  #   group_by(SS.Factor,Method,Covar,name)%>%
  #   summarise(Upper=quantile(value,upq,na.rm=TRUE),
  #             Lower=quantile(value,lowq,na.rm=TRUE),
  #             value=mean(value))%>%mutate(Pi=pis)%>%ungroup()%>%
  #   mutate(Method=ifelse(Method=="SCANGB","SCANG-B",
  #                        ifelse(Method=="SCANGO","SCANG-O","SCANG-S")))
  # 
  # #### DC
  # load(paste0("Output_DC",pis*10,".RData"))
  # nonnulls = sim_map$alt_snps
  # nulls = sim_map$null_snps
  # m_dc = length(sim_map$domain_snp_list)
  # nonnulls_dc = unlist(sim_map$domain_w_sa_alt_snps)
  # nulls_dc = setdiff(seq(m_dc),nonnulls_dc)
  # 
  # out <- out_cov <- NULL
  # 
  # Padj=output%>%dplyr::select(c("pvals","domainID","Test","SS.Factor","Seed","Covar"))%>%
  #   distinct_all()%>%group_by(Test,SS.Factor,Seed,Covar)%>%
  #   mutate(padj=p.adjust(pvals,method="BH"))%>%ungroup()%>%
  #   left_join(output,by=c("domainID","Test","SS.Factor","Seed","Covar"))
  # Padj_cov=output_cov%>%dplyr::select(c("pvals","domainID","Test","SS.Factor","Seed","Covar"))%>%
  #   distinct_all()%>%group_by(Test,SS.Factor,Seed,Covar)%>%
  #   mutate(padj=p.adjust(pvals,method="BH"))%>%ungroup()%>%
  #   left_join(output_cov,by=c("domainID","Test","SS.Factor","Seed","Covar"))
  # Padj <- Padj%>%bind_rows(Padj_cov)
  # out=Padj%>%filter(padj<=0.05)%>%group_by(Test,SS.Factor,Covar,Seed)%>%
  #   summarize("DC_FDP"=uniqueN(intersect(domainID,nulls_dc))/uniqueN(domainID),
  #             "DC_power"=uniqueN(intersect(domainID,nonnulls_dc))/length(nonnulls_dc),
  #             "SNP_FDP"=uniqueN(intersect(snpID,nulls))/uniqueN(snpID),
  #             "SNP_power"=uniqueN(intersect(snpID,nonnulls))/length(nonnulls),
  #             .groups = 'drop')%>%ungroup()
  # 
  # Out_DC <- out%>%pivot_longer(cols=c("DC_FDP","DC_power","SNP_FDP","SNP_power"))%>%
  #   mutate(Method=ifelse(Test=="FET","DC-FL","DC-SS"))%>%
  #   group_by(SS.Factor,Method,Covar,name)%>%
  #   summarise(Upper=quantile(value,upq,na.rm=TRUE),
  #             Lower=quantile(value,lowq,na.rm=TRUE),
  #             value=mean(value))%>%
  #   ungroup()%>%mutate(Pi=pis)
  # 
  #### Method comparison
  Out=out_DATED_mt%>%
    filter(!(name%in%c("DC_FPR","SNP_FPR")))%>%
    mutate(Measure=ifelse(name=="DC_FDP","Avg Domain-FDP",
                          ifelse(name=="DC_power","Avg Domain-Sensitivity",
                                 ifelse(name=="SNP_FDP","Avg RV-FDP",
                                        ifelse(name=="DC_F1","Avg Domain-F1 score",
                                               ifelse(name=="SNP_F1","Avg RV-F1 score",
                                                      ifelse(name=="SNP_power","Avg RV-Sensitivity",
                                                             "Avg Node-FDP")))))))
  
  
  # Out=Out_DATED%>%bind_rows(out_DATED_mt)%>%
  #   bind_rows(Out_sc,Out_DC)%>%
  #   filter(!(name%in%c("DC_FPR","SNP_FPR")))%>%
  #   mutate(Measure=ifelse(name=="DC_FDP","Avg Domain-FDP",
  #                         ifelse(name=="DC_power","Avg Domain-Sensitivity",
  #                                ifelse(name=="SNP_FDP","Avg RV-FDP",
  #                                       ifelse(name=="SNP_power","Avg RV-Sensitivity","Avg GS-FDP")))))

  return(list("out1"=Out,"out_mt"=out_DATED_mt))
}

Out2 <- Simu_analysis(pis=0.2,upq=0.95,lowq=0.05)
Out4 <- Simu_analysis(pis=0.4,upq=0.95,lowq=0.05)
Out6 <- Simu_analysis(pis=0.6,upq=0.95,lowq=0.05)
Out8 <- Simu_analysis(pis=0.8,upq=0.95,lowq=0.05)

lev=paste0(c("DATED-FL (","DATED-SS ("),rep(c(5,7,9,11,19,21),each=2),")")

Out <-  Out6$out1%>% bind_rows(Out2$out1,Out4$out1,Out8$out1)%>%
  filter(substr(Method,1,2)=="DA")%>%
  rename(Pis=Pi)%>%
  mutate(covar=ifelse(Covar,"With Covariates","Without Covariates"),
         Pi=paste0("pi=",Pis),SS.Factor=as.factor(SS.Factor*1000))%>%
  mutate(Measure=factor(Measure,levels=c("Avg Node-FDP","Avg RV-FDP","Avg Domain-FDP",
                                         "Avg RV-Sensitivity","Avg RV-F1 score",
                                         "Avg Domain-Sensitivity","Avg Domain-F1 score")))%>%
  mutate(Method=factor(Method,levels=lev))

#col.value <- c("red","purple","seagreen4","seagreen3","darkolivegreen3","darkolivegreen4","grey","black")
Out <- Out%>%filter(Measure %in% c("Avg Node-FDP","Avg RV-FDP","Avg RV-Sensitivity","Avg RV-F1 score"),
                    SS.Factor!=4000)%>%
  filter(Measure %in% c("Avg Node-FDP","Avg RV-FDP","Avg RV-Sensitivity"))

Outl <- Out%>%filter(Measure%in%c("Avg Node-FDP","Avg RV-FDP"))%>%
  dplyr::select(-c("Method","value","Upper","Lower"))%>%
  distinct_all()%>%
  mutate(FDP=ifelse(Measure=="Avg Node-FDP",0.05,1-Pis))

Plot_compare=ggplot(Out, aes(x = SS.Factor, y = value,fill=Method))+
  facet_nested(covar+Pi~Measure)+
  #scale_fill_manual(values=col.value)+
  scale_y_continuous(breaks=c(0.25,0.75))+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  geom_hline(data=Outl,aes(yintercept=FDP),linetype="dashed")+
  xlab("N0")+ylab("")+#theme_classic()+
  theme(strip.text.y = element_text(size=8, angle=270),
        strip.text.x = element_text(size=8),
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))


Plot_compare_S=ggplot(Out%>%filter(Pis%in%c(0.4,0.8)), aes(x = SS.Factor, y = value,fill=Method))+
  facet_nested(covar+Pi~Measure)+
  #scale_fill_manual(values=col.value)+
  scale_y_continuous(breaks=c(0.25,0.75))+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.2,width=0.2,position=position_dodge(.9))+
  geom_hline(data=Outl%>%filter(Pis%in%c(0.4,0.8)),aes(yintercept=FDP,),linetype="dashed")+
  xlab("N0/N1")+ylab("")+#theme_classic()+
  theme(strip.text.y = element_text(size=8, angle=270),
        strip.text.x = element_text(size=8),
        panel.spacing=unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
        strip.background = element_rect(color = "black", size = 0.1))


pdf.options(reset = TRUE, onefile = FALSE)
pdf("Comparison_ALL_MT.pdf",height=4,width=7.5)
print(Plot_compare)
dev.off()

pdf.options(reset = TRUE, onefile = FALSE)
pdf("Comparison_ALL_MT_S.pdf",height=4,width=7.5)
print(Plot_compare_S)
dev.off()

