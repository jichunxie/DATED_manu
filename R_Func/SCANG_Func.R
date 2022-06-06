DGP2_SCANG_subsamp <- function(pop,N0,N1,covars,sim_map,samps){
  #subsample without struct_map
  #Randomly select N1 cases and N0 controls
  #dgTMatrix objects
  sampcase <- samps$sampcase
  sampctrl <- samps$sampctrl
  
  #Get attribute matrices # 5mins
  geno_case_matrix <- pop$Gmat[sampcase,]
  geno_ctrl_matrix <- pop$Gmat[sampctrl,]
  
  X1ctrl <- pop$nulldata[sampctrl,"X1"]
  X1case <- pop$nulldata[sampcase,"X1"]
  X2ctrl <- pop$nulldata[sampctrl,"X2"]
  X2case <- pop$nulldata[sampcase,"X2"]
  
  #Output
  if(covars){
    list(nulldata=data.frame(phenotype=rep(c(1,0),times=c(N1,N0)),
                             X1=c(X1case,X1ctrl),
                             X2=c(X2case,X2ctrl)),#Note order,
         geno_case_matrix=geno_case_matrix,
         geno_ctrl_matrix=geno_ctrl_matrix)
  } else{
    list(nulldata=data.frame(phenotype=rep(c(1,0),times=c(N1,N0))),
         geno_case_matrix=geno_case_matrix,
         geno_ctrl_matrix=geno_ctrl_matrix)
  }
}

Find_L <- function(snp_loc){
  snp_loc$loc3 <- snp_loc$snp_loc+3000
  snp_loc$loc7 <- snp_loc$snp_loc+7000

  st3 <- st7 <- NULL
  for(chrno in unique(snp_loc$chr)){
    snp_chr <- snp_loc%>%filter(chr==chrno)%>%mutate(snp_order_chr=seq(n()))
    
    res3 <- t(outer(snp_chr$snp_loc,snp_chr$loc3,"<="))
    res3[lower.tri(res3, diag = FALSE)] <- FALSE
    res7 <- t(outer(snp_chr$snp_loc,snp_chr$loc7,"<="))
    res7[lower.tri(res7, diag = FALSE)] <- FALSE
    
    st3 <- c(st3,rowSums(res3))
    st7 <- c(st7,rowSums(res7))
  }
  Lmin <- quantile(st3,0.01)      
  Lmax <- quantile(st7,0.99)
  return(c(Lmin,Lmax))
}

Test_by_SCANG <- function(N1,
                          pop,
                          sampsize.factors,
                          sim_map,
                          params,
                          covars=FALSE,
                          filter=2e-5,
                          f=0.5,
                          alpha,
                          seed,
                          SampleID){
  
  set.seed(seed)
  
  ### Generate Sample Data ###
  Sim.mat <- NULL
  out <- NULL
  for(i in seq_along(sampsize.factors)){
    ## Generate case-control samples
    N0 <- ceiling(N1*sampsize.factors[i])
    
    ## Create data
    t1 <- Sys.time()
    attrDat <- DGP2_SCANG_subsamp(pop,N0 = N0,N1 = N1,covars=covars,sim_map=sim_map,
                                  samps=SampleID[[i]])
    t2 <- Sys.time()
    
    t2-t1
    
    ## Phenotype
    phenotypedata <- attrDat$nulldata
    
    ## Genotype
    if(class(attrDat$geno_case_matrix)%in%c("CsparseMatrix","ngCMatrix")){
      mat1 <- unname(as(attrDat$geno_case_matrix,"dgCMatrix"))
      mat0 <- unname(as(attrDat$geno_ctrl_matrix,"dgCMatrix"))
    } else{
      mat1 <- unname(as(attrDat$geno_case_matrix,"dgTMatrix"))
      mat0 <- unname(as(attrDat$geno_ctrl_matrix,"dgTMatrix"))
    }
    
    genotype0 <- as(rbind(mat1,mat0),"dgCMatrix")
    
    index1 <- which(colSums(genotype0)>0)
    
    ans <- map_df(sim_map$domain_snp_list, ~as.data.frame(.x), .id="DomainID")
    colnames(ans) <- c("DomainID","snp_name")
    snp_loc <- data.table(sim_map$snp_loc)%>%left_join(ans,by="snp_name")
    snp_loc <- snp_loc[index1,]
    snp_loc <- snp_loc %>% mutate("loc_order"=order(snp_loc))
    
    genotype <- genotype0[,index1]
    genotype <- genotype[,snp_loc$loc_order]
    
    snp_loc <- snp_loc%>%arrange(snp_loc)%>%
      mutate(snp_name2=seq(n()))
    Ls <- Find_L(snp_loc)
    #see=snp_loc%>%filter(snp_name%in% sim_map$alt_snps)
    #diff(see$snp_name2)

    # need additional algorithm to find Lmin and Lmax
    # Lmax <- max(lengths(sim_map$domain_snp_list))         
    # Lmin <- min(lengths(sim_map$domain_snp_list))       
    # Ls <- quantile(lengths(sim_map$domain_snp_list),range)
    Lmin <- Ls[1]
    Lmax <- Ls[2]
    
    # annotation_phred=NULL
    # rare_maf_cutoff=1
    # steplength=5
    # alpha=0.05
    # filter=1e-4
    # subseq_num=2000
    

    if(covars){
      phenotype <- as.matrix(phenotypedata$phenotype)
      Covs <- as.matrix(phenotypedata%>%dplyr::select(-c(phenotype)))
      t1 <- Sys.time()
      obj_nullmodel <- fit_null_glm_SCANG(phenotype~Covs,data=phenotypedata,family=binomial)
      res_lm <- tryCatch(SCANG(genotype,obj_nullmodel,Lmin,Lmax,filter=filter,f=f,
                      rare_maf_cutoff=1))
      ptime <- Sys.time()-t1
    }else{
      t1 <- Sys.time()
      obj_nullmodel <- fit_null_glm_SCANG(phenotype~1,data=phenotypedata,family=binomial)
      res_lm <- tryCatch(SCANG(genotype,obj_nullmodel,Lmin,Lmax,filter=filter,f=f,
                      rare_maf_cutoff=1),error=function(e) e)
      ptime <- Sys.time()-t1
    }

    if(inherits(res_lm,"error")){
      CV.power.o <- CV.power.s  <- CV.power.b <- NA
      FDP.power.o <- FDP.power.s <- FDP.power.b <- NA
      Domain.o <- Domain.s <- Domain.b <- Prop.snp.o <- Prop.snp.s <- Prop.snp.b <- NA
      lrejo <- lrejs <- lrejb <- avg.seg.len.o <- avg.seg.len.s <- avg.seg.len.b <- NA
      med.seg.len.o <- med.seg.len.s <- med.seg.len.b <- NA
      FDP.node.o <- FDP.node.s <- FDP.node.b <- NA
      FDP.dc.o <- FDP.dc.s <- FDP.dc.b <- NA
    }else{
      reso <- res_lm$SCANG_O_res
      ress <- res_lm$SCANG_S_res
      resb <- res_lm$SCANG_B_res
      alt_snps <- snp_loc$snp_name2[which(snp_loc$snp_name%in%sim_map$alt_snps)]
      
      avg.seg.len.o <- mean(reso[,3]-reso[,2])
      avg.seg.len.s <- mean(ress[,3]-ress[,2])
      avg.seg.len.b <- mean(resb[,3]-resb[,2])
      
      med.seg.len.o <- median(reso[,3]-reso[,2])
      med.seg.len.s <- median(ress[,3]-ress[,2])
      med.seg.len.b <- median(resb[,3]-resb[,2])
      
      snpo <- outer(reso[,2],alt_snps,"<=")*outer(reso[,3],alt_snps,">=")
      snps <- outer(ress[,2],alt_snps,"<=")*outer(ress[,3],alt_snps,">=")
      snpb <- outer(resb[,2],alt_snps,"<=")*outer(resb[,3],alt_snps,">=")
      
      rejo <- alt_snps[colSums(snpo)>0]
      rejs <- alt_snps[colSums(snps)>0]
      rejb <- alt_snps[colSums(snpb)>0]
      
      Domain.o <- uniqueN(snp_loc[rejo,"DomainID"])
      Domain.s <- uniqueN(snp_loc[rejs,"DomainID"])
      Domain.b <- uniqueN(snp_loc[rejb,"DomainID"])
      
      CV.power.o <- sum(colSums(snpo)>0)/length(sim_map$alt_snps)
      CV.power.s <- sum(colSums(snps)>0)/length(sim_map$alt_snps)
      CV.power.b <- sum(colSums(snpb)>0)/length(sim_map$alt_snps)
      
      FDP.power.o <- mean(rowSums(snpo)==0)
      FDP.power.s <- mean(rowSums(snps)==0)
      FDP.power.b <- mean(rowSums(snpb)==0)
      
      FDP.node.o <- mean(1-rowMeans(snpo))
      FDP.node.s <- mean(1-rowMeans(snps))
      FDP.node.b <- mean(1-rowMeans(snpb))
      
      rejb <- rejs <- rejo <- NULL
      for(j in seq(nrow(reso))){
        rejo <- union(rejo,reso[j,2]:reso[j,3])
      }
      for(j in seq(nrow(ress))){
        rejs <- union(rejs,ress[j,2]:ress[j,3])
      }
      for(j in seq(nrow(resb))){
        rejb <- union(rejb,resb[j,2]:resb[j,3])
      }
      alt.ds <- unique(unlist(sim_map$domain_w_sa_alt_snps))
      rejdsID.o <- as.numeric(unique(snp_loc$DomainID[rejo]))
      rejdsID.s <- as.numeric(unique(snp_loc$DomainID[rejs]))
      rejdsID.b <- as.numeric(unique(snp_loc$DomainID[rejb]))

      FDP.dc.o=uniqueN(intersect(rejdsID.o,alt.ds))/uniqueN(rejdsID.o)
      FDP.dc.s=uniqueN(intersect(rejdsID.s,alt.ds))/uniqueN(rejdsID.s)
      FDP.dc.b=uniqueN(intersect(rejdsID.b,alt.ds))/uniqueN(rejdsID.b)

      Prop.snp.o=length(intersect(alt_snps,rejo))/length(rejo)
      Prop.snp.s=length(intersect(alt_snps,rejs))/length(rejs)
      Prop.snp.b=length(intersect(alt_snps,rejb))/length(rejb)
      
      lrejo <- length(rejo)
      lrejs <- length(rejs)
      lrejb <- length(rejb)
      
    }

    out <- rbind(out,c("SS_Factor"=sampsize.factors[i],
             "CV.power.o"=CV.power.o,
             "CV.power.s"=CV.power.s,
             "CV.power.b"=CV.power.b,
             "FDP.power.o"=FDP.power.o,
             "FDP.power.s"=FDP.power.s,
             "FDP.power.b"=FDP.power.b,
             "FDP.node.o"=FDP.node.o,
             "FDP.node.s"=FDP.node.s,
             "FDP.node.b"=FDP.node.b,
             "Domain.o"=Domain.o,
             "Domain.s"=Domain.s,
             "Domain.b"=Domain.b,
             "Prop.snp.o"=Prop.snp.o,
             "Prop.snp.s"=Prop.snp.s,
             "Prop.snp.b"=Prop.snp.b,
             "rej.len.o"=lrejo,
             "rej.len.s"=lrejs,
             "rej.len.b"=lrejb,
             "avg.seg.len.o"=avg.seg.len.o,
             "avg.seg.len.s"=avg.seg.len.s,
             "avg.seg.len.b"=avg.seg.len.b,
             "med.seg.len.o"=med.seg.len.o,
             "med.seg.len.s"=med.seg.len.s,
             "med.seg.len.b"=med.seg.len.b,
             "FDP.dc.o"=FDP.dc.o,
             "FDP.dc.s"=FDP.dc.s,
             "FDP.dc.b"=FDP.dc.b
             ))
  }
  return(out)
}



Test_by_SCANG_ROC <- function(N1,
                          pop,
                          sampsize.factors,
                          sim_map,
                          params,
                          covars=FALSE,
                          filter=2e-5,
                          f=0.5,
                          seed,
                          SampleID,
                          Rocp=seq(0.05,0.95,0.05)){
  
  set.seed(seed)
  
  ### Generate Sample Data ###
  out <- NULL
  Alt_snps <- NULL
  for(i in seq_along(sampsize.factors)){
    ## Generate case-control samples
    N0 <- ceiling(N1*sampsize.factors[i])
    
    ## Create data
    attrDat <- DGP2_SCANG_subsamp(pop,N0 = N0,N1 = N1,covars=covars,sim_map=sim_map,
                                  samps=SampleID[[i]])

    ## Phenotype
    phenotypedata <- attrDat$nulldata
    
    ## Genotype
    if(class(attrDat$geno_case_matrix)%in%c("CsparseMatrix","ngCMatrix")){
      mat1 <- unname(as(attrDat$geno_case_matrix,"dgCMatrix"))
      mat0 <- unname(as(attrDat$geno_ctrl_matrix,"dgCMatrix"))
    } else{
      mat1 <- unname(as(attrDat$geno_case_matrix,"dgTMatrix"))
      mat0 <- unname(as(attrDat$geno_ctrl_matrix,"dgTMatrix"))
    }
    
    genotype0 <- as(rbind(mat1,mat0),"dgCMatrix")
    
    index1 <- which(colSums(genotype0)>0)
    
    ans <- map_df(sim_map$domain_snp_list, ~as.data.frame(.x), .id="DomainID")
    colnames(ans) <- c("DomainID","snp_name")
    snp_loc <- data.table(sim_map$snp_loc)%>%left_join(ans,by="snp_name")
    snp_loc <- snp_loc[index1,]
    snp_loc <- snp_loc %>% mutate("loc_order"=order(snp_loc))
    
    genotype <- genotype0[,index1]
    genotype <- genotype[,snp_loc$loc_order]
    
    snp_loc <- snp_loc%>%arrange(snp_loc)%>%
      mutate(snp_name2=seq(n()))
    Ls <- Find_L(snp_loc)
    Lmin <- Ls[1]
    Lmax <- Ls[2]
    
    alt_snps <- data.frame("snp"=snp_loc$snp_name2[which(snp_loc$snp_name%in%sim_map$alt_snps)],
                           "seed"=seed,"SS_Factor"=i,"Covar"=covars)
    Alt_snps <- rbind(Alt_snps,alt_snps)
    
    for(alpha in Rocp){
      if(covars){
        phenotype <- as.matrix(phenotypedata$phenotype)
        Covs <- as.matrix(phenotypedata%>%dplyr::select(-c(phenotype)))
        t1 <- Sys.time()
        obj_nullmodel <- fit_null_glm_SCANG(phenotype~Covs,data=phenotypedata,family=binomial)
        res_lm <- tryCatch(SCANG(genotype,obj_nullmodel,Lmin,Lmax,filter=filter,f=f,
                                 rare_maf_cutoff=1,alpha=alpha),error=function(e) e)
        ptime <- Sys.time()-t1
      }else{
        t1 <- Sys.time()
        obj_nullmodel <- fit_null_glm_SCANG(phenotype~1,data=phenotypedata,family=binomial)
        res_lm <- tryCatch(SCANG(genotype,obj_nullmodel,Lmin,Lmax,filter=filter,f=f,
                                 rare_maf_cutoff=1,alpha=alpha),error=function(e) e)
        ptime <- Sys.time()-t1
      }
      
      if(inherits(res_lm,"error")){
        Res <- NULL;
      }else{
        reso <- data.frame(res_lm$SCANG_O_res)%>%mutate("Test"="SCANGO")
        ress <- data.frame(res_lm$SCANG_S_res)%>%mutate("Test"="SCANGS")
        resb <- data.frame(res_lm$SCANG_B_res)%>%mutate("Test"="SCANGB")
        Res <- rbind(reso,ress,resb)[,c("X2","X3","Test")]%>%
          mutate("seed"=seed,"SS_Factor"=i,"alpha"=alpha,"Covar"=covars)
        colnames(Res)[1:3] <- c("Start","End","Test")
      }
      out <- rbind(out,Res)
    }
  }
  return(list("SCANG"=out,"alt_snps"=Alt_snps))
}


Get_Alt_info <- function(N1,
                              pop,
                              sampsize.factors,
                              sim_map,
                              params,
                              covars=FALSE,
                              filter=2e-5,
                              f=0.5,
                              seed,
                              SampleID,
                              Rocp=seq(0.05,0.95,0.05)){
  set.seed(seed)
  ### Generate Sample Data ###
  Snp_loc <- NULL
  for(i in seq_along(sampsize.factors)){
    ## Generate case-control samples
    N0 <- ceiling(N1*sampsize.factors[i])
    
    ## Create data
    attrDat <- DGP2_SCANG_subsamp(pop,N0 = N0,N1 = N1,covars=covars,sim_map=sim_map,
                                  samps=SampleID[[i]])
    
    ## Phenotype
    phenotypedata <- attrDat$nulldata
    
    ## Genotype
    if(class(attrDat$geno_case_matrix)%in%c("CsparseMatrix","ngCMatrix")){
      mat1 <- unname(as(attrDat$geno_case_matrix,"dgCMatrix"))
      mat0 <- unname(as(attrDat$geno_ctrl_matrix,"dgCMatrix"))
    } else{
      mat1 <- unname(as(attrDat$geno_case_matrix,"dgTMatrix"))
      mat0 <- unname(as(attrDat$geno_ctrl_matrix,"dgTMatrix"))
    }
    
    genotype0 <- as(rbind(mat1,mat0),"dgCMatrix")
    
    index1 <- which(colSums(genotype0)>0)
    
    ans <- map_df(sim_map$domain_snp_list, ~as.data.frame(.x), .id="DomainID")
    colnames(ans) <- c("DomainID","snp_name")
    snp_loc <- data.table(sim_map$snp_loc)%>%left_join(ans,by="snp_name")
    snp_loc <- snp_loc[index1,]
    snp_loc <- snp_loc %>% mutate("loc_order"=order(snp_loc))
    
    snp_loc <- snp_loc%>%arrange(snp_loc)%>%
      mutate(snp_name2=seq(n()),
             "seed"=seed,"SS_Factor"=i,"Covar"=covars)
    
    #alt_snps <- data.frame("snp"=snp_loc$snp_name2[which(snp_loc$snp_name%in%sim_map$alt_snps)],
    #                       "seed"=seed,"SS_Factor"=i,"Covar"=covars)
    Snp_loc <- snp_loc%>%bind_rows(Snp_loc)
  
  }
  return(Snp_loc)
}
