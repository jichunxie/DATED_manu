########################################################
#### Original Code by: John Pura
#### Modified by: Xuechan Li
#### DATED package is built based on the functions in this file. The Test_DGS function is equivalent to the Test_Leaf function in DATED package.
########################################################
require(data.table)
require(purrr)
require(dplyr)
require(pracma)
require(Matrix)
require(ggplot2)
require(reshape2)
require(tidyverse)


#Source Saddlepoint Approximation functions
#source("/hpc/home/xl110/SATET/R_Func/SPA_functions.R")

################ MAIN FUNCTION ###################
Test_dc <- function(Gmat_case,Gmat_ctrl,
                    struct_map,
                    sim_map,
                    glm_input=NULL,
                    teststat="FET",
                    midp=TRUE,
                    seed=1){
  
  
  #Some processing
  if(class(Gmat_case)%in%c("CsparseMatrix","ngCMatrix")){
    mat1 <- unname(as(Gmat_case,"dgCMatrix"))
    mat0 <- unname(as(Gmat_ctrl,"dgCMatrix"))
  } else{
    mat1 <- unname(as(Gmat_case,"dgTMatrix"))
    mat0 <- unname(as(Gmat_ctrl,"dgTMatrix"))
  }
  
  mat_all <- rbind(mat1,mat0)
  
  N1 <- nrow(mat1)
  N0 <- nrow(mat0)
  N <- N1 + N0
  
  struct_map <- struct_map%>%mutate(L1=rleid(domainID))
  setkey(struct_map,snpID)
  
  leaf_mat_all <- create_leaf_attribute(mat_all,struct_map)
  
  if(teststat=="FET"){
    
    #leaf_mat_all <- create_leaf_attribute(mat_all,struct_map)
    leaf_mat_all@x <- ifelse(leaf_mat_all@x>0,1,0) #binarize
    
    case_colSums = Matrix::colSums(leaf_mat_all[1:N1,])
    all_colSums = Matrix::colSums(leaf_mat_all)
    
    pvals.1 <- calcFETpval_per_leaf(N1=N1,N0=N0,
                                    case_colSums = case_colSums,
                                    all_colSums = all_colSums,
                                    midp=midp)
    
  } else{
    if(!is.null(glm_input) && ncol(glm_input)>1){
      #With covariates
      #Will consider the intercept term when there are covariates
      score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                             pheno=glm_input[,1],cov=glm_input[,-1],
                                             minmac=1,Cutoff=2)
      
    } else{
      #Without covariates
      score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                             pheno=glm_input[,1],cov=NULL,
                                             minmac=1,Cutoff=2)
    }
    
    score.1 <- score.test$Tstat.sign
    pvals.1 <- score.test$p.value
  }
  
  #If NaN, set to zero
  #(occurs when genotype vector is zero, so denominator of score stat is zero)
  pvals.1[is.na(pvals.1)] <- 1
  out <- data.frame("pvals"=pvals.1,"Test"=teststat,
                    "Seed"=seed,"L1"=seq((pvals.1)))%>%
    left_join(struct_map,by="L1")
  return(out)
}



Test_DGS <- function(snp_dat=NULL,thresh_val=10,
                     struct_map=NULL,
                     Gmat_case=NULL,Gmat_ctrl=NULL,
                     glm_input=NULL,
                     teststat="FET",seed=1){
  if(is.null(struct_map)){
    if(is.null(snp_dat)){
      N0 <- nrow(Gmat_ctrl)
      N1 <- nrow(Gmat_case)
      see=rbind(Gmat_case,Gmat_ctrl)
      Gmat <- data.table("Sample.Type"=c(rep("case",N1),rep("ctrl",N0)),
                         "Sample.Name"=seq(N0+N1),
                         as.matrix(rbind(Gmat_case,Gmat_ctrl)))
      colnames(Gmat)[-(1:2)] <- as.character(seq(ncol(pop$Gmat)))
      domaindf <- reshape2::melt(sim_map$domain_snp_list)
      Gmat2 <- reshape2::melt(Gmat,id.vars=c("Sample.Type","Sample.Name"),
                              value.name="snp")
      Gmat2 <- Gmat2[Gmat2$snp,]
      Gmat2$variable <- as.numeric(Gmat2$variable)
      Gmat2=merge(Gmat2,domaindf,by.x=3,by.y=1)
      Gmat2$snp <- ifelse(is.na(Gmat2$snp),FALSE,TRUE)
      snp_dat <- Gmat2
    }
    struct_map <- construct_leafs(snp_dat=snp_dat,
                                  #geno_case_matrix=Gmat_case,
                                  #geno_ctrl_matrix=Gmat_ctrl,
                                  thresh_val=thresh_val)
  }
  total_leaves <- uniqueN(struct_map$L1)
  D_approx_prev = FD_approx_prev = 0
  
  if(is.null(Gmat_case)|is.null(Gmat_ctrl)){
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
  }else{
    #Some processing
    if(class(Gmat_case)%in%c("CsparseMatrix","ngCMatrix")){
      mat1 <- unname(as(Gmat_case,"dgCMatrix"))
      mat0 <- unname(as(Gmat_ctrl,"dgCMatrix"))
    } else{
      mat1 <- unname(as(Gmat_case,"dgTMatrix"))
      mat0 <- unname(as(Gmat_ctrl,"dgTMatrix"))
    }
    
    mat_all <- rbind(mat1,mat0)
    
    n.snps <- m <- ncol(mat_all)
    N1 <- nrow(mat1)
    N0 <- nrow(mat0)
    N <- N1 + N0
  }
  #domain_snp_list <- unique(struct_map$domainID)
  
  leaf_mat_all <- create_leaf_attribute(mat_all,struct_map)
  leaf_mat_all@x <- ifelse(leaf_mat_all@x>0,1,0)
  
  if(teststat=="FET"){
    #Get marginals - qualifying variants for cases and all
    case_colSums = Matrix::colSums(leaf_mat_all[1:N1,])
    all_colSums = Matrix::colSums(leaf_mat_all)
    
    pvals.1 <- calcFETpval_per_leaf(N1=N1,N0=N0,
                                    case_colSums = case_colSums,
                                    all_colSums = all_colSums,
                                    midp=TRUE)
    pvals.1 <- pmin(pvals.1,1)
   
  } else{
    if(!is.null(glm_input) && ncol(glm_input)>1){
      #With covariates
      score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                             pheno=glm_input[,1],cov=glm_input[,-1],
                                             minmac=1,Cutoff=2)
    } else{
      #Without covariates
      score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                             pheno=glm_input[,1],cov=NULL,
                                             minmac=1,Cutoff=2)
    }
    score.1 <- score.test$Tstat.sign
    pvals.1 <- score.test$p.value
  }
  
  Z.1 <- qnorm(1-pvals.1)
  pvals.1[is.na(pvals.1)] <- 1
  Z.1[is.na(pvals.1)] <- -100
  Z.1 <- pmin(Z.1,100)
  Z.1 <- pmax(Z.1,-100)
  S <- NULL
  #Do not test the leaf with pvals.1=1###########check!!!!
  struct_map0 <- data.table("L1"=seq_along(pvals.1),"pvals"=pvals.1)%>%
    left_join(struct_map,by="L1")%>%filter(pvals!=1)%>%mutate(L0=L1)%>%
    mutate(L1=rleid(L1))%>%mutate("Test"=teststat,"Seed"=seed)

  return(struct_map0)
}


SATET_data <- function(struct_map, 
                       L=5,
                       alpha=0.05,
                       alpha1=NULL,thresh_val=10){
  
  total_leaves <- uniqueN(struct_map$L1)
  
  struct_map <- data.table(struct_map)
  pvals.1 <- struct_map[,c("pvals","L0")]%>%distinct_all()
  pvals.1 <- c(struct_map$pvals)
  
  # get the estimated mixed node structure
  total_leaves <- uniqueN(struct_map$L1)
  hatn1 <- ceiling(sqrt(total_leaves))
  pn1 <- pvals.1[rank(pvals.1,ties.method="first")==hatn1]
  
  setkey(struct_map,snpID)
  struct_map <- struct_map %>% mutate(hatm1=(pvals<=pn1),pvals1=pvals) %>%
    group_by(L1)%>%mutate(wt=1/n())%>%ungroup()
  
  struct_map_ <- struct_map
  
  Sps <- NULL
  p1s <- p0s <- NULL
  smap_res <- NULL
  S<-NULL
  D_approx_prev = FD_approx_prev = 0
  for(l in seq(L)){
    Ll <- paste0("L",l)
    Lm1 <- paste0("L",l-1)
    if(l>1){
      #First remove domains with fewer than 2^(l-1) L1 leaves
      removed_map <- struct_map%>%
        group_by(domainID)%>%filter_at(Lm1,any_vars(uniqueN(.)<2))%>%
        ungroup()
      struct_map <- struct_map%>%
        group_by(domainID)%>%filter_at(Lm1,any_vars(uniqueN(.)>=2))%>%
        ungroup()
      
      setDT(struct_map)
      setkey(struct_map, domainID)
      struct_map <- struct_map[, (Ll) := {
        ## modify rle values
        x <- ceiling(rleid(get(Lm1)) / 2)
        n <- uniqueN(get(Lm1))
        if(n %% 2==1){
          x[x == x[.N]] <- x[.N] - 1
        }
        x
      }, by = .(domainID)][, (Ll) := rleid(domainID,get(Ll))] #reassign groups by domain and Ll
      
      struct_map <- struct_map%>%group_by_at(Ll)%>%#mutate(wt=)%>%
        mutate(pvals=stouffer.p2(pvals1,L1,wt))%>%ungroup()
      Sps <- samp.pvals.leaf_L1(struct_map,Ll,p1s=p1s,p0s=p0s)
    }
    
    m.l = uniqueN(data.frame(struct_map)[,Ll])
    pvals.l=struct_map%>%dplyr::select(c("pvals",all_of(Ll)))%>%distinct_all()
    #Obtain layer specific-threshold
    p.hat = est.p.hat_samp(l=l,
                           D_prev = D_approx_prev,
                           FD_prev = FD_approx_prev,
                           pvals_l=c(pvals.l$pvals),
                           alpha=alpha,
                           alpha1=alpha1,
                           Sps=Sps)
    p.hat.l=p.hat$p.hat
    D_approx_prev = p.hat$D_approx_prev
    FD_approx_prev = p.hat$FD_approx_prev
    if(l==1){
      pvals.1 <- c(pvals.l$pvals)
      S.l1 <- which(pvals.1 <= p.hat$p.hat1)
      if(length(S.l1)==0){
        p1s=min(pvals.1)
        p0s=pvals.1[-which.min(pvals.1)]
      }else{
        p1s=pvals.1[S.l1]
        p0s=pvals.1[-S.l1]
        
        #new (R(1-alpha)), more conservative, 02242022
        p1s1=sort(p1s)[seq(floor(length(p1s)*(1-alpha)))]
        p1s2=sort(p1s)[-seq(floor(length(p1s)*(1-alpha)))]
        p1s=p1s1
        p0s=c(p0s,p1s2)
        #new end
      }
    }
    rej_map <- struct_map%>%filter(pvals<=p.hat.l)%>%mutate(Layer=l)
    smap_res <- smap_res%>%bind_rows(rej_map)
    struct_map <- struct_map%>%filter(pvals>p.hat.l)
    if(nrow(struct_map)==0){next}
  }
  smap_res <- smap_res%>%mutate("thresh_val"=thresh_val)%>%
    dplyr::select(-c("pvals","hatm1"))
  return(smap_res)
}


SATET_qq <- function(struct_map, 
                     L=5,
                     alpha=0.05,
                     alpha1=NULL,
                     sim_map,thresh_val=10,seed=1){
  
  m = sim_map$num_snps
  nonnulls = sim_map$alt_snps
  nonnulls_omni = sim_map$alt_sparse_snps
  nulls = sim_map$null_snps
  nonnulls_leaf=unique(struct_map$L1[struct_map$snpID%in% nonnulls])
  m_dc = length(sim_map$domain_snp_list)
  nonnulls_dc = unlist(sim_map$domain_w_sa_alt_snps)
  nulls_dc = setdiff(seq(m_dc),nonnulls_dc)
  alt_leaves <- c(na.omit(unique((struct_map%>%filter(snpID%in%nonnulls))$L1)))
  null_leaves <- c(na.omit(setdiff(unique(struct_map$L1),alt_leaves)))
  total_leaves <- uniqueN(struct_map$L1)
  
  struct_map <- data.table(struct_map)
  pvals.1 <- struct_map[,c("pvals","L0")]%>%distinct_all()
  pvals.1 <- c(struct_map$pvals)
  
  # get the estimated mixed node structure
  total_leaves <- uniqueN(struct_map$L1)
  hatn1 <- ceiling(sqrt(total_leaves))
  pn1 <- pvals.1[rank(pvals.1,ties.method="first")==hatn1]
  
  setkey(struct_map,snpID)
  struct_map <- struct_map %>% mutate(hatm1=(pvals<=pn1),pvals1=pvals) %>%
    group_by(L1)%>%mutate(wt=1/n())%>%ungroup()%>%
    mutate(alt_snps=ifelse(snpID%in%nonnulls,TRUE,FALSE))%>%
    group_by(Test,L0)%>%mutate(alt_leaf=as.integer(sum(alt_snps)>0))%>%
    ungroup()
  struct_map_ <- struct_map
  
  struct_map <- struct_map_
  Sps <- NULL
  p1s <- p0s <- NULL
  smap_res <- NULL
  S<-NULL
  Pvals <- NULL
  D_approx_prev = FD_approx_prev = 0
  for(l in seq(L)){
    Ll <- paste0("L",l)
    Lm1 <- paste0("L",l-1)
    if(l>1){
      #First remove domains with fewer than 2^(l-1) L1 leaves
      removed_map <- struct_map%>%
        group_by(domainID)%>%filter_at(Lm1,any_vars(uniqueN(.)<2))%>%
        ungroup()
      struct_map <- struct_map%>%
        group_by(domainID)%>%filter_at(Lm1,any_vars(uniqueN(.)>=2))%>%
        ungroup()
      
      setDT(struct_map)
      setkey(struct_map, domainID)
      struct_map <- struct_map[, (Ll) := {
        ## modify rle values
        x <- ceiling(rleid(get(Lm1)) / 2)
        n <- uniqueN(get(Lm1))
        if(n %% 2==1){
          x[x == x[.N]] <- x[.N] - 1
        }
        x
      }, by = .(domainID)][, (Ll) := rleid(domainID,get(Ll))] #reassign groups by domain and Ll
      
      struct_map <- struct_map%>%group_by_at(Ll)%>%#mutate(wt=)%>%
        mutate(pvals=stouffer.p2(pvals1,L1,wt))%>%ungroup()
      Sps <- samp.pvals.leaf_L1(struct_map,Ll,p1s=p1s,p0s=p0s)
    }
    
    m.l = uniqueN(data.frame(struct_map)[,Ll])
    pvals.l=struct_map%>%dplyr::select(c("pvals",all_of(Ll),"alt_leaf"))%>%distinct_all()
    Pval <- pvals.l%>%filter(alt_leaf==0)%>%mutate(Layer=l)%>%dplyr::select(-c(all_of(Ll),"alt_leaf"))
    Pvals <- rbind(Pvals,Pval)
    #Obtain layer specific-threshold
    p.hat = est.p.hat_samp(l=l,
                           D_prev = D_approx_prev,
                           FD_prev = FD_approx_prev,
                           pvals_l=c(pvals.l$pvals),
                           alpha=alpha,
                           alpha1=alpha1,
                           Sps=Sps)
    p.hat.l=p.hat$p.hat
    D_approx_prev = p.hat$D_approx_prev
    FD_approx_prev = p.hat$FD_approx_prev
    if(l==1){
      pvals.1 <- c(pvals.l$pvals)
      S.l1 <- which(pvals.1 <= p.hat$p.hat1)
      if(length(S.l1)==0){
        p1s=min(pvals.1)
        p0s=pvals.1[-which.min(pvals.1)]
      }else{
        p1s=pvals.1[S.l1]
        p0s=pvals.1[-S.l1]
        
        #new (R(1-alpha)), more conservative, 02242022
        p1s1=sort(p1s)[seq(floor(length(p1s)*(1-alpha)))]
        p1s2=sort(p1s)[-seq(floor(length(p1s)*(1-alpha)))]
        p1s=p1s1
        p0s=c(p0s,p1s2)
        #new end
      }
    }
    rej_map <- struct_map%>%filter(pvals<=p.hat.l)%>%mutate(Layer=l)
    smap_res <- smap_res%>%bind_rows(rej_map)
    struct_map <- struct_map%>%filter(pvals>p.hat.l)
    if(nrow(struct_map)==0){next}
  }
  return(Pvals)
}

SATET_multithresh <- function(struct_map, 
                              L=5,
                              alpha=0.05,
                              alpha1=NULL,
                              sim_map,thresh_val=10,seed=1){
  # function to dealing with p-values generated from Test_DGS/SATET_simu
  # rejections only, no FDR,power assessment
  
  m = sim_map$num_snps
  nonnulls = sim_map$alt_snps
  nonnulls_omni = sim_map$alt_sparse_snps
  nulls = sim_map$null_snps
  nonnulls_leaf=unique(struct_map$L1[struct_map$snpID%in% nonnulls])
  m_dc = length(sim_map$domain_snp_list)
  nonnulls_dc = unlist(sim_map$domain_w_sa_alt_snps)
  nulls_dc = setdiff(seq(m_dc),nonnulls_dc)
  alt_leaves <- c(na.omit(unique((struct_map%>%filter(snpID%in%nonnulls))$L1)))
  null_leaves <- c(na.omit(setdiff(unique(struct_map$L1),alt_leaves)))
  total_leaves <- uniqueN(struct_map$L1)
  
  struct_map <- data.table(struct_map)
  pvals.1 <- struct_map[,c("pvals","L0")]%>%distinct_all()
  pvals.1 <- c(struct_map$pvals)
  
  # get the estimated mixed node structure
  total_leaves <- uniqueN(struct_map$L1)
  hatn1 <- ceiling(sqrt(total_leaves))
  pn1 <- pvals.1[rank(pvals.1,ties.method="first")==hatn1]
  
  setkey(struct_map,snpID)
  struct_map <- struct_map %>% mutate(hatm1=(pvals<=pn1),pvals1=pvals) %>%
    group_by(L1)%>%mutate(wt=1/n())%>%ungroup()%>%
    mutate(alt_snps=ifelse(snpID%in%nonnulls,TRUE,FALSE))%>%
    group_by(Test,L0)%>%mutate(alt_leaf=as.integer(sum(alt_snps)>0))%>%
    ungroup()
  struct_map_ <- struct_map
  
  struct_map <- struct_map_
  Sps <- NULL
  p1s <- p0s <- NULL
  smap_res <- NULL
  S<-NULL
  D_approx_prev = FD_approx_prev = 0
  for(l in seq(L)){
    Ll <- paste0("L",l)
    Lm1 <- paste0("L",l-1)
    if(l>1){
      #First remove domains with fewer than 2^(l-1) L1 leaves
      removed_map <- struct_map%>%
        group_by(domainID)%>%filter_at(Lm1,any_vars(uniqueN(.)<2))%>%
        ungroup()
      struct_map <- struct_map%>%
        group_by(domainID)%>%filter_at(Lm1,any_vars(uniqueN(.)>=2))%>%
        ungroup()
      
      setDT(struct_map)
      setkey(struct_map, domainID)
      struct_map <- struct_map[, (Ll) := {
        ## modify rle values
        x <- ceiling(rleid(get(Lm1)) / 2)
        n <- uniqueN(get(Lm1))
        if(n %% 2==1){
          x[x == x[.N]] <- x[.N] - 1
        }
        x
      }, by = .(domainID)][, (Ll) := rleid(domainID,get(Ll))] #reassign groups by domain and Ll
      
      struct_map <- struct_map%>%group_by_at(Ll)%>%#mutate(wt=)%>%
        mutate(pvals=stouffer.p2(pvals1,L1,wt))%>%ungroup()
      Sps <- samp.pvals.leaf_L1(struct_map,Ll,p1s=p1s,p0s=p0s)
    }
    
    m.l = uniqueN(data.frame(struct_map)[,Ll])
    pvals.l=struct_map%>%dplyr::select(c("pvals",all_of(Ll),"alt_leaf"))%>%distinct_all()
    #Obtain layer specific-threshold
    p.hat = est.p.hat_samp(l=l,
                           D_prev = D_approx_prev,
                           FD_prev = FD_approx_prev,
                           pvals_l=c(pvals.l$pvals),
                           alpha=alpha,
                           alpha1=alpha1,
                           Sps=Sps)
    p.hat.l=p.hat$p.hat
    D_approx_prev = p.hat$D_approx_prev
    FD_approx_prev = p.hat$FD_approx_prev
    if(l==1){
      pvals.1 <- c(pvals.l$pvals)
      S.l1 <- which(pvals.1 <= p.hat$p.hat1)
      if(length(S.l1)==0){
        p1s=min(pvals.1)
        p0s=pvals.1[-which.min(pvals.1)]
      }else{
        p1s=pvals.1[S.l1]
        p0s=pvals.1[-S.l1]
        
        #new (R(1-alpha)), more conservative, 02242022
        p1s1=sort(p1s)[seq(floor(length(p1s)*(1-alpha)))]
        p1s2=sort(p1s)[-seq(floor(length(p1s)*(1-alpha)))]
        p1s=p1s1
        p0s=c(p0s,p1s2)
        #new end
      }
    }
    rej_map <- struct_map%>%filter(pvals<=p.hat.l)%>%mutate(Layer=l)
    smap_res <- smap_res%>%bind_rows(rej_map)
    struct_map <- struct_map%>%filter(pvals>p.hat.l)
    if(nrow(struct_map)==0){next}
  }
  smap_res <- smap_res%>%mutate("thresh_val"=thresh_val,"seed"=seed)%>%
    dplyr::select(-c("pvals","hatm1"))
  return(smap_res)
}




SATET_dc <- function(snp_dat,mat_all,
                     glm_input=NULL,
                     teststat="FET",
                     score.cutoff=2,
                     alpha=0.05){
  snp_dat <- snp_dat%>%arrange(domainID,snpID)
  sd <- distinct_all(snp_dat[,c("Sample.Type","Sample.ID")])
  N1 <- table(sd$Sample.Type)[1]
  N0 <- table(sd$Sample.Type)[2]
  N <- N1 + N0
  
  domain_snp_list <- distinct_all(snp_dat[,c("domainID","snpID")])
  
  #Create DT that keeps track of leaf indices at each layer
  struct_map <- data.table(domain_snp_list,key = "snpID")
  
  struct_map <- struct_map%>%arrange(domainID)%>%mutate(L1=rleid(domainID))
  
  leaf_mat_all <- create_leaf_attribute(snp_mat=mat_all,snp_leaf_map=struct_map)
  
  if(teststat=="FET"){
    leaf_mat_all@x <- ifelse(leaf_mat_all@x>0,1,0) #binarize
    
    case_colSums = Matrix::colSums(leaf_mat_all[1:N1,])
    all_colSums = Matrix::colSums(leaf_mat_all)
    
    pvals.1 <- calcFETpval_per_leaf(N1=N1,N0=N0,
                                    case_colSums = case_colSums,
                                    all_colSums = all_colSums,
                                    midp=TRUE)
  }else{
    if(!is.null(glm_input) && ncol(glm_input)>1){
      #With covariates
      #Will consider the intercept term when there are covariates
        score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                               pheno=glm_input[,1],cov=glm_input[,-1],
                                               minmac=1,Cutoff=score.cutoff)
    }else{
      #Without covariates
        score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                               pheno=glm_input[,1],cov=NULL,
                                               minmac=1,Cutoff=score.cutoff)
    }
    score.1 <- score.test$Tstat.sign
    pvals.1 <- score.test$p.value
  }
  #If NaN, set to zero
  #(occurs when genotype vector is zero, so denominator of score stat is zero)
  pvals.1[is.na(pvals.1)] <- 1
  
  S.1 <- which(p.adjust(pvals.1,method="BH")<=alpha)
 
  Genes_domain <- snp_dat%>%
    dplyr::select(c("Gene.Name","subRVIS.Domain.Name","domainID"))%>%
    distinct_all()%>%
    left_join(struct_map,by="domainID")%>%filter(L1%in%S.1)

  return(list(#"p.hats"=p.hat.l,
              "S.list.dc"=S.1,
              "Genes_domain"=unique(Genes_domain$Gene.Name)))
}


Get_Pop <- function(pi1=1,
                    snp_clust_num=20,
                    snp_sparse_clust_num=20,
                    sampsize.factors=seq(4),
                    nrep=5,loadpop,R=200){
  set.seed(123) #Functions for generating population data
  
  # 1. Generate genetic variant
  sim_map <- simGenomicStructure_global_omnigenic_random3(snp_dat=snp_dat,
                                                          snp_clust_num=snp_clust_num,
                                                          pi1=pi1,
                                                          snp_sparse_clust_num=snp_sparse_clust_num)
  pop_name <- paste0(c("pop","pop_cov"),pi1*10,".RData")
  sim_map_name <- paste0(c("sim_map"),pi1*10,".RData")
  save(sim_map,file=paste0(loadpop,sim_map_name))
  
  # 2. Generate population
  #Setup parameters
  n_vars = length(unlist(sim_map$domain_snp_list))
  params = list(beta_rho_factor=c(0.5,1.5)*5, # log OR effect
                gamma=c(0,0), # covariate effects
                maf=10^runif(n_vars,-3,-1.7), # maf 0.0005-0.01
                score.cutoff = 2)
  p0 <- 0.01 # population prevalence, defined as the average disease rate in the simulated population based on their parameters.
  N1 = 1000
  
  Npop=ceiling(nrep*N1/p0)
  
  t1 <- Sys.time()
  # Generate population data (cov=FLASE)
  pop <- DGP2_population(Npop=Npop,sim_map=sim_map,
                         covars=FALSE,params=params,p0=p0) # Estimate intercept term based on large Npop.
  print(Sys.time()-t1)
  
  Y <- pop$nulldata[,"Y"]
  SampleID <- NULL
  for(r in seq(R)){
    Samps <- NULL
    for(i in seq_along(sampsize.factors)){
      N0 <- ceiling(N1*sampsize.factors[i])
      sampcase <- sample(which(Y==1),N1)
      sampctrl <- sample(which(Y==0),N0)
      Samps[[i]] <- list("sampcase"=sampcase,"sampctrl"=sampctrl)
    }
    SampleID[[r]] <- Samps
  }
  save(p0,params,sim_map,pop,SampleID,file=paste0(loadpop,pop_name[1]))
  
  t1 <- Sys.time()
  # Generate population data (cov=TRUE)
  pop <- DGP2_population(Npop=Npop,sim_map=sim_map,
                         covars=TRUE,params=params,p0=p0) # Estimate intercept term based on large Npop.
  print(Sys.time()-t1)
  
  Y <- pop$nulldata[,"Y"]
  SampleID <- NULL
  for(r in seq(R)){
    Samps <- NULL
    for(i in seq_along(sampsize.factors)){
      N0 <- ceiling(N1*sampsize.factors[i])
      sampcase <- sample(which(Y==1),N1)
      sampctrl <- sample(which(Y==0),N0)
      Samps[[i]] <- list("sampcase"=sampcase,"sampctrl"=sampctrl)
    }
    SampleID[[r]] <- Samps
  }
  save(p0,params,sim_map,pop,SampleID,file=paste0(loadpop,pop_name[2]))
}



############# AUXILIARY FUNCTIONS #################
est.p.hat_samp <- function(l,D_prev,FD_prev,pvals_l,alpha,alpha1=NULL,Sps=NULL){
  
  ##print(paste("l:",l))
  if(is.null(alpha1)) {alpha1=alpha}
  m.l = length(pvals_l)
  
  #Threshold very small p-values based on constant
  p.m = ifelse(m.l==0,1,min(1/(m.l*sqrt(log(m.l))),0.05))
  
  filter <- which(pvals_l<=p.m|pvals_l>alpha1)
  
  if(length(filter)>0){
    p.vec = unique(sort(pvals_l[-filter],decreasing = FALSE))
  }else {
    p.vec = unique(sort(pvals_l,decreasing = FALSE))
  }
  if(length(p.vec)==0){p.vec=p.m}
  
  p.indx = 0
  emp.fdr = 0
  
  addi <- 0
  if(!is.null(Sps)){
    propct <- Sps$sps$ct*Sps$sps$altprop
    sps <- Sps$sps%>%dplyr::select(-c("type","altct","nullct","altprop","ct"))
    
    rows=nrow(sps);cols=ncol(sps)
    if(rows==1){
      sps0 <- sapply(sps,function(x){propct%*%outer(x,p.vec,"<=")})
      addi <- mean(sps0)
    }else{
      sps0 <- sapply(sps,function(x){colSums(diag(propct)%*%outer(x,p.vec,"<="))})
      if(length(p.vec)==1){
        addi <- mean(sps0)
      }else{
        addi <- rowMeans(sps0)
      }

    }
  }
  
  fdr.num = FD_prev + m.l*p.vec+addi
  fdr.denom = D_prev + sapply(p.vec,function(x){sum(pvals_l<=x)})
  emp.fdr = fdr.num/pmax(fdr.denom,1)
  p.vec <- c(p.m,p.vec)
  index <- max(which(c(0,emp.fdr)<=alpha))
  index1 <- max(which(c(0,emp.fdr)<=alpha1))
  p.hat <- p.vec[index]
  p.hat1 <- p.vec[index1]
  D_approx_prev=fdr.denom[index]
  FD_approx_prev=min(fdr.num[index],D_approx_prev)
  
  return(list("p.hat"=p.hat,"p.hat1"=p.hat1,
              "FD_approx_prev"=FD_approx_prev,
              "D_approx_prev"=D_approx_prev))
}



####### Functions to compute test-statistics

fisher.exact.test <- function(z,midp=TRUE){
  
  x <- z[1]
  sampTot <- z[2]
  pop1Tot <- z[3]
  pop2Tot <- z[4]
  
  lo <- max(0L, sampTot - pop2Tot)
  hi <- min(sampTot, pop1Tot)
  
  support <- lo : hi
  out <- dhyper(support, pop1Tot, pop2Tot, sampTot)
  
  if(midp){
    #mid p-val with minimum likelihood method
    return(sum(out[out < out[x - lo + 1]]) + sum(out[out==out[x-lo+1]])/2)
  } else{
    #minimum likelihood method
    return(sum(out[out <= out[x - lo + 1]]))
  }
}

compute_counts <- function(thresh, ID, domain_end,snp_end) {
  see_idss <- seen_ids <- NULL
  count <- 0L
  countall <- 0L
  adjust_count <- function(id, domain_end,snp_end) {
    if (!(id %in% seen_ids)) {
      seen_ids <<- c(seen_ids,id)
      count <<- count + 1L
    }
    countall <<- countall+ 1L
    
    if ((snp_end & (uniqueN(seen_ids) >= thresh))|domain_end) {
      count <- count # copy enclosed value locally
      countall <- countall
      seen_ids <<- NULL
      count <<- 0L
      countall <<- 0L
    }
    paste(count,countall,sep=",")
  }
  unlist(Map(adjust_count, ID, domain_end, snp_end))
}


create_leaf_attribute <- function(snp_mat,snp_leaf_map){
  
  snp_mat@Dimnames <- list(NULL,NULL)
  
  #Convert mat to dgTMatrix if not already
  if(class(snp_mat)!="dgTMatrix"){
    snp_mat <- as(snp_mat, "dgTMatrix") 
  }
  
  snp_mat2 <- snp_mat #copy object
  
  #Replace column indices with new set of indices 
  #Make sure initial indices start with zero
  snp_mat2@j <- as.integer(snp_leaf_map[.(snp_mat@j+1)]$L1-1) 
  #Correct dimensions of new matrix
  smij <- distinct_all(data.frame(snp_mat2@i,snp_mat2@j))
  snp_mat2@i <- smij[,1]
  snp_mat2@j <- smij[,2]
  snp_mat2@Dim <- as.integer(c(nrow(snp_mat2),
                               length(unique(snp_leaf_map$L1))))
  
  #Convert to dgCMatrix
  y <- as(snp_mat2,"dgCMatrix")
  return(y)
}

stouffer.p <- function(x){
  pnorm(sum(x)/sqrt(length(x)),lower.tail=FALSE)
}

stouffer.p2 <- function(y,g,wt){
  x <- qnorm(1-y)
  out=pnorm(sum(x*wt)/sqrt(uniqueN(g)),lower.tail=FALSE)
  return(rep(out,length(x)))
}


construct_leafs <- function(snp_dat,
                            thresh_val=10){ ##### data,domainID24 #####
  t1 <- Sys.time()
  #function(thresh, ID, domain_end,snp_end)
  struct_map <- data.table(snp_dat) %>%
    arrange(domainID,snpID)%>%
    mutate_at("Sample.Name",as.character) %>% 
    #Coerce ID from factor to character
    group_by(domainID) %>% 
    mutate(lastObsFlagDomain = as.integer(row_number() == n())) %>%
    group_by(snpID,.add=TRUE) %>%
    mutate(lastObsFlagSnp = as.integer(row_number()==n())) %>%
    ungroup() %>% 
    mutate(num_comb = compute_counts(thresh_val, Sample.Name, 
                                     lastObsFlagDomain,lastObsFlagSnp)) %>%
    mutate(num_unique = colsplit(num_comb,",",c(1:2))[,1]) %>%
    mutate(num_all = colsplit(num_comb,",",c(1:2))[,2]) %>%
    mutate(group = cumsum(c(-1L, diff(num_all)) <= 0L)) %>%
    group_by(domainID) %>%
    mutate(group2 = ifelse(group==max(group) & last(num_unique) < thresh_val,
                           max(max(group)-1L,min(group)),group)) %>%
    ungroup() %>%
    mutate(L1 = rleid(group2)) %>%
    dplyr::select(-c(contains("group"),lastObsFlagDomain,num_unique,num_comb,num_all)) %>%
    mutate_at("Sample.Name",as.factor) %>% #coerce ID back to factor
    data.table() #15-37s, based on N0 size
  
  struct_map <- struct_map %>% dplyr::select(c(snpID,L1,domainID)) %>% distinct(snpID,L1,domainID)
  setkey(struct_map,snpID)
  t2 <- Sys.time()
  
  return(struct_map)
}

samp.pvals.leaf_L1 <- function(struct_map,Ll,Nsamp=100,p1s=NULL,p0s=NULL){
  # change to a more efficient algorithm
  # sampling (check the probability of the entile mixed population)
  delname <- intersect(c("snpID","alt_snps"),colnames(struct_map))
  struct_map2=struct_map%>%group_by_at(Ll)%>%
    mutate(mixed=((prod(hatm1)+prod(!hatm1))==0),
           altprop=sum(!hatm1)/length(hatm1),#altprop here means the estimated proportion of false rejection
           altct=sum(hatm1),
           nullct=sum(!hatm1))%>%ungroup()%>%
    dplyr::select(-delname)%>%distinct_all()
  if(any(sum(struct_map2$mixed)==0,length(p1s)==0)){return(NULL)}else{
    #p0s <- struct_map2$pvals
    struct_map2 <- struct_map2%>%filter(mixed)%>%arrange(hatm1)
                   #%>%dplyr::select(-c("pvals"))%>%distinct_all()
    struct_map3 <- struct_map2%>%group_by(altct,nullct,altprop)%>%
      summarize(ct=n())%>%ungroup()%>%rowid_to_column(var="type")
    sp1s <- sample(p1s,Nsamp*sum(struct_map3$altct),replace=TRUE)
    sp0s <- sample(p0s,Nsamp*sum(struct_map3$nullct),replace=TRUE)
    Sp1s <- data.frame(matrix(qnorm(sp1s,lower.tail=FALSE),ncol=Nsamp))%>%
      mutate(type=rep(struct_map3$type,time=struct_map3$altct))
    Sp0s <- data.frame(matrix(qnorm(sp0s,lower.tail=FALSE),ncol=Nsamp))%>%
      mutate(type=rep(struct_map3$type,time=struct_map3$nullct))
    sps <- rbind(Sp1s,Sp0s)%>%group_by(type)%>%
      summarise_at(paste0("X",seq(Nsamp)), stouffer.p)%>%ungroup()%>%
      left_join(struct_map3,by="type")
    return(list("sps"=sps,"M"=nrow(struct_map2),"altprop"=struct_map2$altprop))
  }
}

############# SIMULATION FUNCTIONS ############
DGP2_population <- function(Npop,
                            sim_map,
                            covars=FALSE,
                            params,
                            p0=0.01){
  n_vars <- length(unlist(sim_map$domain_snp_list))
  
  #Generate covariates
  X1 <- rnorm(Npop,0,1)
  X2 <- rbinom(Npop,1,0.5)
  
  #Generate sparse genotype matrix
  Gmat <- genSpMat(Npop,n_vars,params$maf)
  
  #Alternative SNP locations
  altSNPloc <- sim_map$alt_dense_snps
  alt_omniSNPloc <- sim_map$alt_sparse_snps
  
  #Generate betas's for causal, non-omnigenic variants
  betas <- rep(0,n_vars)
  
  betas[altSNPloc] <- -params$beta_rho_factor[1]*log10(params$maf[altSNPloc])
  betas[alt_omniSNPloc] <- -params$beta_rho_factor[2]*log10(params$maf[alt_omniSNPloc])
  
  #Starting value for alpha0 
  alpha0 <- log(p0/(1-p0)) 
  
  #Covariates and G'beta
  Gb <- as.vector(Gmat[,c(altSNPloc,alt_omniSNPloc)] %*% betas[c(altSNPloc,alt_omniSNPloc)])
  Xg_Gb <- params$gamma[1]*X1 + params$gamma[2]*X2 + Gb
  
  #Starting value for pr
  pr <- plogis(alpha0 + Xg_Gb)
  
  #Adjust intercept until population prevalence is reached (to a certain tolerance)
  if(covars){
    if(mean(pr)>p0){
      #Decrement alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Xg_Gb)
        alpha0 <- alpha0 - 0.05
      }
    } else{
      #Increment alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Xg_Gb)
        alpha0 <- alpha0 + 0.05
      }
    }
  } else{
    if(mean(pr)>p0){
      #Decrement alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Gb)
        alpha0 <- alpha0 - 0.05
      }
    } else{
      #Increment alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Gb)
        alpha0 <- alpha0 + 0.05
      }
    }
  }
  Y <- rbinom(Npop,1,pr)
  return(list(nulldata=data.frame(Y=Y,X1=X1,X2=X2),Gmat=Gmat,Npop=Npop,betas=betas))
}  


DGP2_subsamp <- function(pop,N0,N1,covars,sim_map,samps){

  sampcase <- samps$sampcase
  sampctrl <- samps$sampctrl
  
  #Get attribute matrices
  t1 <- Sys.time()
  geno_case_matrix <- pop$Gmat[sampcase,]
  geno_ctrl_matrix <- pop$Gmat[sampctrl,]
  t2 <- Sys.time()
  #Get Leaf information
  Gmat <- data.table("Sample.Type"=c(rep("case",N1),rep("ctrl",N0)),
                     "Sample.Name"=seq(N0+N1),
                     as.matrix(rbind(geno_case_matrix,geno_ctrl_matrix)))
  colnames(Gmat)[-(1:2)] <- as.character(seq(ncol(pop$Gmat)))
  domaindf <- reshape2::melt(sim_map$domain_snp_list)
  Gmat2 <- reshape2::melt(Gmat,id.vars=c("Sample.Type","Sample.Name"),
                          value.name="snp")%>%
    filter(snp)%>%mutate(variable=as.numeric(variable))%>%
    left_join(domaindf,by=c("variable"="value"))%>%
    rename(domainID = L1,snpID=variable)
  #struct_map <- construct_leafs(Gmat2,thresh_val=thresh_val)

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
         geno_ctrl_matrix=geno_ctrl_matrix,
         Gmat2=Gmat2)
  } else{
    list(nulldata=data.frame(phenotype=rep(c(1,0),times=c(N1,N0))),
         geno_case_matrix=geno_case_matrix,
         geno_ctrl_matrix=geno_ctrl_matrix,
         Gmat2=Gmat2)
  }
}

simGenomicStructure_global_omnigenic_random3 <- function(snp_dat,
                                                         snp_clust_num=20,
                                                         pi1=0.8,
                                                         snp_sparse_clust_num=20
                                                        ){
  #Renumber variants to use less memory
  snp_dat <- snp_dat%>%group_by(subRVIS.Domain.Name)%>%arrange(chrno,loc_adj)
  domain_snpchr_list = split(data.frame(snp_dat[,c("loc_adj","chrno")]),
                          as.character(snp_dat$subRVIS.Domain.Name))
  domain_unique_snpchr_list = unname(lapply(domain_snpchr_list,unique)) 
  
  domain_snploc_list <- unname(lapply(domain_unique_snpchr_list,
                                 function(df){df[,1]}))
  domain_chr_list <- unname(lapply(domain_unique_snpchr_list,
                                 function(df){df[,2]}))
  
  # domain_snp_list = split(snp_dat$loc_adj,
  #                         as.character(snp_dat$subRVIS.Domain.Name))
  # domain_unique_snp_list = unname(lapply(domain_snploc_list,unique))  

  #For each domain keep track of the snp indices - immutable
  domain_snp_list <- relist(seq_along(unlist(domain_snploc_list)),
                            skeleton = domain_snploc_list)
  
  domain_snp_list_ <- domain_snp_list
  
  d.orig <- data.frame("snp_name"=unlist(domain_snp_list),
                  "snp_loc"=unlist(domain_snploc_list),
                  "chr"=unlist(domain_chr_list))
  
  #Get snp IDs for each cluster
  #domain_w_sa_alt_snps <- sample(seq_along(domain_snp_list_)[lengths(domain_unique_leaf_list)>=2],snp_clust_num)
  domain_w_sa_alt_snps <- sample(seq_along(domain_snp_list_)[lengths(domain_snp_list)>=30],snp_clust_num)
  alt_snps_region_len <- floor(pmin(lengths(domain_snp_list_[domain_w_sa_alt_snps])*0.5,20))
  alt_snps_num <- ceiling(alt_snps_region_len*pi1)
  
  cand_sparse <- seq_along(domain_snp_list_)[-domain_w_sa_alt_snps]
  domain_sparse_alt_snps <- sample(cand_sparse,snp_sparse_clust_num)
  alt_sparse_snps <- NULL
  alt_dense_snps <- NULL
  
  for(i in seq_along(domain_sparse_alt_snps)){
    alt_sparse_snps <- c(alt_sparse_snps,
                         sample(domain_snp_list_[[domain_sparse_alt_snps[i]]],1))
  }
  
  for(i in seq_along(domain_w_sa_alt_snps)){
    ds <- domain_snp_list_[[domain_w_sa_alt_snps[i]]]
    alt_dense_snps <- c(alt_dense_snps,
                        sample(domain_snp_list_[[domain_w_sa_alt_snps[i]]][1:alt_snps_region_len[i]],
                               alt_snps_num[i]))
  }
  
  alt_snps <- c(alt_dense_snps,alt_sparse_snps)
  null_snps = setdiff(unlist(domain_snp_list),alt_snps)
  
  return(list("type"="global",
              "num_snps"=length(unique(snp_dat$loc_adj)),
              "domain_snp_list"=domain_snp_list,
              "domain_w_sa_alt_snps"=c(domain_w_sa_alt_snps,domain_sparse_alt_snps),#vector of domains with alt snps.
              "alt_snps"=alt_snps,
              "alt_dense_snps"=alt_dense_snps,
              "alt_sparse_snps"=alt_sparse_snps,
              "null_snps"=null_snps,
              "snp_loc"=d.orig))
}


genSpMat <- function(nrows, ncols, col_probs) {
  r <- lapply(1:ncols, function(x) {
    p <- col_probs[x]
    i <- sample.int(2L, size = nrows, replace = T, prob = c(1 - p, p))
    which(i == 2L)
  })
  rl <- lengths(r)
  nc <- rep(1:ncols, times = rl) # col indexes
  nr <- unlist(r) # row index
  ddims <- c(nrows, ncols)
  #ngTMatrix output
  sparseMatrix(i = nr, j = nc, dims = ddims,giveCsparse = FALSE)
}


SATET_simu <- function(N1,pop,
                       sampsize.factors,
                       covars=FALSE,
                       SampleID,
                       thresh_val,seed=1,pi1=1){
  # combined test DGS to generate p-values
  ### Generate Sample Data ###
  Res <- NULL
  Res_dc <- NULL
  for(i in seq_along(sampsize.factors)){
    ## Generate case-control samples
    N0 <- ceiling(N1*sampsize.factors[i])
    
    ## Create data
    attrDat <- DGP2_subsamp(pop,N0 = N0,N1 = N1,covars=covars,sim_map=sim_map,
                            samps=SampleID[[i]])
    struct_map <- construct_leafs(attrDat$Gmat2,thresh_val=thresh_val)
    
    res_FET_dc =Test_dc(Gmat_case=attrDat$geno_case_matrix,
                        Gmat_ctrl=attrDat$geno_ctrl_matrix,
                        sim_map=sim_map,
                        struct_map=struct_map,
                        glm_input = attrDat$nulldata,
                        teststat = "FET",
                        seed=seed)
    
    res_score_wSPA_dc =Test_dc(Gmat_case=attrDat$geno_case_matrix,
                               Gmat_ctrl=attrDat$geno_ctrl_matrix,
                               sim_map=sim_map,
                               struct_map=struct_map,
                               glm_input = attrDat$nulldata,
                               teststat = "score",
                               seed=seed)
    
    res_FET =Test_DGS(snp_dat=NULL,thresh_val=thresh_val,
                      struct_map=struct_map,
                      Gmat_case=attrDat$geno_case_matrix,
                      Gmat_ctrl=attrDat$geno_ctrl_matrix,
                      glm_input = attrDat$nulldata,
                      teststat = "FET",
                      seed=seed)
    
    res_score_wSPA =Test_DGS(snp_dat=NULL,thresh_val=thresh_val,
                             struct_map=struct_map,
                             Gmat_case=attrDat$geno_case_matrix,
                             Gmat_ctrl=attrDat$geno_ctrl_matrix,
                             glm_input = attrDat$nulldata,
                             teststat = "score",
                             seed=seed)
    
    res_tmp <- data.frame(res_FET)%>%bind_rows(data.frame(res_score_wSPA))%>%
      mutate(SS.Factor=sampsize.factors[i])
    
    res_tmp_dc <- data.frame(res_FET_dc)%>%
      bind_rows(data.frame(res_score_wSPA_dc))%>%
      mutate(SS.Factor=sampsize.factors[i])
    
    Res <- rbind(Res,res_tmp)
    Res_dc <- rbind(Res_dc,res_tmp_dc)
    
  }
  Res <- Res%>%mutate("Covar"=covars,"pi1"=pi1)
  Res_dc <- Res_dc%>%mutate("Covar"=covars,"pi1"=pi1)
  
  return(list("DGS"=Res,"DC"=Res_dc))
}

SATET_simu_dc <- function(N1,pop,
                       sampsize.factors,
                       covars=FALSE,
                       SampleID,
                       thresh_val,seed=1,pi1=1){
  ### Generate Sample Data ###
  Res <- NULL
  Res_dc <- NULL
  for(i in seq_along(sampsize.factors)){
    ## Generate case-control samples
    N0 <- ceiling(N1*sampsize.factors[i])
    
    ## Create data
    attrDat <- DGP2_subsamp(pop,N0 = N0,N1 = N1,covars=covars,sim_map=sim_map,
                            samps=SampleID[[i]])
    struct_map <- construct_leafs(attrDat$Gmat2,thresh_val=thresh_val)
    
    res_FET_dc =Test_dc(Gmat_case=attrDat$geno_case_matrix,
                        Gmat_ctrl=attrDat$geno_ctrl_matrix,
                        sim_map=sim_map,
                        struct_map=struct_map,
                        glm_input = attrDat$nulldata,
                        teststat = "FET",
                        seed=seed)
    
    res_score_wSPA_dc =Test_dc(Gmat_case=attrDat$geno_case_matrix,
                               Gmat_ctrl=attrDat$geno_ctrl_matrix,
                               sim_map=sim_map,
                               struct_map=struct_map,
                               glm_input = attrDat$nulldata,
                               teststat = "score",
                               seed=seed)
    
    res_tmp_dc <- data.frame(res_FET_dc)%>%
      bind_rows(data.frame(res_score_wSPA_dc))%>%
      mutate(SS.Factor=sampsize.factors[i])
    
    Res_dc <- rbind(Res_dc,res_tmp_dc)
    
  }
  Res_dc <- Res_dc%>%mutate("Covar"=covars,"pi1"=pi1)
  
  return("DC"=Res_dc)
}


qqPlot <- function(pval) {
  require(ggplot2)
  pval <- pval[!is.na(pval)]
  n <- length(pval)
  x <- 1:n
  dat <- data.frame(obs=sort(pval),
                    exp=x/n,
                    upper=qbeta(0.025, x, rev(x)),
                    lower=qbeta(0.975, x, rev(x)))
  ggplot(dat, aes(-log10(exp), -log10(obs))) +
    geom_line(aes(-log10(exp), -log10(upper)), color="gray") +
    geom_line(aes(-log10(exp), -log10(lower)), color="gray") +
    geom_point(size=0.3) +
    geom_abline(intercept=0, slope=1, color="red") +
    xlab(expression(paste(-log[10], "(expected P)"))) +
    ylab(expression(paste(-log[10], "(observed P)"))) +
    theme_bw()
  
}



