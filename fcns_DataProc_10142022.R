##update for nonlinear scenarios

dataProcLOD <- function(analys5a){
  
require(tidyverse)
  
analys_LOD <- analys5a %>%
  filter(!(is.na(tshnew_log)|is.na(total_t3_log)|is.na(total_t4_log)|is.na(free_t4_log))) %>%
  select(starts_with("PPB_LW_LOD")) %>%
  summarise_all(quantile,probs=0.25,na.rm=TRUE) %>%
  pivot_longer(cols=everything(),names_to=c("chemical"),values_to="LOD_25") %>%
  mutate(chemical=substr(chemical,12,length(chemical)),LOD_25=signif(log(LOD_25),4)) %>%
  filter(chemical %in% chemlist)

  return(analys_LOD)
}


nmiProc <- function(regs){
  
  require(mice)
  require(miceadds)
  
qhat <- NMIextract(results=regs,fun=coef)
u <- rep(list(vector(mode="list",length=10)),10)
u <- lapply(regs,FUN=function(u44){
  lapply(u44,FUN=function(u55){
    vcov(u55)
  })
})

pval_s <- miceadds::NMIcombine(qhat=qhat,u=u) 


#return(pval)
return(pval_s)
}

meanData <- function(analysMean,missType){
  
  if (missType %in% c('NM','CC','SI')){
  # measured_miss_log1 <- mean(analysMean$measured_miss_log)
  bmi_impM <- analysMean %>% group_by(sexC) %>% 
              summarize(bmi_imp_cen = mean(bmi_imp_cen)) %>%
              mutate(bmi_imp_cen = round(bmi_imp_cen,3))
  #bmi_impM <- round(bmi_impM,1)
  sexCM <- 0
  sexCF <- 1
  age_THM <- analysMean %>% group_by(sexC) %>%
              summarize(age_TH_cen=mean(age_TH_cen)) %>%
              mutate(age_TH_cen=round(age_TH_cen,3))
  #age_THM <- round(age_THM,3)
  measured_imp_log1 <- analysMean$measured_imp_log
  # 
  } else if (missType == 'NMI') {
    #calculate the mean for the test cases..
    peop <- nrow(analysMean[[1]][[1]])
    measured_imp_log1 <- vector(mode="numeric",length=peop)
    bmi_imp1 <- vector(mode="numeric",length=peop)
    for (kq in 1:peop){
      for (kr in 1:10){
        for (ks in 1:10){
          measured_imp_log1[[kq]] <- measured_imp_log1[[kq]] + analysMean[[kr]][[ks]]$measured_imp_log[[kq]]
          bmi_imp1[[kq]] <- bmi_imp1[[kq]] + analysMean[[kr]][[ks]]$bmi_imp_cen[[kq]]
        }
      }
    }
    measured_imp_log1 <- measured_imp_log1/100
    bmi_imp1 <- (bmi_imp1/100)
    # mil1 <- data.frame(measured_imp_log1,sexC=dataInns[[1]][[1]]$sexC)
    bmi_impM <- data.frame(bmi_imp1,sexC=analysMean[[1]][[1]]$sexC) %>%
      group_by(sexC) %>%
      summarise(bmi_imp_cen=mean(bmi_imp1))
    
    #bmi_impM <- mean(bmi_imp1)
    #bmi_impM <- round(bmi_impM,1)
    sexCM <- 0
    sexCF <- 1
    age_THM <-  analysMean[[1]][[1]] %>%
      select(sexC,age_TH_cen) %>%
      group_by(sexC) %>%
      summarise(age_TH_cen=mean(age_TH_cen))
    meas_avg <- analysMean[[1]][[1]] %>%
      select(sexC,LOD) %>%
      mutate(LOD_log=log(LOD)) %>%
      group_by(sexC) %>%
      summarize(LOD_log = mean(LOD_log))
    
  }
  #################
    
  testdataM <- data.frame(measured_imp_log=measured_imp_log1,bmi_imp_cen=bmi_impM[1,2],sexC=sexCM,
                          age_TH_cen=age_THM[1,2]) %>%
              mutate(sexC=as.factor(sexC))
  testdataF <- data.frame(measured_imp_log=measured_imp_log1,bmi_imp_cen=bmi_impM[2,2],sexC=sexCF,
                          age_TH_cen=age_THM[2,2]) %>%
              mutate(sexC=as.factor(sexC))
  
  
  tdTibs <- list(testdataM,testdataF)
  return(tdTibs)
}

predInd <- function(testDat,glmOut){
  
  testdataM <- testDat[[1]]
  testdataF <- testDat[[2]]
  #outs2 <- dataInns[[1]][[1]][,outs]
    
    fitsM <- predict(glmOut,newdata=testdataM,type='response',se=T)
    fitsF <- predict(glmOut,newdata=testdataF,type='response',se=T)
    
    fit <- fitsM$fit*0.54 + fitsF$fit*0.46
    se.fit <- sqrt((fitsM$se.fit*0.54)^2+(fitsF$se.fit*0.46)^2)
    fits <- data.frame(fit,se.fit)
    
  return(fits)
}
  
predComb <- function(typeM,testdataF,fits,dataInns,ins,outs){
    
  if (typeM == 'NMI'){
    
    peop <- nrow(fits[[1]][[1]])
    if (is.null(peop) == TRUE){
      peop <- length(fits[[1]][[1]]$fit)
    }
    fits1 <- vector(mode="numeric",length=peop)
    se.fits1 <- vector(mode="numeric",length=peop)
    for (kq in 1:peop){
      for (kr in 1:10){
        for (ks in 1:10){
          fits1[[kq]] <- fits1[[kq]] + fits[[kr]][[ks]]$fit[[kq]]
          se.fits1[[kq]] <- se.fits1[[kq]] + (fits[[kr]][[ks]]$se.fit[[kq]])^2
        }
      }
    }
  #calculate Qj
    Qj <- vector(mode="numeric",length=peop*10)
    for (kqj in 1:peop){
    for (krj in 1:10){
      ind = kqj+krj-1
      for (ksj in 1:10){
      Qj[[ind]] = Qj[[ind]] + fits[[krr]][[kss]]$fit[[kqq]]
     } 
    }
    } 
    Qj = 1/10*Qj
  #Between nest variance
    B <- vector(mode="numeric",length=peop)
    for (kqq in 1:peop){
    for (krr in 1:10){
      indx = kqq+krr-1
      B[[kqq]] <- B[[kqq]] + (Qj[[indx]]-fits1[[kqq]]/100)^2
    }
    }  
    B = 1/9*B

    #within nest variance
    W <- vector(mode="numeric",length=peop)
    for (kqk in 1:peop){
      for (krk in 1:10){
        indx2 = kqk+krr-1
        for (ksk in 1:10){
          W[[kqk]] <- W[[kqk]] + (fits[[krk]][[ksk]]$fit[[kqk]] - Qj[[indx2]])^2
          }
        }
      }
    W <- 1/10*1/9*W
    
    se.fits1 <- sqrt(se.fits1/100 + (1-1/10)*W + (1+1/10)*B)
    fits1 <- fits1/100
    fits2 <- data.frame(fit=fits1,se.fit=se.fits1) 
  } else {
    fits2 <- fits
  }

    predicts <- data.frame(testdataF,fits=fits2$fit,se.fits=fits2$se.fit,outs2=outs,ins2=ins) %>% 
                mutate(lower=fits-1.96*se.fits,
                               upper=fits+1.96*se.fits)
  
  return(predicts)
}

##create average analysis
avgGAM <- function(regs,pval_s){
  
  pval <- pval_s %>% summary()
  pval <- pval %>% as.data.frame(.) %>%
    mutate(covs=rownames(.))
  ##pval_s is regsOuta
  ##load up a "dummy" gam to find overall s-statistic and fit statistics
  ovr_gam <- regs[[3]][[3]]
  ## the coefficients
  coefs <- pval$est
  names(coefs) <- pval$covs
  ovr_gam$coefficients <- coefs
  ## the variances
  ovr_gam$Vp <- pval_s$Tm
  
  #average degrees of freedom and R matrix
  dfr <- 0
  rM <- matrix(0,nrow=nrow(pval),ncol=nrow(pval))
  edf <- vector(mode="numeric",length=nrow(pval))
  edf1 <- vector(mode="numeric",length=nrow(pval))
  for (ss in 1:10){
    for (sj in 1:10){
      dfr <- dfr + regs[[ss]][[sj]]$df.residual
      rM <- rM + regs[[ss]][[sj]]$R
      edf <- edf + regs[[ss]][[sj]]$edf
      edf1 <- edf1 + regs[[ss]][[sj]]$edf1
    }  
  }
  dfr <- dfr/100
  ovr_gam$df.residual <- dfr
  rM <- rM/100
  ovr_gam$R <- rM
  #average R matricies
  edf <- edf/100
  ovr_gam$edf <- edf
  edf1 <- edf1/100
  ovr_gam$edf1 <- edf1
  
  return(ovr_gam)
  
}
