
gamMissPar <- function(sdyType,missType,chemlst,thylst,msPt2){
  
  msPt <- msPt2[[1]]
  msBmi <- msPt2[[2]]
  sd2 <- msPt2[[3]]
  dirLoc <- c('C:/Users/Michelle Mellers/Documents/ThyAnalysR')
  

  require(tidyverse)
  require(parallel)
  
  # load the functions
  
  source(file.path(dirLoc,"funcs_reg_8_24_21.R"))
 f <- file.path(dirLoc,"nonLinChk",c("fcns_DataProc_10142022.R",
                                      "fcns_viz_10142022.R",
                                      "dataAnaly_fcn_08052022.R",
                                       "fcn_missCreate_08052022.R",
                                       "fcns_basic_08012022.R",
                                       "simStudy_fcns_08052022.R"))
  lapply(f, source)

  
  nlen <- sum(!is.na(chemlst))*sum(!is.na(thylst))*5*length(missType)*length(sdyType)
  ## declare output datasets
  rout <- data.frame(missType=character(nlen),out=character(nlen),ins=character(nlen),
                     covs=character(nlen),est=double(nlen),se=double(nlen),t=double(nlen),
                     df=double(nlen),pv=double(nlen),r2=double(nlen),sigma=double(nlen),
                     sigmaT=double(nlen),
                     miss_chem=double(nlen),miss_bmi=double(nlen))
  tt <- 1
  
  sz <- 1
  nlh <- length(missType)
  o_all3 <- data.frame(outs=character(nlh),rt=double(nlh),medians=double(nlh),mins=double(nlh),maxs=double(nlh),
                       f25=double(nlh),f75=double(nlh))
  if (sdyType=='Sim'){
  nlen2 <- length(chemlst)*length(thylst)*length(missType)
  o_rmserA <-data.frame(simNum=double(nlen2),rmse=double(nlen2),type=character(nlen2),outs=character(nlen2))
  nl2 <- 1
  nlen3 <- 741*length(missType)*length(thylst)
  o_biasA <- data.frame(outs=character(nlen3),ID=double(nlen3),trialm=double(nlen3),bias=double(nlen3),
                        type=character(nlen3),measured_true=double(nlen3),fitTrue=double(nlen3))
  nl3 <- 1
  nls4 <- 741*length(chemlst)*length(thylst)*length(missType)
  o_allTA <- data.frame(outs=character(nls4),simNum=double(nls4),ID=double(nls4),measTrue=double(nls4),
                       measured_imp_log=double(nls4),se=double(nls4),fit=double(nls4),
                       fitTrue=double(nls4),fitTrueSecDev=double(nls4),msTp=character(nls4))
  }
  
  for (rr in 1:length(missType)){
  
  msTp <- missType[[rr]]
 # sdTp <- sdyType[[zz]]
  
  
 # analys_LOD <- dataProcLOD(analys5a)
  if(sdyType == 'Sim'){
    nls <- 741*length(chemlst)*length(thylst)
    o_all <- data.frame(outs=character(nls),simNum=double(nls),ID=double(nls),measTrue=double(nls),
                        measured_imp_log=double(nls),se=double(nls),fit=double(nls),
                        fitTrue=double(nls),fitTrueSecDev=double(nls))
    oo <- 1
    o_allT <- data.frame(outs=character(nls),simNum=double(nls),ID=double(nls),measTrue=double(nls),
                         measured_imp_log=double(nls),se=double(nls),fit=double(nls),
                         fitTrue=double(nls),fitTrueSecDev=double(nls))
    o2 <- 1
  }
  
  for (sst in 1:length(thylst)){
  for (ss in 1:length(chemlst)){
    
    if (sdyType == 'Sim'){
      LODT <- 0
      fcnT <- thylst[[sst]]
      foo <-  formData(msPt,sd2,LODT,fcnT)
      expos_impL <- genM(foo,msBmi)
      dataInnsC <- expos_impL
      if (msTp == 'NM'){
        dataInnsC <- dataInnsC %>% 
          select(-bmi,-measured_imp_log) %>%
          mutate(measured_imp_log=measured_true,bmi=bmiComp,bmi_imp=bmiComp)
        expos_impL <- expos_impL %>%  select(-bmi,-measured_imp_log) %>%
          mutate(measured_imp_log=measured_true,bmi=bmiComp,bmi_imp=bmiComp)
      }
    }  
   
      ins <- 'Sim'
      thySim <- simFcnNms(thylst[[sst]])
      outs <- 'yNoisy'
    
  if (msTp %in% c('NM','CC','SI')){
   # ins = chemlist[[ss]]
    
      dataInns <- expos_impL

    
  if (msTp=='CC'){
    dataInns <- createCC(dataInns)
  } else if (msTp == 'SI') {
    dataInns <- createSI(dataInns)
  } 
  dataInns <- dataCh(dataInns)
  
  testDataMean <- meanData(dataInns,msTp)
  knoter <- NA_real_
  regs <- gamFnc(dataInns,outs,knoter)
  regsOut <- regs
  fits <- predInd(testDataMean,regs)
  if (sdyType=='Sim'){
    if (msTp %in% c('SI','NM')){
  fits22 <- predict(regs,newdata=dataInns,type='response',se=T)
    }else{
      dataInnsC2 <- dataInnsC %>% select(-measured_imp_log) %>% mutate(measured_imp_log=measured_true) %>%
        rename(bmi_imp_cen=bmi_cenB,age_TH_cen=age_TH_cenB)
      fits22 <- predict(regs,newdata=dataInnsC2,type='response',se=T)
    }}
  } else if (msTp == 'NMI'){
  

      dataInns <- createNMI(expos_impL)
  
  testDataMean <- meanData(dataInns,msTp)
  regs <- rep(list(vector(mode="list",length=10)),10)
  fits <- rep(list(vector(mode="list",length=10)),10)
  sigmas <- rep(list(vector(mode="list",length=10)),10)
  fits22 <- rep(list(vector(mode="list",length=10)),10)
  
  knoter <- gamFcnKn(dataInns,ins,outs)
  for (t in 1:10){  
    for (s in 1:10){
         regs[[t]][[s]] <- gamFnc(dataInns[[t]][[s]],outs,knoter)
      fits[[t]][[s]] <- predInd(testDataMean,regs[[t]][[s]])
      if (sdyType=='Sim'){
        fits22[[t]][[s]] <- predict(regs[[t]][[s]],newdata=dataInns[[t]][[s]],type='response',se=T)
      }
      }
  }
  
  regsOuta <- nmiProc(regs)
  regsOut <- avgGAM(regs,regsOuta)
  }
  

    ins2 <- as.character(chemlst[[ss]])
    outs2 <- thySim

  sigR <- sigR2G(regs,msTp)
  rout <- addRoutG(rout,regsOut,sigR,msTp,ins2,outs2,tt)
  
    rout$miss_chem[tt:(tt+4)] <- msPt*100
    rout$sigmaT[tt:(tt+4)] <- sd2


      rout$miss_bmi[tt:(tt+4)] <- signif(sum(is.na(dataInnsC$bmi))/nrow(dataInnsC)*100,2)
  
  
  tt <- tt+5
  

  predicts <- predComb(msTp,testDataMean[[2]],fits,dataInns,ins,outs)
  
    if (msTp=='NMI'){
      dI2 <- dataInns[[1]][[1]]$ID
    } else {
      dI2 <- dataInns$ID
    }
    predicts2 <- data.frame(predicts,ID=dI2)
    predicts3 <- foo %>% select(yB,yBd2,ID,measured_true) %>% merge(.,predicts2,by="ID") %>%  
      mutate(fit=fits,se=se.fits,simNum=ss,fitTrue=yB,measTrue=measured_true,outs=thySim,fitTrueSecDev=yBd2)
    o_all <- addOG(o_all,predicts3,oo)
    oo <-oo+nrow(predicts3)
    
    ##add predicts here
    if (msTp == 'CC'){
      predictor3 <- predComb(msTp,dataInnsC2,fits22,dataInnsC2,ins,outs)
    } else {
      predictor <- predComb(msTp,testDataMean[[2]],fits22,dataInns,ins,outs)
      predictor2 <- data.frame(predictor,ID=dI2)
      predictor3 <- foo %>% select(yB,yBd2,ID,measured_true) %>% merge(.,predictor2,by="ID")
    }
    predictor4 <- predictor3 %>% mutate(fit=fits,simNum=ss,se=se.fits,fitTrue=yB,fitTrueSecDev=yBd2,measTrue=measured_true
                                        ,outs=thylst[[sst]])
    o_allT <- addOG(o_allT,predictor4,o2)
    o2 <- o2+nrow(predictor4)
    
    #####
   ##missing data type
  
  } ##end chem list 
  } ##end thyroid list
  
  ##put
  if (sdyType=="Sim"){
    
    o_rmser <- o_allT %>% filter(ID!=0) %>% group_by(simNum,outs) %>%  mutate(diffs=(fit-as.numeric(fitTrue))^2) %>%
      summarize(sums=sum(diffs)) %>%
      mutate(rmse=signif(sqrt(sums/741),digits=3)) %>% select(-sums) %>% 
      mutate(type=msTp)
  o_rmserA$simNum[nl2:(nl2-1+nrow(o_rmser))] <- o_rmser$simNum
  o_rmserA$rmse[nl2:(nl2-1+nrow(o_rmser))] <- o_rmser$rmse
  o_rmserA$type[nl2:(nl2-1+nrow(o_rmser))] <- o_rmser$type
  o_rmserA$outs[nl2:(nl2-1+nrow(o_rmser))] <- o_rmser$outs
  nl2 <- nl2+nrow(o_rmser)
  
  for (sz in 1:length(thylst)){
  o_allm <- o_allT %>% filter(ID != 0,outs==thylst[[sz]]) %>% mutate(se.type2= se^2,maxSim=max(simNum)) %>% group_by(ID) %>% 
    summarize(se2=sum(se.type2),fitTrueSecDev=mean(fitTrueSecDev),
                                                   fit=mean(fit),maxSim=max(maxSim),fitTrue=mean(as.numeric(fitTrue)),measured_true=mean(as.numeric(measTrue))) %>%
    mutate(se=sqrt(se2/maxSim)) %>% select(-maxSim,-se2)
    outsNm <- simFcnNms(thylst[[sz]])
    vizMiss(msTp,o_allm,ins,outs,outsNm)
}

    o_allTA$outs[nl3:(nl3-1+nrow(o_allT))] <- o_allT$outs
    o_allTA$simNum[nl3:(nl3-1+nrow(o_allT))] <- o_allT$simNum
    o_allTA$ID[nl3:(nl3-1+nrow(o_allT))] <- o_allT$ID
    o_allTA$measTrue[nl3:(nl3-1+nrow(o_allT))] <- o_allT$measTrue
    o_allTA$measured_imp_log[nl3:(nl3-1+nrow(o_allT))] <- o_allT$measured_imp_log
    o_allTA$se[nl3:(nl3-1+nrow(o_allT))] <- o_allT$se
    o_allTA$fit[nl3:(nl3-1+nrow(o_allT))] <- o_allT$fit
    o_allTA$fitTrue[nl3:(nl3-1+nrow(o_allT))] <- o_allT$fitTrue
    o_allTA$fitTrueSecDev[nl3:(nl3-1+nrow(o_allT))] <- o_allT$fitTrueSecDev
    o_allTA$msTp[nl3:(nl3-1+nrow(o_allT))] <- msTp
     nl3 <- nl3+nrow(o_allT)
    # 

  }
  
  } ##end missing type
  
  routBB <- list(rout,o_rmserA,o_allTA)
  return(routBB)
}
  