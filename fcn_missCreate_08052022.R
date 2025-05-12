
createCC <- function(dtObs){

  dtObs2 <- dtObs[complete.cases(dtObs),]
  dtObs2 <- dtObs2 %>% mutate(bmi_imp=bmi)
  return(dtObs2)
}

createSI <- function(dtObs){
  
  dtObs2 <- dtObs %>%
    mutate(measured_imp_log = case_when(is.na(measured_imp_log) == TRUE ~ lod_log-log(sqrt(2)),
                            TRUE~measured_imp_log),
           bmi_imp = case_when((is.na(bmi)==TRUE)~mean(bmi,na.rm=TRUE),
                           TRUE~bmi))
  
  return(dtObs2)
}

createNMI <- function(dtObs){
  dirLoc <- file.path("C:","Users","Michelle Mellers","Documents","ThyAnalysR","Functions")
  #source(file.path(dirLoc,"simStudy_fcns_051822.R"))
  #nlen <- ntrials*length(sdm2)*741*length(LODs)
  
  # if  fcnT == 0, use linear
  # if fcnT == 1, use polynomial
  # if fcnT == 2, use exponential
  require(mice)
  require(miceadds)
  require(miWQS)
  
  ##set variables to impute
  predM <- matrix(rbind(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)),nrow=16,
                  dimnames=list(c("id","y","e","yNoisy","lod_log","offs","measured_imp_log","measured_true","bmi_cenB",
                                  "age_TH_cenB","sexC","raceC","bmi","age_TH","yB","bmiComp"),
                                c("id","y","e","yNoisy","lod_log","offs","measured_imp_log","measured_true","bmi_cenB",
                                  "age_TH_cenB","sexC","raceC","bmi","age_TH","yB","bmiComp")))
  meth <- c("","","","","","","","","","","","","pmm","","","")
  names(meth) <- c("id","y","e","yNoisy","lod_log","offs","measured_imp_log","measured_true","bmi_cenB",
                   "age_TH_cenB","sexC","raceC","bmi","age_TH","yB","bmiComp")
  fooImp <- mice(dtObs,method=meth,m=10,predictorMatrix=predM)
  fooImp2 <- complete(fooImp,"repeated",include=TRUE)
  
  fooImp3 <- pivot_longer(fooImp2,cols=everything(),names_to=c(".value","set"),names_pattern="([A-Za-z_]*).([0-9]*)")
  
  fooImp3 <- dtObs %>% select(id,bmi) %>% rename(bmiM = bmi) %>% merge(.,fooImp3,by="id",all=TRUE) %>% rename(bmi_imp=bmi)
  
  ####
  nset <- 10
  for (kk in 1:nset){
    fooImp3a <- fooImp3 %>% filter(set==kk) %>% select(-yB) 
    rownames(fooImp3a) <- NULL
    lods <- exp(fooImp3a$lod_log[1])
    xlode <- exp(fooImp3a$measured_imp_log)
    covs <- fooImp3a %>% select(age_TH,raceC)
    fooImp4 <- impute.boot(X=xlode,DL=lods,Z=covs, K=10L ,verbose=FALSE)
    #fooImp4 <- impute.Lubin(chemcol=fooImp3a$xlod,dlcol=lods,K=10L,verbose=FALSE)
    fooImp3a <- fooImp3a %>% mutate(idT=row.names(.))
    fooImp5 <- fooImp4[[1]] %>% as.data.frame() %>% mutate(idT=row.names(.)) %>%
      merge(.,fooImp3a,by="idT") %>% 
      mutate(set=as.numeric(set)) %>%
      arrange(id,set) %>%
      pivot_longer(cols=starts_with("1.Imp"),names_to=c("meas_imp","imp_cong"),names_pattern="([A-Za-z_]+).([0-9]+)") %>%
      select(-meas_imp,-idT) %>% rename(imp_bmi = set) %>%
      filter(imp_bmi != 0) %>% mutate(measured_imp_log=log(value))
    if (kk==1){fooImp6<- fooImp5
    }else{
      fooImp6 <- rbind(fooImp6,fooImp5)
    }
  }
  
  fooImp6 <- fooImp6 %>% mutate(sexC=as.factor(sexC)) %>%
    rename(bmi=bmiM)
  fooImp6 <- dataCh(fooImp6)
  #fooImp6 <- fooImp6 %>% mutate(meas_imp <- log(meas_imp))
  expos_impL <- rep(list(vector(mode="list",length=10)),10)
  for (ll in 1:10){
    for (mm in 1:10){
      expos_impL[[ll]][[mm]] <- fooImp6 %>% 
        filter((imp_bmi == ll)&(imp_cong == mm)) %>%
        data.frame()
    }
  }
  
  return(expos_impL)
}  
