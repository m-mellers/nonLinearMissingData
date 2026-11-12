
##formatting tables here
formatTabs <- function(o_rmserA,o_allT,thylst,dirLoc,msPt2){
  msPt <- msPt2[[1]]
  msBmi <- msPt2[[2]]
  sd2 <- msPt2[[3]]
  
source(file.path(dirLoc,'nonLinChk','fcns_viz_10142022.R'),local=TRUE)
require(tidyverse)
  
 #3/28: added abs here for the absolute of the fit-fitTrue 
   o_biasA <-  o_allT %>% filter(ID!=0) %>% group_by(msTp,ID,outs) %>% mutate(fitTrue = as.numeric(fitTrue)) %>%
     mutate(diffs=(fit-fitTrue)) %>%
     summarize(sums=sum(diffs),trialm=max(simNum),measured_true=mean(as.numeric(measTrue)),fitTrue=mean(fitTrue),
               fitTrueSecDev=mean(fitTrueSecDev),sum2=mean(fit-fitTrue)) %>%
     mutate(bias=signif(sums/trialm,digits=3),fitTrueSecDev=signif(fitTrueSecDev,digits=3)) %>% 
     select(-sums) %>% rename(type=msTp)

  
  o_rmserA <- o_rmserA %>% filter(type %in% c('NM','CC','SI','NMI')) %>% 
    mutate(type=case_when(type=='SI'~"S",
                          TRUE~type)) %>%  mutate(type=ordered(type,levels=c('NM','CC','S','NMI')))
  o_biasA <- o_biasA %>% filter(type %in% c('NM','CC','SI','NMI')) %>% 
    mutate(type=case_when(type=='SI'~"S",
                          TRUE~type)) %>% mutate(type=ordered(type,levels=c('NM','CC','S','NMI'))) 
  for (zz in 1:length(thylst)){
    o_rmserA1 <- o_rmserA %>% filter(outs==thylst[[zz]])
    o_biasA1 <- o_biasA %>% filter(outs==thylst[[zz]])
    ggRMSE(o_rmserA1)
    ggBias(o_biasA1)
    ggTrue(o_biasA1)
  }
  o_rmserB <- o_rmserA %>% group_by(outs,type) %>% summarize(val.min=min(rmse),val.25=quantile(rmse,0.25),
                                                             val.50=quantile(rmse,0.50),val.75=quantile(rmse,0.75),
                                                             val.max=max(rmse),val.sd=sd(rmse)) %>% 
    mutate(stat="RMSE")
  o_biasB <- o_biasA %>% group_by(outs,type) %>% mutate(bias=abs(bias)) %>%
    summarize(val.min=min(bias),val.25=quantile(bias,0.25),
              val.50=quantile(bias,0.50),val.75=quantile(bias,0.75),
              val.max=max(bias)) %>% 
    mutate(stat="Bias")
  routMiss <- rout %>% filter(covs=="s(measured_imp_log)") %>% 
    group_by(out) %>%
    summarize(missC=mean(miss_chem),missB=mean(miss_bmi),sigmaTrue=mean(sigmaT)) %>%
    rename(Model=out) %>% mutate(missC=signif(missC,3),missB=signif(missB,3)) %>% ungroup()
  o_biasB2 <- o_biasA %>% ungroup() %>% 
    mutate(fitTrueSecDevA=abs(fitTrueSecDev)) %>% filter(fitTrueSecDevA != 0) %>%
    summarize(val.min=min(fitTrueSecDevA),val.50=quantile(fitTrueSecDevA,0.50),
              val.75=quantile(fitTrueSecDevA,0.75),
              val.max=max(fitTrueSecDevA)) %>%
    merge(.,o_biasA) %>% mutate(fitTrueSecDevA=abs(fitTrueSecDev)) %>%
    mutate(quant=case_when(fitTrueSecDevA == 0 ~"Zero",
                           (fitTrueSecDevA > 0)&(fitTrueSecDevA<=val.50)~"<50",
                           (fitTrueSecDevA>val.50)&(fitTrueSecDevA<=val.75)~"50to75",
                           fitTrueSecDevA>val.75~">75"),
           vals=paste(signif(val.50,3),',',signif(val.75,3),',',signif(val.max,3))) %>%
    select(-val.50,-val.min,-val.max,-val.75)  %>%
    group_by(outs,type,quant,vals) %>% summarize(bias=mean(abs(bias))) %>%
    mutate(bias=signif(bias,2)) %>% ungroup(outs) %>%
    pivot_wider(names_from=c("quant"),values_from=bias) %>%
    mutate(Model = case_when(outs==0 ~ "Zero",
                             outs==1 ~ "One",
                             outs==2 ~ "Two",
                             outs==3 ~ "Three",
                             outs==4 ~ "Four",
                             outs==5 ~ "Five",
                             outs==6 ~ "Six",
                             outs==7 ~ "Seven",
                             outs==8 ~ "Eight",
                             TRUE ~ "Unknown")) %>% 
    arrange(outs) %>%
    select(-outs) %>%
    select(any_of(c("Model","type","Zero","<50","50to75",">75","vals"))) %>% 
    merge(.,routMiss,by=c("Model"))
  o_rb <- rbind(o_rmserB,o_biasB) %>% pivot_longer(cols=starts_with("val"),values_to="Values",names_to="Summary") %>%
    mutate(Values=signif(Values,3)) %>%
    pivot_wider(names_from=c("type"),values_from=Values) %>% 
    ungroup() %>%
    mutate(Model = case_when(outs==0 ~ "Zero",
                             outs==1 ~ "One",
                             outs==2 ~ "Two",
                             outs==3 ~ "Three",
                             outs==4 ~ "Four",
                             outs==5 ~ "Five",
                             outs==6 ~ "Six",
                             outs==7 ~ "Seven",
                             outs==8 ~ "Eight",
                             TRUE ~ "Unknown")) %>%
    arrange(outs) %>%
    select(Model,stat,Summary,everything()) %>% 
    merge(.,routMiss,by=c("Model")) %>% 
    arrange(stat) 
  
  
  dater <- Sys.Date()
  dater <- gsub("-","",paste(gsub("2023","",dater),gsub("(-[0-9]+-[0-9]+)","",dater),sep=""))
  write.csv(o_rb,file.path(dirLoc,'nonLinChk','1_Tables',"Sim",
                           paste('rsmeBias_',msPt*100,msBmi*100,sd2*100,"_",dater,'.csv',sep='')),row.names=FALSE)
  write.csv(o_biasB2,file.path(dirLoc,'nonLinChk','1_Tables',"Sim",
                               paste('Bias2d_',msPt*100,msBmi*100,sd2*100,"_",dater,'.csv',sep='')),row.names=FALSE)
      write.csv(o_allT,file.path(dirLoc,'nonLinChk','1_Tables',"Sim",
                            paste("rbiasA_",msPt*100,msBmi*100,sd2*100,"_",dater,'.csv',sep="")),row.names=FALSE)  
}
