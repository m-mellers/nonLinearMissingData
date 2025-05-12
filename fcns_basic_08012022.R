
labsy <- function(outs){

    y_lab <- "LN Thyroid Level"
  
  return(y_lab)
}
#
simFcnNms <- function(simFcn){
  
  if (simFcn==0){
    outer9 <- "Zero"
  } else if (simFcn==1){
    outer9 <- "One"
  } else if (simFcn==2){
    outer9 <- "Two"
  } else if (simFcn==3){
    outer9 <- "Three"
  } else if (simFcn==4){
    outer9 <- "Four"
  } else if (simFcn==5){
    outer9 <- "Five"
  } else if (simFcn==6){
    outer9 <- "Six"
  } else if (simFcn==7){
    outer9 <- "Seven"
  } else {
    outer9 <- "Eight"
  } 
}


outsy <- function(outs){
 
    outs_y <- outs
  
  return(outs_y)
} 

addRout <- function(rout,glmout,sigRR,missTyper,ins,outs,bb){
  
  
  if (missTyper %in% c('CC','SI')){

    glmout2 <- summary(glmout)
    lens <- glmout2$df[[1]]
    rout$missType[bb:(bb+lens-1)] <- missTyper
    rout$out[bb:(bb+lens-1)] <- outs
    rout$ins[bb:(bb+lens-1)] <- ins
    rout$covs[bb:(bb+lens-1)] <- names(glmout2$coefficients[,1])
    rout$est[bb:(bb+lens-1)] <- glmout2$coefficients[,1]
    rout$se[bb:(bb+lens-1)] <- glmout2$coefficients[,2]
    rout$t[bb:(bb+lens-1)] <- glmout2$coefficients[,3]
    rout$df[bb:(bb+lens-1)] <- glmout2$df[1]
    rout$pv[bb:(bb+lens-1)] <- glmout2$coefficients[,4]
    rout$r2[bb:(bb+lens-1)] <- sigRR$r2
    rout$sigma[bb:(bb+lens-1)] <- sigRR$sigma
    rout$lo.95[bb:(bb+lens-1)] <- 999
    rout$hi.95[bb:(bb+lens-1)] <- 999
    rout$fmi[bb:(bb+lens-1)] <- 999
    rout$fmi_bet[bb:(bb+lens-1)] <- 999
    rout$fmi_within[bb:(bb+lens-1)] <- 999
    
  } else {
  
    lens <- nrow(glmout)
  rout$missType[bb:(bb+lens-1)] <- missTyper
  rout$out[bb:(bb+lens-1)] <- outs
  rout$ins[bb:(bb+lens-1)] <- ins
  rout$covs[bb:(bb+lens-1)] <- glmout$covs
  rout$est[bb:(bb+lens-1)] <- glmout$est
  rout$se[bb:(bb+lens-1)] <- glmout$se
  rout$t[bb:(bb+lens-1)] <- glmout$t
  rout$df[bb:(bb+lens-1)] <- glmout$df
  rout$pv[bb:(bb+lens-1)] <- glmout$`Pr(>|t|)`
  rout$r2[bb:(bb+lens-1)] <- sigRR$r2
  rout$sigma[bb:(bb+lens-1)] <- sigRR$sigma
  rout$lo.95[bb:(bb+lens-1)] <- glmout$`lo 95`
  rout$hi.95[bb:(bb+lens-1)] <- glmout$`hi 95`
  rout$fmi[bb:(bb+lens-1)] <- glmout$fmi
  rout$fmi_bet[bb:(bb+lens-1)] <- glmout$fmi_Betw
  rout$fmi_within[bb:(bb+lens-1)] <- glmout$fmi_Within
  
  }
  
  return(rout)
}



addRoutG <- function(rout,gamout,sigRR,missTyper,ins,outs,bb){
  
  gamout2 <- summary(gamout)
  #require(tidyverse)
  #require(data.table)
    #paremetric
    temp_gamp <- as.data.frame(gamout2$p.table) %>%
   # temp_gamp <- temp_gamp %>%
      mutate(covs=rownames(.))
    lens_p <- nrow(temp_gamp)
    
    #smoother
    temp_gams <- as.data.frame(gamout2$s.table) %>%
 #   temp_gams <- temp_gams %>%
      mutate(covs=rownames(.))
    lens_s <- nrow(temp_gams)
    
    rout$missType[bb:(bb+lens_p-1)] <- missTyper
    rout$out[bb:(bb+lens_p-1)] <- outs
    rout$ins[bb:(bb+lens_p-1)] <- ins
    rout$covs[bb:(bb+lens_p-1)] <- temp_gamp$covs
    rout$est[bb:(bb+lens_p-1)] <- temp_gamp$Estimate
    rout$se[bb:(bb+lens_p-1)] <- temp_gamp$`Std. Error`
    rout$t[bb:(bb+lens_p-1)] <- temp_gamp$`t value`
    rout$df[bb:(bb+lens_p-1)] <- temp_gams$edf
    rout$pv[bb:(bb+lens_p-1)] <- temp_gamp$`Pr(>|t|)`
    rout$r2[bb:(bb+lens_p-1)] <- sigRR$r2
    rout$sigma[bb:(bb+lens_p-1)] <- sigRR$sigma
    
    bb <- bb + lens_p
    
    rout$missType[bb:(bb+lens_s-1)] <- missTyper
    rout$out[bb:(bb+lens_s-1)] <- outs
    rout$ins[bb:(bb+lens_s-1)] <- ins
    rout$covs[bb:(bb+lens_s-1)] <- temp_gams$covs
    rout$est[bb:(bb+lens_s-1)] <- NA
    rout$se[bb:(bb+lens_s-1)] <- NA
    rout$t[bb:(bb+lens_s-1)] <- temp_gams$`F`
    rout$df[bb:(bb+lens_s-1)] <- temp_gams$edf
    rout$pv[bb:(bb+lens_s-1)] <- temp_gams$`p-value`
    rout$r2[bb:(bb+lens_s-1)] <- sigRR$r2
    rout$sigma[bb:(bb+lens_s-1)] <- sigRR$sigma

  return(rout)
}



reLabThy <- function(out){

    out2 <- out
  
  return(out2)
}

missLab <- function(missType){
  if (missType == 'CC'){
    missType2 <- 'Complete Case'
  } else if (missType == 'SI'){
    missType2 <- 'Substitution'
  } else if (missType == 'NMI'){
    missType2 <- 'Nested Multiple Imputation'
  } else if (missType == 'NM'){
    missType2 <- "No Missing"
  }
  return(missType2)
}


aotAdds <- function(file1,file2){
  
  bb1 <- nrow(file1)
  bb2 <- nrow(file2)
  nlen <- bb1+bb2
  
  routs <- data.frame(missType=character(nlen),out=character(nlen),ins=character(nlen),
                      covs=character(nlen),est=double(nlen),se=double(nlen),t=double(nlen),
                      df=double(nlen),pv=double(nlen),lo.95=double(nlen),hi.95=double(nlen),
                      fmi=double(nlen),fmi_bet=double(nlen),fmi_within=double(nlen))
  
  bb <- 1
  routs[bb:(bb+bb1-1)] <- file1
  bb <- bb+bb1
  routs[bb:(bb+bb2-1)] <- file2
  
  return(routs)
  
}

sigR2 <- function(regss,msT){
  
  sigR <- data.frame(sigma=0,r2=0)
  if (msT == 'NMI'){
    
    for (s in 1:length(regss)){
      for (t in 1:length(regss[[1]])){
        tmp <- summary(regss[[s]][[t]])
        sigR$r2 <- sigR$r2 + tmp$adj.r.squared
        sigR$sigma <- sigR$sigma + tmp$sigma^2
      }}
    sigR <- sigR %>% mutate(r2=r2/100, sigma=sqrt(sigma/100))
  } else {
    tmp <- summary(regss)
    sigR <- sigR %>% mutate(sigma=tmp$sigma,r2=tmp$adj.r.squared)
  }
  
  sigR <- sigR %>% mutate(sigma=signif(sigma,3),r2=signif(r2,3))
  return(sigR)
}


sigR2G <- function(regss,msT){
  
  sigR <- data.frame(sigma=0,r2=0)
  if (msT == 'NMI'){
    
    for (s in 1:length(regss)){
      for (t in 1:length(regss[[1]])){
        tmp <- summary(regss[[s]][[t]])
        sigR$r2 <- sigR$r2 + tmp$r.sq
        sigR$sigma <- sigR$sigma + sigma(regss[[s]][[t]])^2
      }}
    sigR <- sigR %>% mutate(r2=r2/100, sigma=sqrt(sigma/100))
  } else {
    tmp <- summary(regss)
    sigR <- sigR %>% mutate(sigma=sigma(regss),r2=tmp$r.sq)
  }
  
  sigR <- sigR %>% mutate(sigma=signif(sigma,3),r2=signif(r2,3))
  return(sigR)
}


gamlabsy <- function(outs){

    y_lab <- "LN Thyroid Level (Sim)"
  
  return(y_lab)
}
#

gamoutsy <- function(outs){

    outs_y <- outs
  
  return(outs_y)
}  


taber <- function(outgm,outdat,ins,out,s,rs){
  
  tmp <- summary(outgm)
  nleng <- length(tmp$p.coeff)-1
  outdat$chem[s:(nleng+s)] <- ins
  outdat$thyroid[s:(nleng+s)] <- out
  outdat$coeff[s:(nleng+s)] <- names(tmp$p.coeff)
  outdat$p.coeff[s:(nleng+s)] <- tmp$p.coeff
  outdat$se[s:(nleng+s)] <- tmp$se[1:(nleng+1)]
  outdat$p.pv[s:(nleng+s)] <- tmp$p.pv
  outdat$dev.expl[s:(nleng+s)] <- tmp$dev.expl
  outdat$dispersion[s:(nleng+s)] <- tmp$dispersion
  outdat$r.sq[s:(nleng+s)] <- tmp$r.sq
  outdat$s.pv[s:(nleng+s)] <- tmp$s.pv
  outdat$refdf[s:(nleng+s)] <- tmp$s.table[2]
  outdat$edf[s:(nleng+s)] <- tmp$edf
  outdat$typeSpl[s:(nleng+s)] <- rs
  
  return(outdat) 
  
}

addO <- function(oold,onew,sp){
  
  nl <- nrow(onew)-1
  oold$measured_imp_log[sp:(nl+sp)] <- onew$measured_imp_log
  oold$se[sp:(nl+sp)] <- onew$se
  oold$fit[sp:(nl+sp)] <- onew$fit
  oold$ins[sp:(nl+sp)] <- onew$ins
  oold$outs[sp:(nl+sp)] <- onew$outs
  oold$rt[sp:(nl+sp)] <- onew$rt
  oold$typ[sp:(nl+sp)] <- onew$typ

    return(oold)
}

addOG <- function(oold,onew,sp){
  
  nl <- nrow(onew)-1
  oold$measured_imp_log[sp:(nl+sp)] <- onew$measured_imp_log
  oold$se[sp:(nl+sp)] <- onew$se
  oold$fit[sp:(nl+sp)] <- onew$fit
  oold$measTrue[sp:(nl+sp)] <- onew$measTrue
  oold$fitTrue[sp:(nl+sp)] <- onew$fitTrue
  oold$fitTrueSecDev[sp:(nl+sp)] <- onew$fitTrueSecDev
  oold$outs[sp:(nl+sp)] <- onew$outs
  oold$simNum[sp:(nl+sp)] <- onew$simNum
  oold$ID[sp:(nl+sp)] <- onew$ID

  return(oold)
}


addOall <- function(oa3,oa2,sz){
  
  nl2 <- nrow(oa2)-1
  oa3$outs[sz:(nl2+sz)] <- oa2$outs
  oa3$rt[sz:(nl2+sz)] <- oa2$rt
  oa3$medians[sz:(nl2+sz)] <- oa2$medians
  oa3$mins[sz:(nl2+sz)] <- oa2$mins
  oa3$maxs[sz:(nl2+sz)] <- oa2$maxs
  oa3$f25[sz:(nl2+sz)] <- oa2$f25
  oa3$f75[sz:(nl2+sz)] <- oa2$f75
  return(oa3)
}

classNms <- function(chemss){
  
    groups <- "Sim"

  return(groups)
}
