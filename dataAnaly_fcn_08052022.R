

###data processing function
dataCh <- function(dataIn){
dataIn2 <- dataIn %>%
        mutate(bmi_imp_cen = (bmi_imp-25)/10,
       age_TH_cen = (age_TH-30)/10)
return(dataIn2)
}

###add glm function here
glmFnc <- function(dataInns,outs){ 
  regs <- lm(eval(parse(text=outs)) ~ measured_imp_log+bmi_imp_cen + sexC + age_TH_cen,data=dataInns)   
return(regs)
}

gamLinFnc <- function(dataInns,outs){
  require(mgcv)
  
  #out_gam1 <- gam(eval(parse(text=outs)) ~ s(measured_imp_log,bs='bs',k=7) + sexC + age_TH_cen + bmi_imp_cen,method="REML",data=dataInns)  
  #print(summary(out_gam1))
  out_gam1 <- gam(eval(parse(text=outs)) ~ s(measured_imp_log) + sexC + age_TH_cen + bmi_imp_cen,data=dataInns,method="REML")
  #print(summary(out_gam2))
  out_gam2 <- gam(eval(parse(text=outs)) ~ measured_imp_log + s(measured_imp_log, m=c(2,0)) + sexC + age_TH_cen + bmi_imp_cen, method="REML",data=dataInns)
  
  out_gam <- list(out_gam1,out_gam2)
  return(out_gam)
}

gamFnc <-function(dataInns2,outs2,knoter){
  require(mgcv)
  if (is.na(knoter)==TRUE){
    regsGAM <- gam(eval(parse(text=outs2)) ~ s(measured_imp_log,bs="bs",k=7) + sexC + age_TH_cen + bmi_imp_cen, data=dataInns2,method="REML")
  } else{
      regsGAM <- gam(eval(parse(text=outs2)) ~ s(measured_imp_log,bs="bs",k=7) + sexC + age_TH_cen + bmi_imp_cen, data=dataInns2,method="REML",knots=knoter)
  }
      return(regsGAM)
}

gamFcnKn <- function(dataInns3,ins,outs3){
  require(mgcv)
  tt <- 1
  ss <- 1
  minP <- 99
  for (ss1 in 1:10){
    for (ss2 in 1:10){
      minVal <- min(dataInns3[[ss1]][[ss2]]$measured_imp_log)
      if (minP > minVal){
      minP <- minVal
      tt <- ss1
      ss <- ss2
      }
    }
  }
  tmp_reg <- gam(eval(parse(text=outs3)) ~ s(measured_imp_log, bs='bs',k=7) + sexC + age_TH_cen + bmi_imp_cen, data=dataInns3[[tt]][[ss]],method="REML")
  regs_knots <- list(measured_imp_log=c(tmp_reg$smooth[[1]]$knots))
  return(regs_knots)
}
