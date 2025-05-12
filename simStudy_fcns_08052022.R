
### list of functions for the simulation study

formData <- function(missPct,sd2,LODloc,fcnT){
#missPct is the percent missing
#LODloc is the location along the outcome of the limit of detection
  

require(simstudy)
require(data.table)
nsub=741
#sd2 = 1.5
#sd2 = 0.18
foot=data.frame(x=numeric(nsub),y=numeric(nsub),e=numeric(nsub),yNoisy=numeric(nsub),
               id=numeric(nsub))
q <- 1
#foot$x=seq(0,5,by=(5)/740) # x is a grid x...a known regular grid
#foot$x <- runif(n=741,min=0,max=5)
#foot$x <- rnorm(n=741,mean=3,sd=1.5) #%>%0 as.data.frame(.) %>% rename(tester=".")
dd <- data.frame(id=1:741,meas=double(741))
#original variance is 0.88
#def <- defData(varname="meas",dist="normal",formula=2.5,variance=0.88)
def <- defData(varname="meas",dist="normal",formula=2.5,variance=1.2)

#dd <- genData(nsub,def)
dd3 <- genData(nsub,def)
dd_tmp <- dd3 %>% filter((meas>=0)&(meas<=5))
dd$meas[q:(nrow(dd_tmp))] <- dd_tmp$meas
q <- q+nrow(dd_tmp)
for (var in 1:500){
dd2 <- nrow(dd3[((dd3$meas<0)|(dd3$meas>5)),])
if (dd2 == 0){
  break
}
dd3 <- genData(dd2,def)
dd_tmp <- dd3 %>% filter((meas>=0)&(meas<=5))
if (nrow(dd_tmp) > 0){
dd$meas[q:(q+nrow(dd_tmp)-1)] <- dd_tmp$meas
q <- q+nrow(dd_tmp)
}
}

foot$x <- round(dd$meas,digits=4)

if (missPct != 0){
foot <- foot %>% arrange(x) %>%
  mutate(lod = LODloc,
         ms = round(missPct*741),
         rn = row_number()) %>%
  rename(xa=x) %>%
  mutate(xloda = case_when(rn <= ms ~ NA_real_,
                           TRUE ~ xa))
foot$offs <- foot[(foot$rn==(foot$ms+1)),names(foot)=='xa']      
foot <- foot %>% mutate(xlod=xloda-offs+lod,
                        x=xa-offs+lod) %>% 
  select(-xloda,-ms,-rn,-xa)

} else {
 foot <- foot %>% mutate(offs = LODloc,xlod=x+LODloc,lod=LODloc,x=x+LODloc)
}

foot <- foot %>% as.data.table()
def2 <- defDataAdd(varname="bmi_cenB",dist="normal",formula=0,variance=0.089)
#def2 <- defDataAdd(varname="bmi_cenB",dist="normal",formula=0,variance=0.089)
def2 <- defDataAdd(def2,varname="age_TH_cenB",dist="normal",formula="-0.517 + 0.559*x+0.245*bmi_cenB",variance=0.293)
#def2 <- defDataAdd(def2,varname="age_TH_cenB",dist="normal",formula="-0.517 + 0.559*x+0.245*bmi_cenB",variance=2.93)
def2 <- defDataAdd(def2,varname="sexC",dist="binary",formula="-0.126-0.288*x-1.98*bmi_cenB",link="logit")
def2 <- defDataAdd(def2,varname="raceC",dist="binary",formula="-2.62+0.619*x-0.361*age_TH_cenB+0.915*sexC",link="logit")
foot <- addColumns(def2,foot)
foot <- foot %>% rowid_to_column(.,"ID")
foot <- foot %>%
          mutate(bmi = bmi_cenB*10+25,
                 age_TH = age_TH_cenB*10+30)
  
footM <- foot %>% group_by(sexC) %>%
                  summarize(bmi_cenBM=mean(bmi_cenB),age_TH_cenBM=mean(age_TH_cenB))
               

#tbmi <- analys5a %>% summarise(means=mean(bmi,na.rm=TRUE),sd=sd(bmi,na.rm=TRUE))
#pm <- tst %>% as.data.frame() %>% ggplot(.,aes(x=tester)) + geom_density()
#add a random x

##limit of detection comes in through the x-value rather than the y-value
#foot$y=3-0.00007*exp(foot$x)

offer <- 0.1
if (fcnT == 0) {
  
  
  b0 <- 2.13
  #b1 is for chemicals
 # b1 <- -0.0232
  b1 <- -0.05
  #b1 <- -0.232
 # b1 <- -2.32
 b2 <- 0.0752
  #b2 <- 0.752
  b3 <- -0.0243
 # b3 <- -0.243
 # b3 <- 0
  #b4 is for bmi
  b4 <- -0.0401
 # b4 <- -0.401
#  b4 <- -4.01
  
  foot$y <-  b0+b1*foot$x+b2*foot$sexC+b3*foot$age_TH_cenB+b4*foot$bmi_cenB
  foot$yBM <- b0+b1*foot$x+as.numeric(b2*foot$sexC+b3*footM[which(footM$sexC==0),names(footM)=="age_TH_cenBM"]+b4*footM[which(footM$sexC==0),names(footM)=="bmi_cenBM"])
  foot$yBF <-  b0+b1*foot$x+as.numeric(b2*foot$sexC+b3*footM[which(footM$sexC==1),names(footM)=="age_TH_cenBM"]+b4*footM[which(footM$sexC==1),names(footM)=="bmi_cenBM"])
  pctW <- signif(sum(foot$sexC/nsub),2)
  foot$yB <- (1-pctW)*foot$yBM+pctW*foot$yBF
  foot$yBd1 <- b1
  foot$yBd2 <- 0
  
} else if (fcnT == 1) {
  
foot <- foot %>% mutate(xz=case_when(x<offer~offer,
                                     TRUE~x))
  b0=2.11
  b1=0.080
  b2=-0.030
  b3=-0.033
  #bmi and chem
  b4=-0.001
  b5=-0.004
  b6=-0.001
  b7=-0.001
  #b4=-0.01
  #b5=-0.04
  #b6=-0.01
  #b7=-0.01
  foot$y <- b0+b1*foot$sexC+b2*foot$age_TH_cenB+b3*foot$bmi_cenB+b4*foot$xz+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4
  #foot$y = 0.01*(foot$x+1)-0.01*(foot$x)^2+0.01*foot$bmi_cenB+1
  foot$yBM <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==0),names(footM)=="age_TH_cenBM"]+
                                        b3*footM[which(footM$sexC==0),names(footM)=="bmi_cenBM"])+b4*foot$xz+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4
  foot$yBF <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==1),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==1),names(footM)=="bmi_cenBM"])+b4*foot$xz+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4                                     
  pctW <- signif(sum(foot$sexC/nsub),2)
  foot$yB <- (1-pctW)*foot$yBM+pctW*foot$yBF
  foot <- foot %>% mutate(yBd1 = b4+2*b5*xz+3*b6*xz^2+4*b7*xz^3,
                          yBd2 = case_when(xz==offer~0,
                                           TRUE~2*b5+6*b6*xz+12*b7*xz^2)) %>%
    select(-xz)

 # foot$yBd2 <- if (foot$xz==0.75){
   
  #foot <- foot %>% select(-xz)
} else if (fcnT ==2){
  
  offer2 <- offer+0.2
  foot <- foot %>% mutate(xz=case_when(x<offer2~offer2,
                                       TRUE~x))
  b0=2.11
  b1=0.080
  b2=-0.030
  b3=-0.033
  #bmi and chem
  b4=-0.001
  b5=-0.004
  b6=-0.001
  b7=-0.001
  # b4=-0.05
  #  b5=-0.04
  # b6=-0.01
  #  b7=-0.01
  foot$y <- b0+b1*foot$sexC+b2*foot$age_TH_cenB+b3*foot$bmi_cenB+b4*foot$x+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4
  #foot$y = 0.01*(foot$x+1)-0.01*(foot$x)^2+0.01*foot$bmi_cenB+1
  foot$yBM <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==0),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==0),names(footM)=="bmi_cenBM"])+b4*foot$x+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4
  foot$yBF <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==1),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==1),names(footM)=="bmi_cenBM"])+b4*foot$x+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4                                     
  pctW <- signif(sum(foot$sexC/nsub),2)
  foot$yB <- (1-pctW)*foot$yBM+pctW*foot$yBF
  
  foot <- foot %>% mutate(yBd1=b4+2*b5*foot$xz+3*b6*xz^2+4*b7*xz^3,
                          yBd2=case_when(xz==offer2~0,
                                         TRUE~2*b5+6*b6*xz+12*b7*xz^2)) %>% select(-xz) 
  
} else if (fcnT == 3){
  
  offer3 <- offer + 0.3
  foot <- foot %>% mutate(xz=case_when(x<offer3~offer3,
                                       TRUE~x))
  b0=2.11
  b1=0.080
  b2=-0.030
  b3=-0.033
  #bmi and chem
  b4=-0.001
  b5=-0.004
  b6=-0.001
  b7=-0.001
  # b4=-0.05
  #  b5=-0.04
  # b6=-0.01
  #  b7=-0.01
  foot$y <- b0+b1*foot$sexC+b2*foot$age_TH_cenB+b3*foot$bmi_cenB+b4*foot$x+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4
  #foot$y = 0.01*(foot$x+1)-0.01*(foot$x)^2+0.01*foot$bmi_cenB+1
  foot$yBM <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==0),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==0),names(footM)=="bmi_cenBM"])+b4*foot$x+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4
  foot$yBF <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==1),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==1),names(footM)=="bmi_cenBM"])+b4*foot$x+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4                                     
  pctW <- signif(sum(foot$sexC/nsub),2)
  foot$yB <- (1-pctW)*foot$yBM+pctW*foot$yBF
  
  foot <- foot %>% mutate(yBd1=b4+2*b5*foot$xz+3*b6*xz^2+4*b7*xz^3,
                          yBd2=case_when(xz==offer3~0,
                                         TRUE~2*b5+6*b6*xz+12*b7*xz^2)) %>% select(-xz) 
  
} else if (fcnT == 4){  

  foot <- foot %>% mutate(xz=case_when(x<offer~offer,
                                       TRUE~x))
  b0=2.11
  b1=0.078
  b2=-0.022
  b3=-0.032
  ##bmi and chems
  #b4=-0.001
  #b5=-0.004
  #b6=-0.0005
  b4=-0.01
  b5=-0.04
  b6=-0.005

  
  foot$y <- b0+b1*foot$sexC+b2*foot$age_TH_cenB+b3*foot$bmi_cenB+b4*(foot$xz)+b5*foot$xz^2+b6*foot$xz^3
  #foot$y = 0.01*(foot$x+1)-0.01*(foot$x)^2+0.01*foot$bmi_cenB+1
  foot$yBM <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==0),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==0),names(footM)=="bmi_cenBM"])+b4*foot$xz+b5*foot$xz^2+b6*foot$xz^3
  foot$yBF <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==1),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==1),names(footM)=="bmi_cenBM"])+b4*foot$xz+b5*foot$xz^2+b6*foot$xz^3                                     
  pctW <- signif(sum(foot$sexC/nsub),2)
  foot$yB <- (1-pctW)*foot$yBM+pctW*foot$yBF
 
  foot <- foot %>% mutate(yBd1=b4+2*b5*xz+3*b6*xz^2,
                          yBd2=case_when(xz==offer~0,
                                         TRUE~2*b5+6*b6*xz)) %>%
    select(-xz)
  
  
} else if (fcnT == 5){

  foot <- foot %>% mutate(xz=case_when(x<offer~offer,
                                       TRUE~x))
  b0=2.11
  b1=0.077
  b2=-0.027
  b3=-0.036
  ##bmi and chems
  #b4=-0.001
  #b5=-0.004
  #b6=-0.001
  b4=-0.01
  b5=-0.04
  b6=-0.01
  a1=0
  #a2=-2.0
  a2=-1.0
  a3=0
  foot$y <- b0+b1*foot$sexC+b2*foot$age_TH_cenB+b3*foot$bmi_cenB+b4*(foot$xz+a1)+b5*(foot$xz+a2)^2+b6*(foot$xz+a3)^3
  #foot$y = 0.01*(foot$x+1)-0.01*(foot$x)^2+0.01*foot$bmi_cenB+1
  foot$yBM <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==0),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==0),names(footM)=="bmi_cenBM"])+b4*(foot$xz+a1)+b5*(foot$xz+a2)^2+b6*(foot$xz+a3)^3
  foot$yBF <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==1),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==1),names(footM)=="bmi_cenBM"])+b4*(foot$xz+a1)+b5*(foot$xz+a2)^2+b6*(foot$xz+a3)^3                                     
  pctW <- signif(sum(foot$sexC/nsub),2)
  foot$yB <- (1-pctW)*foot$yBM+pctW*foot$yBF
  foot <- foot %>% mutate(yBd1 = b4+2*b5*(xz+a2)+3*b6*(xz+a3)^2,
                               yBd2 = case_when(xz==offer~0,
                                                TRUE~2*b5+6*b6*(xz+a3))) %>% 
    select(-xz)
  
} else if (fcnT == 6){  

  foot <- foot %>% mutate(xz=x)
  b0=2.11
  b1=0.076
  b2=-0.018
  b3=-0.034
  ##bmi and chems
 b4=-0.03
  b5=-0.005
  b6=0
  #b4=-0.3
  #b5=-0.05
  #b6=0
  a1=0
  a2=-2.0
  a3=0
  foot$y <- b0+b1*foot$sexC+b2*foot$age_TH_cenB+b3*foot$bmi_cenB+b4*(foot$xz+a1)+b5*(foot$xz+a2)^2+b6*(foot$xz+a3)^3
  #foot$y = 0.01*(foot$x+1)-0.01*(foot$x)^2+0.01*foot$bmi_cenB+1
  foot$yBM <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==0),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==0),names(footM)=="bmi_cenBM"])+b4*(foot$xz+a1)+b5*(foot$xz+a2)^2+b6*(foot$xz+a3)^3
  foot$yBF <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==1),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==1),names(footM)=="bmi_cenBM"])+b4*(foot$xz+a1)+b5*(foot$xz+a2)^2+b6*(foot$xz+a3)^3                                     
  pctW <- signif(sum(foot$sexC/nsub),2)
  foot$yB <- (1-pctW)*foot$yBM+pctW*foot$yBF

  
  foot <- foot %>% mutate(yBd1=b4+b5*(xz+a2)+3*b6*(xz+a3),
                          yBd2=case_when(xz==offer~0,
                                         TRUE~2*b5+6*b6*(xz+a3))) %>%
      select(-xz)
  
} else if (fcnT == 7){  

  foot <- foot %>% mutate(xz=case_when(x<offer~offer,
                                       TRUE~x))
  b0=2.11
  b1=0.080
  b2=-0.030
  b3=-0.033
  #bmi and chem
  b4=-0.001
  b5=-0.004
  b6=-0.001
  b7=-0.001
 # b4=-0.05
#  b5=-0.04
 # b6=-0.01
#  b7=-0.01
  foot$y <- b0+b1*foot$sexC+b2*foot$age_TH_cenB+b3*foot$bmi_cenB+b4*foot$x+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4
  #foot$y = 0.01*(foot$x+1)-0.01*(foot$x)^2+0.01*foot$bmi_cenB+1
  foot$yBM <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==0),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==0),names(footM)=="bmi_cenBM"])+b4*foot$x+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4
  foot$yBF <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==1),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==1),names(footM)=="bmi_cenBM"])+b4*foot$x+b5*foot$xz^2+b6*foot$xz^3+b7*foot$xz^4                                     
  pctW <- signif(sum(foot$sexC/nsub),2)
  foot$yB <- (1-pctW)*foot$yBM+pctW*foot$yBF
 
  foot <- foot %>% mutate(yBd1=b4+2*b5*foot$xz+3*b6*xz^2+4*b7*xz^3,
                          yBd2=case_when(xz==offer~0,
                                         TRUE~2*b5+6*b6*xz+12*b7*xz^2)) %>% select(-xz) 
  
} else if (fcnT == 8) {
  
  #foot$y=2-0.0025*exp(foot$x)+0.01*foot$bmi_cenB
  b0=2.11
  b1=0.080
  b2=-0.030
  #bmi and chems
  b3=-0.033
  b4=-0.001
  #b3=-0.33
  #b4=-0.01
  foot$y <- b0+b1*foot$sexC+b2*foot$age_TH_cenB+b3*foot$bmi_cenB+b4*exp(foot$x)
  foot$yBM <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==0),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==0),names(footM)=="bmi_cenBM"])+b4*exp(foot$x)
  foot$yBF <- b0+as.numeric(b1*foot$sexC+b2*footM[which(footM$sexC==1),names(footM)=="age_TH_cenBM"]+
                              b3*footM[which(footM$sexC==1),names(footM)=="bmi_cenBM"])+b4*exp(foot$x)                                      
  pctW <- signif(sum(foot$sexC/nsub),2)
  foot$yB <- (1-pctW)*foot$yBM+pctW*foot$yBF
  foot$yBd1 <- b4*exp(foot$x)
  foot$yBd2 <- b4*exp(foot$x)
  
} else {
  foot$y <- 0
  foot$yBM <- 0
  foot$yBF <- 0
  foot$yB <- 0
  foot$yBd1 <- 0
  foot$yBd2 <- 0
}

foot$e=rnorm(nsub,mean=0,sd=sd2)
foot <- foot %>% mutate(yNoisy=y+e)
foot$id <- 1:nrow(foot)

foot <- foot %>% select(-yBM,-yBF) %>% rename(measured_true=x,measured_imp_log=xlod,lod_log=lod) %>%
  mutate(LOD=lod_log,sexC=as.factor(sexC)) %>% as.data.table()
return(foot)
}

genM <- function(foos,msBmi){

  foos <- foos %>% mutate(bmiComp=bmi)
  fms <- paste(msBmi,"+age_TH*0.0285-0.015*raceC",sep="")
#  fms <- "age_TH*0.0285-0.015*raceC"
defM <- defMiss(varname="bmi",formula=fms,logit.link=FALSE)
#defM <- defMiss(varname="bmi",formula="0.45+age_TH_cenB*0.25-0.15*raceC",logit.link=FALSE)
  missM <- genMiss(foos,defM,idvars="id")
dtObs <- genObs(foos,missM,idvars="id")
return(dtObs)
}



rmsers <- function(dats){
  rmse <- dats %>% group_by(trial,sdm,lod) %>%
    filter(is.na(yGAM)!=TRUE) %>%
    mutate(diffs=(yGAM-yB)^2) %>%
    summarize(sums=sum(diffs)) %>%
    mutate(rmse=signif(sqrt(sums/741),digits=3)) %>% select(-sums)
 # rmse <- as.numeric(rmse)
  return(rmse)
}

bias2 <- function(daty){
  biasy <- daty %>% group_by(lod,sdm,x) %>%
          mutate(diffs=(yGAM-yB)) %>%
          summarize(sums=sum(diffs),trialm=max(trial)) %>%
          mutate(bias=signif(sums/trialm,digits=3))
return(biasy)
}

addVars <- function(tmp_exp,dataInn,m){
  
  nleng <- nrow(tmp_exp)-1
  
 # dataInn$X[m:(nleng+m)] <- tmp_exp$X
  dataInn$id[m:(nleng+m)] <- tmp_exp$id
  dataInn$meas_imp1[m:(nleng+m)] <- tmp_exp$meas_imp1
  dataInn$measured_true[m:(nleng+m)] <- tmp_exp$measured_true
  dataInn$y[m:(nleng+m)] <- tmp_exp$y
  dataInn$e[m:(nleng+m)] <- tmp_exp$e
  dataInn$yNoisy[m:(nleng+m)] <- tmp_exp$yNoisy
  dataInn$bmi[m:(nleng+m)] <- tmp_exp$bmi
  dataInn$yB[m:(nleng+m)] <- tmp_exp$yB
  dataInn$sex[m:(nleng+m)] <- tmp_exp$sex
  dataInn$race[m:(nleng+m)] <- tmp_exp$race
  dataInn$age[m:(nleng+m)] <- tmp_exp$age
  dataInn$lod[m:(nleng+m)] <- tmp_exp$lod
  dataInn$measured_imp_log[m:(nleng+m)] <- tmp_exp$measured_imp_log
  dataInn$yGAM[m:(nleng+m)] <- tmp_exp$yGAM
  dataInn$se.fit[m:(nleng+m)] <- tmp_exp$se.fit
  dataInn$trial[m:(nleng+m)] <- tmp_exp$trial
  dataInn$sdm[m:(nleng+m)] <- tmp_exp$sdm
  
  return(dataInn)
}

