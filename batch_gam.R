
fileDir <- "C:/Users/Michelle Mellers/Documents/ThyAnalysR"
if (file.exists(file.path(fileDir,"nonLinChk","2_Dump"))==FALSE){
        dir.create(file.path(fileDir,'nonLinChk','2_Dump'))
}
setwd(file.path(fileDir,"nonLinChk","2_Dump"))

########################################################################

library(tictoc)
library(futile.logger)
library(tryCatchLog)
library(parallel)
library(tidyverse)

source(file.path(fileDir,'nonLinChk','gamMissPar_03032023.R'),local=TRUE)
source(file.path(fileDir,'nonLinChk','dirCreate_gam_10142022.R'),local=TRUE)
source(file.path(fileDir,'nonLinChk','formatTabs_03032023.R'),local=TRUE)

## set-up tryCatchLog
options(keep.source=TRUE)
options("tryCatchLog.write.error.dump.file"=TRUE)
flog.appender(appender.file("my_analy.log"))
flog.threshold(INFO)


#library(tidyverse)

######################################################################
# 
# ##list of covariates I can possibly use for analysis
 covl <- c("raceC", "sexC","year_diffC2","service_pbdeC", "hor_state_reg3", "age_PBDE_C2", 
           "age_th_C2", "year_pbdeC2", "year_th_C","bmiC", "bmiC_imp","officer", "deploy", "byC")
# 
# #######################################################################
# 
 dirCreate()
# 

 ##########################################################################
## SIMULATION

 stdyType <- "Sim"
 #chemlist <- 1:60
 chemlist <- 1:156
# thylist <- c(0,1,6,7,8)
 ## change points at different spots
 thylist <- c(0,2,3,6,7,8)
 #thylist <- c(2,3,7)
 #0 linear; 6 polynomial; 8 exponential
 #thylist <- c(0,6,8)
 #thylist <- 5
 #thylist <- c(0,3)
 missType <- c("NM","CC","SI","NMI")


 #the no missing (NM) only works for the sim case
 
 ##########################################################################
 ## Sim 1, low miss for chem, high miss BMI
 ##(<miss chem %>,<miss bmi offset>,<sd of error>)
#msPt <- c(0.1,0.05,0.2)
#rout<-tryCatchLog(gamMiss(stdyType,missType,chemlist,thylist,msPt))
#writeoutGAM(rout,stdyType,msPt)
#rm(rout)

 
 tic("low chem high bmi")
 numCores <- detectCores()
 cl <- makeCluster(numCores)
 clusterExport(cl,"fileDir")
 clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','gamMissPar_03032023.R')))
 clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','errorSetup_03032023.R')))
 clusterEvalQ(cl,errorLoger())
 clusterExport(cl,"stdyType")
 clusterExport(cl,"chemlist")
 clusterExport(cl,"thylist")
 clusterExport(cl,"missType")
 
 
 msPt <- c(0.1,0.01,0.3)
 
 clusterExport(cl,"msPt")
 #chemlistPar <- ceiling(max(chemlist2)/numCores)
 router <- clusterEvalQ(cl,rout<-tryCatchLog(gamMissPar(stdyType,missType,chemlist,thylist,msPt)))
 stopCluster(cl)
 #clusterExport(cl,rout)
 #rout <- parLapply(cl,1:chemlistPar,f)
 rout <- router[[1]][[1]]
 o_rmserA <- router[[1]][[2]]
 o_allT <- router[[1]][[3]]
 maxC <- max(chemlist)
 for (i in(2:numCores)){
   routtmp <- router[[i]][[1]] %>% mutate(ins=as.character(as.numeric(ins)+maxC*(i-1)))
   rout <- rbind(rout,routtmp)
   o_rmersAtmp <- router[[i]][[2]] %>% mutate(simNum=simNum+maxC*(i-1))
   o_rmserA <- rbind(o_rmserA,o_rmersAtmp)
   o_allTmp <- router[[i]][[3]] %>% mutate(simNum=simNum+maxC*(i-1))
   o_allT <- rbind(o_allT,o_allTmp)
 }
 
 formatTabs(o_rmserA,o_allT,thylist,fileDir,msPt)
 writeoutGAM(rout,stdyType,msPt)
 toc()
 
 
##########################################################################
tic("GAM_Miss")
## Sim 2, low miss for chem, low miss BMI, avg sigma
##(<miss chem %>,<miss bmi offset>,<sd of error>)
 
 numCores <- detectCores()
 cl <- makeCluster(numCores)
 clusterExport(cl,"fileDir")
 clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','gamMissPar_03032023.R')))
 clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','errorSetup_03032023.R')))
 clusterEvalQ(cl,errorLoger())
 clusterExport(cl,"stdyType")
 clusterExport(cl,"chemlist")
 clusterExport(cl,"thylist")
 clusterExport(cl,"missType") 
 
 
msPt <- c(0.1,0.31,0.3)
clusterExport(cl,"msPt")

#chemlistPar <- ceiling(max(chemlist2)/numCores)
router <- clusterEvalQ(cl,rout<-tryCatchLog(gamMissPar(stdyType,missType,chemlist,thylist,msPt)))
stopCluster(cl)
#clusterExport(cl,rout)
#rout <- parLapply(cl,1:chemlistPar,f)
rout <- router[[1]][[1]]
o_rmserA <- router[[1]][[2]]
o_allT <- router[[1]][[3]]
maxC <- max(chemlist)
for (i in(2:numCores)){
  routtmp <- router[[i]][[1]] %>% mutate(ins=as.character(as.numeric(ins)+maxC*(i-1)))
  rout <- rbind(rout,routtmp)
  o_rmersAtmp <- router[[i]][[2]] %>% mutate(simNum=simNum+maxC*(i-1))
  o_rmserA <- rbind(o_rmserA,o_rmersAtmp)
   o_allTmp <- router[[i]][[3]] %>% mutate(simNum=simNum+maxC*(i-1))
  o_allT <- rbind(o_allT,o_allTmp)
 }

formatTabs(o_rmserA,o_allT,thylist,fileDir,msPt)
writeoutGAM(rout,stdyType,msPt)
#rm(rout)
#stopCluster(cl)

rm(routtmp,router,rout,o_rmserA,o_rmersAtmp,o_allTmp,o_allT)
toc()
##########################################################################

## Sim 3, low miss for chem, low miss BMI, high sigma
##(<miss chem %>,<miss bmi offset>,<sd of error>)
#msPt <- c(0.1,0.42,0.5)

#rout<-tryCatchLog(gamMiss(stdyType,missType,chemlist,thylist,msPt))
#writeoutGAM(rout,stdyType,msPt)
#rm(rout)

########################################################################## 
## Sim 4, high miss for chem, high miss BMI, avg sigma
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterExport(cl,"fileDir")
clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','gamMissPar_03032023.R')))
clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','errorSetup_03032023.R')))
clusterEvalQ(cl,errorLoger())
clusterExport(cl,"stdyType")
clusterExport(cl,"chemlist")
clusterExport(cl,"thylist")
clusterExport(cl,"missType")



##(<miss chem %>,<miss bmi offset>,<sd of error>
msPt <- c(0.45,0.01,0.3)
clusterExport(cl,"msPt")
#chemlistPar <- ceiling(max(chemlist2)/numCores)
router <- clusterEvalQ(cl,rout<-tryCatchLog(gamMissPar(stdyType,missType,chemlist,thylist,msPt)))
stopCluster(cl)
#clusterExport(cl,rout)
#rout <- parLapply(cl,1:chemlistPar,f)
rout <- router[[1]][[1]]
o_rmserA <- router[[1]][[2]]
o_allT <- router[[1]][[3]]
maxC <- max(chemlist)
for (i in(2:numCores)){
  routtmp <- router[[i]][[1]] %>% mutate(ins=as.character(as.numeric(ins)+maxC*(i-1)))
  rout <- rbind(rout,routtmp)
  o_rmersAtmp <- router[[i]][[2]] %>% mutate(simNum=simNum+maxC*(i-1))
  o_rmserA <- rbind(o_rmserA,o_rmersAtmp)
  o_allTmp <- router[[i]][[3]] %>% mutate(simNum=simNum+maxC*(i-1))
  o_allT <- rbind(o_allT,o_allTmp)
}

formatTabs(o_rmserA,o_allT,thylist,fileDir,msPt)
writeoutGAM(rout,stdyType,msPt)
#rm(rout)

#stopCluster(cl)

toc()
########################################################################## 

## Sim 5, high miss for chem, high miss BMI, high sigma
##(<miss chem %>,<miss bmi offset>,<sd of error>)
tic("high everything")
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterExport(cl,"fileDir")
clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','gamMissPar_03032023.R')))
clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','errorSetup_03032023.R')))
clusterEvalQ(cl,errorLoger())
clusterExport(cl,"stdyType")
clusterExport(cl,"chemlist")
clusterExport(cl,"thylist")
clusterExport(cl,"missType")


msPt <- c(0.45,0.01,0.6)

clusterExport(cl,"msPt")
#chemlistPar <- ceiling(max(chemlist2)/numCores)
router <- clusterEvalQ(cl,rout<-tryCatchLog(gamMissPar(stdyType,missType,chemlist,thylist,msPt)))
stopCluster(cl)
#clusterExport(cl,rout)
#rout <- parLapply(cl,1:chemlistPar,f)
rout <- router[[1]][[1]]
o_rmserA <- router[[1]][[2]]
o_allT <- router[[1]][[3]]
maxC <- max(chemlist)
for (i in(2:numCores)){
  routtmp <- router[[i]][[1]] %>% mutate(ins=as.character(as.numeric(ins)+maxC*(i-1)))
  rout <- rbind(rout,routtmp)
  o_rmersAtmp <- router[[i]][[2]] %>% mutate(simNum=simNum+maxC*(i-1))
  o_rmserA <- rbind(o_rmserA,o_rmersAtmp)
  o_allTmp <- router[[i]][[3]] %>% mutate(simNum=simNum+maxC*(i-1))
  o_allT <- rbind(o_allT,o_allTmp)
}

formatTabs(o_rmserA,o_allT,thylist,fileDir,msPt)
writeoutGAM(rout,stdyType,msPt)
toc()
#rout<-tryCatchLog(gamMiss(stdyType,missType,chemlist,thylist,msPt))
#writeoutGAM(rout,stdyType,msPt)
#rm(rout)

##########################################################################


## Sim 6, avg miss for chem, avg miss BMI, avg sigma
##(<miss chem %>,<miss bmi offset>,<sd of error>)
#msPt <- c(0.25,0.34,0.2)

#rout<-tryCatchLog(gamMiss(stdyType,missType,chemlist,thylist,msPt))
#writeoutGAM(rout,stdyType,msPt)
#rm(rout)

tic("average everything")
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterExport(cl,"fileDir")
clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','gamMissPar_03032023.R')))
clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','errorSetup_03032023.R')))
clusterEvalQ(cl,errorLoger())
clusterExport(cl,"stdyType")
clusterExport(cl,"chemlist")
clusterExport(cl,"thylist")
clusterExport(cl,"missType")


msPt <- c(0.25,0.24,0.3)

clusterExport(cl,"msPt")
#chemlistPar <- ceiling(max(chemlist2)/numCores)
router <- clusterEvalQ(cl,rout<-tryCatchLog(gamMissPar(stdyType,missType,chemlist,thylist,msPt)))
stopCluster(cl)
#clusterExport(cl,rout)
#rout <- parLapply(cl,1:chemlistPar,f)
rout <- router[[1]][[1]]
o_rmserA <- router[[1]][[2]]
o_allT <- router[[1]][[3]]
maxC <- max(chemlist)
for (i in(2:numCores)){
  routtmp <- router[[i]][[1]] %>% mutate(ins=as.character(as.numeric(ins)+maxC*(i-1)))
  rout <- rbind(rout,routtmp)
  o_rmersAtmp <- router[[i]][[2]] %>% mutate(simNum=simNum+maxC*(i-1))
  o_rmserA <- rbind(o_rmserA,o_rmersAtmp)
  o_allTmp <- router[[i]][[3]] %>% mutate(simNum=simNum+maxC*(i-1))
  o_allT <- rbind(o_allT,o_allTmp)
}

formatTabs(o_rmserA,o_allT,thylist,fileDir,msPt)
writeoutGAM(rout,stdyType,msPt)
toc()


##########################################################################

## Sim 7, high miss for chem, high miss BMI, higher sigma
##(<miss chem %>,<miss bmi offset>,<sd of error>)
#msPt <- c(0.45,0.22,0.75)

#rout<-tryCatchLog(gamMiss(stdyType,missType,chemlist,thylist,msPt))
#writeoutGAM(rout,stdyType,msPt)
#rm(rout)
##########################################################################


## Sim 8, low miss for chem, low miss BMI, higher sigma
##(<miss chem %>,<miss bmi offset>,<sd of error>)
#msPt <- c(0.1,0.42,0.75)

#rout<-tryCatchLog(gamMiss(stdyType,missType,chemlist,thylist,msPt))
#writeoutGAM(rout,stdyType,msPt)
#rm(rout)


##########################################################################

## Sim 9, high miss for chem, low miss BMI, avg sigma
##(<miss chem %>,<miss bmi offset>,<sd of error>)
#msPt <- c(0.45,0.59,0.2)

#rout<-tryCatchLog(gamMiss(stdyType,missType,chemlist,thylist,msPt))
#writeoutGAM(rout,stdyType,msPt)
#rm(rout)


tic("high chem low bmi")
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterExport(cl,"fileDir")
clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','gamMissPar_03032023.R')))
clusterEvalQ(cl,source(file.path(fileDir,'nonLinChk','errorSetup_03032023.R')))
clusterEvalQ(cl,errorLoger())
clusterExport(cl,"stdyType")
clusterExport(cl,"chemlist")
clusterExport(cl,"thylist")
clusterExport(cl,"missType")


msPt <- c(0.45,0.31,0.3)

clusterExport(cl,"msPt")
#chemlistPar <- ceiling(max(chemlist2)/numCores)
router <- clusterEvalQ(cl,rout<-tryCatchLog(gamMissPar(stdyType,missType,chemlist,thylist,msPt)))
stopCluster(cl)
#clusterExport(cl,rout)
#rout <- parLapply(cl,1:chemlistPar,f)
rout <- router[[1]][[1]]
o_rmserA <- router[[1]][[2]]
o_allT <- router[[1]][[3]]
maxC <- max(chemlist)
for (i in(2:numCores)){
  routtmp <- router[[i]][[1]] %>% mutate(ins=as.character(as.numeric(ins)+maxC*(i-1)))
  rout <- rbind(rout,routtmp)
  o_rmersAtmp <- router[[i]][[2]] %>% mutate(simNum=simNum+maxC*(i-1))
  o_rmserA <- rbind(o_rmserA,o_rmersAtmp)
  o_allTmp <- router[[i]][[3]] %>% mutate(simNum=simNum+maxC*(i-1))
  o_allT <- rbind(o_allT,o_allTmp)
}

formatTabs(o_rmserA,o_allT,thylist,fileDir,msPt)
writeoutGAM(rout,stdyType,msPt)
toc()


##########################################################################
