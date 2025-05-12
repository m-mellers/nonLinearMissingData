
dirCreate <- function(){
ylst <- c("Zero","One","Two","Three","Four","Five","Six","Seven","Eight")
#xlst <- "Sim"

fileDir <- "C:/Users/Michelle Mellers/Documents/ThyAnalysR"
source(file.path(fileDir,"nonLinChk","fcns_basic_08012022.R"))
#dir.create(file.path(fileDir,'nonLinChk','2_Dump'))
## groups


## table file structure
if (file.exists(file.path(fileDir,'nonLinChk',"1_Tables"))==FALSE){
  dir.create(file.path(fileDir,'nonLinChk',"1_Tables"))
}
if (file.exists(file.path(fileDir,'nonLinChk',"1_Tables","Sim"))==FALSE){
dir.create(file.path(fileDir,'nonLinChk',"1_Tables","Sim"))
}
if (file.exists(file.path(fileDir,'nonLinChk','1_Sim'))==FALSE){
dir.create(file.path(fileDir,'nonLinChk',"1_Sim"))
}  
#dir.create(file.path(fileDir,'MissingDataGLM','Sim','Sim'))
for (k in 1:length(ylst)){
  ylsts <- ylst[[k]]
  if (file.exists(file.path(fileDir,'nonLinChk',"1_Sim",ylsts))==FALSE){
  dir.create(file.path(fileDir,"nonLinChk","1_Sim",ylsts))
  }
}
}

upDateDirGAM <- function(lst,typeReq){
  chemM2 <- lst[[1]]*100
  bmiM2 <- lst[[2]]*100
  sigma2 <- lst[[3]]*100
  fileDir <- "C:/Users/Michelle Mellers/Documents/ThyAnalysR/nonLinChk"
  dater <- Sys.Date()
  dater <- gsub("-","",paste(gsub("2023","",dater),gsub("(-[0-9]+-[0-9]+)","",dater),sep=""))
  newDir <- paste(typeReq,chemM2,bmiM2,sigma2,"_",dater,sep="")
    if (file.exists(file.path(fileDir,'1_Sim',newDir))==FALSE){
  dir.create(file.path(fileDir,"1_Sim",newDir))
    }
  fls <- list.files(path=file.path(fileDir,"1_Sim",typeReq))
  
  if (length(fls) != 0){
  for (s in 1:length(fls)){
    from <- file.path(fileDir,"1_Sim",typeReq,fls[[s]])
    to <- file.path(fileDir,"1_Sim",newDir,fls[[s]])
    file.rename(from=from,to=to)
  }
  }
}
  
writeoutGAM <- function(rout,sdyType,lst){
  fileDirs <- "C:/Users/Michelle Mellers/Documents/ThyAnalysR"
 if (sdyType=='Sim'){
    #rout$out = simFcnNms(as.numeric(rout$out))
    syst <- unique(rout$out)
    for (mn in 1:length(syst)){
      typReq <- syst[[mn]]
    upDateDirGAM(lst,typReq) 
  }
    chemM <- lst[[1]]*100
    bmiM <- lst[[2]]*100
    sigma2 <- lst[[3]]*100
    dater <- Sys.Date()
    dater <- gsub("-","",paste(gsub("2023","",dater),gsub("(-[0-9]+-[0-9]+)","",dater),sep=""))
    write.csv(rout,file.path(fileDirs,'nonLinChk','1_Tables','Sim',paste('taballsim_',chemM,bmiM,sigma2,'_',dater,'.csv',sep="")),row.names=FALSE)
    
             
  }
}