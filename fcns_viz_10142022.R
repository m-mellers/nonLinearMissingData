###10/14 update for non-linear scenarios


vizMiss <- function(msTyp,predicts,ins,outs,msTp3){
  
y_lab <- labsy(outs)
outs_y <- outsy(outs)

iners <- gsub("_","-",ins)   
x_lab <- paste("LN ",iners," (ng/g)",sep="")
msTyp2 <- missLab(msTyp)
##ggplot
#average dataset
 
    pglm <- ggplot(predicts,aes(measured_imp_log,fits)) + 
      geom_ribbon(aes(ymin=lower,ymax=upper),fill='grey70',alpha=0.6) +
      geom_line(size=1) + 
      labs(x =x_lab,y=y_lab,
           title=paste("Missing Data Type: ",msTyp2,"\nChemical: ",iners,", Thyroid: ",outs_y,sep="")) +
      theme_minimal() +
      theme(axis.text.x=element_text(size=12,face='bold'),plot.title=element_text(size=14,hjust=0.5,face='bold'), axis.text.y=element_text(size=12,face='bold'),axis.title.x=element_text(size=12,face='bold'),
            axis.title.y=element_text(size=12,face='bold'),legend.position="bottom",
            legend.title=element_text(size=12,face='bold'),legend.text=element_text(size=10,face='bold')) 
    pglm <- 
    pathLab <- file.path('C:','Users','Michelle Mellers','Documents','ThyAnalysR')
    
ggsave(filename=file.path(pathLab,'nonLinChk',paste("1_",classNms(ins),sep=""),ins,outs_y,paste(msTyp,'_',ins,'_',outs_y,'.png',sep="")),width=4,height=4,plot=pglm)
}


tabProc <- function(rout){

rout2 <- rout %>% filter(covs=='s(measured_imp_log)') %>%
  mutate(pts=paste(signif(est,3),' (',signif(se,2),')',sep=''),pv=as.character(signif(pv,2)),
         r2=as.character(signif(r2,2)),sigma=as.character(signif(sigma,2)),
         out=case_when(out=="tshnew_log"~"TSH",
                       out=="total_t3_log"~"TT3",
                       out=="total_t4_log"~"TT4",
                       out=="free_t4_log"~"FT4"),
         pct=paste(signif((2^(est)-1)*100,3),"%",sep="")) %>%
  select(missType,miss_chem,miss_bmi,out,ins,pts,pv,pct,r2,sigma) %>%
  pivot_longer(cols=c(pts,pv,pct,r2,sigma)) %>%
  pivot_wider(names_from=missType,values_from=c(value)) %>%
  mutate(chem_num = case_when(substr(ins,1,3)=="PBB"~gsub("PBB","",ins),
                              substr(ins,1,3)=="PBD"~gsub("PBDE","",ins),
                              substr(ins,1,3)=="PCB"~gsub("PCB","",ins),
                              TRUE~ins),
         type = case_when(substr(ins,1,3)=="PBD"~"PBDE",
                          (substr(ins,1,2)=="PP")|(substr(ins,1,1)=="T") ~"OP",
                          TRUE~substr(ins,1,3)),
         chem_num=sub("_","-",chem_num)) %>%
  select(-ins) %>%
  relocate(out,type,chem_num,miss_chem,miss_bmi,name,CC,SI,NMI) %>%
  group_by(out) %>%
  arrange(chem_num) %>%
  mutate(name=case_when(name=='pts'~"Estimate (SE)",
                        name=="pv"~"p-value",
                        name=="pct"~"Percent Change",
                        name=="r2"~"R2",
                        TRUE~name)) %>%
  filter(!(name %in% c("Estimate (SE)","Percent Change")))
rout2[is.na(rout2)] <- '--'

return(rout2)
}

tabProcSim <- function(rout){
  
  Real <- data.frame(REAL=c(-0.0232,paste(signif((2^(-0.0232)-1)*100,3),'%',sep=''),'--','--','--',rout$sigmaT[1]))
  ## mean estimate and standard deviation of the beta-estimates compared to "true" value
  ## mean estimate and standard deviation of the standard-deviation compared to "true" value
  ##by group
  ## linear model with no missing and no error comparison
  rout2 <- rout %>% filter(covs=="s(measured_imp_log)") %>%
              mutate(out=as.numeric(ins)) %>%
              group_by(missType) %>%
              select(missType,ins,out,miss_bmi,miss_chem,est,se,pv,r2,sigma)# %>%
              summarize(nsim=max(ins),
                        miss_bmi=median(miss_bmi),
                        miss_chem=median(miss_chem),
                      est.25=quantile(est,0.25),
                        est.50=quantile(est,probs=0.5),
                        est.75=quantile(est,0.75),
                        se.25=quantile(se,0.25),
                        se.50=quantile(se,0.5),
                      se.75=quantile(se,0.75),
                      pv.25=quantile(pv,0.25),
                      pv.50=quantile(pv,0.5),
                      pv.75=quantile(pv,0.75),
                      r2.25=quantile(r2,0.25),
                      r2.50=quantile(r2,0.50),
                      r2.75=quantile(r2,0.75),
                      sigma.25=quantile(sigma,0.25),
                      sigma.50=quantile(sigma,0.50),
                      sigma.75=quantile(sigma,0.75)) %>%
    ungroup() %>%
    mutate(est=paste(signif(est.50,2),"(",signif(est.25,2),",",signif(est.75,2),")",sep=''),
           pct=paste(signif((2^(est.50)-1)*100,3),"%(",signif((2^(est.25)-1)*100,3),"%,",signif((2^(est.75)-1)*100,3),"%)",sep=""),
             se=paste(signif(se.50,2),"(",signif(se.25,2),",",signif(se.75,2),")",sep=''),
      pv=paste(signif(pv.50,2),"(",signif(pv.25,2),",",signif(pv.75,2),")",sep=''),
      r2=paste(signif(r2.50,2),"(",signif(r2.25,2),",",signif(r2.75,2),")",sep=''),
        sigma=paste(signif(sigma.50,2),"(",signif(sigma.25,2),",",signif(sigma.75,2),")",sep='')) %>%
    select(missType,nsim,miss_bmi,miss_chem,est,pct,se,pv,r2,sigma) %>%
    pivot_longer(cols=c("est","pct","se","pv","r2","sigma"),names_to=c("namer"),values_to="vals") %>%
    pivot_wider(names_from=missType,values_from=vals) %>%
    cbind(.,Real) %>%
    relocate(nsim,miss_bmi,miss_chem,namer,REAL,CC,SI,NMI) %>%
    mutate(namer=case_when(namer=='est'~"Estimate",
                           namer=='pct'~"Percent Change",
                           namer=='se'~"Std. Error",
                           namer=='pv'~"p-value",
                           namer=="r2"~"R2",
                           TRUE~namer))
  
  return(rout2)
}


ggRMSE <- function(o_rmse2){
  
  y_lab <- labsy("yNoisy")
  outs_y <- simFcnNms(o_rmse2$outs[1])
  
  iners <- "Sim"  
  x_lab <- paste("LN ",iners," (ng/g)",sep="")
  
 # pathLab <- file.path('C:','Users','Michelle Mellers','Documents','ThyAnalysR')
  ps2 <- ggplot(data=o_rmse2,aes(rmse,fill=type)) + geom_density(alpha=0.6) + 
    scale_fill_manual(values=c('grey30','grey50','grey70','grey90')) + 
    labs(x =x_lab,y=y_lab,
         title=paste("RMSE for Each Missing Data Type,\nModel:",outs_y,sep=""),fill="Type") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=12,face='bold'),plot.title=element_text(size=14,hjust=0.5,face='bold'), axis.text.y=element_text(size=12,face='bold'),axis.title.x=element_text(size=12,face='bold'),
          axis.title.y=element_text(size=12,face='bold'),legend.position="bottom",
          legend.title=element_text(size=12,face='bold'),legend.text=element_text(size=10,face='bold')) 
  pathLab <- file.path('C:','Users','Michelle Mellers','Documents','ThyAnalysR')
  
  ggsave(filename=file.path(pathLab,'nonLinChk','1_Sim',outs_y,paste('rmse_',iners,'_',outs_y,'.png',sep="")),width=4,height=4,plot=ps2)
}


ggBias <- function(o_bias2){
  
  y_lab <- labsy("yNoisy")
  outs_y <- simFcnNms(o_bias2$outs[1])
 # outs_y <- outsy(outsy1)
  
  iners <- "Sim"  
  x_lab <- paste("LN ",iners," (ng/g)",sep="")
  ##ggplot

 # pathLab <- file.path('C:','Users','Michelle Mellers','Documents','ThyAnalysR')
  ps3 <- ggplot(data=o_bias2,aes(x=measured_true,y=bias,group=type)) + 
    geom_smooth(aes(linetype=type,color=type),method='loess',se=FALSE,span=0.2) + 
   # geom_smooth(aes(x=measured_true,y=fitTrueSecDev,color='2ndDer',linetype='2ndDer'))+
    scale_linetype_discrete(labels=c('NM','CC','S','NMI')) +
    scale_color_manual(labels=c('NM','CC','S','NMI'),
                       values=c('grey25','grey40','grey55','grey70')) +
    labs(x =x_lab, y='Bias',
         title=paste("Bias for Each Missing Data Type,\nModel:",outs_y,sep=""),
         linetype="Type",color="Type") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=12,face='bold'),plot.title=element_text(size=14,hjust=0.5,face='bold'), axis.text.y=element_text(size=12,face='bold'),axis.title.x=element_text(size=12,face='bold'),
          axis.title.y=element_text(size=12,face='bold'),legend.position="bottom",
          legend.title=element_text(size=10,face='bold'),legend.text=element_text(size=8,face='bold')) 
  #ps4 <- ps3 + geom_hline(yintercept=0,color='black',linetype='solid',size=2)  
  pathLab <- file.path('C:','Users','Michelle Mellers','Documents','ThyAnalysR')
  
  ggsave(filename=file.path(pathLab,'nonLinChk','1_Sim',outs_y,paste('bias_',iners,'_',outs_y,'.png',sep="")),width=4,height=4,plot=ps3)
}

ggTrue <- function(o_bias22){
    
  require(patchwork)
    y_lab <- labsy("yNoisy")
    outs_y <- simFcnNms(o_bias22$outs[1])
    # outs_y <- outsy(outsy1)
    o_bias3 <- o_bias22 %>% select(measured_true,fitTrue,fitTrueSecDev) %>%
      pivot_longer(cols=c("fitTrue","fitTrueSecDev"),names_to="groups",values_to="fitTrueAll")
    iners <- "Sim"  
    x_lab <- paste("LN ",iners," (ng/g)",sep="")
    ##ggplot
    # pathLab <- file.path('C:','Users','Michelle Mellers','Documents','ThyAnalysR')
   ps33a <-  o_bias3 %>% filter(groups=='fitTrue') %>%
      ggplot(data=.,aes(x=measured_true,y=fitTrueAll)) + 
     labs(x =x_lab, y='f(x) (Simulated)') + 
      geom_smooth(color='grey20',method='loess',se=FALSE,span=0.2) + 
      theme_minimal() +
      theme(axis.text.x=element_text(size=12,face='bold'),plot.title=element_text(size=14,hjust=0.5,face='bold'), axis.text.y=element_text(size=12,face='bold'),axis.title.x=element_text(size=12,face='bold'),
            axis.title.y=element_text(size=12,face='bold'),legend.position="bottom",
            legend.title=element_text(size=12,face='bold'),legend.text=element_text(size=10,face='bold')) 
   ps33b <-  o_bias3 %>% filter(groups=='fitTrueSecDev') %>%
     ggplot(data=.,aes(x=measured_true,y=fitTrueAll)) + 
     geom_line(color='grey20')+
     #geom_smooth(color='grey20',method='loess',se=FALSE,span=0.2) + 
     theme_minimal() +
     labs(x =x_lab, y='f(x) (Simulated)') +
     theme(axis.text.x=element_text(size=12,face='bold'),plot.title=element_text(size=14,hjust=0.5,face='bold'), axis.text.y=element_text(size=12,face='bold'),axis.title.x=element_text(size=12,face='bold'),
           axis.title.y=element_text(size=12,face='bold'),legend.position="bottom",
           legend.title=element_text(size=12,face='bold'),legend.text=element_text(size=10,face='bold')) 
    pathLab <- file.path('C:','Users','Michelle Mellers','Documents','ThyAnalysR')
  ps33 <- ps33a + ps33b + plot_annotation( title=paste("Plot of Simulated Functions, Model:",outs_y,sep="")) 
    
    ggsave(filename=file.path(pathLab,'nonLinChk','1_Sim',outs_y,paste('fitTrue_',iners,'_',outs_y,'.png',sep="")),width=4,height=4,plot=ps33)
}
