
packages <- c('tidyverse','car',
              'psych','Matrix','purrr','broom','rjson',
              'rstatix','gridExtra','lme4','lmerTest','data.table','entropy','effectsize',
              'ggridges','lattice','grid','gridExtra','gtable','ggpubr','cowplot','rstatix',
              'gghalves')

have = packages %in% rownames(installed.packages())
if ( any(!have) ) { install.packages(packages[!have]) }
invisible(lapply(packages, require, character.only = TRUE))
summarise <- dplyr::summarise
dodge <- position_dodge(width = 0.9)


saveFunction <- function(funToSave,body=TRUE,names=TRUE,environ=FALSE)
{
  fullFunction <- deparse(funToSave)
  body=rlang::fn_body(funToSave)
  names<- rlang::fn_fmls_names(funToSave)
  env<-ifelse(environ==TRUE,rlang::fn_env(funToSave),"NULL")
  returnFunInfo <- list(fullFunction=fullFunction,body=body,formals=names,env=env)
}

engagementOutlier <- function(sbjWiseData,varName, CT="median", cutoff=2.5,printInfo=TRUE){
  
  # compute mean,median, and strdv for given variable
  meanDv=mean(sbjWiseData[,varName])
  medDv=median(sbjWiseData[,varName])
  sdDv=sd(sbjWiseData[,varName])
  # do subsetting
  if (CT=="mean") { outliers <- subset(sbjWiseData,(sbjWiseData[,varName]>=(meanDv+(cutoff*sdDv))));threshVal=meanDv+(cutoff*sdDv)}
  else{outliers <- subset(sbjWiseData,(sbjWiseData[,varName]>=(medDv+(cutoff*sdDv))));threshVal=medDv+(cutoff*sdDv)}
  
  outliers <- arrange(outliers,sbjCode)
  if (printInfo==TRUE){
    descriptiveStat <- paste(varName,"has mean",round(meanDv,1),"and Median: ",medDv,"with SD of: ", round(sdDv,1))
    cutOff <- paste("excluding ",nrow(outliers), "subjects with ",CT, "greater than or equal to: ", round(threshVal,2) )
    print(descriptiveStat)
    print(cutOff)
    print(outliers)
  }
  return(outliers)
}




barPlotFg <- function(dat2Plot,xVar,yVar,groupVar,colorVar,fg,title="",xl="",yl="",lt="",leg=TRUE,lim.y=NULL)
{
  
  gg <-  ggplot(dat2Plot,aes_string(x=deparse(substitute(xVar)),y=deparse(substitute(yVar)),
                                    group=deparse(substitute(groupVar)),fill=deparse(substitute(colorVar))))+
    geom_bar(stat="summary",position=dodge,fun=mean)+ 
    
    #facet_grid(fg~fg2,labeller=label_parsed)+
    #facet_grid(noquote(substitute(fg)))+
    #facet_grid(fg,labeller=label_parsed)+
    #facet_grid(.~fg,labeller=label_parsed)+
    stat_summary(fun.data=mean_se,geom="errorbar",position=dodge,width=.5)
  print(fg)
  if(!(fg =='null')){
    gg <- gg+facet_grid(substitute(fg))
  }
    gg <- gg+ylab(yl) +xlab(xl)+theme(plot.title = element_text(hjust = 0.5))+
    guides(fill=guide_legend(title=lt))+theme(legend.title.align=.25)
    
    if(leg==FALSE){
      gg <- gg+theme(legend.position ="none")
    }
    
    if(!(is.null(lim.y))){
      gg <- gg + coord_cartesian(ylim=lim.y)
      
      
    }
  
   gg<- gg+ ggtitle(title)
   return(gg)
   
  }
 
  
  
#rm(list=ls())

#library(raster)
#library(ggplotify)

#library(plot.matrix)



plotWeightMat <- function(wm)
{
  wm  <- reshape2::melt(wm) 
  colnames(wm) <- c("output","input","weight")
  trainWeights <- ggplot(wm)+ geom_tile(aes(x=input, y=output, fill=weight))+  # fill + legend, gray border
    coord_equal() # squares
  trainWeights <- trainWeights + theme_bw() # remove some chart junk
  trainWeights <- trainWeights + theme(panel.grid=element_blank())
  trainWeights <- trainWeights + theme(panel.border=element_blank()) +ggtitle("Weights")
  trainWeights
  
}


plotIndv_tt <- function(sbjDat,nInclude,fullDat,inType="optimOut"){
  funGmat <- list()
  
  for (i in 1:nInclude)
  {
    sbj=sbjDat[i]
    sbj1=sbjDat[i]
    
    if(inType=="optimOut"){
      fitDat <- fullDat %>% filter(sbjCode==sbj1)
      fitDat<- fitDat$optimPred[[1]]
      lr=round(fitDat$lr[1],4); c = round(fitDat$c[1],4); trainRmsd=round(fitDat$rmsd[1],1);
      condit=fitDat$condit[1]
      wb = round(fitDat$weightBias[1],4)
      gtitle=paste("Sbj: ",sbj,",",fitDat$condit[1],", "," lr:",lr," c:",c," wb: ",wb," train rmsd:",trainRmsd,sep=" ")
    }
    if(inType=="df"){
      fitDat=fullDat
      lr=round(fullDat$lr[i],4); c = round(fullDat$c[i],4); condit=fitDat$condit[1]; rmsd=round(fitDat$rmsd[1],1)
      gtitle=paste("Sbj: ",sbj,",",fitDat$condit,", "," lr:",lr," c:",c,"  rmsd:",rmsd,sep="")
      
    }
    title = ggdraw()+draw_label(gtitle,fontface = 'bold',x=0,hjust=0)+theme(plot.margin = margin(0, 0, 0, 7))
    
    
    # gTrain <- fitDat %>% filter(expMode=="train")%>% ggplot(aes(nGoodTrial,vxb,color=vb))+
    # stat_summary_bin(geom="line",fun=mean,bins=3)+stat_summary_bin(aes(y=almPred.vx,color=vb),geom="point",fun=mean,bins=3
    gTrain <- fitDat %>% filter(expMode=="train")%>% ggplot(aes(vb,vxb,fill=vb))+
      stat_summary(geom="bar",fun=mean)+stat_summary(geom="errorbar",fun.data=mean_se,width=.2)+
      stat_summary(geom="point",aes(x=vb,y=almPred.vx,fill=vb),fun=mean)+ggtitle("training average")
    
    
    s.wm <- fitDat$weightMat[1][[1]][[1]]
    lt.wm <- fitDat$weightMat[fitDat$lastTrain[1]][[1]][[1]]
    f.wm <- fitDat$weightMat[fitDat$lastTrial[1]][[1]][[1]]
    
    
    dat_long0 <- reshape2::melt(s.wm) %>%dplyr::mutate(stage="Start_Training")
    dat_long1 <- reshape2::melt(lt.wm) %>%dplyr::mutate(stage="End_Training")
    dat_long2 <- reshape2::melt(f.wm) %>% dplyr::mutate(stage="Final_Testing")
    wm <- do.call(rbind,list(dat_long0,dat_long1,dat_long2)) 
    colnames(wm) <- c("output","input","weight","stage")
    wm$stage <- factor(wm$stage,levels=c("Start_Training","End_Training","Final_Testing"))
    
    
    trainWeights <- ggplot(wm)+ geom_tile(aes(x=input, y=output, fill=weight))+  # fill + legend, gray border
      coord_equal() +facet_wrap(~stage) # squares
    trainWeights <- trainWeights + theme_bw() # remove some chart junk
    trainWeights <- trainWeights + theme(panel.grid=element_blank())
    trainWeights <- trainWeights + theme(panel.border=element_blank()) +ggtitle("Weights")
    
    gTest <- fitDat %>% filter(expMode=="test-Nf" | expMode=="test-train-nf") %>% ggplot(aes(x=vb,y=vxb,fill=vb))+
      stat_summary(geom="bar",fun=mean)+stat_summary(geom="errorbar",fun.data=mean_se,width=.2)+
      stat_summary(geom="point",aes(x=vb,y=examPred.vx,fill=vb),fun=mean)+ggtitle("testing -no feedback")
    
    # fitDat %>% filter(expMode=="test-feedback") %>% ggplot(aes(x=nGoodTrial,y=vxb,color=vb))+
    #     stat_summary_bin(geom="line",fun=mean,bins=3)+
    #   stat_summary_bin(aes(y=almPred.vx,color=vb),geom="point",fun=mean,bins=3)
    
    gFinal <- fitDat %>% filter(expMode=="test-feedback") %>% ggplot(aes(x=vb,y=vxb,fill=vb))+
      stat_summary(geom="bar",fun=mean)+stat_summary(geom="errorbar",fun.data=mean_se,width=.2)+
      stat_summary(geom="point",aes(x=vb,y=almPred.vx,fill=vb),fun=mean)+ggtitle("final test w/feedback")
    
    
    p=plot_grid(title,NULL,gTrain,trainWeights,gTest,gFinal,ncol=2,rel_heights=c(.1,1,1))
    
    
    funGmat[[(length(funGmat)+1)]]=p
    
  }
  
  return(funGmat)
}





scaleMad <- function(grpM,grpV,sbjM){
  (sbjM-grpM) / grpV
}

scaleVar <- function(grpM,grpV,sbjM){
  (sbjM-grpM) / grpV
}


scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}



linearModel <- function(x){
  nls(vxMed ~ m * lowBound+b, data = x, start = list(m =0,b = 0), trace = F)
}

linearModelAvg <- function(x){
  nls(vxAvg ~ m * lowBound+b, data = x, start = list(m =0,b = 0), trace = F)
}

linearModelTrial <- function(x){
  nls(vx ~ m * lowBound+b, data = x, start = list(m =0,b = 0), trace = F)
}




Round2 <- function(x,y) {
  if((y - x %% y) <= x %% y) { x + (y - x %% y)}
  else { x - (x %% y)}
}

##### External Functions #######

