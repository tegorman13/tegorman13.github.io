packages <- c('ez','plyr','dplyr','sciplot','ggplot2','magrittr','car','psych','tidyr',
              'Matrix','cluster','emdbook','hexbin','caret','glmnet',
              'afex','purrr','broom','ggpubr','rstatix','stargazer','data.table','cowplot','magick',
              'stringr')
have = packages %in% rownames(installed.packages())
if ( any(!have) ) {print("installing missing packages"); install.packages(packages[!have]) }
invisible(lapply(packages, require, character.only = TRUE))
#install.packages("data.table", type="source", repos="https://Rdatatable.gitlab.io/data.table")
select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter
summarise <- dplyr::summarise
#rm(list=ls())
dodge <- position_dodge(width = 0.9)





startCluster <- function(progressVar=100,ncores=NULL)
{
  if(is.null(ncores))
  {
  ncores=detectCores(all.tests = FALSE, logical = TRUE)
  }
  cl       <<- makeSOCKcluster(ncores)
  pb       <<- txtProgressBar(max = progressVar, style = 3)
  progress <- function(n){setTxtProgressBar(pb, n)}
  sopts     <<- list(progress = progress)
  registerDoSNOW(cl)
  
print(paste0("started cluster with ",ncores," cores"))
}






















# grab subset of subjects who are >= the cutoff # of strdvs for the given variable. Default central tendency is median. 
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







sbjInfoIgas <- function(dat=d,sbjId){
  
  sbj <- dat %>% filter(sbjCode == sbjId)
  
  if (is.null(sbj) || nrow(sbj)==0){
    print("sbjCode does not exist within dataframe provided"); }
  
  else {
    print(paste("data for subject:",sbjId))
    sbjStats <- sbj %>% filter(mode==1,mouseType==0,trial<=320) %>% group_by(stage,positionX) %>% 
      summarise(n=n(),
                drops = sum(ballDrop == TRUE),
                unr = sum(trialType == 55),
                br = sum(trialType == 44),
                outBounds = sum(trialType == 99),
                wrngSide = sum(trialType == 22),
                wrngDir = sum(wrongDirection == TRUE),
                unusableThrow=drops+unr+br,
                AbsDev = median(AbsDistFromCenter),
                signedDev = mean(predDistFromCent*-1)
      )
    print(sbjStats)
  }
}



sbjInDfInfo <- function(sbjDat,df=d){
  s=as.character(sbjDat[,"sbjCode"])
lapply(s, function(i) sbjInfoIgas(dat=df,sbjId=i))
}

grpInfoIgas <- function(dat=d,grpCode){

  grp<- dat %>% filter(condit2 == grpCode)
  
  if (is.null(grp) || nrow(grp)==0){
    print("sbjCode does not exist within dataframe provided"); }
  
  else {
    print(paste("data for group:",grpCode))
    grpStats <- grp %>% filter(mode==1,mouseType==0,trial<=320) %>% group_by(stage,positionX) %>% 
      summarise(n=n(),
                drops = mean(ballDrop == TRUE),
                unr = mean(trialType == 55),
                br = mean(trialType == 44),
                outBounds = mean(trialType == 99),
                wrngSide = mean(trialType == 22),
                wrngDir = mean(wrongDirection == TRUE),
                unusableThrow=drops+unr+br,
                AbsDev = median(AbsDistFromCenter),
                signedDev = mean(predDistFromCent*-1)
                
      )
    
    (as.matrix(grpStats))
    return 
  }
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







###################### Just Model Functions Below #############################

cvSim <- function(cv,cc,expt=NULL) {
  exp=ifelse(!is.null(expt),ifelse(expt==2,2,1),2)
  
  
  # code for experiment 2
  if(exp==2) {
    
    difMatFull= data.frame()
    for (i in 1:nrow(goodSubjects2)) # iterate over each subject
    {
      s = filter(allTrain,sbjCode==goodSubjects2[,1][i]) # select individual subject
      t = filter(transfer,sbjCode==goodSubjects2[,1][i])
      c= ifelse((goodSubjects2[,1][i] %in% goodSubjectsVaried$sbjCode),cv,cc)
      
      
      t <- t %>% mutate(sim400=simComp1(s,sol400,c),sim500=simComp1(s,sol500,c),sim625=simComp1(s,sol625,c),sim675=simComp1(s,sol675,c),
                        sim800=simComp1(s,sol800,c),sim900=simComp1(s,sol900,c))
      
      difMatFull= rbind(as.matrix(t),difMatFull)
    }
  }
  # code for experiment1
  if(exp==1) {
    difMatFull= data.frame()
    for (i in 1:nrow(goodSubjects2)) # iterate over each subject
    {
      s = filter(allTrain,sbjCode==goodSubjects[,1][i]) # select individual subject
      t = filter(transfer,sbjCode==goodSubjects[,1][i])
      c= ifelse((goodSubjects2[,1][i] %in% goodSubjectsVaried$sbjCode),cv,cc)
      
      t <- t %>% mutate(sim610=simComp1(s,sol610,c),sim760=simComp1(s,sol760,c),sim835=simComp1(s,sol835,c),sim910=simComp1(s,sol910,c))
      difMatFull= rbind(as.matrix(t),difMatFull)
    }
    
  } 
  
  mFull=sortMatE2(difMatFull)
  
  return(mFull)
}



cSim <- function(c,expt=NULL) {
  
  exp=ifelse(!is.null(expt),ifelse(expt==2,2,1),2)
  # code for experiment 2
  if(exp==2) {
    difMatFull= data.frame()
    for (i in 1:length(goodSubjects)) # iterate over each subject
    {
      s = filter(allTrain,sbjCode==goodSubjects[i])  # select individual subject
      t = filter(transfer,sbjCode==goodSubjects[i])
      t <- t %>% mutate(sim400=simComp1(s,sol400,c),sim500=simComp1(s,sol500,c),sim625=simComp1(s,sol625,c),sim675=simComp1(s,sol675,c),
                        sim800=simComp1(s,sol800,c),sim900=simComp1(s,sol900,c))
      
      difMatFull= rbind(as.matrix(t),difMatFull)
    }
  }
  # code for experiment1
  else {
    for (i in 1:nrow(goodSubjects)) # iterate over each subject
    {
      s = filter(allTrain,sbjCode==goodSubjects[,1][i]) # select individual subject
      t = filter(transfer,sbjCode==goodSubjects[,1][i])
      t <- t %>% mutate(sim610=simComp1(s,sol610,c),sim760=simComp1(s,sol760,c),sim835=simComp1(s,sol835,c),sim910=simComp1(s,sol910,c))
      difMatFull= rbind(as.matrix(t),difMatFull)
    }
  } 
  mFull=sortMatE2(difMatFull)
  return(mFull)
}





distComp1 <- function(dat,sol){
  dif = sqrt(((outer(dat$initialVelocityX,sol$initialVelocityX,"-"))^2)+((outer(dat$initialVelocityY,sol$initialVelocityY,"-"))^2))
  return(dif)
}




simComp1 <- function(dat,sol,c){
  dif = sqrt(((outer(dat$initialVelocityX,sol$initialVelocityX,"-"))^2)+((outer(dat$initialVelocityY,sol$initialVelocityY,"-"))^2))
  sim = sum(exp(-c*dif))/nrow(dat)
  
  return(sim)
}


simComp1p <- function(dat,sol,c,p){
  dif = sqrt(((outer(dat$initialVelocityX,sol$initialVelocityX,"-"))^2)+((outer(dat$initialVelocityY,sol$initialVelocityY,"-"))^2))
  sim = sum(exp(-c*dif)^p)/nrow(dat)
  
  return(sim)
}

simCompTrial <- function(x,y,sol,c){
  dif = sqrt(((outer(x,sol$initialVelocityX,"-"))^2)+((outer(y,sol$initialVelocityY,"-"))^2))
  sim = sum(exp(-c*dif))/nrow(sol)
  
  return(sim)
}


simCompTrialOld <- function(x,y,sol,c,p=1){
  dif = sqrt(((outer(x,sol$initialVelocityX,"-"))^2)+((outer(y,sol$initialVelocityY,"-"))^2))
  sim = sum(exp(-c*dif)^p)/nrow(sol)
  return(sim)
}

sortMatE2 <- function(df){
  difMatFull=df
  
  m400 <- difMatFull  %>% rename(sim=sim400,dev="400",scaleDev=sp400)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="400",scalePos="sp400")
  m500 <- difMatFull  %>% rename(sim=sim500,dev="500",scaleDev=sp500)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="500",scalePos="sp500")
  m625 <-difMatFull  %>% rename(sim=sim625,dev="625",scaleDev=sp625)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="625",scalePos="sp625")
  m675 <- difMatFull  %>% rename(sim=sim675,dev="675",scaleDev=sp675)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="675",scalePos="sp675")
  m800 <- difMatFull  %>% rename(sim=sim800,dev="800",scaleDev=sp800)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="800",scalePos="sp800")
  m900 <- difMatFull  %>% rename(sim=sim900,dev="900",scaleDev=sp900)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="900",scalePos="sp900")
  mFull <- rbind(m400,m500,m625,m675,m800,m900)
  
  # mFull$sim = as.numeric(levels(mFull$sim)[mFull$sim])
  # mFull$dev= as.numeric(levels(mFull$dev)[mFull$dev])
  # mFull$scaleDev= as.numeric(levels(mFull$scaleDev)[mFull$scaleDev])
  
  mFull$sim = as.numeric(mFull$sim)
  mFull$dev = as.numeric(mFull$dev)
  mFull$scaleDev = as.numeric(mFull$scaleDev)
  mFull$pos = as.factor(mFull$pos)
  return(mFull)
}

# 
# sortMatE2 <- function(df){
#   difMatFull=df
# 
#   m400 <- difMatFull  %>% rename(sim=sim400,dev=p400,scaleDev=sp400)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="p400",scalePos="sp400")
#   m500 <- difMatFull  %>% rename(sim=sim500,dev=p500,scaleDev=sp500)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="p500",scalePos="sp500")
#   m625 <-difMatFull  %>% rename(sim=sim625,dev=p625,scaleDev=sp625)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="p625",scalePos="sp625")
#   m675 <- difMatFull  %>% rename(sim=sim675,dev=p675,scaleDev=sp675)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="p675",scalePos="sp675")
#   m800 <- difMatFull  %>% rename(sim=sim800,dev=p800,scaleDev=sp800)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="p800",scalePos="sp800")
#   m900 <- difMatFull  %>% rename(sim=sim900,dev=p900,scaleDev=sp900)  %>% select(sbjCode,conditType,sim,dev,scaleDev) %>% mutate(pos="p900",scalePos="sp900")
#   mFull <- rbind(m400,m500,m625,m675,m800,m900)
#   mFull$sim = as.numeric(levels(mFull$sim)[mFull$sim])
#   mFull$dev= as.numeric(levels(mFull$dev)[mFull$dev])
#   mFull$scaleDev= as.numeric(levels(mFull$scaleDev)[mFull$scaleDev])
#   mFull$pos = as.factor(mFull$pos)
# 
#   return(mFull)
# 
# }

# 
# sortMatE2 <- function(df){
#   difMatFull=df
# 
#   m400 <- difMatFull  %>% rename(sim=sim400,dev=p400)  %>% select(sbjCode,conditType,sim,dev) %>% mutate(pos="p400")
#   m500 <- difMatFull  %>% rename(sim=sim500,dev=p500)  %>% select(sbjCode,conditType,sim,dev) %>% mutate(pos="p500")
#   m625 <-difMatFull  %>% rename(sim=sim625,dev=p625)  %>% select(sbjCode,conditType,sim,dev) %>% mutate(pos="p625")
#   m675 <- difMatFull  %>% rename(sim=sim675,dev=p675)  %>% select(sbjCode,conditType,sim,dev) %>% mutate(pos="p675")
#   m800 <- difMatFull  %>% rename(sim=sim800,dev=p800)  %>% select(sbjCode,conditType,sim,dev) %>% mutate(pos="p800")
#   m900 <- difMatFull  %>% rename(sim=sim900,dev=p900)  %>% select(sbjCode,conditType,sim,dev) %>% mutate(pos="p900")
#   mFull <- rbind(m400,m500,m625,m675,m800,m900)
#   mFull$sim = as.numeric(levels(mFull$sim)[mFull$sim])
#   mFull$dev= as.numeric(levels(mFull$dev)[mFull$dev])
#   mFull$pos = as.factor(mFull$pos)
#   return(mFull)
# }

sortMatE1 <- function(df){
  difMatFull=df
  m610 <- difMatFull  %>% rename(sim=sim610,dev=p610)  %>% select(sbjCode,conditType,sim,dev) %>% mutate(pos="p610")
  m760 <- difMatFull  %>% rename(sim=sim760,dev=p760)  %>% select(sbjCode,conditType,sim,dev) %>% mutate(pos="p760")
  m835 <-difMatFull  %>% rename(sim=sim835,dev=p835)  %>% select(sbjCode,conditType,sim,dev) %>% mutate(pos="p835")
  m910 <- difMatFull  %>% rename(sim=sim910,dev=p910)  %>% select(sbjCode,conditType,sim,dev) %>% mutate(pos="p910")
  
  mFull <- rbind(m610,m760,m835,m910)
  mFull$sim = as.numeric(levels(mFull$sim)[mFull$sim])
  mFull$dev= as.numeric(levels(mFull$dev)[mFull$dev])
  mFull$pos = as.factor(mFull$pos)
  
  return(mFull)
}





# 
# testDataBySbj %>% filter(outly==FALSE) %>%ggplot(aes(sdDev, madDev)) + 
#   geom_smooth(method = "loess", colour = "red", fill = "red") + 
#   geom_smooth(method = "lm", colour = "blue", fill = "blue") + 
#   geom_point() + facet_grid(condit2 ~ positionX, margins=TRUE)+xlim(c(0,800))
##

