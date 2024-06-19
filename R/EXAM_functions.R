packages <- c('plyr','dplyr','tidyr','data.table','magrittr','ggplot2','Matrix')
lapply(packages, require, character.only = TRUE)








exam_htw <- function(df,learnRate,scaling,wb,netArc,opts) # calls some cpp functions
{

  #mute=TRUE
  weightMat = netArc$weightMat
  inputNode=netArc$inputNode
  outputNode=netArc$outputNode
  input=df$input
  fbVec=df$feedback
  output=df$vxb
  lastT=df$lastTrial[1]
  lastTrain=df$lastTrain[1]
  # weightMat = random_weights(length(outputNode),length(inputNode),0.0025,0.001)
  
  #weightMat <- diag(x=wb,length(outputNode))[,1:length(inputNode)]
  weightMat[1,1]=wb
  #weightMat <- matrix(0,nrow=length(outputNode),ncol=length(inputNode)) # weights initialized to 0 (as in Delosh 1997)
  
  
  endVec = c(1,df$lastTrain,df$lastTrial)
  trainVec=df$trainVec[[1]]
  extrapVec=c(1,trainVec)
  
  #trainVec=c()
  for (i in 1:length(input))
    # for (i in 1:109)
  {
    # input.activations=exp(-scaling*(input[i]-inputNode)^2) # input activity = psychological distance between node and stimulus value
    # input.activations <- input.activations/sum(input.activations)
    # output.activations<- Matrix::crossprod(t(weightMat),input.activations)
    # output.probs <- output.activations /( sum(output.activations))
    # alm.output <- round(sum(outputNode*output.probs),3)
    
    
    id <- pointVecDist(input[i],inputNode)
    input.activations<- inputActivation(scaling,id)
    output.activations <- aocalc(weightMat,input.activations)
    alm.output <- alm(output.activations,outputNode)
    
    
    # alm.output=inputToALM(input[i],scaling,inputNode,outputNode,weightMat)
    
    
    df[i,'almPred.vx']<- alm.output
    df[i,'Model_Error']<- output[i]-alm.output # set here, overwritten if exam runs
    
    
    if (fbVec[i] %in% opts$fbOn){
      trainVec<- c(trainVec,which.max(input.activations)[! which.max(input.activations) %in% trainVec])
      
      # if(opts$mode=="trainModel"){feedbackDist = computeFeedbackDist(input[i],outputNode)}
      # if(opts$mode=="fitSbj") {feedbackDist=outputNode-output[i]}
      
      #feedbackDist=outputNode-output[i]
      #feedbackSignal <- exp(-scaling*((feedbackDist)^2)) #feedback signal includes gaussian similarity gradient (mirroring input nodes)
      feedbackDist <- pointVecDist(output[i],outputNode)
      feedbackSignal <- inputActivation(scaling,feedbackDist)
      ws = learnRate*(feedbackSignal-output.activations)
      #deltaW=ws %*% t(input.activations) # calculate the amount that each weight needs to change (outer product)
      
      
      #deltaW= Matrix::tcrossprod(ws,input.activations)
      oe <- feedbackSignal-output.activations
      deltaW=deltawcalc(learnRate,oe,input.activations)
      
      # weightMat <- weightMat + deltaW
      
      weightMat <- wupdate(weightMat,deltaW)
      # plotWeightMat(weightMat)
      
      #df[[i,'weightMat']]<- list(weightMat)
      
      if(df$modeStart[i]==TRUE || i==lastT || i==lastTrain){ # save weightMat at these indices
        df[[i,'weightMat']]<- list(weightMat)
      }
      
    }
    
    if ( opts$exam && fbVec[i] %in% opts$fbOff) # Exam generalization?
    {
      
      #mean.output
      maxIn=which.max(input.activations)
      nearestTrain <- trainVec[which.min(abs(maxIn-trainVec))] # training item most similar to test item
      
      #xUnder = ifelse(min(trainVec)==nearestTrain,nearestTrain,trainVec[which(trainVec==nearestTrain)-1])
      xUnder = ifelse(min(extrapVec)==nearestTrain,nearestTrain,extrapVec[which(extrapVec==nearestTrain)-1])
      
      
      #xOver = ifelse(max(trainVec)==nearestTrain,nearestTrain,trainVec[which(trainVec==nearestTrain)+1])
      xOver = ifelse(max(extrapVec)==nearestTrain,nearestTrain,extrapVec[which(extrapVec==nearestTrain)+1])
      
      
      mOver1=exp(-scaling*(inputNode[xOver]-inputNode)^2)
      #mOver2 = weightMat %*%(mOver1/sum(mOver1))
      mOver2=Matrix::crossprod(t(weightMat),(mOver1/sum(mOver1)))
      
      mOver3= sum(outputNode*(mOver2/(sum(mOver2))))
      
      mUnder1=exp(-scaling*(inputNode[xUnder]-inputNode)^2)
      #mUnder2 = weightMat %*%(mUnder1/sum(mUnder1))
      mUnder2 = Matrix::crossprod(t(weightMat),(mUnder1/sum(mUnder1)))
      mUnder3=sum(outputNode*(mUnder2/(sum(mUnder2))))
      
      exam.output= round(alm.output + ((mOver3-mUnder3)/(xOver-xUnder))*(maxIn-nearestTrain),3)
      df[i,'examPred.vx']<-exam.output
      
      df[i,'Model_Error']<-output[i]-exam.output
      
    }
    
    if (opts$mute==FALSE){
      print(paste0("trial:",i,".stage= ",df$feedback[i],".input=",input[i],".Correct output: ",df$vb[i],
                   " obsv.out=",df$vx[i]," alm pred=",df$almPred.vx[i],
                   " .exam pred=",df$examPred.vx[i]," deviance=",df$Model_Error[i]))
      
    }
    
  }
  
  run.rmsd=rmsd(df$Model_Error); #print(run.rmsd)
  if(opts$optimModel==TRUE) {return(run.rmsd)}
  if (opts$returnAll==TRUE){
    df <- df %>% mutate(lr=learnRate,c=scaling,weightBias=wb,rmsd=run.rmsd)
    df <- df %>% relocate(Model_Error,almPred.vx,examPred.vx,.after=vxb)
    return(df)}
}










fitExamIndv <- function(sbj,parStart,fd,in.opts)
{
  #print(sbj);print(par); 
  dat=fd[sbjCode==sbj,]
  dat$Model_Error=0;dat$almPred.vx=0;dat$examPred.vx=0;dat$weightMat<-0
  #tdat=dat[expMode=="train" | expMode=="train-Nf",]
  #tdat=dat
  tdat=dat[expMode %in% in.opts$nnTrainModes,]
  inputNode = matrix(seq(1,7,1) ) # 
  outputNode = matrix(seq(50,1600,50))
  weightMat <- matrix(.0000001,nrow=length(outputNode),ncol=length(inputNode)) # weights initialized to 0 (as in Delosh 1997)
  
  arc <- list(weightMat=weightMat,inputNode=inputNode,outputNode=outputNode)
  
  
  if(in.opts$optimMethod=="NM")
    indvTime = system.time(
      out.indv<-optim(par=parStart,fitHtw,df=tdat,netArc=arc,opts=in.opts) 
    )
  
  if(in.opts$optimMethod=="BFGS")
  indvTime = system.time(
  out.indv<-optim(par=parStart,fitHtw,df=tdat,netArc=arc,opts=in.opts,method="L-BFGS-B",lower=in.opts$BFGS.low,control=list(trace=0)) 
  )
  
  fit.lr=out.indv$par[1]; #.54
  fit.c=out.indv$par[2] #1.29
  fit.w=out.indv$par[3] 
  
  
  fbValues=list(1,"1","train","test-feedback")
  noFbValues=list(0,"0","test-Nf","test-train-nf","train-Nf")
  newOpts <- list(mode="fitSbj",optimModel=FALSE,returnAll=TRUE,mute=TRUE,exam=in.opts$exam,fbOn=fbValues,fbOff=noFbValues)
  optimPred=list(exam_htw(dat,learnRate = fit.lr,scaling = fit.c,wb=fit.w,netArc=arc,opts=newOpts)) # generate predictions using fit parameters
  optimResult=list(out.indv)
  runTime=round(indvTime[3],3)
  finfo=list(sbj,out.indv$par[1],out.indv$par[2],out.indv$par[3],
             out.indv$value,out.indv$convergence,optimPred,optimResult,runTime) 
  names(finfo)=c("sbj","lr","c","wb","rmsd","converge","optimPred","optimOutput","fitTime")
  
  
  print(paste("finish sbj:",sbj,"train.rmsd=",round(finfo$rmsd,2),"totalRMSD=",round(optimPred[[1]]$rmsd[1],2),
              "lr=",round(finfo$lr,3),"c=",round(finfo$c,3),"w=",finfo$w,". time=",runTime))
  return(finfo)
  
}










fitExamIndvParallel <- function(sbj,parStart,fd,in.opts)
{
  #print(sbj);print(par); 
  
  
  
  
  
  
  
  
  dat=fd[sbjCode==sbj,]
  dat$Model_Error=0;dat$almPred.vx=0;dat$examPred.vx=0;dat$weightMat<-0
  #tdat=dat[expMode=="train" | expMode=="train-Nf",]
  #tdat=dat
  tdat=dat[expMode %in% in.opts$nnTrainModes,]
  inputNode = matrix(seq(1,7,1) ) # 
  outputNode = matrix(seq(50,1600,50))
  weightMat <- matrix(.0000001,nrow=length(outputNode),ncol=length(inputNode)) # weights initialized to 0 (as in Delosh 1997)
  
  arc <- list(weightMat=weightMat,inputNode=inputNode,outputNode=outputNode)
  
  
  if(in.opts$optimMethod=="NM")
    indvTime = system.time(
      out.indv<-optim(par=parStart,fitHtw,df=tdat,netArc=arc,opts=in.opts) 
    )
  
  if(in.opts$optimMethod=="BFGS")
    indvTime = system.time(
      out.indv<-optim(par=parStart,fitHtw,df=tdat,netArc=arc,opts=in.opts,method="L-BFGS-B",lower=in.opts$BFGS.low,control=list(trace=0)) 
    )
  
  fit.lr=out.indv$par[1]; #.54
  fit.c=out.indv$par[2] #1.29
  fit.w=out.indv$par[3] 
  
  
  fbValues=list(1,"1","train","test-feedback")
  noFbValues=list(0,"0","test-Nf","test-train-nf","train-Nf")
  newOpts <- list(mode="fitSbj",optimModel=FALSE,returnAll=TRUE,mute=TRUE,exam=in.opts$exam,fbOn=fbValues,fbOff=noFbValues)
  optimPred=list(exam_htw(dat,learnRate = fit.lr,scaling = fit.c,wb=fit.w,netArc=arc,opts=newOpts)) # generate predictions using fit parameters
  
  finfo=list(sbj,out.indv$par[1],out.indv$par[2],out.indv$par[3],
             out.indv$value,out.indv$convergence,optimPred,out.indv,indvTime) 
  names(finfo)=c("sbj","lr","c","wb","rmsd","converge","optimPred","optimOutput","fitTime")
  
  runTime=round(indvTime[3],3)
  print(paste("finish sbj:",sbj,"train.rmsd=",round(finfo$rmsd,2),"totalRMSD=",round(optimPred[[1]]$rmsd[1],2),
              "lr=",round(finfo$lr,3),"c=",round(finfo$c,3),"w=",finfo$w,". time=",runTime))
  return(finfo)
  
}











fitHtw <- function(par,df,netArc, opts)
{
  #print(o.n)
  learnRate=par[1]
  scaling=par[2]
  wb=ifelse(length(par)>2,par[3],.0000001)
  
  if(opts$optimMethod=="NM"){
  if(scaling <= opts$parmLowLimit$cLow) {return(fitOptions$parmLowLimit$returnExtreme)}
  if(learnRate <= fitOptions$parmLowLimit$lrLow) {return(fitOptions$parmLowLimit$returnExtreme)}
  if(wb <=0 || wb>1) {return(10000)}
  }
  
  #t=df[expMode=="train",]
  exam_htw(df,learnRate,scaling,wb,netArc=netArc, opts=opts)
}







# plot(raster(weightMat))
# if(i%%50==0){
# par(mfrow = c(2,2))
# plot(raster(feedbackSignal),main=c("signal",i))
# plot(raster(ws),main="ws")
# plot(raster(deltaW),main=c("deltaW",i))
# plot(raster(weightMat),main=c("wm",i))
#dev.off()














# distance should be 0 for any correct value bin, and deviation from closest side for any other
computeFeedbackDist <- function(band,oN)
{
  # binSize=50
  # rangeMin=50
  # rangeMax=1500
  bandBounds=lookupBandBounds(band)
  corMin=bandBounds[1];corMax=bandBounds[2]
  inBand = oN[which(oN>=corMin & oN<=corMax)]
  lower=oN<corMin
  higher=oN>corMax
  highNode=oN[which(higher)]-corMax
  lowNode=oN[which(lower)]-corMin
  bandNode=rep(0,length(inBand))
  sigZ <- c(lowNode,bandNode,highNode)
  
  return(sigZ)
}


computeFeedbackDistSbj <- function(inputVxb,oN)
{
  # binSize=50
  # rangeMin=50
  # rangeMax=1500
  
  inBand = inputVxb
  lower=oN<inBand
  higher=oN>inBand
  highNode=oN[which(higher)]-corMax
  lowNode=oN[which(lower)]-corMin
  bandNode=rep(0,length(inBand))
  sigZ <- c(lowNode,bandNode,highNode)
  
  return(sigZ)
}




lookupBandBounds <- function(band)
{
  
  # if(band==1){return(c(100,300))}
  # if(band==2){return(c(350,550))}
  # if(band==3){return(c(600,800))}
  # if(band==4){return(c(800,1000))}
  # if(band==5){return(c(1000,1200))}
  # if(band==6){return(c(1200,1400))}
  
  if(band==2){return(c(100,300))}
  if(band==3){return(c(350,550))}
  if(band==4){return(c(600,800))}
  if(band==5){return(c(800,1000))}
  if(band==6){return(c(1000,1200))}
  if(band==7){return(c(1200,1400))}
  
  
}




unSignedDev <- function(band,output)
{
  bandBounds=lookupBandBounds(band)
  corMin=bandBounds[1];corMax=bandBounds[2];
  if(output>=corMin && output<=corMax) {return (0)}
  if(output<corMin) {return (corMin-output)}
  if(output>corMax) {return (output-corMax)}
  
}


random_weights <- function(nrow,ncol,mu,sigma){
  random_vector = mu + sigma*runif(nrow*ncol,0,1)
  temp_weights <- matrix(random_vector, nrow = nrow, ncol = ncol)
  return(temp_weights)
}


rmsd <- function(devVec)
{
  sqrt(mean(devVec^2))
  
}

Round2 <- function(x,y) {
  if((y - x %% y) <= x %% y) { x + (y - x %% y)}
  else { x - (x %% y)}
}
