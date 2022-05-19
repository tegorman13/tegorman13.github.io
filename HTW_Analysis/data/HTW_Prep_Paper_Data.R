#rm(list=ls())

source('Functions/htw_functions.R')
source('Functions/htw_plotting.R')
select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter
summarise <- dplyr::summarise
ifelse = base::ifelse

d <- readRDS("data/htw-04-07-22.rds")
options(contrasts = c("contr.sum", "contr.poly"))
defaultContrasts = options()$contrasts
theme_set(theme_classic())

d <- d %>% filter(catOrder=="orig",feedbackType=="continuous",goodThrow=="TRUE") %>% droplevels() 
d$trainStage=factor(d$trainStage,ordered = TRUE)

d <- d %>% mutate(distCapped=ifelse(dist>1500,1500,dist),
                  vxCapped =ifelse(vx>2000,2000,vx),
                  rectWidth=.4,
                  vbn=as.numeric(throwCategory),
                  vbLag=vbn-rectWidth,vbLead=vbn+rectWidth,
                  condit=recode_factor(condit,constant="Constant",varied="Varied"))

dt <- d %>% filter(expMode=="train")
dtf <- d %>% filter(expMode=="test-feedback")
dtest <- d %>% filter(expMode=="test-Nf" | expMode=="test-train-nf") %>% 
  droplevels() %>%   
  mutate(testMode=factor(ifelse(expMode=="test-Nf","Novel Velocity Band","Trained Velocity Band \n(Constant trained from only 800-1000)")))



dtest$vbLabel <- factor(dtest$throwCategory, 
                        labels= c("100-300\nNovel", "350-550\nNovel", "600-800\nNovel",
                                  "800-1000\nTrained\nBoth", "1000-1200\nTrained\nVaried","1200-1400\nTrained\nVaried"))


dtrainTest <- d %>% filter(expMode=="train" | expMode=="test-Nf" | expMode=="test-train-nf") %>% 
  mutate(phase=factor(ifelse(expMode=="train","train","test"),levels=c("train","test"))) %>%
  group_by(sbjCode,throwCategory,phase) %>% mutate(ct=row_number()) %>% relocate(ct,.after="trial")


# separate dataframe for showing solution bands, makes alpha overlap work
vbRect<- dtest %>% group_by(throwCategory,vbLabel) %>% 
  summarise(vxMean=mean(vx,trim=.05),lowBound=first(bandInt),highBound=first(highBound),vx=median(vx),
            vxCapped=median(vxCapped),devMean=mean(dist)) %>%
  mutate(rectWidth=.4,vbn=as.numeric(throwCategory),vbLag=vbn-rectWidth,vbLead=vbn+rectWidth)



vbRect2<- dtest %>% group_by(throwCategory,vbLabel,condit) %>% 
  summarise(vxMean=mean(vx,trim=.05),lowBound=first(bandInt),highBound=first(highBound),vx=median(vx),
            vxCapped=median(vxCapped)) %>%
  mutate(rectWidth=.4,vbn=as.numeric(condit),vbLag=vbn-rectWidth,vbLead=vbn+rectWidth)



bandLines4 <- list(geom_segment(data=vbRect,aes(x=vbLag,xend=vbLead,y=highBound,yend=highBound),alpha=1,linetype="dashed"),
                   geom_segment(data=vbRect,aes(x=vbLag,xend=vbLead,y=lowBound,yend=lowBound),alpha=1,linetype="dashed"),
                   geom_text(data=vbRect,aes(x=vbLag-.03,y=lowBound+100,label=throwCategory),angle=90,size=2.5,fontface="bold") )    


bandLines5 <- list(geom_segment(data=vbRect2,aes(x=vbLag,xend=vbLead,y=highBound,yend=highBound),alpha=1,linetype="dashed"),
                   geom_segment(data=vbRect2,aes(x=vbLag,xend=vbLead,y=lowBound,yend=lowBound),alpha=1,linetype="dashed"),
                   geom_text(data=vbRect2,aes(x=vbLag-.03,y=lowBound+100,label=throwCategory),angle=90,size=2.5,fontface="bold") )  

bandLines6 <- list(geom_segment(data=vbRect,aes(x=vbLag,xend=vbLead,y=highBound,yend=highBound),alpha=1,linetype="dashed"),
                   geom_segment(data=vbRect,aes(x=vbLag,xend=vbLead,y=lowBound,yend=lowBound),alpha=1,linetype="dashed"),
                   geom_text(data=vbRect,aes(x=vbLag-.001,y=lowBound+100,label=throwCategory),angle=90,size=3.7,fontface="bold") )   
