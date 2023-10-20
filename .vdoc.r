#
#
#
#
#
#
#
#

packages <- c('tidyverse','rmarkdown',"RColorBrewer",'gghalves')
invisible(lapply(packages, require, character.only = TRUE))


e2<- readRDS('data/igas_e2_cleanedData-final.rds')%>% mutate(initialVelocityX=X_Velocity,initialVelocityY=Y_Velocity)
solSpace <- e2 %>% filter(trialType==11)
#solSpace %>% ggplot(aes(x=X_Velocity,y=Y_Velocity)) + geom_point(aes(colour=ThrowPosition),alpha=0.58) + ggtitle("") 

solSpace$Result = ifelse(solSpace$ThrowPosition==400,"400",solSpace$ThrowPosition)
solSpace$Result = ifelse(solSpace$ThrowPosition==500,"500",solSpace$Result)
solSpace$Result= ifelse(solSpace$ThrowPosition==625,"625",solSpace$Result)
solSpace$Result = ifelse(solSpace$ThrowPosition==675,"675",solSpace$Result)
solSpace$Result = ifelse(solSpace$ThrowPosition==800,"800",solSpace$Result)
solSpace$Result = ifelse(solSpace$ThrowPosition==900,"900",solSpace$Result)
theme_set(theme_classic())

d <- readRDS("data/htw-04-07-22.rds")
d <- d %>% mutate(distCapped=ifelse(dist>1500,1500,dist),
                  vxCapped =ifelse(vx>2000,2000,vx),
                  rectWidth=.4,
                  vbn=as.numeric(throwCategory),
                  vbLag=vbn-rectWidth,vbLead=vbn+rectWidth,
                  condit=recode_factor(condit,constant="Constant",varied="Varied"))
dtest <- d %>% filter(expMode=="test-Nf" | expMode=="test-train-nf") %>% 
  droplevels() %>%   
  mutate(testMode=factor(ifelse(expMode=="test-Nf","Novel Velocity Band","Trained Velocity Band \n(Constant trained from only 800-1000)")))


dtest$vbLabel <- factor(dtest$throwCategory, 
                        labels= c("100-300\nNovel", "350-550\nNovel", "600-800\nNovel",
                                  "800-1000\nTrained\nBoth", "1000-1200\nTrained\nVaried","1200-1400\nTrained\nVaried"))


vbRect<- dtest %>% group_by(throwCategory,vbLabel) %>% 
  summarise(vxMean=mean(vx,trim=.05),lowBound=first(bandInt),highBound=first(highBound),vx=median(vx),
            vxCapped=median(vxCapped),devMean=mean(dist)) %>%
  mutate(rectWidth=.4,vbn=as.numeric(throwCategory),vbLag=vbn-rectWidth,vbLead=vbn+rectWidth)


bandLines4 <- list(geom_segment(data=vbRect,aes(x=vbLag,xend=vbLead,y=highBound,yend=highBound),alpha=1,linetype="dashed"),
                   geom_segment(data=vbRect,aes(x=vbLag,xend=vbLead,y=lowBound,yend=lowBound),alpha=1,linetype="dashed"),
                   geom_text(data=vbRect,aes(x=vbLag-.03,y=lowBound+100,label=throwCategory),angle=90,size=2.5,fontface="bold") )    


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#, out.width= "85%", out.extra='style="float:right; padding:2px"'
ss=solSpace %>% ggplot(aes(x=X_Velocity*-1,y=Y_Velocity*-1)) + 
  geom_point(aes(colour=Result),alpha=0.6) + scale_color_manual(values =brewer.pal(n=6,name="Set1"))+
  theme(axis.text=element_blank(),axis.title=element_blank(),axis.line=element_blank(),
        axis.ticks=element_blank(), legend.position = "none" )
ss

#
#
#
#
#
#
#
#
#

dtest %>% group_by(sbjCode,vbLabel,condit,throwCategory) %>%
  summarise(vxMean=mean(vxCapped),lowBound=first(bandInt),highBound=first(highBound),
            vbLag=first(vbLag),vbLead=first(vbLead),.groups = 'keep') %>%
  ggplot(aes(x=vbLabel,y=vxMean,fill=throwCategory))+
  geom_half_violin(color=NA)+ # remove border color
  geom_half_boxplot(position=position_nudge(x=-0.05),side="r",outlier.shape = NA,center=TRUE,
                    errorbar.draw = FALSE,width=.25)+
  geom_half_point(transformation = position_jitter(width = 0.05, height = 0.05),size=.3,aes(color=throwCategory))+
  geom_rect(data=vbRect,aes(xmin=vbLag,xmax=vbLead,ymin=lowBound,ymax=highBound,fill=throwCategory),alpha=.3)+
  bandLines4+
  scale_y_continuous(expand=expansion(add=100),breaks=round(seq(0,2000,by=200),2))+
  scale_fill_discrete(name="Velocity Band")+
  scale_color_discrete(guide="none")+  # remove extra legend
  theme(axis.text=element_blank(),axis.title=element_blank(),
        legend.position = "none" )

#
#
#
#
#
#
#
#
cowplot::ggdraw()+cowplot::draw_image("Assets/LIM_InGame.png",hjust=0)+theme(plot.margin = margin(0, 0, 0, 0))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
