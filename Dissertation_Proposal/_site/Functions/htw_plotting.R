
library(knitr)
library(kableExtra)
dodge <- position_dodge(width = 0.9)

bandLines <- list(geom_hline(aes(yintercept=highBound,color=throwCategory),alpha=.5),
                  geom_hline(aes(yintercept=lowBound,color=throwCategory),alpha=.5))

bandLines2 <- list(geom_hline(aes(yintercept=highBound),alpha=1,linetype="dashed"),
                   geom_hline(aes(yintercept=lowBound),alpha=1,linetype="dashed"))

bandLines3 <- list(geom_segment(aes(x=vbLag,xend=vbLead,y=highBound,yend=highBound),alpha=1,linetype="dashed"),
                   geom_segment(aes(x=vbLag,xend=vbLead,y=lowBound,yend=lowBound),alpha=1,linetype="dashed"))


ggBars <- list(stat_summary(geom="bar",fun=mean,position=dodge,alpha=.8), 
               stat_summary(geom="errorbar",fun.data=mean_se,width=.3,position=dodge))
barJitter=list(stat_summary(geom="bar",fun=mean,position=dodge), 
               stat_summary(geom="errorbar",fun.data=mean_se,width=.3,position=dodge),
               geom_jitter(alpha=.2))

barJitterBox=list(geom_boxplot(notch=TRUE,alpha=.4,width=.35,outlier.shape=NA,lwd=.04,position=dodge),
                  stat_summary(geom="bar",alpha=.8,width=.9,color="black",fun=mean,position=dodge), 
                  stat_summary(geom="errorbar",fun.data=mean_se,width=.25,color="black",alpha=.8,position=dodge),
                  geom_jitter(size=.3,alpha=.2,position=position_jitterdodge()))


barBox=list(geom_boxplot(notch=TRUE,alpha=.4,width=.35,outlier.shape=NA,lwd=.04,position=dodge),
            stat_summary(geom="bar",alpha=.8,width=.9,color="black",fun=mean,position=dodge), 
            stat_summary(geom="errorbar",fun.data=mean_se,width=.25,color="black",alpha=.8,position=dodge))


boxJitter=list(geom_boxplot(notch=TRUE,alpha=.4,width=.35,outlier.shape=NA,lwd=.04,position=dodge),
               geom_jitter(size=.3,alpha=.2,position=position_jitterdodge()))


lineBars <- list(stat_summary(geom="line",fun=mean),
                 stat_summary(geom="point",fun=mean),
                 stat_summary(geom="errorbar",fun.data=mean_se,width=.1))

lineBarsJitter <- list(stat_summary(geom="line",fun=mean),
                       stat_summary(geom="point",fun=mean),stat_summary(geom="errorbar",fun.data=mean_se,width=.25),
                       geom_jitter(size=.3,alpha=.2,position=position_jitterdodge()))



hgrid <- list(theme(panel.grid.major = element_line(color = "black",
                                                    size = 0.1,
                                                    linetype = 1)))