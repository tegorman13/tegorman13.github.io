library(plyr)
library(dplyr) # great library for massaging data
library(sciplot) # for the bar chart and line graph with error bars
library(ggplot2)
library(magrittr) # to allow pipelines
library (car) # for special recode
library(psych)
library(tidyr)


freqLevels <- 3
nCategory <- 6

uniItem <- 6

nHigh=18
nMatched=6
nLow=2


freqLevels <- 3
nHighItem=3
nMatchedItem=5
nLowItem=9
nItemPerCat <- nHighItem+nMatchedItem+nLowItem

total=((nHigh*nHighItem)+(nMatched*nMatchedItem)+(nLow*nLowItem))*nCategory
total/3



zipfItemDist = c(rep(nHigh,nHighItem),rep(nMatched,nMatchedItem),rep(nLow,nLowItem))
repZ=rep(zipfItemDist,times=nCategory)
repUnf = rep(nMatched,nCategory*nItemPerCat)


# generate items and category labels
items = 1:(nCategory*nItemPerCat)
c=1:(length(items)/nItemPerCat)
c2=rep(c,times=nItemPerCat)
cat=sort.int(c2)
d=as.data.frame(cbind(items,cat,repZ,repUnf)) |> mutate(cat=as.factor(cat))



item.sim()objects6 = as.factor(c("object1","object2","object3","object4","object5","object6"))


zipfItemDist
#  18 18 18  6  6  6  6  6  2  2  2  2  2  

repUnf
#  6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 


library(tidyverse)
library(patchwork)


# use ggplot to plot zipfian and repUnf distributions

d %>% 
  ggplot(aes(x=items,y=repZ,fill=cat)) +
  geom_bar(stat="identity") +
  #geom_bar(aes(x=items,y=repUnf),stat="identity",fill="red") +
  facet_wrap(~cat, scale="free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Zipfian Training",x="Items",y="Frequency") +
  scale_x_continuous(breaks=seq(1,108,1)) +
  scale_y_continuous(breaks=seq(0,18,1))


d %>% 
  ggplot(aes(x=items,y=repUnf,fill=cat)) +
  geom_bar(stat="identity") +
  #geom_bar(aes(x=items,y=repUnf),stat="identity",fill="red") +
  facet_wrap(~cat, scale="free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Uniform Training",x="Items",y="Frequency") +
  scale_x_continuous(breaks=seq(1,108,1)) +
  scale_y_continuous(breaks=seq(0,18,1))

# (18*3)+(6*5)+(9*2)


skew_train <- d %>% 
filter(cat==1) %>%
  ggplot(aes(x=items,y=repZ,fill=cat)) +
  geom_bar(stat="identity") +
  #geom_bar(aes(x=items,y=repUnf),stat="identity",fill="red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size=11),
  axis.title.x=element_text(size=14),
  axis.title.y=element_text(size=14),
  plot.title=element_text(size=16),
  legend.position = "none") +
  labs(title="Skewed Training",x="Item",y="Frequency") +
  scale_x_continuous(breaks=seq(1,17,1)) +
  scale_y_continuous(breaks=seq(0,18,1)) + 
  scale_fill_manual(values = c("#00A08A","#ECCBAE", "#FF0000")) +
  ylim(0,18)


unif_train <- d %>% 
filter(cat==1) %>%
  ggplot(aes(x=items,y=repUnf,fill=cat)) +
  geom_bar(stat="identity") +
  #geom_bar(aes(x=items,y=repUnf),stat="identity",fill="red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size=11),
  axis.title.x=element_text(size=14),
  axis.title.y=element_text(size=14),
  plot.title=element_text(size=16),
  legend.position = "none") +
  labs(title="Uniform Training",x="Item",y="Frequency") +
  scale_x_continuous(breaks=seq(1,17,1)) +
  scale_y_continuous(breaks=seq(0,18,1)) + 
  scale_fill_manual(values = c("#FF0000"))+
  ylim(0,18)

csl_plot <- skew_train / unif_train


#ggsave(filename=paste0(here::here(),"/Assets/csl_train.png"),plot=csl_plot,bg='transparent',dpi=400)
