#Figure 7
#!!!All data used for these scripts should be found on sheets in Supplementary Files 1 and 2!!!

#SET working directory to where files are on your computer
setwd("/Users/amywebster/Documents/Baugh_PaperRevisions/AMA1_revisions/")
library(nlme)
library(ggplot2)
library(reshape2)
#Figure 7A-C, data in Supplementary File 2
cec4_TotalBroodSize<-read.csv("N2_RB2301_bs.csv",header = T)
head(cec4_TotalBroodSize)

cec4_TotalBroodSize_melt<-melt(cec4_TotalBroodSize,id.vars = c("strain","rep","days_starved","worm","total","sterile","firstday_egglaying"),variable.name="day",value.name="progeny_per_day")
head(cec4_TotalBroodSize_melt)


head(cec4_TotalBroodSize)
cec4_TotalBroodSize_notSterile<-subset(cec4_TotalBroodSize,sterile=="no")
NumberSterile<-aggregate(cec4_TotalBroodSize_notSterile$sterile ~ cec4_TotalBroodSize_notSterile$days_starved + cec4_TotalBroodSize_notSterile$strain + cec4_TotalBroodSize_notSterile$rep, data = cec4_TotalBroodSize_notSterile, FUN = length)
NumberSterile$SterileNumber<- (18 - NumberSterile$`cec4_TotalBroodSize_notSterile$sterile`)
NumberSterile$PropSterile<- (NumberSterile$SterileNumber / 18)
head(NumberSterile)
colnames(NumberSterile)<-c("days_starved","strain","rep","NotSterileNumber","SterileNumber","PropSterile")


ELOTime<-aggregate(cec4_TotalBroodSize_notSterile$firstday_egglaying ~ cec4_TotalBroodSize_notSterile$days_starved + cec4_TotalBroodSize_notSterile$strain + cec4_TotalBroodSize_notSterile$rep, data = cec4_TotalBroodSize_notSterile, FUN = length)

ELO_day1<-subset(cec4_TotalBroodSize,firstday_egglaying=="day1")
ELOTime2<-aggregate(ELO_day1$firstday_egglaying ~ ELO_day1$days_starved + ELO_day1$strain + ELO_day1$rep, data = ELO_day1, FUN = length)
ELOTime2

ELOTime$day1_ELO<-ELOTime2$`ELO_day1$firstday_egglaying`
ELOTime$delayed<-ELOTime$`cec4_TotalBroodSize_notSterile$firstday_egglaying` - ELOTime$day1_ELO
ELOTime$PropDelayed<-ELOTime$delayed / ELOTime$`cec4_TotalBroodSize_notSterile$firstday_egglaying`
ELOTime

colnames(ELOTime)<-c("days_starved","strain","rep","NotDelayed_number","Day1_number","Delayed_number","Prop_Delayed")
ELOTime<-subset(ELOTime,days_starved!="d0")



d1d8_only<-subset(cec4_TotalBroodSize,days_starved!="d0")
head(d1d8_only)
cec4_d1d8_model<-lme(total~strain*days_starved,random=~1|rep,data = d1d8_only)
summary(cec4_d1d8_model)

d1d0_only<-subset(cec4_TotalBroodSize,days_starved!="d8")
head(d1d0_only)
cec4_d1d0_model<-lme(total~strain*days_starved,random=~1|rep,data = d1d0_only)
summary(cec4_d1d0_model)

#Mixed model for the proportion sterile and proportion delayed 
head(NumberSterile)
d1d8_SterileModel<-lme(PropSterile~strain*days_starved,random = ~1|rep, data=subset(NumberSterile,days_starved=="d1"|days_starved=="d8"))
summary(d1d8_SterileModel) #p=0.0004

d0d1_SterileModel<-lme(PropSterile~strain*days_starved,random = ~1|rep, data=subset(NumberSterile,days_starved=="d0"|days_starved=="d1"))
summary(d0d1_SterileModel) #p=1

head(ELOTime)
d1d8_DelayModel<-lme(Prop_Delayed~strain*days_starved,random = ~1|rep, data=subset(ELOTime,days_starved=="d1"|days_starved=="d8"))
summary(d1d8_DelayModel) #p=0.0001


#Put together the key graphs

TotalBroodSize<-ggplot(cec4_TotalBroodSize,aes(x=strain,y=total,color=strain))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(alpha=0.5)+
  facet_grid(.~days_starved)+
  scale_color_manual(values=c("#000000","#c1272d"))+
  theme_classic(base_size = 12)+theme(aspect.ratio = 1)+labs(tag = "A")+
  labs(x="Days of L1 starvation",y="Total progeny upon recovery")+
  scale_x_discrete(labels = c('N2',expression(italic('cec-4(ok3124)'))))+
  theme(legend.position="none")
  TotalBroodSize

TotalBrood_Sterile<-ggplot(NumberSterile,aes(x=strain,y=PropSterile,color=strain))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(alpha=0.5)+
  facet_grid(.~days_starved)+
  scale_color_manual(values=c("#000000","#c1272d"))+
  theme_classic(base_size = 12)+theme(aspect.ratio = 1)+labs(tag = "B")+
  labs(x="Days of L1 starvation",y="Proportion sterile")+
  scale_x_discrete(labels = c('N2',expression(italic('cec-4(ok3124)'))))+
  theme(legend.position="none")


TotalBrood_Delayed<-ggplot(ELOTime,aes(x=strain,y=Prop_Delayed,color=strain))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(alpha=0.5)+
  facet_grid(.~days_starved)+
  scale_color_manual(values=c("#000000","#c1272d"))+
  theme_classic(base_size = 12)+theme(aspect.ratio = 1)+labs(tag = "C")+
  labs(x="Days of L1 starvation",y="Proportion delayed (among fertile)")+
  scale_x_discrete(labels = c('N2',expression(italic('cec-4(ok3124)'))))+
  theme(legend.position="none")


library(gridExtra)
grid.arrange(TotalBroodSize,TotalBrood_Sterile,TotalBrood_Delayed,nrow=3)

