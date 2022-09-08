#Figure 4 
#!!!All data used for these scripts should be found on sheets in Supplementary Files 1 and 2!!!

#Effects on depletion of AMA-1 during L1 arrest in soma and germline on survival and recovery
colorBlindGrey8   <- c("#56B4E9", "#E69F00", "#999999", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)
install.packages("ggpubr")
library(ggpubr)
library(ggplot2)
library(nlme)

#plots for Figure 4. the hatching efficiency plot in 4A was originally in the supplement, so  it is part of the SupplementaryFigures.R script
#SET working directory to where files are on your computer
setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/FollowUp/")

#Figure 4E (data in Supp File 2)
AID_1hr_wormsizer<-read.csv("AID_1hrExposure_wormsizer.csv",header=T)
head(AID_1hr_wormsizer)
AID_1hr_wormsizer_d1_ama1<-subset(AID_1hr_wormsizer,day=="d1" & condition!="daf-2/CA1200 1mM aux" & condition!="daf-2/CA1200 EtOH" & condition!="CA1202 1mM aux")
AID_1hr_wormsizer_d1_ama1$condition <- factor(AID_1hr_wormsizer_d1_ama1$condition,levels = c("ama-1/CA1200 1mM aux", "ama-1/CA1200 EtOH", "CA1200 1mM aux", "ama-1/CA1352 1mM aux","ama-1/CA1352 EtOH", "CA1352 1mM aux"))
levels(AID_1hr_wormsizer_d1_ama1$condition)

plot_recovery<-ggplot(AID_1hr_wormsizer_d1_ama1,aes(x=condition,y=length,color=condition))+
  geom_boxplot()+
  geom_dotplot(inherit.aes = TRUE,binaxis = "y",alpha=0.5,binwidth=5,stackdir = "center")+
  scale_color_manual(values=c("#56B4E9", "#E69F00","#999999","#56B4E9", "#E69F00","#999999"))+
  theme_classic(base_size = 10)+labs(y="Body length (microns)")+
  ggtitle("48 hr growth after 1 hr exposure to 1 mM auxin")+
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
  labs(tag = "B")
plot_recovery

head(AID_1hr_wormsizer_d1_ama1)
AMA1_recovery_model2<-lme(length~condition,random=~1|replicate,data = subset(AID_1hr_wormsizer_d1_ama1,condition=="CA1200 1mM aux"|condition=="ama-1/CA1200 EtOH"))
summary(AMA1_recovery_model2)
AMA1_recovery_model1<-lme(length~condition,random=~1|replicate,data = subset(AID_1hr_wormsizer_d1_ama1,condition=="ama-1/CA1200 1mM aux"|condition=="ama-1/CA1200 EtOH"))
summary(AMA1_recovery_model1)

AMA1_recovery_model3<-lme(length~condition,random=~1|replicate,data = subset(AID_1hr_wormsizer_d1_ama1,condition=="ama-1/CA1352 1mM aux"|condition=="ama-1/CA1352 EtOH"))
summary(AMA1_recovery_model3) #3e-4
AMA1_recovery_model4<-lme(length~condition,random=~1|replicate,data = subset(AID_1hr_wormsizer_d1_ama1,condition=="CA1352 1mM aux"|condition=="ama-1/CA1352 EtOH"))
summary(AMA1_recovery_model4) #ns
AMA1_recovery_model5<-lme(length~condition,random=~1|replicate,data = subset(AID_1hr_wormsizer_d1_ama1,condition=="CA1352 1mM aux"|condition=="ama-1/CA1352 1mM aux"))
summary(AMA1_recovery_model5)

#Figure 4G (data in Supp File 2)
ama1_germline_count<-read.csv("ama1_CA1352_germcellCount.csv",header = T)
ama1_germline_count_rmvRep2<-subset(ama1_germline_count,replicate!="rep2")
head(ama1_germline_count_rmvRep2)
plot_gonad_recovery<-ggplot(ama1_germline_count_rmvRep2,aes(genotype_treatment,germcell_number,color=genotype_treatment))+
  geom_boxplot()+theme_classic(base_size = 10)+facet_grid(.~time_point)+
  labs(y="Gonad cell number")+ggtitle("GERMLINE AMA-1 DEPLETION Cells in gonad after 24 hr recovery")+
  scale_color_manual(values=colorBlindGrey8)+
  geom_dotplot(inherit.aes = TRUE,binaxis = "y",alpha=0.5,binwidth=0.5,stackdir = "center")+
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 90))+
  labs(tag="C")
plot_gonad_recovery



#Figure 4H (data in Supp File 2)
ama1_aux_timecourse<-read.csv("ama1_timecourse_aux.csv",header = T)
head(ama1_aux_timecourse)
ama1_means<-aggregate(ama1_aux_timecourse$prop_alive,list(condition=ama1_aux_timecourse$condition,time=ama1_aux_timecourse$aux_addition_hr_hatching),mean)

plot_SS_soma<-ggplot(ama1_aux_timecourse,aes(x=aux_addition_hr_hatching,y=prop_alive,color=condition))+
  geom_point(data=ama1_means,aes(x=time,y=x,color=condition),size=2.5)+
  geom_line(data=ama1_means,aes(x=time,y=x,color=condition))+
  geom_point(size=1,alpha=0.7)+
  theme_classic(base_size = 10)+
  scale_color_manual(values=colorBlindGrey8)+
  labs(x="Time of auxin addition relative to hatching",y="Proportion alive at day 12")+
  ggtitle("Survival during L1 arrest upon AMA-1 degradation")+
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 0))+
  labs(tag="E")

plot_SS_soma


#Figure 4I (data in Supp File 2)
#Survival with germline AMA-1 depleted 
ama1_Jul2020_results<-read.csv("ama1_0.1mM_aux_SS.csv",header = T)
ama1_Jul2020_results_d12<-subset(ama1_Jul2020_results,day_scored=="d12")
ama1_Jul2020_results_d12_germline<-subset(ama1_Jul2020_results_d12,tissue=="germline")
ama1_Jul2020_results_d12_germline$hr_aux_addition<-c(12,12,12,84,84,84,132,132,132,12,12,12,36,36,36,84,84,84,132,132,132,12,12,12,84,84,84,132,132,132)
ama1_means_d12germline<-aggregate(ama1_Jul2020_results_d12_germline$prop_alive,list(condition=ama1_Jul2020_results_d12_germline$condition,time=ama1_Jul2020_results_d12_germline$hr_aux_addition),mean)
head(ama1_Jul2020_results_d12_germline)
head(ama1_means_d12germline)
plot_SS_gline<-ggplot(subset(ama1_Jul2020_results_d12_germline,hr_aux_addition!="36"),aes(x=hr_aux_addition,y=prop_alive,color=condition))+
  geom_point(data=subset(ama1_means_d12germline,time!="36"),aes(x=time,y=x,color=condition),size=2.5)+
  geom_line(data=subset(ama1_means_d12germline,time!="36"),aes(x=time,y=x,color=condition))+
  geom_point(size=1,alpha=0.7)+
  theme_classic(base_size = 10)+ylim(0,1)+
  scale_color_manual(values=colorBlindGrey8)+
  labs(x="Time of auxin addition relative to hatching",y="Proportion alive at day 12")+
  ggtitle("Survival during L1 arrest upon AMA-1 degradation")+
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 0))+
  labs(tag="F")

plot_SS_gline



grid.arrange(plot_recovery,plot_gonad_recovery,plot_SS_soma,plot_SS_gline,nrow=2)

grid.arrange(
  grobs = list(plot_recovery,plot_gonad_recovery,plot_SS_soma,plot_SS_gline),
  widths = c(3,3,3,3),
  layout_matrix = rbind(c(1,NA),
                        c(NA,2),
                        c(3,4)))




