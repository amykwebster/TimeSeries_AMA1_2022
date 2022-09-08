 library(ggplot2)
library(forcats)
library(dplyr)
library(edgeR)
library(gridExtra)
library(reshape2)

#Figure 1
#!!!All data used for these scripts should be found on sheets in Supplementary Files 1 and 2!!!
#This also includes the one-way ANOVA-like test for determining differentially expresed genes

#SET working directory to where files are on your computer
#setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/RNAseq_timeseries/")


#hatch curves for RNAseq data
#Figure 1B (data in Supp File 2)
##RNAseq_hatch<-read.csv("HatchCurves_RNAseq.csv",header = T) READ IN FILE 
RNAseq_hatch1<-subset(RNAseq_hatch,replicate=="mean")
RNAseq_hatch2<-subset(RNAseq_hatch,replicate=="rep1" | replicate=="rep5" | replicate=="rep6" | replicate=="rep8")
plot1<-(ggplot(RNAseq_hatch2,aes(x=hour,y=proportion_L1s,group=replicate))+geom_point(size=1.5,alpha=0.5)+
  geom_line(alpha=0.3)+
  geom_line(inherit.aes=FALSE,data=RNAseq_hatch1,aes(x=hour,y=proportion_L1s),size=1)+
  theme_classic(base_size = 10)+theme(aspect.ratio = 1)+labs(x="Hours of L1 arrest", y= "Proportion hatched")
  +labs(tag="B"))
plot1


#survival at day 8 and day 12
#Figure 1C (data in Supp File 2)
survivald8d12<-read.csv("Survival_d8d12.csv",header = T)
survival1<-subset(survivald8d12,rep=="mean")
survival2<-subset(survivald8d12,rep=="rep1"|rep=="rep5"|rep=="rep6"|rep=="rep8")
survival1$day2<- c("8 d","12 d")
survival2$day2<-c("12 d","12 d","12 d","12 d","8 d","8 d","8 d","8 d")
survival2$day3<-c(12,12,12,12,8,8,8,8)

plot2<- survival2 %>%
  mutate(day2= fct_relevel(day2,"8 d","12 d")) %>%
  ggplot(aes(x=day2,y=prop_alive,group=rep))+
  geom_line(alpha=0.3)+
  #geom_bar(inherit.aes = FALSE,data=survival1,aes(x=day2,y=prop_alive),stat="identity",width = 0.5)+
  geom_point(size=1.5,alpha=0.5)+theme_classic(base_size = 10)+theme(aspect.ratio = 1)+
  labs(x="Days of L1 arrest",y="Proportion alive")+ylim(0,1)+labs(tag="C")

plot2

#PCA of the time series
#Figure 1D, pull in the count data from Supp File 1 (Counts_TimeSeries tab) and filter to 16699 protein coding genes with CPM>1 in at least 4 libraries (Those without NA values on the TimeSeries_DiffExpResults tab of Supp File 1)
#Timeseries_Counts<-read.csv("counts_StarvationTimeSeries.csv",header = TRUE,row.names="sequence")
#Maxwell_genes<-read.csv("Maxwell_Kaplan_RNAseqdata.csv",header = T,row.names = "gene_id")
#Timeseries_filtered<-Timeseries_Counts[rownames(Maxwell_genes),]

setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/RNAseq_timeseries")
Timeseries_Groups<-c("-2hr","-2hr","-2hr","-2hr","0hr","0hr","0hr","0hr","2hr","2hr","2hr","2hr","4hr","4hr","4hr","4hr","6hr","6hr","6hr","6hr","9hr","9hr","9hr","9hr","12hr","12hr","12hr","12hr","1d","1d","1d","1d","2d","2d","2d","2d","4d","4d","4d","4d","8d","8d","8d","8d","12d","12d","12d","12d")
d2<-DGEList(counts=Timeseries_filtered,group=factor(Timeseries_Groups))
dim(d2)
cpm_d2<-cpm(d2, normalized.lib.sizes = TRUE)
keep_filter<-rowSums(cpm_d2>1)>=4
cpm_d2.2<-cpm_d2[keep_filter,]
dim(cpm_d2.2) #16699
cpm_d2_data<-data.frame(cpm_d2.2)
cpm_d2_data$mean<-rowMeans(cpm_d2_data)
cpm_d2_data2<-cpm_d2_data[,1:48]/cpm_d2_data$mean
tail(cpm_d2_data2,20)
cpm_d2_data2<-log2(cpm_d2_data2+1)
head(cpm_d2_data2,30)
pca = prcomp(t(cpm_d2_data2))
summary(pca)
pca$x

pca_genes<-pca$x
pca_genes_dataframe<-as.data.frame(pca_genes)
#change the conditions
Timeseries_Groups1<-c("a_-2hr","a_-2hr","a_-2hr","a_-2hr","b_0hr","b_0hr","b_0hr","b_0hr","c_2hr","c_2hr","c_2hr","c_2hr","d_4hr","d_4hr","d_4hr","d_4hr","e_6hr","e_6hr","e_6hr","e_6hr","f_9hr","f_9hr","f_9hr","f_9hr","g_12hr","g_12hr","g_12hr","g_12hr","h_1d","h_1d","h_1d","h_1d","i_2d","i_2d","i_2d","i_2d","j_4d","j_4d","j_4d","j_4d","k_8d","k_8d","k_8d","k_8d","l_12d","l_12d","l_12d","l_12d")
Timeseries_Groups2<-c("-2hr","-2hr","-2hr","-2hr","0hr","0hr","0hr","0hr","2hr","2hr","2hr","2hr","4hr","4hr","4hr","4hr","6hr","6hr","6hr","6hr","9hr","9hr","9hr","9hr","12hr","12hr","12hr","12hr","1d","1d","1d","1d","2d","2d","2d","2d","4d","4d","4d","4d","8d","8d","8d","8d","12d","12d","12d","12d")

conditions<-Timeseries_Groups2
pca_genes_dataframe<-data.frame(conditions,pca_genes_dataframe)
head(pca_genes_dataframe)

#Figure 1D
#graph PCA of filtered genes: 16699 genes that have cpm>1 in at least 4 libraries
plot3<- pca_genes_dataframe %>%
  mutate(conditions= fct_relevel(conditions,"-2hr","0hr","2hr","4hr","6hr","9hr","12hr","1d","2d","4d","8d","12d")) %>%
ggplot(aes(x=PC1,y=PC2,colour=conditions))+
  geom_point(size=1.5)+
  labs(x="PC1 (42.7% of variance)",y="PC2 (14.2% of variance)")+
  theme_classic(base_size = 10)+labs(tag="D")+
  theme(legend.position = "none")+theme(aspect.ratio = 0.8)

plot3

#Differential expression analysis for the time series data (48 libraries, 4 bio replicates for 12 time point)
#differential expression - one-way ANOVA
Timeseries_filtered2<-Timeseries_filtered[rownames(cpm_d2.2),]
y<-DGEList(counts=Timeseries_filtered2)
head(Timeseries_filtered2)
anova_condition<-data.frame(row.names = c("-2hr_1","-2hr_5","-2hr_6","-2hr_8","0hr_1","0hr_5","0hr_6","0hr_8","2hr_1","2hr_5","2hr_6","2hr_8","4hr_1","4hr_5","4hr_6","4hr_8","6hr_1","6hr_5","6hr_6","6hr_8","9hr_1","9hr_5","9hr_6","9hr_8","12hr_1","12hr_5","12hr_6","12hr_8","1d_1","1d_5","1d_6","1d_8","2d_1","2d_5","2d_6","2d_8","4d_1","4d_5","4d_6","4d_8","8d_1","8d_5","8d_6","8d_8","12d_1","12d_5","12d_6","12d_8"))
anova_Group<-Timeseries_Groups
anova_Group
cbind(anova_condition,Group=anova_Group)
design<-model.matrix(~anova_Group)
design

y<-calcNormFactors(y)
y<-estimateDisp(y,design)
plotBCV(y)
fit<-glmFit(y,design)
colnames(fit)

anov <- glmLRT(fit, coef=2:12)
anov_sort<-topTags(anov,n=nrow(anov$table))$table
head(anov_sort)
#anov_top1<-rownames(anov_sort[1:500,])
#tail(anov_sort[1:500,])
anov_top2<-rownames(anov_sort)[anov_sort$FDR<=0.05] #differentially expressed
anov_top3<-rownames(anov_sort)[anov_sort$FDR<=0.000000000000000000000000000001] #used in clustering
#The full data from the time series analysis is under the TimeSeries_DiffExpResults tab of Supplementary Figure 1


#Figure 1F (Figure 1E comes next)
#Pull in the TimeSeries_DiffExpResults tab of Supplementary Figure 1
#make a plot showing FDR and number of genes DE at that cutoff
StarvationTimeSeriesResults<-read.csv("StarvationTimeSeries_v2_resultspage.csv",header = T)
StarveTimeSeriesResults_FDRonly<-StarvationTimeSeriesResults[,c(1,66)]

new_df<-c()
for (k in seq(-300,100,0.1)){
  n<-10^k
  StarveSubsetFDR<-subset(StarveTimeSeriesResults_FDRonly,FDR < n)
  #print(nrow(StarveSubsetFDR))
  new_df<-rbind(new_df,c(n,nrow(StarveSubsetFDR)))
  
}

new_df<-as.data.frame(new_df)
colnames(new_df)<-c("FDR_cutoff","Num_of_DE_genes")
plot4<-ggplot(new_df,aes(x=log10(FDR_cutoff),y=Num_of_DE_genes))+
  geom_point(size=0.5)+
  geom_hline(yintercept = 14304,color="blue",size=1)+
  geom_hline(yintercept = 6027,color="red",size=1)+
  theme_classic(base_size = 10)+theme(aspect.ratio = 1)+
  labs(x="log10(FDR)",y="Differentially expressed genes")+labs(tag="F")

plot4


#Figure 1E
#Timeseries_Counts<-read.csv("counts_StarvationTimeSeries.csv",header = TRUE,row.names="sequence")
#Maxwell_genes<-read.csv("Maxwell_Kaplan_RNAseqdata.csv",header = T,row.names = "gene_id")
Timeseries_filtered<-Timeseries_Counts[rownames(Maxwell_genes),]

Timeseries_Groups<-c("-2hr","-2hr","-2hr","-2hr","0hr","0hr","0hr","0hr","2hr","2hr","2hr","2hr","4hr","4hr","4hr","4hr","6hr","6hr","6hr","6hr","9hr","9hr","9hr","9hr","12hr","12hr","12hr","12hr","1d","1d","1d","1d","2d","2d","2d","2d","4d","4d","4d","4d","8d","8d","8d","8d","12d","12d","12d","12d")
d2<-DGEList(counts=Timeseries_filtered,group=factor(Timeseries_Groups))
cpm_d2<-cpm(d2, normalized.lib.sizes = TRUE)
keep_filter<-rowSums(cpm_d2>1)>=4
cpm_d2.2<-cpm_d2[keep_filter,]
cpm_d2_data<-data.frame(cpm_d2.2)
cpm_d2_data$mean<-rowMeans(cpm_d2_data)
cpm_d2_data2<-cpm_d2_data[,1:48]/cpm_d2_data$mean
cpm_d2_data2<-log2(cpm_d2_data2+1)

head(cpm_d2_data2)
genes_of_interest_dataframe<-as.data.frame(cpm_d2_data2)
#remove NA rows if gene is not included 
row.has.na <- apply(genes_of_interest_dataframe, 1, function(x){any(is.na(x))})
genes_of_interest_dataframe<-genes_of_interest_dataframe[!row.has.na,]
#change the rownames so include gene symbol and not just sequence name, if there is one
genes_of_interest_dataframe$sequence<-rownames(genes_of_interest_dataframe)
genes_of_interest_dataframe<-melt(genes_of_interest_dataframe)
genes_of_interest_dataframe<-genes_of_interest_dataframe[order(genes_of_interest_dataframe$sequence),]
genes_of_interest_dataframe$groups<-Timeseries_Groups
gene_averages<-aggregate(genes_of_interest_dataframe$value,list(variable=genes_of_interest_dataframe$groups,genes=genes_of_interest_dataframe$sequence),mean)
gene_averages<-unique(gene_averages)

gene_averages_2<-dcast(gene_averages,genes ~ variable, value='x')
rownames(gene_averages_2)<-gene_averages_2$genes
gene_averages_2<-gene_averages_2[,-1]
head(gene_averages_2)
euclidean_genes<-round(dist(t(gene_averages_2)),digits = 1)
df_euclidean<-melt(as.matrix(euclidean_genes))

agg_euclid<-aggregate(df_euclidean$value,list(df_euclidean$Var1),mean)


#now try to figure out the average pairwise TIME difference between each of these
timepoints<-c(-2,0,2,4,6,9,12,24,48,96,192,288)
timepoints2<-c(-2,0,2,4,6,9,12,24,48,96,192,288)
timepoints3<-as.data.frame(cbind(timepoints,timepoints2))

new_tp<-c()
for (i in 1:12){
  new_subtract<-c()
  for (j in 1:12){
    a<-abs(timepoints3[i,1]-timepoints3[j,2])
    new_subtract<-cbind(new_subtract,a)
  }
  new_tp<-rbind(new_tp,new_subtract)
}
new_tp
new_tp<-as.data.frame(new_tp)
new_tp$mean<-rowMeans(new_tp)
new_tp$timepoints<-c("-2hr","0hr","2hr","4hr","6hr","9hr","12hr","1d","2d","4d","8d","12d")
euc_dist_time<-merge(new_tp,agg_euclid,by.x="timepoints",by.y = "Group.1")
euc_dist_time$rate<-euc_dist_time$x/euc_dist_time$mean
euc_dist_time2<-euc_dist_time[,c("timepoints","mean","x","rate")]
euc_dist_time2



plot5<- euc_dist_time2 %>%
  mutate(timepoints= fct_relevel(timepoints,"-2hr","0hr","2hr","4hr","6hr","9hr","12hr","1d","2d","4d","8d","12d")) %>%
  ggplot(aes(x=timepoints,y=x))+
  #geom_point(aes(x=timepoints,y=mean,color="Time (hrs)"),size=3,alpha=0.8)+
  geom_point(size=3,alpha=0.5)+theme_classic(base_size = 10)+theme(aspect.ratio = 1)+
  labs(x="Time point",y="Average pairwise difference")+labs(tag="F")

plot5

plot6<- euc_dist_time2 %>%
  mutate(timepoints= fct_relevel(timepoints,"-2hr","0hr","2hr","4hr","6hr","9hr","12hr","1d","2d","4d","8d","12d")) %>%
  ggplot(aes(x=timepoints,y=rate))+
  #geom_bar(inherit.aes = FALSE,data=survival1,aes(x=day2,y=prop_alive),stat="identity",width = 0.5)+
  geom_point(size=3,alpha=0.5)+theme_classic(base_size = 10)+theme(aspect.ratio = 1)+
  labs(x="Days of L1 arrest",y="Avg pairwise Euclidean distance / avg pairwise time (hr)")+labs(tag="G")

plot6

#pull the df_euclidean dataframe from before plot5 to only pull out the adjacent time points
df_euclidean1<-df_euclidean[c(13,35,48,52,65,74,90,103,117,128,142),]
df_euclidean1
df_euclidean1$hours_diff<-c(2,96,3,12,24,2,48,2,2,96,3)
df_euclidean1$norm_euclid<- df_euclidean1$value / df_euclidean1$hours_diff
df_euclidean1$hours<-c(0,288,12,24,48,2,96,4,6,192,9)

plot7<- df_euclidean1 %>%
  mutate(Var2= fct_relevel(Var2,"0hr","2hr","4hr","6hr","9hr","12hr","1d","2d","4d","8d","12d")) %>%
  ggplot(aes(x=Var2,y=norm_euclid))+
  geom_point(size=3,alpha=0.5)+theme_classic(base_size = 10)+theme(aspect.ratio = 1)+
  labs(x="Time point",y="Normalized adjacent euclidian distance")+labs(tag="F")

plot7

log(2.71)
plot8<- ggplot(df_euclidean1,aes(x=hours,y=log2(norm_euclid)))+
  geom_point(size=1.5,alpha=0.5)+theme_classic(base_size = 10)+theme(aspect.ratio = 1)+
  labs(x="Hours of L1 arrest",y="Rate of change (Euclidean distance / time)")+labs(tag="F")

plot8
head(cpm_d2_data2)

new_df2<-c()
for (i in 1:11){
  euclidean_genes2<-round(dist(t(cpm_d2_data2[,(1+4*(i-1)):(1+4*(i-1)+7)])),digits = 1)
  euclidean_genes2
  df_euclidean2<-melt(as.matrix(euclidean_genes2))
  df_euc_subset<-df_euclidean2[c(5,6,7,8,13,14,15,16,21,22,23,24,29,30,31,32),]
  new_df2<-rbind(new_df2,df_euc_subset)
}

new_df2
time_points<-c(0,2,4,6,9,12,24,48,96,192,288)
time_diff<-c(2,2,2,2,3,3,12,24,48,96,96)
new_df2$timepoint<-rep(time_points,each=16)
new_df2$time_diff<-rep(time_diff,each=16)
new_df2$norm_dist<-new_df2$value / new_df2$time_diff
library(Hmisc)

#This is the plot that is Figure 1E
plot9<-ggplot(new_df2,aes(x=timepoint,y=norm_dist))+
  #stat_summary(geom="ribbon",fun.data=mean_cl_normal, fun.args=list(conf.int=0.95),fill="gray")+  
  geom_jitter(alpha=0.2,size=1.5)+
  geom_abline(slope = 0,intercept = 0)+
  theme_classic(base_size = 10)+#ylim(0,1)+
  theme(aspect.ratio = 1)+labs(tag = "E")+
  labs(x="Hours of L1 arrest",y="Pairwise Euclidean distance / time for adjacent time points")

plot9


#All plots for Figure 1
grid.arrange(plot1,plot2,plot3,plot9,plot4,nrow=2)

grid.arrange(
  grobs = list(plot1,plot2,plot3,plot9,plot4),
  widths = c(1,1,1,1,1),
  layout_matrix = rbind(c(NA,NA,1,2),
                        c(NA,3,4, 5))
)

